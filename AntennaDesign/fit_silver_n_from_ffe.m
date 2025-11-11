function [n_hat, A_hat, stats, fitTbl] = fit_silver_n_from_ffe(ffeFile, opts)
% fit_silver_n_from_ffe
%   Fit Silver's cos^n feed model to a FEKO .ffe far-field pattern
%   WITHOUT normalizing the FEKO data.
%
% Model: G_model(theta) = A * cosd(theta).^n,    0 <= theta <= 90 deg
%
% [n_hat, A_hat, stats, fitTbl] = fit_silver_n_from_ffe(ffeFile, opts)
%
% INPUTS
%   ffeFile : path to FEKO .ffe (Format 8) file
%   opts    : struct with optional fields:
%       .phiMethod        : 'phi0' (default) | 'phi90' | 'avg' | 'max'
%       .theta_max_deg    : use data up to this angle (default: 90)
%       .theta_min_deg    : (default: 0)
%       .n_init           : initial n for optimizer (default: 10)
%       .n_bounds         : [nmin nmax] bounds for search (default: [0 200])
%       .plot             : true/false to plot fit (default: true)
%
% OUTPUTS
%   n_hat  : best-fit exponent
%   A_hat  : best-fit scale (same units as the .ffe directivity/gain)
%   stats  : struct with fields SSE, RMSE, R2, G0_from_Balanis (=2*(n_hat+1))
%   fitTbl : table with columns [theta_deg, G_data, G_fit]
%
% NOTES
%   * The code detects if "Directivity(Total)" is in dBi; if so it converts
%     to linear; otherwise it uses values as-is. No normalization is applied.
%   * You can compare A_hat against the closed-form Silver normalization
%     constant G0^(n)=2(n+1) from Balanis (15-58) if desired.

    if nargin < 2, opts = struct(); end
    phiMethod     = getOpt(opts,'phiMethod','phi0');
    theta_max     = getOpt(opts,'theta_max_deg',90);
    theta_min     = getOpt(opts,'theta_min_deg',0);
    n_init        = getOpt(opts,'n_init',10);
    n_bounds      = getOpt(opts,'n_bounds',[0 200]);
    doPlot        = getOpt(opts,'plot',true);

    % ---- Read the .ffe and aggregate to a 1-D G(theta) ----
    [theta_all, phi_all, Dtot_all] = readFFE_local(ffeFile);
    [theta_deg, Gtheta] = make_cut(theta_all, phi_all, Dtot_all, phiMethod);

    % Keep only theta in [theta_min, min(theta_max, 90)]
    theta_cap = min(90, theta_max);
    use = (theta_deg >= theta_min) & (theta_deg <= theta_cap);
    theta_deg = theta_deg(use);
    Gtheta    = Gtheta(use);

    % Guard
    assert(numel(theta_deg) >= 5, 'Too few samples in requested range.');

    % ---- Fit n (A solved in closed form for every n) ----
    obj = @(n) sse_for_n(n, theta_deg, Gtheta);
    % Constrain n to [nmin,nmax] via fminbnd
    nmin = n_bounds(1); nmax = n_bounds(2);
    n_hat = fminbnd(obj, nmin, nmax);

    % Recover the corresponding A(n_hat)
    x = cosd(theta_deg).^n_hat;
    A_hat = (x' * Gtheta) / (x' * x);

    % Fitted curve and stats
    G_fit = A_hat * x;
    resid = Gtheta - G_fit;
    SSE   = sum(resid.^2);
    RMSE  = sqrt(mean(resid.^2));
    SSyy  = sum( (Gtheta - mean(Gtheta)).^2 );
    R2    = 1 - SSE/max(SSyy,eps);

    stats = struct( ...
        'SSE',  SSE, ...
        'RMSE', RMSE, ...
        'R2',   R2, ...
        'G0_from_Balanis', 2*(n_hat+1) );   % (15-58), for reference

    fitTbl = table(theta_deg(:), Gtheta(:), G_fit(:), ...
        'VariableNames', {'theta_deg','G_data','G_fit'});

    if doPlot
        figure; plot(theta_deg, Gtheta, 'o','MarkerSize',4); hold on; grid on;
        plot(theta_deg, G_fit, 'LineWidth',1.8);
        xlabel('\theta (deg)'); ylabel('G(\theta) (linear units, no norm)');
        title(sprintf('Silver fit:  G(\\theta)=A cos^n(\\theta),  n=%.3f,  A=%.3g,  R^2=%.4f', ...
              n_hat, A_hat, R2));
        legend('FEKO data','Silver fit','Location','best');
    end
end

%==================== helpers ====================

function sse = sse_for_n(n, theta_deg, G)
    % For a given n, solve A in least squares and return SSE
    x = cosd(theta_deg).^n;
    % If any cos becomes negative (shouldn't in [0,90]), zero it
    x(x<0) = 0;
    A = (x' * G) / (x' * x);
    Ghat = A * x;
    r = G - Ghat;
    sse = sum(r.^2);
end

function val = getOpt(s, field, default)
    if isfield(s, field) && ~isempty(s.(field)), val = s.(field); else, val = default; end
end

function [theta_deg, phi_deg, Dtot] = readFFE_local(fp)
    fid = fopen(fp,'r'); assert(fid>0, 'Cannot open %s', fp);
    gotHdr = false; cols = {}; data = [];
    while true
        t = fgetl(fid);
        if ~ischar(t), break; end
        if ~gotHdr && startsWith(strtrim(t),'#') && contains(t,'"')
            toks = regexp(t,'"(.*?)"','tokens');
            cols = cellfun(@(c)c{1}, toks, 'UniformOutput', false);
            gotHdr = true; continue
        end
        if gotHdr && ~startsWith(strtrim(t),'#') && ~isempty(strtrim(t))
            v = textscan(t,'%f','Delimiter',' \t','MultipleDelimsAsOne',true);
            if ~isempty(v{1}), data(end+1,1:numel(v{1})) = v{1}.'; end %#ok<AGROW>
        end
    end
    fclose(fid);
    assert(~isempty(cols) && ~isempty(data), 'Failed to parse .ffe');

    ix = containers.Map(); for k=1:numel(cols), ix(cols{k})=k; end
    needed = {'Theta','Phi','Directivity(Total)'};
    for k=1:numel(needed)
        assert(isKey(ix,needed{k}), 'Missing column "%s"', needed{k});
    end
    theta_deg = data(:, ix('Theta'));
    phi_deg   = data(:, ix('Phi'));
    Dtot      = data(:, ix('Directivity(Total)'));
end

function [thetas, Gtheta] = make_cut(theta_all, phi_all, D_all, method)
    thetas = unique(theta_all);
    Gtheta = nan(size(thetas));
    tol = 1e-6;

    for k = 1:numel(thetas)
        mth = (theta_all == thetas(k));
        switch lower(method)
            case 'phi0'
                vals = D_all(mth & abs(phi_all) <= tol);
            case 'phi90'
                vals = D_all(mth & abs(phi_all-90) <= tol);
            case 'avg'
                vals = D_all(mth);
            case 'max'
                vals = D_all(mth);
            otherwise
                error('Unknown phiMethod: %s', method);
        end

        if isempty(vals)
            Gtheta(k) = NaN;
        else
            % If it looks like dBi, convert to linear; otherwise keep.
            if max(vals) <= 100
                linVals = 10.^(vals/10);
            else
                linVals = vals;
            end
            switch lower(method)
                case 'avg', Gtheta(k) = mean(linVals,'omitnan');
                case 'max', Gtheta(k) = max(linVals);
                otherwise,  Gtheta(k) = linVals(1); % single cut
            end
        end
    end

    good = ~isnan(Gtheta);
    thetas = thetas(good);
    Gtheta = Gtheta(good);
end


%% Run this
%[n_hat, A_hat, stats, fitTbl] = fit_silver_n_from_ffe( ...
%    'PyramidalCustomHornV0_SingleFreq_StableVersion_FarField1.ffe', ...
%    struct('phiMethod','phi0','theta_max_deg',90,'plot',true));

%fprintf('Best n = %.3f,  A = %.6g\n', n_hat, A_hat);
%fprintf('R^2 = %.4f,  SSE = %.3g,  RMSE = %.3g\n', stats.R2, stats.SSE, stats.RMSE);
%fprintf('Silver normalization constant G0^(n)=2(n+1)= %.3f (for reference)\n', stats.G0_from_Balanis);

