function [theta0_opt_deg, eta_opt, sweepTbl] = OptimizeBasedOnGivenFEKORadPattern(ffeFile, opts)
% maximize_eta_ap_ffe_noNorm
%   Sweep theta0 to maximize eta_ap using FEKO .ffe pattern WITHOUT NORMALIZATION.
%
% Inputs:
%   ffeFile : path to FEKO far-field (.ffe, Format 8)
%   opts    : optional struct
%       .phiMethod        : 'phi0' (default) | 'phi90' | 'avg' | 'max'
%       .theta0_min_deg   : default 5
%       .theta0_max_deg   : default = min(85, max theta in file)
%       .theta0_step_deg  : default 0.1
%       .plot             : true/false (default true)
%
% Outputs:
%   theta0_opt_deg : theta0 (deg) at which eta_ap is maximal over the sweep
%   eta_opt        : maximal eta_ap value
%   sweepTbl       : table with [theta0_deg, eta_ap]

    if nargin < 2, opts = struct(); end
    phiMethod       = getOpt(opts,'phiMethod','phi0');
    th0_min         = getOpt(opts,'theta0_min_deg',5);
    th0_step        = getOpt(opts,'theta0_step_deg',0.1);
    doPlot          = getOpt(opts,'plot',true);

    % --- Read .ffe ---
    [theta_deg_all, phi_deg_all, Dtot_all] = readFFE_local(ffeFile);

    % Aggregate to 1-D Gf(theta) (NO normalization)
    [thetas_deg, Gf_theta] = build_Gtheta_noNorm(theta_deg_all, phi_deg_all, Dtot_all, phiMethod);

    % Determine sweep range
    th0_max_default = min(85, max(thetas_deg));
    th0_max         = getOpt(opts,'theta0_max_deg', th0_max_default);
    theta0_vec      = (th0_min:th0_step:th0_max).';

    % Sweep theta0
    eta_vals = arrayfun(@(t0) eta_from_Gtheta_noNorm(thetas_deg, Gf_theta, t0), theta0_vec);

    % Pick optimum
    [eta_opt, idx]  = max(eta_vals);
    theta0_opt_deg  = theta0_vec(idx);
    sweepTbl        = table(theta0_vec, eta_vals, 'VariableNames', {'theta0_deg','eta_ap'});

    % Plots
    if doPlot
        figure; plot(thetas_deg, Gf_theta, 'LineWidth', 1.6); grid on;
        xlabel('\theta (deg)'); ylabel('G_f(\theta) (no normalization)');
        title(sprintf('Feed pattern from %s  [%s]', ffeFile, upper(phiMethod)));

        figure; plot(theta0_vec, eta_vals, 'LineWidth', 1.8); grid on; hold on;
        plot(theta0_opt_deg, eta_opt, 'o', 'MarkerSize', 7);
        xlabel('\theta_0 (deg)'); ylabel('\eta_{ap} (no normalization)');
        title('Aperture “efficiency” vs. \theta_0 (no normalization)');
        legend('sweep','maximum','Location','best');
    end
end

%==================== helpers ====================

function val = getOpt(s, field, default)
    if isfield(s, field) && ~isempty(s.(field)), val = s.(field); else, val = default; end
end

function [theta_deg, phi_deg, Dtot] = readFFE_local(ffeFile)
    fid = fopen(ffeFile,'r');
    assert(fid>0, 'Cannot open %s', ffeFile);

    headerFound = false; colNames = {}; data = [];

    while true
        line = fgetl(fid);
        if ~ischar(line), break; end
        if ~headerFound && startsWith(strtrim(line),'#') && contains(line,'"')
            toks = regexp(line,'"(.*?)"','tokens');
            colNames = cellfun(@(c)c{1}, toks, 'UniformOutput', false);
            headerFound = true;
            continue
        end
        if headerFound && ~startsWith(strtrim(line),'#') && ~isempty(strtrim(line))
            nums = textscan(line,'%f','Delimiter',' \t','MultipleDelimsAsOne',true);
            if ~isempty(nums{1})
                data(end+1,1:numel(nums{1})) = nums{1}.'; %#ok<AGROW>
            end
        end
    end
    fclose(fid);
    assert(~isempty(colNames) && ~isempty(data), 'Failed to parse .ffe');

    imap = containers.Map(); for k=1:numel(colNames), imap(colNames{k})=k; end
    needed = {'Theta','Phi','Directivity(Total)'};
    for k=1:numel(needed)
        assert(isKey(imap,needed{k}), 'Missing column "%s"', needed{k});
    end

    theta_deg = data(:, imap('Theta'));
    phi_deg   = data(:, imap('Phi'));
    Dtot      = data(:, imap('Directivity(Total)'));
end

function [thetas, Gtheta] = build_Gtheta_noNorm(theta_all, phi_all, Dtot_all, method)
    thetas = unique(theta_all);
    Gtheta = nan(size(thetas));

    for k = 1:numel(thetas)
        mask_th = (theta_all == thetas(k));
        switch lower(method)
            case 'phi0'
                tol = 1e-6;
                mask = mask_th & (abs(phi_all) <= tol);
                vals = Dtot_all(mask);
            case 'phi90'
                tol = 1e-6;
                mask = mask_th & (abs(phi_all - 90) <= tol);
                vals = Dtot_all(mask);
            case 'avg'
                vals = Dtot_all(mask_th);
            case 'max'
                vals = Dtot_all(mask_th);
            otherwise
                error('Unknown phiMethod: %s', method);
        end

        if isempty(vals)
            Gtheta(k) = NaN;
        else
            % If it looks like dBi, convert to linear; else keep as-is.
            if max(vals) <= 100
                linVals = 10.^(vals/10);   % dBi -> linear
            else
                linVals = vals;            % already linear
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

function eta = eta_from_Gtheta_noNorm(theta_deg, Gtheta, theta0_deg)
    % Interpolate Gf over [0, theta0] and integrate
    theta0_deg = min(theta0_deg, max(theta_deg));
    thFine_deg = (0:0.01:theta0_deg).';
    GFine = interp1(theta_deg, Gtheta, thFine_deg, 'pchip', 'extrap');
    GFine(GFine<0) = 0;

    thFine_rad = deg2rad(thFine_deg);
    integrand  = sqrt(GFine) .* tan(thFine_rad/2);
    I          = trapz(thFine_rad, integrand);
    eta        = (cot(deg2rad(theta0_deg)/2)^2) * (I^2);
end

%% Use like this
%[theta0_opt_deg, eta_opt, sweepTbl] = OptimizeBasedOnGivenFEKORadPattern( ...
%    'PyramidalCustomHornV0_SingleFreq_StableVersion_FarField1.ffe', ...
%    struct('phiMethod','phi0','theta0_min_deg',5,'theta0_max_deg',70,'theta0_step_deg',0.05,'plot',true));

%fprintf('Max at theta0 = %.3f deg, eta_ap = %.6f (no normalization)\n', theta0_opt_deg, eta_opt);
