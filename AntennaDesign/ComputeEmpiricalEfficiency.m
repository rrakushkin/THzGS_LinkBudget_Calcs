function eta_ap = ComputeEmpiricalEfficiency(ffeFile, theta0_deg)
% eta_ap_from_ffe  Compute aperture efficiency using Balanis (15-55)
%   eta_ap = cot^2(theta0/2) * | ∫_0^{theta0} sqrt(Gf(θ)) * tan(θ/2) dθ |^2
%
% INPUTS
%   ffeFile     : path to FEKO .ffe export (Format 8)
%   theta0_deg  : edge angle in degrees (e.g., 46.1)
%
% ASSUMPTIONS
%   - File contains columns labeled (typical FEKO format 8):
%       "Theta" "Phi" "Re(Etheta)" "Im(Etheta)" "Re(Ephi)" "Im(Ephi)"
%       "|Etheta|" "|Ephi|" "Directivity(Theta)" "Directivity(Phi)" "Directivity(Total)"
%   - We will use the φ = 0° cut and the "Directivity(Total)" column.
%   - If directivity appears in dBi, we convert to linear. If already linear,
%     the heuristic below keeps it unchanged.
%
% OUTPUT
%   eta_ap      : aperture efficiency in [0,1]

    arguments
        ffeFile        (1,:) char
        theta0_deg     (1,1) double {mustBePositive,mustBeReal}
    end

    %% --- Read file ---
    fid = fopen(ffeFile,'r');
    assert(fid>0, 'Could not open file: %s', ffeFile);

    headerDone = false;
    colNames   = {};
    rows = [];

    while true
        t = fgetl(fid);
        if ~ischar(t), break; end

        if startsWith(strtrim(t),'#No. of Header Lines')
            % nothing to do
        elseif startsWith(strtrim(t),'#') && contains(t,'"')
            % This is the single header line with column names
            % Example:
            % # "Theta" "Phi" ... "Directivity(Total)"
            rawNames = regexp(t,'"(.*?)"','tokens');
            colNames = cellfun(@(c)c{1}, rawNames, 'UniformOutput', false);
            headerDone = true;
        elseif headerDone
            % Data line (space-separated in scientific notation)
            nums = textscan(t,'%f','Delimiter',' \t','MultipleDelimsAsOne',true);
            v = nums{1};
            if ~isempty(v)
                rows(end+1,1:numel(v)) = v.'; %#ok<AGROW>
            end
        end
    end
    fclose(fid);

    if isempty(colNames) || isempty(rows)
        error('Failed to parse .ffe file (no column names or data).');
    end

    % Build an index map
    idx = containers.Map();
    for k = 1:numel(colNames), idx(colNames{k}) = k; end

    needCols = {"Theta","Phi","Directivity(Total)"};
    for k = 1:numel(needCols)
        assert(isKey(idx,needCols{k}), 'Missing column "%s" in file.', needCols{k});
    end

    theta_deg_all = rows(:, idx('Theta'));
    phi_deg_all   = rows(:, idx('Phi'));
    Dtot_all      = rows(:, idx('Directivity(Total)'));

    %% --- Take φ = 0° cut and 0 ≤ θ ≤ θ0 ---
    % (Use a small tolerance to catch 359.999.../0.000... variants.)
    tol = 1e-6;
    cut = (abs(phi_deg_all) <= tol) & (theta_deg_all >= 0) & (theta_deg_all <= theta0_deg);
    theta_deg = theta_deg_all(cut);
    Dtot      = Dtot_all(cut);

    % Sort and deduplicate on theta
    [theta_deg, u] = unique(theta_deg,'stable');
    Dtot = Dtot(u);

    if isempty(theta_deg)
        error('No samples found for φ=0° within 0..θ0.');
    end

    %% --- Convert Directivity(Total) to linear if needed, then normalize
    % Heuristic: if values look like dBi (typ. < 100), convert; if already linear
    % (typ. <= few tens), this still works correctly because linear directivity at boresight
    % for a feed is usually << 100; but we also try to detect > 100 as already-linear anomaly.
    if max(Dtot) <= 100  % likely dBi
        G_lin = 10.^(Dtot/10);
    else
        G_lin = Dtot;     % already linear
    end

    % Normalize to unity at boresight using the maximum within [0, θ0]
    G_norm = G_lin ;%/ max(G_lin);

    %% --- Interpolate to a fine grid and integrate (15-55)
    dth = 0.01; % degrees
    thFine_deg = (0:dth:min(theta0_deg, max(theta_deg))).';
    GFine = interp1(theta_deg, G_norm, thFine_deg, 'pchip', 'extrap');
    GFine(GFine < 0) = 0;  % defensive clamp

    thFine_rad = deg2rad(thFine_deg);
    integrand  = sqrt(GFine) .* tan(thFine_rad/2);

    I = trapz(thFine_rad, integrand);   % ∫_0^{θ0} sqrt(Gf) tan(θ/2) dθ
    eta_ap = (cot(deg2rad(theta0_deg)/2)^2) * (I^2);

    %% --- Optional: quick report
    fprintf('File: %s\n', ffeFile);
    fprintf('phi = 0° cut, theta0 = %.4f deg, samples used: %d\n', theta0_deg, numel(thFine_deg));
    fprintf('eta_ap (Balanis 15-55) = %.6f  (%.2f %%)\n', eta_ap, 100*eta_ap);
end
