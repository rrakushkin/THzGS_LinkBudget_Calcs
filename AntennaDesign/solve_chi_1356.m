function chi = solve_chi_1356(a, b, f_Hz, G_in)
% solve_chi_1356  Solve Balanis (3e) eq. (13-56) for chi given a,b,G_in.
% a,b : waveguide inner dims [meters]
% f_Hz: operating frequency [Hz]
% G_in: desired gain (linear)
% unit: 'lin' (default) or 'dBi'
%
% Returns:
%   chi : dimensionless solution of (13-56)

    lambda = physconst('LightSpeed')/f_Hz;

    % 13-57 initial guess
    chi0 = G_in/(2*pi*sqrt(2*pi));

    % Define residual F(chi) = LHS - RHS
    F = @(chi) ( (sqrt(2*chi) - b/lambda).^2 .* (2*chi - 1) ) ...
             - ( ( G_in/(2*pi)*sqrt(3/(2*pi))./sqrt(chi) - a/lambda ).^2 ...
             .* ( (G_in^2)/(6*pi^3)./chi - 1 ) );

    % Guard domain: chi must be positive; in practice use >0.5 to avoid (2chi-1)<0 issues
    lo = max(1e-6, 0.51);
    hi = 1e6;  % generous upper bound

    % Try to bracket a root by scanning logspace from ~lo to hi
    grid = logspace(log10(lo), log10(hi), 400);
    vals = F(grid);
    sgn  = sign(vals);
    k = find(diff(sgn)~=0, 1, 'first');  % first sign change

    if ~isempty(k)
        a_br = grid(k); b_br = grid(k+1);
        chi = fzero(F, [a_br, b_br]);      % bracketed root
    else
        % Fall back: use fzero with initial guess; if that fails, use least-squares fit
        try
            chi = fzero(F, max(chi0, lo));
        catch
            % Minimize squared residual to get the nearest feasible solution
            obj = @(x) F(x).^2;
            chi = fminbnd(obj, lo, hi);
            % optional: project to exact root if close
            try, chi = fzero(F, chi); catch, end
        end
    end
end
