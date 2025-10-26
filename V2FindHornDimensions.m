function O = V2FindHornDimensions(G0_in, lambda, a_throat, b_throat, bounds_mm)
% Choose the (a1,b1) pair (with a1*b1 fixed by G0, lambda) that minimizes
% the total phase-error measure J = Pe^2 + Ph^2, under:
%   rho_e = a1/sqrt(3),  rho_h = b1/sqrt(2).
%
% Inputs:
%   G0_in     : desired gain (dBi or linear). If >100, treated as dBi.
%   lambda    : wavelength [m]
%   a_throat  : waveguide broad dimension at throat [m]
%   b_throat  : waveguide narrow dimension at throat [m]
%   bounds_mm : [a1_min_mm, a1_max_mm] search bounds for a1 (aperture)
%
% Output O: struct with the chosen geometry and phase metrics.

    if G0_in > 100
        G0 = 10^(G0_in/10);
    else
        G0 = G0_in;
    end
    C = (G0 * lambda^2) / (2*pi);

    a1_min = max(a_throat*1.05, bounds_mm(1)*1e-3);  % keep > throat
    a1_max = max(a1_min*1.2, bounds_mm(2)*1e-3);

    % Objective in terms of a1 only (b1 = C/a1)
    function J = objective(a1)
        b1   = C ./ a1;
        if any(b1 <= b_throat) || any(a1 <= a_throat)
            J = inf; return;
        end
        rho_e = a1./sqrt(3);
        rho_h = b1./sqrt(2);
        Pe    = (a1.^2 - a_throat^2) ./ (8*rho_e);
        Ph    = (b1.^2 - b_throat^2) ./ (8*rho_h);
        J     = Pe.^2 + Ph.^2;
    end

    % 1D bounded search (fminbnd)
    opts = optimset('TolX',1e-9,'Display','off');
    [a1_opt, ~] = fminbnd(@objective, a1_min, a1_max, opts);

    % Recompute full set at optimum
    b1_opt = C / a1_opt;
    rho_e  = a1_opt/sqrt(3);
    rho_h  = b1_opt/sqrt(2);
    Pe     = (a1_opt^2 - a_throat^2)/(8*rho_e);
    Ph     = (b1_opt^2 - b_throat^2)/(8*rho_h);
    Hh_e   = Pe/rho_e;
    Hh_h   = Ph/rho_h;

    O = struct();
    O.G0_lin     = G0;
    O.lambda     = lambda;
    O.C          = C;
    O.a_throat   = a_throat;
    O.b_throat   = b_throat;

    O.a1         = a1_opt;
    O.b1         = b1_opt;
    O.rho_e      = rho_e;
    O.rho_h      = rho_h;
    O.Pe         = Pe;
    O.Ph         = Ph;
    O.H_over_h_e = Hh_e;
    O.H_over_h_h = Hh_h;
    O.H_over_h   = 0.5*(Hh_e + Hh_h);
    O.mismatch   = Hh_e - Hh_h;

    % quick console dump
    fprintf('\n--- Min-phase solution ---\n');
    fprintf('a1 = %.3f mm, b1 = %.3f mm\n', 1e3*O.a1, 1e3*O.b1);
    fprintf('rho_e = %.3f mm, rho_h = %.3f mm\n', 1e3*O.rho_e, 1e3*O.rho_h);
    fprintf('Pe = %.4g m, Ph = %.4g m\n', O.Pe, O.Ph);
    fprintf('H/h: E-plane = %.5f, H-plane = %.5f (avg = %.5f)\n', O.H_over_h_e, O.H_over_h_h, O.H_over_h);
end
