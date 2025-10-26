function S = V1FindHornDimensions(G0_in, lambda, a_throat, b_throat, a1_range)
% Enumerate all (a1,b1) pairs satisfying G0 = (2*pi/lambda^2)*a1*b1,
% then compute rho_e,rho_h, Pe,Ph, and H/h from each plane.
%
% Inputs:
%   G0_in     : desired gain (dBi or linear). If > 100, treated as dBi.
%   lambda    : wavelength [m]
%   a_throat  : waveguide broad dimension at throat [m]
%   b_throat  : waveguide narrow dimension at throat [m]
%   a1_range  : vector of aperture a1 candidates [m] (solution space sweep)
%
% Output struct S contains vectors for each quantity across the sweep.

    % --- normalize gain to linear
    if G0_in > 100
        G0 = 10^(G0_in/10);
    else
        G0 = G0_in;
    end

    % --- constant from Step 1
    C = (G0 * lambda^2) / (2*pi);      % so that a1.*b1 = C

    % --- solve b1 for each a1
    a1 = a1_range(:);
    b1 = C ./ a1;

    % keep only physically valid (a1 > a_throat, b1 > b_throat)
    valid = (a1 > a_throat) & (b1 > b_throat);
    a1 = a1(valid);
    b1 = b1(valid);

    % --- Step 2: rho_e, rho_h from your notes
    rho_e = a1./sqrt(3);      % a1 = sqrt(3)*rho_e  -> rho_e = a1/sqrt(3)
    rho_h = b1./sqrt(2);      % b1 = sqrt(2)*rho_h  -> rho_h = b1/sqrt(2)

    % --- Step 3: phase errors and H/h ratios (two ways)
    Pe = (a1.^2 - a_throat^2) ./ (8*rho_e);
    Ph = (b1.^2 - b_throat^2) ./ (8*rho_h);

    H_over_h_e = Pe ./ rho_e;   % = H/h from E-plane expression
    H_over_h_h = Ph ./ rho_h;   % = H/h from H-plane expression
    mismatch   = H_over_h_e - H_over_h_h;

    % --- package results
    S = struct();
    S.G0_lin     = G0;
    S.lambda     = lambda;
    S.C          = C;
    S.a_throat   = a_throat;
    S.b_throat   = b_throat;

    S.a1         = a1;
    S.b1         = b1;
    S.rho_e      = rho_e;
    S.rho_h      = rho_h;
    S.Pe         = Pe;
    S.Ph         = Ph;
    S.H_over_h_e = H_over_h_e;
    S.H_over_h_h = H_over_h_h;
    S.H_over_h   = 0.5*(H_over_h_e + H_over_h_h);  % averaged
    S.mismatch   = mismatch;

    % --- quick plot of solution curve (a1,b1)
    figure; plot(a1*1e3, b1*1e3, 'LineWidth',1.5);
    xlabel('a_1 (mm)'); ylabel('b_1 (mm)');
    title('Solution space: a_1 b_1 = C (fixed gain)');
    grid on;

    % --- optional: plot H/h from each plane to see consistency
    figure; 
    plot(a1*1e3, S.H_over_h_e, 'LineWidth',1.5); hold on;
    plot(a1*1e3, S.H_over_h_h, '--', 'LineWidth',1.5);
    legend('H/h from E-plane','H/h from H-plane','Location','best');
    xlabel('a_1 (mm)'); ylabel('H/h (dimensionless)'); grid on;
    title('Consistency of H/h from both planes');
end
