% Approach: Antenna as a whole is driven by the directivity and efficiency.
% First, the directivity is fixed by fixing a diameter dimension. Then a
% reasonable focal length is set. This fixes the f/d ratio of the antenna.
% For a given f/d ratio, there is one unique pyramidal horn geometry which will
% yield the smallest taper and spillover losses. This script calculates this horn geometry.
%% Primary Driving Formulas From Text
% ε_ap = cot^2(θ0/2) * | ∫_0^{θ0} sqrt(Gf(θ')) * tan(θ'/2) dθ' |^2
% With Gf(θ') = (n+1) * (cos θ')^n  (Antenna Theory Analysis and Design, Chapter 8)
% theta0 is the subtended angle (subtended angle measures from the vertix to the edge of dish)
%------------------------------------------------------------------------------------------------------------%

close all; clear; clc;
%% General Antenna Parameters (F/D and or theta0)   
%---set the subtended angle or the f/d ratio, then set f_m or d_m ---% 46.1
[theta0_deg, f_over_d, f_m, d_m] = BasicAntennaGeometry('theta0_deg', 46.1, 'd_m', 1);
theta0 = deg2rad(theta0_deg);     % convert to radians
prefactor = cot(theta0/2)^2;      % cot^2(theta0/2)

%% Define Integrand and Aperture Efficiency
% sqrt(Gf) = sqrt(n+1) * (cos θ')^(n/2)
integrand = @(tp, n) sqrt(2*(n+1)) .* (cos(tp)).^(n/2) .* tan(tp/2);

I_of_n = @(n) arrayfun(@(nn) ...
    integral(@(tp) integrand(tp, nn), 0, theta0, ...
             'RelTol', 1e-10, 'AbsTol', 1e-12), n);
eff_ap = @(n) prefactor .* abs(I_of_n(n)).^2;

%% Sweep n and Evaluate
n_min = 0; n_max = 10; Npts = 601;
n_vec = linspace(n_min, n_max, Npts);
eff_vec = eff_ap(n_vec);

%% Find Maximum Efficiency and Corresponding n
[eff_max, index_max] = max(eff_vec);
n_opt = n_vec(index_max); %n_opt is the n which yielded the maximum efficiency

%% Half-power angle and HPBW at n_opt
% Half-power condition: (cos θ_HP)^n = 1/2  -> θ_HP = acos(2^(-1/n))
theta_HP_rad = acos(2.^(-1 ./ n_opt)); %yields radians result (also half cone only)
HPBW_rad     = 2 * theta_HP_rad; 
theta_HP_deg = rad2deg(theta_HP_rad);
HPBW_deg_opt = rad2deg(HPBW_rad);

%% Display Results
fprintf('-------------------------------------------\n');
fprintf('Aperture Efficiency Analysis\n');
fprintf('-------------------------------------------\n');
fprintf('θ₀ = %.4f degrees\n', theta0_deg);
fprintf('Maximum ε_ap = %.6f\n', eff_max);
fprintf('Optimal n = %.4f\n', n_opt);
fprintf('θ_HP (at n_opt) = %.4f°\n', theta_HP_deg);
fprintf('HPBW  (at n_opt) = %.4f°\n', HPBW_deg_opt);
%% Plot ε_ap vs n
figure;
plot(n_vec, eff_vec, 'LineWidth', 1.8); hold on;
plot(n_opt, eff_max, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
grid on;
xlabel('n'); ylabel('\epsilon_{ap}');
title(sprintf('\\epsilon_{ap} vs n  (\\theta_0 = %.4f^\\circ)', theta0_deg));
text(n_opt + 1, eff_max, sprintf('n = %.2f', n_opt), 'Color', 'r');
xlim([n_min n_max]);

%% Plot G_f(θ) at n = n_opt over [0, θ0]
theta = linspace(0, theta0, 800);           % radians
Gf_opt = 2*(n_opt + 1) .* (cos(theta)).^n_opt;

figure;
plot(rad2deg(theta), Gf_opt, 'LineWidth', 1.8); hold on;
grid on;
xlabel('\theta (degrees)', 'Interpreter', 'tex');
ylabel('G_f(\theta) = (n+1)(cos\theta)^n', 'Interpreter', 'tex');
title(sprintf('Feed Pattern at n = %.3f  over [0, \\theta_0=%.4f^\\circ]', ...
      n_opt, theta0_deg), 'Interpreter', 'tex');
xlim([0 theta0_deg]);

% Mark half-power level and θ_HP on the feed plot (if within plotted span)
Gmax = (n_opt + 1);                 % at θ = 0
yline(Gmax/2, '--k', 'Half-Power Level');
if theta_HP_deg <= theta0_deg
    xline(theta_HP_deg, '--r', sprintf('\\theta_{HP} = %.2f^\\circ', theta_HP_deg), ...
          'LabelOrientation', 'horizontal', 'Interpreter', 'tex');
end

%% Compute and Plot HPBW vs n (global curve) and mark n_opt
n_vals = linspace(0.5, 60, 200);
HPBW_rad_all = 2 * acos(2.^(-1 ./ n_vals));
HPBW_deg_all = rad2deg(HPBW_rad_all);

figure;
plot(n_vals, HPBW_deg_all, 'LineWidth', 1.8); hold on;
plot(n_opt, HPBW_deg_opt, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
grid on;
xlabel('n');
ylabel('HPBW (degrees)');
title('HPBW of Feed Pattern (cos^n model)');
text(n_opt + 1, HPBW_deg_opt, sprintf('n_{opt}=%.2f, HPBW=%.2f^\\circ', ...
     n_opt, HPBW_deg_opt), 'Color', 'r', 'Interpreter', 'tex');
%% Find Horn Dimensions for Determined HPBW
a_m = 1.092e-3;
b_m = 0.546e-3;
freq_Hz = 225e9;
HornGain = 2*(n_opt+1)
[chiVal, rho_eVal, rho_hVal, a1Val, b1Val, P_eVal, P_hVal] = FindHornDimensions(HornGain, a_m, b_m, freq_Hz)

h1 = sqrt(rho_eVal^2 - (b1Val/2)^2)
h2 = sqrt(rho_hVal^2 - (a1Val/2)^2)

H1 = h1*P_eVal/rho_eVal
H2 = h2*P_hVal/rho_hVal

%sol = solve([H/h==P_hVal/rho_1, H/h==P_eVal/rho_2],[H,h]);


fprintf('\n---------------------------------------------');
fprintf('\n---------------------------------------------\n');
fprintf('Horn Dimensions\n');
fprintf('---------------------------------------------\n');
fprintf('a = WG (E-plane)\n');
fprintf('b = WG (H-plane)\n');
fprintf('a1 = Aperture (E-plane)\n');
fprintf('b1 = Aperture (H-plane)\n');
fprintf('rho_e = Cone Slant Length from peak (E-plane)\n');
fprintf('rho_h = Cone Slant Length from peak (H-plane)\n');
fprintf('P_e = Cone Slant Length from a (E-plane)\n');
fprintf('P_h = Cone Slant Length from b (H-plane)\n');
fprintf('***Note P_e and P_h should be equal***\n')
fprintf('---------------------------------------------\n');
fprintf('a: %.4f mm\n', a_m*1e3);
fprintf('b: %.4f mm\n', b_m*1e3);
fprintf('a1: %.4f mm\n', a1Val*1e3);
fprintf('b1: %.4f mm\n', b1Val*1e3);
fprintf('rho_e %.4f mm\n', rho_eVal*1e3);
fprintf('rho_h: %.4f mm\n', rho_hVal*1e3);
fprintf('P_e %.4f mm\n', P_eVal*1e3);
fprintf('P_h: %.4f mm\n', P_hVal*1e3);
fprintf('---------------------------------------------\n');