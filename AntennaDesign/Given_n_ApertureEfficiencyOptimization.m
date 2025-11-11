%% Max-eta_ap for Silver feed (script)
% G_f(θ) = 2*(n+1)*cos^n(θ),  0 ≤ θ ≤ π/2
% η_ap(θ0) = cot^2(θ0/2) * |∫_0^{θ0} sqrt(G_f(θ)) * tan(θ/2) dθ|^2
% All computations in radians.

close all; clear; clc;

%% --- Inputs ---
n                  = 12.02;             % Silver exponent
Npts               = 6001;           % resolution of sweep/grid (odd is nice)
theta0_min_deg     = 1e-6;           % avoid exactly 0°
theta0_max_deg     = 90;             % Silver model valid to 90°
theta0_vec_deg     = linspace(theta0_min_deg, theta0_max_deg, Npts).';  % column

% Convert to radians once
theta0_vec = deg2rad(theta0_vec_deg);

%% --- Define integrand √Gf(θ)*tan(θ/2) ---
integrand = @(th) sqrt(2*(n+1)) .* (cos(th)).^(n/2) .* tan(th/2);

%% --- Efficient cumulative integration over θ (radians) ---
% Build a fine θ grid identical to the sweep grid so the cumulative integral
% evaluated at θ0_i is simply I(i).
I_cum = cumtrapz(theta0_vec, integrand(theta0_vec));  % ∫_0^θ integrand dθ

%% --- η_ap(θ0) on same grid ---
eta_ap_vec = (cot(theta0_vec/2)).^2 .* (I_cum.^2);

%% --- Find maximum ---
[eta_ap_max, idx_max] = max(eta_ap_vec);
theta0_opt_rad        = theta0_vec(idx_max);
theta0_opt_deg        = rad2deg(theta0_opt_rad);

%% --- Report ---
fprintf('-------------------------------------------\n');
fprintf('Aperture Efficiency (Balanis 15-55)\n');
fprintf('Silver model: n = %.6g\n', n);
fprintf('Max η_ap = %.8f\n', eta_ap_max);
fprintf('θ0* = %.6f deg\n', theta0_opt_deg);
fprintf('-------------------------------------------\n');

%% --- Plot ---
figure; plot(theta0_vec_deg, eta_ap_vec, 'LineWidth', 1.6); grid on; hold on;
plot(theta0_opt_deg, eta_ap_max, 'o', 'MarkerSize', 7);
xlabel('\theta_0 (deg)'); ylabel('\eta_{ap}');
title(sprintf('Silver feed, n = %.3g — maximum at \\theta_0 = %.3f^\\circ', n, theta0_opt_deg));
legend('\eta_{ap}(\theta_0)','maximum','Location','best');

%% --- (Optional) value at a specific θ0 (example)
theta0_check_deg = 46.1;
th = deg2rad(linspace(0, theta0_check_deg, 1 + ceil(theta0_check_deg/0.01)));
I  = trapz(th, integrand(th));
eta_at_46p1 = (cot(deg2rad(theta0_check_deg)/2)^2) * (I^2);
fprintf('η_ap at θ0 = %.4f deg: %.8f\n', eta_at_46p1, theta0_check_deg);
