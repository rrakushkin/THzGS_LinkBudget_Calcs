%% Ruze surface-roughness requirement vs dish diameter
% Plots the RMS surface roughness sigma [um] needed to achieve a target
% effective gain at a given frequency, using Ruze's formula.
% G_eff = eta0 * eta_s * (pi*D/lambda)^2, where eta_s = exp(-(4*pi*sigma/lambda)^2)

clear; clc; close all;

%% -------------------- User inputs --------------------
f_GHz       = 225;      % Operating frequency [GHz]
G_eff_dBi   = 68;       % Target effective gain [dBi] (includes all losses)
eta0        = 0.70;     % "Other" efficiency excluding surface (0..1)

D_min       = 0.5;      % [m]
D_max       = 2.5;      % [m]
num_pts     = 600;      % number of samples along D

%% -------------------- Constants ---------------------
c      = physconst('LightSpeed');
lambda = c/(f_GHz*1e9);              % wavelength [m]

%% -------------------- Compute sigma(D) ----------------
D = linspace(D_min, D_max, num_pts); % dish diameters [m]
G_eff_lin = 10^(G_eff_dBi/10);       % linear gain target

% Max gain (with perfect surface, eta_s=1) for each D:
G_max_nosurface = eta0 * (pi*D/lambda).^2;

% Ratio used inside the log; must be <= 1 for real sigma
ratio = G_eff_lin ./ G_max_nosurface;

% Initialize sigma (RMS) [m]; invalid where ratio>1 (target > theoretical max)
sigma_m = nan(size(D));

valid = ratio <= 1 & ratio > 0;
sigma_m(valid) = (lambda/(4*pi)) .* sqrt(-log(ratio(valid)));

% Convert to micrometers for plotting
sigma_um = sigma_m * 1e6;

%% -------------------- Plot ---------------------------
figure('Color','w','Position',[80 80 900 480]); 
plot(D, sigma_um, 'LineWidth', 2); hold on; grid on;

% Shade infeasible region where target gain exceeds perfect-surface limit
infeas = ratio > 1;
if any(infeas)
    yl = ylim;
    x_infeas = D; y1 = yl(1)*ones(size(D)); y2 = yl(2)*ones(size(D));
    area(x_infeas(infeas), y2(infeas), yl(1), 'FaceAlpha', 0.10, 'LineStyle', 'none');
    ylim(yl);
end

xlabel('Dish diameter, D [m]');
ylabel('Required RMS surface roughness, \sigma [\mum]');
title(sprintf('Ruze Requirement: \\sigma vs D  (f = %.0f GHz, G_{eff}=%.1f dBi, \\eta_0=%.2f)', ...
              f_GHz, G_eff_dBi, eta0));

% Add helpful annotations
text(0.01+min(D), 0.92*max(sigma_um(~isnan(sigma_um))), ...
    ['Ruze:  \eta_s = e^{-(4\pi\sigma/\lambda)^2}', newline, ...
     'G_{eff} = \eta_0 \eta_s (\pi D/\lambda)^2'], ...
     'FontSize', 10, 'BackgroundColor', [1 1 1 0.85], 'EdgeColor',[0.7 0.7 0.7]);

if any(infeas)
    txtx = D(find(infeas,1,'first'));
    text(txtx, 0.1*max(sigma_um(~isnan(sigma_um))), ...
        'Target gain exceeds perfect-surface limit (no solution)', ...
        'Color',[0.4 0 0], 'FontWeight','bold');
end

%% -------------------- Optional: print a quick table --------------------
% Show a few example points for convenience
probe_D = [0.75 1.0 1.5 2.0];  % meters
for d = probe_D
    Gmax_no_surf = eta0*(pi*d/lambda)^2;
    r = G_eff_lin/Gmax_no_surf;
    if r<=1 && r>0
        sig = (lambda/(4*pi))*sqrt(-log(r))*1e6;
        fprintf('D = %.2f m  ->  sigma = %.1f um (Ruze)\n', d, sig);
    else
        fprintf('D = %.2f m  ->  infeasible for this G_eff (needs eta_0 or D larger)\n', d);
    end
end

%% -------------------- Notes --------------------
% - eta0 should capture everything EXCEPT surface roughness (e.g., illumination,
%   spillover, blockage, ohmic, polarization, phase, pointing, etc.).
% - If your target G_eff is very aggressive for a given D and eta0, the curve
%   will show an infeasible region (no real sigma). Increase D or eta0, or
%   reduce the target G_eff.
% - As the target G_eff gets easier relative to D and eta0, the required
%   sigma grows (curve rises) because rougher surfaces can still meet the target.


% Compute derivative dσ/dD
R = G_eff_lin ./ (eta0 * (pi*D/lambda).^2);
dsigma_dD = (lambda ./ (4*pi*D .* sqrt(-log(R))));   % [m per m]

% Convert to um per m
dsigma_dD_um = dsigma_dD * 1e6;

% Plot derivative
figure('Color','w','Position',[80 80 900 480]);
plot(D, dsigma_dD_um,'r','LineWidth',2); grid on;
xlabel('Dish diameter, D [m]');
ylabel('d\sigma/dD [\mum per m]');
title(sprintf('Derivative of Ruze Requirement  (f=%.0f GHz, G_{eff}=%.1f dBi, \\eta_0=%.2f)', ...
               f_GHz, G_eff_dBi, eta0));

%% -------------------- Notes --------------------
% - eta0 should capture everything EXCEPT surface roughness (e.g., illumination,
%   spillover, blockage, ohmic, polarization, phase, pointing, etc.).
% - If your target G_eff is very aggressive for a given D and eta0, the curve
%   will show an infeasible region (no real sigma). Increase D or eta0, or
%   reduce the target G_eff.
% - As the target G_eff gets easier relative to D and eta0, the required
%   sigma grows (curve rises) because rougher surfaces can still meet the target.
