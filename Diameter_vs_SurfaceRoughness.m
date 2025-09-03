%% Ruze surface-roughness requirement vs dish diameter
clear; clc; close all;

% --- Inputs ---
f_GHz     = 225;      % [GHz]
G_eff_dBi = 68;       % [dBi]
eta0      = 0.70;
D_min     = 0.5;      % [m]
D_max     = 3.0;      % [m]
num_pts   = 600;

% --- Constants ---
c      = physconst('LightSpeed');
lambda = c/(f_GHz*1e9);

% --- Compute sigma(D) ---
D            = linspace(D_min, D_max, num_pts);
G_eff_lin    = 10^(G_eff_dBi/10);
G_max_ideal  = eta0 * (pi*D/lambda).^2;
ratio        = G_eff_lin ./ G_max_ideal;

sigma_m        = nan(size(D));
valid          = (ratio <= 1) & (ratio > 0);
sigma_m(valid) = (lambda/(4*pi)) .* sqrt(-log(ratio(valid)));
sigma_um       = sigma_m * 1e6;


%% --- Plot with green line and shading under curve ---
figure('Color','w','Position',[80 80 980 520]);
hold on; grid on;

% Mask valid region
xvalid = D(valid);
yvalid = sigma_um(valid);

% Shade under curve (down to y=0) with green, 25% transparent
fill([xvalid fliplr(xvalid)], [yvalid zeros(size(yvalid))], ...
     [0.0 0.6 0.0], 'FaceAlpha',0.75, 'EdgeColor','none');

% Plot curve on top in matching green
plot(xvalid, yvalid, 'Color',[0.0 0.5 0.0], 'LineWidth', 2);

xlabel('Dish diameter, $D$ [m]', 'Interpreter','latex');
ylabel('Required RMS surface roughness, $\sigma$ [$\mu$m]', 'Interpreter','latex');
title(sprintf(['Ruze Requirement: $\\sigma$ vs $D$ ', ...
               '$(f=%.0f\\,\\mathrm{GHz},\\; G_{\\mathrm{eff}}=%.1f\\,\\mathrm{dBi},\\; \\eta_0=%.2f)$'], ...
               f_GHz, G_eff_dBi, eta0), 'Interpreter','latex');
ylim([10,150]);

% Formula box
xl = xlim; yl = ylim;
xpos = xl(1) + 0.03*(xl(2)-xl(1));
ypos = yl(1) + 0.90*(yl(2)-yl(1));
formulaStr = ['$\eta_s = e^{-(4\pi\sigma/\lambda)^2}$', newline, ...
              '$G_{\mathrm{eff}} = \eta_0\,\eta_s\,(\pi D/\lambda)^2$'];
text(xpos, ypos, formulaStr, ...
     'Interpreter','latex', 'FontSize',12, ...
     'BackgroundColor',[1 1 1], 'EdgeColor',[0.6 0.6 0.6]);
