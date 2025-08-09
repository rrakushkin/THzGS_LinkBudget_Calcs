clear; clc; close all;

% Parameters
f = 225e9;                 % [Hz] frequency
c = physconst('lightspeed');
lambda = c / f;            % [m]
D = 1.5;                   % [m] dish diameter
G_max = 10^(65/10);        % Linear gain (assume boresight 65 dBi)

% Angle sweep (theta off-boresight)
theta_deg = linspace(-5, 5, 1000);      % [deg]
theta_rad = deg2rad(theta_deg);

% Compute kD*sin(theta)/2
u = pi * D / lambda * sin(theta_rad);   % u = kD*sin(theta)/2

% Airy pattern (normalized gain)
G_norm = (2 * besselj(1, u) ./ u).^2;
G_norm(u == 0) = 1;   % avoid divide by 0 at theta = 0

% Convert to dB
G_dB = 10 * log10(G_norm);
G_dB = G_dB - max(G_dB);  % Normalize to 0 dB

% Plot
figure;
plot(theta_deg, G_dB, 'b', 'LineWidth', 2);
xlabel('Angle from Boresight [deg]');
ylabel('Normalized Gain [dB]');
title('Dish Antenna Radiation Pattern');
grid on; xlim([-2 2]);

% Compute HPBW numerically
idx_hp = find(G_dB >= -3);         % indices above -3 dB
theta_hp = theta_deg(idx_hp);      % corresponding angles
HPBW = abs(theta_hp(end) - theta_hp(1));

% Annotate HPBW
yline(-3, '--r', '3 dB Line');
xline(theta_hp(1), '--k');
xline(theta_hp(end), '--k');
text(0, -3.5, sprintf('HPBW = %.2f deg', HPBW), ...
     'HorizontalAlignment', 'center', 'FontSize', 12);

