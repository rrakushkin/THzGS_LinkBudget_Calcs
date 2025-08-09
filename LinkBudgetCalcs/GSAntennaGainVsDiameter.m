% MATLAB script to plot Antenna Gain vs. Antenna Diameter for 225 GHz

clc; clear; close all;

% Given parameters
f = 225e9; % Frequency in Hz
c = 3e8;   % Speed of light in m/s
lambda = c / f; % Wavelength in meters

% Antenna efficiencies to plot
efficiencies = [0.70, 0.60, 0.50, 0.40, 0.25];

% Antenna diameter range (0.01m to 1m)
D = linspace(0.01, 5, 1000); % Diameter values in meters

% Figure setup
figure;
hold on;
grid on;
colors = lines(length(efficiencies)); % Generate different colors for each curve

% Compute and plot gain for each efficiency
for i = 1:length(efficiencies)
    eta = efficiencies(i);
    G = eta .* (pi * D / lambda).^2; % Peak antenna gain (linear scale)
    G_dB = 10 * log10(G); % Convert to dB
    plot(D, G_dB, 'LineWidth', 2, 'Color', colors(i, :), 'DisplayName', sprintf('\\eta = %.0f%%', eta * 100));
end

% Labels and title
xlabel('Antenna Diameter (m)');
ylabel('Antenna Gain (dB)');
title('Antenna Gain vs. Antenna Diameter at 225 GHz');
legend('Location', 'southeast');
hold off;
