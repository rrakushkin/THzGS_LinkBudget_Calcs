%Delete this script later!

% Raw data
Diameter = [1.5, 1.5, 1.6, 1.6, 1.6, 1.7, 1.7, 1.9, 1.9, 2.0, 2.0]';
FD_ratio = [0.5, 0.7, 0.5, 0.6, 0.7, 0.5, 0.6, 0.4, 0.5, 0.48, 0.5]';
HPBW = [0.109, 0.072, 0.102, 0.078, 0.067, 0.096, 0.073, NaN, 0.085, 0.088, 0.0818]';
Gain = [62.33, 65.76, 62.89, 64.86, 66.32, 63.42, 65.38, 61.539, 64.38, 64.34, 64.82]';

% Remove NaNs for interpolation
validIdx = ~isnan(HPBW);
x = Diameter(validIdx);
y = FD_ratio(validIdx);
zGain = Gain(validIdx);
zHPBW = HPBW(validIdx);

% Interpolants
F_gain = scatteredInterpolant(x, y, zGain, 'natural', 'none');
F_hpbw = scatteredInterpolant(x, y, zHPBW, 'natural', 'none');

% Grid for interpolation
[xx, yy] = meshgrid(linspace(min(x), max(x), 100), linspace(min(y), max(y), 100));
zzGain = F_gain(xx, yy);
zzHPBW = F_hpbw(xx, yy);

% Plotting
figure;
surf(xx, yy, zzGain, zzHPBW, 'EdgeColor', 'none');
xlabel('Diameter (m)');
ylabel('F/D Ratio');
zlabel('Interpolated Gain (dBi)');
title('Interpolated Gain vs Diameter and F/D Ratio');
colormap('jet');
colorbar;
ylabel(colorbar, 'HPBW (degrees)');
view(135, 30);
grid on;
