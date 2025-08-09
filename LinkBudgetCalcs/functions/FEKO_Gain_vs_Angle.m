function [hpbw,hpbw_points] = FEKO_Gain_vs_Angle(filename)

% Load the data manually (assumes 2 columns: theta [deg] and gain [dBi])
data = readmatrix(filename);

theta = data(:,1);        % angles in degrees from col 1
gain_dBi = data(:,2);     % gain in dBi from col 2

% Step 1: Find peak gain
[max_gain, idx_max] = max(gain_dBi);
half_power_level = max_gain + (10 * log10(0.5)); %10log10(0.5) yieldsv negative ~3.01...etc (aka half power in dB)

% Step 2: Find indices where gain crosses the half-power level
above = gain_dBi > half_power_level; %above is boolean matrix with 1 where gain_dBi element is in fact > half_power_level and 0s elsewhere
crossings = find(diff(above) ~= 0); %gives the indices between which a change occurs; so 0->1 or 1->0 is seen

% Step 3: Interpolate exact -3 dB crossing points
hpbw_points = [];
for i = 1:length(crossings) %interpolate for each crossing using linear formula
    idx1 = crossings(i);
    x1 = theta(idx1);    y1 = gain_dBi(idx1);
    x2 = theta(idx1+1);  y2 = gain_dBi(idx1+1);
    
    % Linear interpolation to find theta at -3 dB
    x_cross = x1 + (half_power_level - y1) * (x2 - x1) / (y2 - y1);
    hpbw_points(end+1) = x_cross;
end

% Step 4: Calculate HPBW

if length(hpbw_points) >= 2

    hpbw = abs(hpbw_points(2) - hpbw_points(1));

    fprintf('Total HPBW around boresight: %.2f degrees\n', hpbw);
    fprintf('3dB points at: %.2f° and %.2f°\n', hpbw_points(1), hpbw_points(2));
else
    fprintf('Could not determine HPBW — insufficient crossings found.\n');
end



% Step 5: Plot for visualization
figure;
plot(theta, gain_dBi, 'b-', 'LineWidth', 1.5); hold on;
yline(half_power_level, 'r--', 'LineWidth', 1);
xline(hpbw_points(1), 'g--');
xline(hpbw_points(2), 'g--');
xlabel('Theta (degrees)');
ylabel('Gain (dBi)');
title(sprintf('Gain Pattern with HPBW = %.2f°', hpbw));
legend('Gain', '-3 dB Line', 'HPBW Points');
grid on;

end