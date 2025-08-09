clear, clc%, close all;

addpath(genpath('./functions'));


% Constants
k = physconst('boltzman');
Re = physconst('earthRadius');

%  Parameters
LinkParams.fc = 225e9;
LinkParams.Bw = [10e6, 100e6, 1e9, 2e9];
LinkParams.PtxDBm = 23;
LinkParams.Gsat = 45;
LinkParams.Ggs = 0;
LinkParams.nfdb = 8;
LinkParams.targetBER = 1e-5;
LinkParams.rolloff = 0.3;
LinkParams.LinkMargin = 3;
LinkParams.OrbitalAltitude = 400e3;
Data = [];
snrData = [];
BwData = [];

linestyles = {'-', '--', '-.', ':'};
markers = {'o', 'diamond', '^', 'v'};
colors = {'#D95319','#0072BD', '#77AC30', '#A2142F'}; % {Orange, Blue, Green, Red}
handles = [];
markEvery = 12;
offset = 6;
linewidth=2;

elevAngle = linspace(0, 90,100);
slantPathDistance = sqrt((Re+LinkParams.OrbitalAltitude)^2 ...
                    -Re^2.*cosd(elevAngle).^2) ...
                    - Re.*sind(elevAngle);

% Link Budget
labs = absLossSlant(LinkParams.OrbitalAltitude * 1e-3, LinkParams.fc * 1e-9, elevAngle, 0.1, 0, "globalAnnual");
labs = labs(1,:);
p_rx_dbm = linkBudget(LinkParams.PtxDBm, LinkParams.Gsat, LinkParams.Ggs, LinkParams.fc, slantPathDistance, labs);
for i = 1:length(LinkParams.Bw)
    bw = LinkParams.Bw(i);
    snr_db = receivedSNR(p_rx_dbm, LinkParams.nfdb, bw);
    snrData = [snrData; snr_db];
end

figure

%% Link Margin
% BPSK
M = 2;
b = log2(M);
% Find EbNo for the desired BER
BER_function = @(EbNo_dB) berawgn(EbNo_dB, "psk", M, "nondiff") - LinkParams.targetBER;
initial_guess = 10; % Initial guess in dB
EbNo_min_dB = fzero(BER_function, initial_guess);

% Received EbNo
EbNo_dB = snrData + 10*log10((1+LinkParams.rolloff)/b);
% Link Margin
linkMargin = EbNo_dB - EbNo_min_dB;

% Plotting
for i = 1:length(LinkParams.Bw)
    h = plot(elevAngle, linkMargin(i,:), 'Color', colors{1}, 'LineWidth', linewidth, 'LineStyle', linestyles{i});
    hold on;
    handles = [handles h];
end

% QPSK
M = 4;
b = log2(M);
% Find EbNo for the desired BER
BER_function = @(EbNo_dB) berawgn(EbNo_dB, "psk", M, "nondiff") - LinkParams.targetBER;
initial_guess = 10; % Initial guess in dB
EbNo_min_dB = fzero(BER_function, initial_guess);

% Received EbNo
EbNo_dB = snrData + 10*log10((1+LinkParams.rolloff)/b);
% Link Margin
linkMargin = EbNo_dB - EbNo_min_dB;

% Plotting
for i = 1:length(LinkParams.Bw)
    h = plot(elevAngle, linkMargin(i,:), 'Color', colors{2}, 'LineWidth', linewidth, 'LineStyle', linestyles{i});
    handles = [handles h];
end

% 8PSK
M = 8;
b = log2(M);
% Find EbNo for the desired BER
BER_function = @(EbNo_dB) berawgn(EbNo_dB, "psk", M, "nondiff") - LinkParams.targetBER;
initial_guess = 10; % Initial guess in dB
EbNo_min_dB = fzero(BER_function, initial_guess);

% Received EbNo
EbNo_dB = snrData + 10*log10((1+LinkParams.rolloff)/b);
% Link Margin
linkMargin = EbNo_dB - EbNo_min_dB;

% Plotting
for i = 1:length(LinkParams.Bw)
    h = plot(elevAngle, linkMargin(i,:), 'Color', colors{3}, 'LineWidth', linewidth, 'LineStyle', linestyles{i});
    handles = [handles h];
end

% 16PSK
M = 16;
b = log2(M);
% Find EbNo for the desired BER
BER_function = @(EbNo_dB) berawgn(EbNo_dB, "psk", M, "nondiff") - LinkParams.targetBER;
initial_guess = 10; % Initial guess in dB
EbNo_min_dB = fzero(BER_function, initial_guess);

% Received EbNo
EbNo_dB = snrData + 10*log10((1+LinkParams.rolloff)/b);
% Link Margin
linkMargin = EbNo_dB - EbNo_min_dB;

% Plotting
for i = 1:length(LinkParams.Bw)
    h = plot(elevAngle, linkMargin(i,:), 'Color', colors{4}, 'LineWidth', linewidth, 'LineStyle', linestyles{i});
    handles = [handles h];
end

hlgd = gridlegend(handles', ...
     {'BPSK', 'QPSK', '8PSK', '16PSK'}, ...
     {'10MHz', '100MHz', '1GHz', '2GHz'}, ...
     'Location', 'northeast', ...
     'FontSize',14);

hline = yline(3, 'Color', 'k', 'LineStyle','--', 'LineWidth', 3);
vline = xline(15, 'Color', 'k', 'LineStyle','--', 'LineWidth', 3);

grid on;
xlabel('Elevation Angle [deg]')
ylabel('Link Margin [dB]')
xlim([0, 90])
ylim([-40, 100])
% yscale('log')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Power Received vs Elevation Angle Plot

figure;

plot(elevAngle, p_rx_dbm, 'b', 'LineWidth', 2);
grid on;
xlabel('Elevation Angle [deg]');
ylabel('Received Power [dBm]');
title(sprintf('Received Power vs Elevation at %.0f GHz', LinkParams.fc / 1e9));
xlim([0, 90]);
ylim([-200, 0]);

% Add reference lines if desired
yline(-100, '--k', 'Typical Rx Threshold', 'LabelVerticalAlignment','bottom');
xline(15, '--r', 'Min Elevation', 'LabelVerticalAlignment','top');


% Find the index closest to 10 degrees
[~, idx10] = min(abs(elevAngle - 10));

% Print the actual elevation and corresponding P_rx
fprintf("Actual Elevation: %.12f deg\n", elevAngle(idx10));
fprintf("P_rx at %.2f deg: %.4f dBm\n", elevAngle(idx10), p_rx_dbm(idx10));
