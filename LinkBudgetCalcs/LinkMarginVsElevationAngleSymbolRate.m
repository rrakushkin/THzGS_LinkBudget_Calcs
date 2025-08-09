clear, clc, close all;

addpath('./functions')

% Constants
k = physconst('boltzman');
Re = physconst('earthRadius');

%  Parameters
LinkParams.fc = 225e9;
LinkParams.Rs = [10e6, 100e6, 1e9];
LinkParams.PtxDBm = 24;
LinkParams.Gsat = 44;
LinkParams.Ggs = 65;
LinkParams.nfdb = 7;
LinkParams.targetBER = 1e-4;
LinkParams.rolloff = 0.3;
LinkParams.LinkMargin = 3;
LinkParams.OrbitalAltitude = 416e3;
Data = [];
snrData = [];
BwData = [];

linestyles = {'--', '-', ':'};
markers = {'o', 'diamond', '^', 'v'};
colors = {'#D95319','#0072BD', '#77AC30', '#A2142F'}; % {Orange, Blue, Green, Red}
handles = [];
markEvery = 12;
offset = 6;
linewidth=2;

elevAngle = linspace(0, 90);
slantPathDistance = sqrt((Re+LinkParams.OrbitalAltitude)^2 ...
                    -Re^2.*cosd(elevAngle).^2) ...
                    - Re.*sind(elevAngle);

% Link Budget
labs = absLossSlant(LinkParams.OrbitalAltitude * 1e-3, LinkParams.fc * 1e-9, elevAngle, 0.1, 0, "globalAnnual");
labs = labs(1,:);
p_rx_dbm = linkBudget(LinkParams.PtxDBm, LinkParams.Gsat, LinkParams.Ggs, LinkParams.fc, slantPathDistance, 0)-labs;
for i = 1:length(LinkParams.Rs)
    bw = LinkParams.Rs(i)*(1+LinkParams.rolloff);
    snr_db = receivedSNR(p_rx_dbm, LinkParams.nfdb, bw);
    snrData = [snrData; snr_db];
end

figure

%% Link Margin
% For the different received SNR, we compute the corresponding link margin
% based on the Modulation chosen
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

% Finding Elevation angle mask based on link margin
[~, ind] = min(abs(linkMargin(2,:)-3));
ElevMask = elevAngle(ind);

% Plotting
for i = 1:length(LinkParams.Rs)
    h = plot(elevAngle, linkMargin(i,:), 'Color', colors{1}, 'LineWidth', linewidth, 'LineStyle', linestyles{i});
    hold on;
    handles = [handles h];
end

% for e = 0:10:90
%     [~, i] = min(abs(e-elevAngle));
%     sprintf('Elev:{%d}\n', e)
%     sprintf('Link Margin = {%d}\n', linkMargin(2,i))
% end

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
for i = 1:length(LinkParams.Rs)
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
for i = 1:length(LinkParams.Rs)
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
for i = 1:length(LinkParams.Rs)
    h = plot(elevAngle, linkMargin(i,:), 'Color', colors{4}, 'LineWidth', linewidth, 'LineStyle', linestyles{i});
    handles = [handles h];
end

M = reshape(handles, 3, 4); % Convert to 4x4 matrix
M_T = M'; % Transpose the matrix
handles = M_T(:); % Convert back to a 1D array
hlgd = gridlegend(handles, ...
    {'BPSK (x1bps/Bd)', 'QPSK (x2bps/Bd)', '8PSK (x3bps/Bd)', '16PSK (x4bps/Bd)'}, ...
    {'10MBd', '100MBd', '1GBd'}, ...
     'Location', 'southeast', ...
     'FontSize',14);

hline = yline(3, 'Color', 'k', 'LineStyle','--', 'LineWidth', 3);
vline = xline(15, 'Color', 'k', 'LineStyle','--', 'LineWidth', 3);

grid on;
xlabel('Elevation Angle [deg]')
ylabel('Link Margin [dB]')
xlim([0, 90])
ylim([-60,50])
xticks(0:10:90)
% yscale('log')

%% Contact time and Max contact over 24 hours
% Define time span: 24 hours from now
startTime = datetime(2025, 1, 1);
stopTime = startTime + hours(24);

% Create satellite scenario
sc = satelliteScenario(startTime, stopTime, 1); % 60s sample time

% Create satellite scenario
% Define orbital elements for Satellite
semiMajorAxis = Re + 416e3;   % [km]
eccentricity = 0; 
inclination = 51.64;       % [deg]
rightAscension = 0;                         % [deg]
argumentOfPerigee = 0;                      % [deg]
trueAnomaly = 0;                            % [deg]

% Create CubeSat object using these parameters
sat = satellite(sc, ...
                  semiMajorAxis, ...
                  eccentricity, ...
                  inclination, ...
                  rightAscension, ...
                  argumentOfPerigee, ...
                  trueAnomaly);

% Create ground station in Boston
gs = groundStation(sc, "Name", "Boston", ...
                   "Latitude", 42.3601, ...
                   "Longitude", -71.0589, ...
                   "Altitude", 0, ...
                   "MinElevationAngle", ElevMask);

% Set up access analysis
access = access(sat, gs);
accessIntervals = accessIntervals(access);

% Extract durations of each contact
durations = seconds(accessIntervals.EndTime - accessIntervals.StartTime);

% Compute results
totalContactTime = sum(durations);             % In seconds
longestContactTime = max(durations);           % In seconds

% Display results
fprintf('Elevation Mask for 100MBd and BPSK, BER<1e-4=%.2fdeg\n', ElevMask)
fprintf('Total contact time over 24h: %.2f minutes\n', totalContactTime / 60);
fprintf('Longest single contact: %.2f minutes\n', longestContactTime / 60);

% Visualize the scenario (optional)
% play(sc);