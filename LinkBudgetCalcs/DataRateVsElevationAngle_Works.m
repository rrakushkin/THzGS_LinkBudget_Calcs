clear, clc, close all;
addpath(genpath('./functions'));

% Constants
k = physconst('boltzman')
Re = 6371

% ---- Constants ----
K_boltz = physconst('boltzman');
T0      = 290;
c       = physconst('LightSpeed');

% ---- Link budget params ----
sat_tx        = 23;        % [dBm]
geff_satAnt   = 44;        % [dBi]
geff_gsAnt    = 65;        % [dBi]
directivity_gsAnt =  [60,62,64,65,66,68,70,72,74];   % [dBi]
numElem       = 16;       % satellite array (num per length)
distElem      = 0.0045;    % [m]
D             = 1.5;       % [m]
%surface_rms   = 100;       % [um]
freq          = 225;       % [GHz]
freq_Hz       = freq * 1e9;
gs_pol_type   = 'linear';
pol_angle     = 45;        % [deg]
gs_ptg_error  = 0.01;      % [deg]
sat_ptg_error = 0.10;       % [deg]
targetBER     = 1e-4;
rolloff       = 0.3;
maxRs = 2e9
NF = 8
LinkMargin = 3
%Data = [];
%BwData = [];
%handles = [];


% ---- Geometry ----
HOSL = 37e-3;              % [km] (GS height above sea level)
alt  = 460;                % [km] 
Elev = linspace(0,90,100); % [deg]

% ---- Atmospheres ----
hstep    = 0.1;            % [km]
gs_lat   = 42.3378054237531;
atmType = "InterpSummer"; %["Summer 45","Winter 45","Annual 15","InterpWinter","InterpSummer", "globalAnnual"];

elevAngle = linspace(0, 90);
Re_m = Re*1000
alt_m = alt*1000
slantPathDistance = (sqrt((Re_m+alt_m)^2 ...
                    -Re_m^2.*cosd(elevAngle).^2) ...
                    - Re_m.*sind(elevAngle))

% Link Budget
% Atmosphere loss (function assumes altitude is in km already and f in GHz)
labs3d = absLossSlant(alt, freq, elevAngle, 0.1, 0, atmType, gs_lat);
labs   = squeeze(labs3d(1,1,:)).';   % 1×N row vector over elevation angles


p_rx_dbm = zeros(numel(elevAngle), numel(directivity_gsAnt));
%HPBW estimate from directivity
Dlin        = 10.^(directivity_gsAnt/10);
GS_HPBW_deg = 2*sqrt(41253 ./ Dlin);  
for k = 1:numel(elevAngle)
    d_k   = slantPathDistance(k);  % [m]
    LabsK = labs(k);               % [dB]
    for i = 1:numel(directivity_gsAnt)
        Gg_dBi   = directivity_gsAnt(i);     % GS directivity (dBi)
        HPBW_deg = GS_HPBW_deg(i);           % its HPBW (deg)
        p_rx_dbm(k,i) = linkBudget1( ...
            sat_tx, geff_satAnt, Gg_dBi, HPBW_deg, ...
            freq_Hz, numElem, distElem, ...
            d_k, LabsK, ...
            sat_ptg_error, gs_ptg_error, gs_pol_type, pol_angle);
    end
end
%p_rx_dbm = linkBudget1(sat_tx,geff_satAnt, geff_gsAnt, 0.08, freq_Hz, numElem, distElem, slantPathDistance, labs, sat_ptg_error, gs_ptg_error, gs_pol_type, pol_angle)

% System Noise Temperature
Ta = 300;                           % [K] Antenna noise temperature when pointing at the Earth (Worst case)
T0 = 290;                           % [K] Reference noise temperature
%Te = (10^(NF/10)-1)*T0;% [K] Receiver equivalent noise temperature = (F-1)T0
%T_sys = (Ta + Te);                  % [K] System noise temperature
[T1,P1,e1]     = atmProfile(HOSL,"Annual 15");
[T2S,P2S,e2S]  = atmProfile(HOSL,"Summer 45");
[gs_s,a,c] = InterpAtm({T1,P1,e1},{T2S,P2S,e2S},gs_lat);
T_gs_s = max(gs_s, 300)
T_sys = (10^(NF/10)-1)*T0 + T_gs_s

%% Data rate
% BPSK
M = 2;
b = log2(M);
% Find EbNo for the desired BER
BER_function = @(EbNo_dB) berawgn(EbNo_dB, "psk", M, "nondiff") - targetBER;
initial_guess = 10; % Initial guess in dB
EbNo_dB = fzero(BER_function, initial_guess)
% "Converting" EbNo to SNR based on BW, baudrate, and bit-rate
SNR_min_db = EbNo_dB + 10*log10((1+rolloff)/b);
% Find allowable bandwidth given noise figure (30 for converting dBm to dB)
N = p_rx_dbm - 30 - SNR_min_db - LinkMargin
Bwa = 10.^(N/10)/k/T_sys
Bwb = maxRs*(1+rolloff)
Bw = min(Bwa,Bwb)
size(Bwa)
dataRate = Bw * b / (1 + rolloff);


%Plot
% Ensure elevation is a column so columns map to scenarios
elevAngle = elevAngle(:);
Ne = numel(elevAngle);
Ng = numel(directivity_gsAnt);

% Sanity: dataRate and Bw should be Ne×Ng
% If they are not, fix upstream or use reshape accordingly.

% Fill Data/BwData matrices (columns = directivity cases)
Data   = zeros(Ne, Ng);
BwData = zeros(Ne, Ng);
for j = 1:Ng
    Data(:,  j) = dataRate(:, j);   % << use existing matrices
    BwData(:, j) = Bw(:, j);
end

% -------- Plot data rate --------
figure; hold on
for j = 1:Ng
    plot(elevAngle, Data(:, j) * 1e-9, 'DisplayName', sprintf('%.0f dBi', directivity_gsAnt(j)),'LineWidth', 1.5);  % Gbps
end
grid on
xlabel('Elevation Angle [deg]')
ylabel('Achievable Bit Rate R_b [Gbps]')
xlim([0, 90]); xticks(0:10:90)
legend('show','Location','northwest');


% -------- Plot bandwidth --------
figure; hold on
for j = 1:Ng
    plot(elevAngle, BwData(:, j) * 1e-9, 'DisplayName', sprintf('%.0f dBi', directivity_gsAnt(j)),'LineWidth', 1.5); % GHz
end
grid on
xlabel('Elevation angle [deg]')
ylabel('Allowable Bandwidth B [GHz]')
xlim([0, 90])
legend('show','Location','northwest');
