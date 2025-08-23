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
sat_tx        = 25;        % [dBm]
geff_satAnt   = 44;        % [dBi]
%geff_gsAnt    = 65;        % [dBi]
directivity_gsAnt = 60:2:70;   % [dBi]
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
targetBER     = 10e-4;
rolloff       = 0.3;
maxRs = 2e9
NF = 7
LinkMargin = 3
Data = [];
BwData = [];
handles = [];


% ---- Geometry ----
HOSL = 37e-3;              % [km] (GS height above sea level)
alt  = 420;                % [km] 
Elev = linspace(0,90,100); % [deg]

% ---- Atmospheres ----
hstep    = 0.1;            % [km]
gs_lat   = 42.3378054237531;
atmType = "InterpSummer"; %["Summer 45","Winter 45","Annual 15","InterpWinter","InterpSummer"];

elevAngle = linspace(0, 90);
slantPathDistance = (sqrt((Re+alt)^2 ...
                    -Re^2.*cosd(elevAngle).^2) ...
                    - Re.*sind(elevAngle))*1e3;

% Link Budget
% Atmosphere loss (function assumes altitude is in km already and f in GHz)
labs3d = absLossSlant(alt, freq, elevAngle, 0.1, 0, atmType, gs_lat);
labs   = squeeze(labs3d(1,1,:)).';   % 1×N row vector over elevation angles

p_rx_dbm = linkBudget1(23,44, 65, 0.05, 225e9, 16, 0.0045, slantPathDistance, labs, 0.1, 0.01, 'linear', 45)

% System Noise Temperature
Ta = 300;                           % [K] Antenna noise temperature when pointing at the Earth (Worst case)
T0 = 290;                           % [K] Reference noise temperature
Te = (10^(NF/10)-1)*T0;% [K] Receiver equivalent noise temperature = (F-1)T0
%T_sys = (Ta + Te);                  % [K] System noise temperature
[T1,P1,e1]     = atmProfile(HOSL,"Annual 15");
[T2S,P2S,e2S]  = atmProfile(HOSL,"Summer 45");
[gs_s,a,c] = InterpAtm({T1,P1,e1},{T2S,P2S,e2S},gs_lat);
T_sys = (10^(NF/10)-1)*T0 + gs_s
%% Data rate
% BPSK
M = 2;
b = log2(M);
% Find EbNo for the desired BER
BER_function = @(EbNo_dB) berawgn(EbNo_dB, "psk", M, "nondiff") - targetBER;
initial_guess = 10; % Initial guess in dB
EbNo_dB = fzero(BER_function, initial_guess);
% "Converting" EbNo to SNR based on BW, baudrate, and bit-rate
SNR_min_db = EbNo_dB + 10*log10(b/(1+rolloff));
% Find allowable bandwidth given noise figure (30 for converting dBm to dB)
N = p_rx_dbm - 30 - SNR_min_db - LinkMargin;
Bwa = 10.^(N/10)/k/T_sys
Bwb = maxRs*(1+rolloff)
Bw = min(Bwa,Bwb)
dataRate = Bw * b / (1 + rolloff);

figure
h = plot(elevAngle, dataRate * 1e-9);
Data = [Data; dataRate];
BwData = [BwData; Bw];
handles = [handles h];

for e = 0:10:90
    [~, i] = min(abs(e-elevAngle));
    %sprintf('Elev:{%d}\n', e)
    %sprintf('Rb = {%d}\n', dataRate(i).*1e-6)
end


grid on
xlabel('Elevation Angle [deg]')
ylabel('Achievable Bit Rate R_b [Gbps]')
xlim([0, 90])
xticks(0:10:90)


%% Bandwidth 

% BPSK
figure
h = plot(elevAngle, BwData(1,:) * 1e-9);

grid on
xlabel('Elevation angle [deg]')
ylabel('Allowable Bandwidth B [GHz]')
xlim([0, 90])
% ylim([1e-2, 100])
% yscale('log')

