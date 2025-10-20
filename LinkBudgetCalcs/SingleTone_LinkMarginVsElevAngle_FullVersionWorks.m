clear, clc, close all;
%addpath(genpath('./functions'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES

% ---- Constants ----
K_boltz = physconst('boltzman');
T0      = 290;
c       = physconst('LightSpeed');

% ---- Link budget params ----
sat_tx        = 23;        % [dBm]
geff_satAnt   = 44;        % [dBi]
geff_gsAnt = [60,62,64,65,66,68,70,72,74];   % [dBi]
freq          = 225;       % [GHz]
freq_Hz       = freq * 1e9;
targetBER     = 1e-4;       %= 10^(-4)
rolloff       = 0.3;
BW       = 130e6;             %[Hz] 130 MHz = 100 Mbps
NF       = 8;               %[dB]

% ---- Geometry ----
HOSL = 18e-3;              % [km] (Conservative Value; GS height above sea level)
alt  = 460;                % [km] 
Elev = linspace(0,90,100); % [deg]

% ---- Atmospheres ----
hstep    = 0.1;            % [km]
gs_lat   = 42.3378054237531; % Egan Roof
atmType = "InterpSummer"; %["Summer 45","Winter 45","Annual 15","InterpWinter","InterpSummer"];

% ---- GS Noise Profiles ----
[T1,P1,e1]     = atmProfile(HOSL,"Annual 15");
[T2S,P2S,e2S]  = atmProfile(HOSL,"Summer 45");
[Tgs_s,~,~]     = InterpAtm({T1,P1,e1},{T2S,P2S,e2S},gs_lat);
%[T2W,P2W,e2W]  = atmProfile(HOSL,"Winter 45");
%[Tgs_w,~,~]     = InterpAtm({T1,P1,e1},{T2W,P2W,e2W},gs_lat);
T_gs_s = max(Tgs_s, 300);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRECOMPUTATIONS

% BPSK required Eb/N0 for target BER
M = 2; b = log2(M);
BER_function = @(x) berawgn(x,"psk",M,"nondiff") - targetBER;
EbNo_min_dB  = fzero(BER_function,10);   % good initial guess

% System noise temperatures for seasonal cases
Tsys_s = (10^(NF/10)-1)*T0 + T_gs_s;
%Tsys_w = (10^(NF/10)-1)*T0 + gs_w;

% Slant range (km->m)
Re          = 6371;
GS_pos      = Re + HOSL;
GS_pos_m      = GS_pos*1e3;
alt_m = alt*1e3;
slant_dist_m  = sqrt(GS_pos_m.^2 .* sind(Elev).^2 + 2*GS_pos_m*alt_m + alt_m.^2) ...
              - GS_pos_m .* sind(Elev);                     % [m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORMATTING
figure
hold on
xlabel('Elevation Angle [deg]');
ylabel('Link Margin (dB)');
title(sprintf('Link Margin vs Elevation Angle at %.0f GHz', freq));
xlim([10,90]); ylim([-40,100]); 
grid on;
legend;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP (per atmosphere)
colors = [
    0.1216 0.4667 0.7059;   % Blue
    1.0000 0.4980 0.0549;   % Orange
    0.1725 0.6275 0.1725;   % Green
    0.8392 0.1529 0.1569;   % Red
    0.5804 0.4039 0.7412;   % Purple
    0.5490 0.3373 0.2941;   % Brown
    0.8902 0.4667 0.7608;   % Pink
    0.4980 0.4980 0.4980;   % Gray
    0.7373 0.7412 0.1333;   % Olive
];
for k = 1:numel(geff_gsAnt)
    % Atmospheric absorption along slant path
    l_abs = absLossSlant(alt, freq, Elev, hstep, HOSL, atmType, gs_lat); % [dB]
    % Received power vs elevation (dBm)
    p_rx_dBm = zeros(size(Elev));
    for j = 1:numel(Elev)
        p_rx_dBm(j) = linkBudget_Simplified(sat_tx, geff_satAnt, geff_gsAnt(k), freq_Hz, slant_dist_m(j), l_abs(1,1,j)); % result is in dBm
    end
   
     % Noise floor (dBm), SNR, Eb/N0, Link Margin
     Tsys = Tsys_s;
     Pnf_dBm    = 10*log10(K_boltz*Tsys*BW) + 30;        %dB
     snr_dB     = p_rx_dBm  - Pnf_dBm;                   %dB

    name1 = sprintf('%.2f dBi',geff_gsAnt);

    plot(Elev, snr_dB, 'LineWidth', 1.5,'Color', colors(k,:), ...
         'DisplayName', sprintf('%.0f dBi', geff_gsAnt(k)));
     if k == numel(geff_gsAnt)
         hline = yline(3, 'Color', 'black', 'LineStyle','--', ...
             'LineWidth', 3, 'DisplayName', '3 dB');
     end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
