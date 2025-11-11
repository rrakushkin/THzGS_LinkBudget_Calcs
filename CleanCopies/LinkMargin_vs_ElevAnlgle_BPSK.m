clear, close all;
%addpath(genpath('./functions'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES
% ---- Constants ----
K_boltz = physconst('boltzman');
T0      = 290;
c       = physconst('LightSpeed');
%% ---- Geographic Specifications----
alt_m  = 550e3;            % [m] 
inclination_angle = 51.6;   %[degrees]
Ant_HOSL = 18+12.192+3; % height above sea level (m)--this is ground level altitude + building height +antenna mast height
GS_pos_m      = (6371000 + Ant_HOSL); %m (6371000 is Earth average radius in meters)
gs_lat   = 42.3378054237531; % Egan Roof
gs_long = -71.08885435; % Egan Roof
Elev = linspace(0,90,91);
%% ---- RF Hardware Parameters ----
freq_Hz       = 2.25 * 1e11;
M = 2; %BPSK
targetBitRate = 1e8;
targetBER     = 1e-4;       
rolloff       = 0.3;
BW       = 1.30e8;            %[Hz] 130 MHz = 100 Mbps
atmType = "InterpSummer"; %["Summer 45","Winter 45","Annual 15","InterpWinter","InterpSummer", or "AndrewsBostonProfile"];
%% ---- Hardware Performance----
sat_tx        = 24.98;   % [dBm]
geff_satAnt   = 44;      % [dBi]
NF       = 7;            %[dB]
GS_feed_network_loss = 3;
GS_radome = 2;
geff_gsAnt = [58,60,62,64,66,68,70,72];   % [dBi]
%% ---- Pointing Losses ----
gs_ptg_error_loss = 1.2; %[dB] NOTE!!!! For 0.05 deg ptg error this varies from 1.2 to 10 dB depending on directivity of dish. Using 1.2 to illustrate case for why we need high directivity
sat_ptg_error_loss = 1; %[dB]
polarization_miss_match_loss = 0.033; %5 deg mismatch is used for this. 5 degrees is a HUGE overestimate, but taking a smaller value doesn't make much of a difference anyways %[dB]
l_ptg = gs_ptg_error_loss+sat_ptg_error_loss+polarization_miss_match_loss;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRECOMPUTATIONS

% ---- Atmosphere ----
hstep    = 0.1;            % [km]
slant_dist_m  = sqrt(GS_pos_m.^2 .* sind(Elev).^2 + 2*GS_pos_m*alt_m + alt_m.^2) ...
              - GS_pos_m .* sind(Elev)      ;           % [m]

% ---- GS Noise Profiles ----
% 18+12.192+3 is height above sea level (MSL)
[T1,P1,e1]     = atmProfile(18+12.192+3,"Annual 15");
[T2S,P2S,e2S]  = atmProfile(18+12.192+3,"Summer 45");
[Tgs_s,~,~]     = InterpAtm({T1,P1,e1},{T2S,P2S,e2S},gs_lat);
[T2W,P2W,e2W]  = atmProfile(18+12.192+3,"Winter 45");
[Tgs_w,~,~]     = InterpAtm({T1,P1,e1},{T2W,P2W,e2W},gs_lat);
T_gs_s = max(Tgs_s, 300);

% ---- System noise temperatures for seasonal cases ----
Tsys_s = (10^(NF/10)-1)*T0 + T_gs_s;
Tsys_w = (10^(NF/10)-1)*T0 + Tgs_w;

% ---- BPSK required Eb/N0 for target BER ----
b = log2(M);
BER_function = @(x) berawgn(x,"psk",M,"nondiff") - targetBER;
EbNo_min_dB  = fzero(BER_function,10)  ; % good initial guess

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORMATTING
figure
hold on
xlabel('Elevation Angle [deg]');
ylabel('Link Margin (dB)');
title(sprintf('Link Margin vs Elevation Angle at %.0f GHz', freq_Hz*1e-9));
xlim([10,90]); ylim([-40,100]); 
grid on;
legend;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP (per atmosphere)

for k = 1:numel(geff_gsAnt)
    % Atmospheric absorption along slant path
    l_abs = absLossSlant(alt_m/1000, freq_Hz*1e-9, Elev, hstep, (18+12.192+3)/1000, atmType, gs_lat); % [dB]
    % Received power vs elevation (dBm)
    p_rx_dBm = zeros(size(Elev));
    for j = 1:numel(Elev)
        p_rx_dBm(j) = linkBudget_Simplified(sat_tx, geff_satAnt, geff_gsAnt(k), freq_Hz, slant_dist_m(j), l_abs(1,1,j), l_ptg); % result is in dBm
    end
   
     % Noise floor (dBm), SNR, Eb/N0, Link Margin
     Tsys = Tsys_s;
     Pnf_dBm    = 10*log10(K_boltz*Tsys*BW) + 30;        %dB
     snr_dB     = p_rx_dBm  - Pnf_dBm   ;                %dB

     %Commented out for SNR to understand pure tone behavior
     EbNo_dB    = snr_dB + 10*log10((1+rolloff)/b) ;     %dB
     linkMargin = EbNo_dB - EbNo_min_dB  ;              %dB

    % ---------- Figure 1: plot link margin for each gain value ----------
    name1 = sprintf('%.1f dBi',geff_gsAnt);

    plot(Elev, linkMargin, 'LineWidth', 1.5,'Color', colors(k,:), ...
        'DisplayName', sprintf('%.2f dBi', geff_gsAnt(k)));
     if k == numel(geff_gsAnt)
         hline = yline(3, 'Color', 'black', 'LineStyle','--', ...
             'LineWidth', 3, 'DisplayName', '3 dB');
     end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%