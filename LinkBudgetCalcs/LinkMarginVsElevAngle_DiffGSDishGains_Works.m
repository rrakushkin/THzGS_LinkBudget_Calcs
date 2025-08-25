clear, close all;
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
%geff_gsAnt    = 65;        % [dBi]
directivity_gsAnt = [60,62,64,65,66,68,70,72,74];   % [dBi]
numElem       = 16;       % satellite array (num per length)
distElem      = 0.0045;    % [m]
D             = 1.5;       % [m]
%surface_rms   = 100;       % [um]
freq          = 225;       % [GHz]
freq_Hz       = freq * 1e9;
gs_pol_type   = 'linear';
pol_angle     = 45;        % [deg]
gs_ptg_error  = 0.01;      % [deg]
sat_ptg_error = 0.1;       % [deg]
targetBER     = 1e-4;
rolloff       = 0.3;

% ---- Geometry ----
HOSL = 37e-3;              % [km] (GS height above sea level)
alt  = 460;                % [km] 
Elev = linspace(0,90,100); % [deg]

% ---- Atmospheres ----
hstep    = 0.1;            % [km]
gs_lat   = 42.3378054237531;
atmType = "InterpSummer"; %["Summer 45","Winter 45","Annual 15","InterpWinter","InterpSummer"];

% ---- GS Noise Profiles ----
[T1,P1,e1]     = atmProfile(HOSL,"Annual 15");
[T2S,P2S,e2S]  = atmProfile(HOSL,"Summer 45");
[gs_s,a,c]     = InterpAtm({T1,P1,e1},{T2S,P2S,e2S},gs_lat);
T_gs_s = max(gs_s, 300)
%[T2W,P2W,e2W]  = atmProfile(HOSL,"Winter 45");
%[gs_Tw]     = InterpAtm({T1,P1,e1},{T2W,P2W,e2W},gs_lat);
BW       = 2e9;   %[Hz] 
NF       = 8;       % [dB]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRECOMPUTATIONS

% Antenna directivities (for reference; geff_* are held fixed above)
lambda   = c ./ freq_Hz;
%g_satAnt = (4.*pi.*(distElem.^2)./(lambda.^2)).*numElem.*numElem;
%g_gsAnt  = (pi.*D./lambda).^2;

% Surface efficiency (not applied because geff_gsAnt provided)
%surface_efficiency = exp(-(4*pi*(surface_rms*1e-6)./lambda).^2);

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

% HPBW estimate from directivity
Dlin     = 10.^(directivity_gsAnt/10);
GS_HPBW  = sqrt(41253 ./ Dlin)   % degrees

% colors for plotting
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
for k = 1:numel(directivity_gsAnt)
    % Atmospheric absorption along slant path
    l_abs = absLossSlant(alt, freq, Elev, hstep, HOSL, atmType, gs_lat); % [dB]
    % Received power vs elevation (dBm)
    p_rx_dBm = zeros(size(Elev));
    for j = 1:numel(Elev)
        p_rx_dBm(j) = linkBudget1(sat_tx, geff_satAnt, directivity_gsAnt(k), GS_HPBW(k), freq_Hz, ...
                             numElem, distElem, ...
                             slant_dist_m(j), l_abs(1,1,j), ...
                             sat_ptg_error, gs_ptg_error, gs_pol_type, ...
                             pol_angle); % result is in dBm
    end
   
     % Noise floor (dBm), SNR, Eb/N0, Link Margin
     Tsys = Tsys_s;
     Pnf_dBm    = 10*log10(K_boltz*Tsys*BW) + 30;        %dB
     snr_dB     = p_rx_dBm  - Pnf_dBm;                    %dB

     %Commented out for SNR to understand pure tone behavior
     EbNo_dB    = snr_dB + 10*log10((1+rolloff)/b) ;     %dB
     linkMargin = EbNo_dB - EbNo_min_dB;                %dB

    % ---------- Figure 1: plot link margin for each gain value ----------
    name1 = sprintf('%.2f dBi',directivity_gsAnt);

    plot(Elev, linkMargin, 'LineWidth', 1.5,'Color', colors(k,:), ...
         'DisplayName', sprintf('%.0f dBi', directivity_gsAnt(k)));
     if directivity_gsAnt(k) == 74
         hline = yline(3, 'Color', 'black', 'LineStyle','--', ...
             'LineWidth', 3, 'DisplayName', '3 dB');
     end
%end 
    % --- Assign colors: gray for all except 60 dBi case ---
    %if directivity_gsAnt(k) == 60
        %thisColor = [162, 20, 47] / 255;   % #A2142F (highlight)
    %else
        %grayShades = [0.8 0.8 0.8;0.7 0.7 0.7;0.6 0.6 0.6;0.55 0.55 0.55; 0.5 0.5 0.5;...
         %0.45 0.45 0.45; 0.4 0.4 0.4; 0.35 0.35 0.35; 0.3 0.3 0.3; 0.2 0.2 0.2];
        %thisColor = grayShades(mod(k-1, size(grayShades,1))+1, :); % cycle through grays
    %end

    % Plot
    %plot(Elev, linkMargin, 'LineWidth', 1.5, ...
         %'Color', thisColor, ...
         %'DisplayName', sprintf('%.0f dBi', directivity_gsAnt(k)));
    %if directivity_gsAnt(k) == 74
    %hline = yline(3, 'Color', '#A2142F', 'LineStyle','--', ...
                  %'LineWidth', 3, 'DisplayName', '3 dB');
    %end
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


