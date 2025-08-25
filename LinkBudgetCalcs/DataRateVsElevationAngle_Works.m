clear; clc; close all;
addpath(genpath('./functions'));

%% -------------------- Constants --------------------
kB = physconst('Boltzmann');     % [J/K]
T0 = 290;                        % [K]
c  = physconst('LightSpeed');    % [m/s]

%% -------------------- Link budget params ------------
sat_tx        = 23;              % [dBm]
geff_satAnt   = 44;              % [dBi]  (TX antenna gain)
directivity_gsAnt = [60 62 64 65 66 68 70 72 74];  % [dBi] (sweep)
numElem       = 16;              % sat array (used inside linkBudget1)
distElem      = 0.0045;          % [m]
freq_GHz      = 225;             % [GHz] (for absLossSlant)
freq_Hz       = freq_GHz*1e9;    % [Hz]   (for linkBudget1)
gs_pol_type   = 'linear';
pol_angle     = 45;              % [deg] polarization mismatch
gs_ptg_error  = 0.01;            % [deg]
sat_ptg_error = 0.10;            % [deg]
targetBER     = 1e-4;
rolloff       = 0.3;
maxRs         = 2e9;             % [Hz] symbol-rate ceiling
NF_dB         = 8;               % [dB] receiver noise figure
LinkMargin    = 3;               % [dB]

%% -------------------- Geometry ----------------------
HOSL_km = 37e-3;                 % [km] GS height above sea level
alt_km  = 460;                   % [km]
elevDeg = linspace(0,90,100);    % [deg]
elevDeg = elevDeg(:);            % make column (Ne×1)

Re_km     = 6371;
GS_pos_km = Re_km + HOSL_km;
GS_pos_m  = GS_pos_km*1e3;
alt_m     = alt_km*1e3;

% Slant range [m]
slant_m = sqrt(GS_pos_m.^2 .* sind(elevDeg).^2 + 2*GS_pos_m*alt_m + alt_m.^2) ...
        - GS_pos_m .* sind(elevDeg);     % (Ne×1)

Ne = numel(elevDeg);
Ng = numel(directivity_gsAnt);

%% -------------------- Atmosphere --------------------
gs_lat  = 42.3378054237531;
atmType = "globalAnnual";

% alt in km, f in GHz (per your function header)
labs3d = absLossSlant(alt_km, freq_GHz, elevDeg.', 0.1, 0, atmType, gs_lat);
labs   = squeeze(labs3d(1,1,:));     % (Ne×1) atmospheric absorption [dB]
labs
%% -------------------- GS HPBW from directivity ------
% Full HPBW (degrees) from directivity (linear): theta ≈ sqrt(41253 / Dlin)
Dlin_gs     = 10.^(directivity_gsAnt/10);    % (1×Ng)
GS_HPBW_deg = sqrt(41253 ./ Dlin_gs);        % (1×Ng) full HPBW in degrees
GS_HPBW_deg
%% -------------------- Received power (Ne×Ng) --------
p_rx_dbm = zeros(Ne, Ng);
for k = 1:Ne
    d_k   = slant_m(k);   % [m]
    LabsK = labs(k);      % [dB]
    for i = 1:Ng
        p_rx_dbm(k,i) = linkBudget1( ...
            sat_tx, geff_satAnt, directivity_gsAnt(i), GS_HPBW_deg(i), ...
            freq_Hz, numElem, distElem, ...
            d_k, LabsK, ...
            sat_ptg_error, gs_ptg_error, gs_pol_type, pol_angle );
    end
end
p_rx_dbm
%% -------------------- System noise temperature -------
% Interpolate GS ambient temperature; clamp to >=300 K for worst-case;
% then add receiver NF contribution.
[T1,P1,e1]    = atmProfile(HOSL_km,"Annual 15");
[T2S,P2S,e2S] = atmProfile(HOSL_km,"Summer 45");
[gs_T,~,~]    = InterpAtm({T1,P1,e1},{T2S,P2S,e2S}, gs_lat);

T_gs_s = max(max(gs_T), 300);          % scalar worst-case ambient (K)
Tsys   = (10^(NF_dB/10)-1)*T0 + T_gs_s;  % scalar system noise temp (K)

%% -------------------- Required Eb/N0 (BPSK) ----------
M = 2;  b = log2(M);
EbNo_req_dB = fzero(@(x) berawgn(x,"psk",M,"nondiff") - targetBER, 10);
SNR_min_dB  = EbNo_req_dB + 10*log10(b/(1+rolloff));

%% -------------------- Allowed BW & Data rate ----------
% Thermal-noise-limited bandwidth, elementwise
N_dB = p_rx_dbm - 30 - SNR_min_dB - LinkMargin;       % (Ne×Ng) dBW
Bwa  = 10.^(N_dB/10) ./ kB ./ Tsys;                   % (Ne×Ng) Hz
Bwb  = maxRs*(1+rolloff);                              % scalar Hz
Bw   = min(Bwa, Bwb);                                  % (Ne×Ng) Hz (min vs scalar)

dataRate = Bw * b / (1 + rolloff);                     % (Ne×Ng) bps

%% -------------------- Sanity prints ------------------
fprintf('Alt=%g km, NF=%g dB, BER=%.1e, Atm=%s\n', alt_km, NF_dB, targetBER, atmType);
fprintf('Tsys=%.1f K; Eb/N0 req=%.2f dB; SNR_min=%.2f dB\n', Tsys, EbNo_req_dB, SNR_min_dB);
fprintf('P_rx median=%.2f dBm; BW median=%.3f GHz; Rb median=%.3f Gbps\n', ...
        median(p_rx_dbm,'all'), median(Bw,'all')/1e9, median(dataRate,'all')/1e9);

%% -------------------- Plot: Data rate ----------------
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
figure; hold on;
%grayShades = linspace(0.8,0.2,Ng);   % Ng gray levels, light→dark
%for i = 1:Ng
%    if directivity_gsAnt(i) == 65
%        thisColor = [162 20 47]/255;   % #A2142F
%    else
%        thisColor = [grayShades(i) grayShades(i) grayShades(i)];
%    end
%    plot(elevDeg, dataRate(:,i)*1e-9, 'LineWidth', 1.5, ...
%        'Color', thisColor, ...
%         'DisplayName', sprintf('%.0f dBi', directivity_gsAnt(i)));
%end
for i = 1:Ng
    plot(elevDeg, dataRate(:,i)*1e-9, 'LineWidth', 1.5,'Color', colors(i,:), ...
         'DisplayName', sprintf('%.0f dBi', directivity_gsAnt(i)));
end
grid on;
xlabel('Elevation Angle [deg]');
ylabel('Achievable Bit Rate R_b [Gbps]');
xlim([0 90]); xticks(0:10:90);
% denser y-ticks (auto based on range)
yr = ylim; yticks(linspace(yr(1), yr(2), 11));
legend('show','Location','northwest');
title(sprintf('Rb vs Elevation @ %.0f GHz', freq_GHz));
%% -------------------- Plot: Bandwidth ----------------
figure; hold on;

%for i = 1:Ng
%    if directivity_gsAnt(i) == 65
%        thisColor = [162 20 47]/255;   % #A2142F
%    else
%        thisColor = [grayShades(i) grayShades(i) grayShades(i)];
%    end
%    plot(elevDeg, Bw(:,i)*1e-9, 'LineWidth', 1.5, ...
%         'Color', thisColor, ...
%         'DisplayName', sprintf('%.0f dBi', directivity_gsAnt(i)));
%end
for i = 1:Ng
    plot(elevDeg, Bw(:,i)*1e-9, 'LineWidth', 1.5,'Color', colors(i,:), ...
         'DisplayName', sprintf('%.0f dBi', directivity_gsAnt(i)));
end
grid on;
xlabel('Elevation angle [deg]');
ylabel('Allowable Bandwidth B [GHz]');
xlim([0 90]); xticks(0:10:90);
yr = ylim; yticks(linspace(yr(1), yr(2), 11));
legend('show','Location','northwest');
title('Allowed Bandwidth vs Elevation');