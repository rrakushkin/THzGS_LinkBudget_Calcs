clear, clc, close all;
addpath(genpath('./functions'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES %

% ---- Constant(s) ----
K_boltz = physconst('boltzman');
T0 = 290 ;
c = physconst('LightSpeed');

% ---- Link budget parameters ----
sat_tx = 23;                % [dBm]
sat_gain = 44;              % [dBi]
numElem = 16;               % number of horns in horn array antenna per dimension (16x16)
distElem = 0.0045;        % [m] CHECK WITH ALBERT!!!
%geff_gsAnt = 65;           % [dBi] optionally can use set value and comment out calculation of this below 
freq = 225;                 % [GHz]
freq_Hz = freq * 1e9;       % [Hz]
gs_pol_type = 'circular';   % polarization of GS antenna (either circular or linear <<lowercase)
pol_angle = 0;              % linspace(0,180,360); %[degrees] this is used only when gs_pol_type is linear 
D = 1.5;                    % [m] diameter of GS dish
gs_ptg_error = 0.01;        % [degree] of axis angle of ground station
sat_ptg_error = 0.1;        % [degree] of axis angle of satellite (Albert said to go with 0.1)
surface_rms = 50;           % [um] rms surface roughness 


% ---- GS Location Geometry parameters ----
HOSL = 37e-3;               % [km] height above sea level 
alt = 420;                  % [km]
Elev = linspace(0,90,100);  % [deg] default 100 points between min,max numbers. Written 100 for completion/later ease of toggling 

% ---- Atmospheric parameters ----
hstep = 0.1;  % [km]
%atmTypes = ["globalAnnual", "Annual 15", "Summer45", ...
            %"Winter45", "Summer60", "Winter60"]; % These are simple seasonal reference atmospheres for low (15° N), mid (45° N),and high (60° N) northern hemisphere latitude regimes.
                                                      % We are targetting ISS orbit so between mid and high 

gs_lat = 42.3378054237531; % GS on Egan Roof Coordinates: (42.3378054237531, -71.08874165317037) (lat,long)
atmTypes = ["Summer 45", "Winter 45", "Annual 15", "InterpWinter", "InterpSummer"];

% ---- GS Noise Parameters ----
[T1, P1, e1] = atmProfile(HOSL, "Annual 15");
[T2S, P2S, e2S] = atmProfile(HOSL, "Summer 45");
[T2W, P2W, e2W] = atmProfile(HOSL, "Winter 45");
[gs_s, ~, ~] = InterpAtm({T1, P1, e1}, {T2S, P2S, e2S}, gs_lat);
[gs_w, ~, ~] = InterpAtm({T1, P1, e1}, {T2W, P2W, e2W}, gs_lat);

%gs_temp = max(gs_s,gs_w);             % [Kelvin] temperature (results in gs_temp = 295.2696)
BW = 2e9;            % [Hz] bandwidth
NF = 7;                % [dB] Noise Figure
SNR_req = 3;           % [dB] minimum SNR necessary to reliably detect signal (aka link margin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% geff_* are vectors over freq (because lambda depends on freq)
% Ensure column vectors for clean indexing
freq = linspace(220, 240, 21);          % [GHz]
freq_Hz = freq(:) * 1e9;                % [Hz], column
lambda = c ./ freq_Hz;                  % [m], column

g_satAnt = (4.*pi.*(distElem.^2)./(lambda.^2)).*numElem.*numElem;
g_gsAnt  = (pi.*D./lambda).^2;

surface_efficiency = exp(-1 * (4 * pi * (surface_rms * 1e-6) ./ lambda).^2);
otherSat_efficiencies = 0.99;
otherGS_efficiencies  = 0.7;

geff_satAnt = 10.*log10(g_satAnt.*otherSat_efficiencies);                 % [dBi], Nx1
geff_gsAnt  = 10.*log10(g_gsAnt .* surface_efficiency.*otherGS_efficiencies); % [dBi], Nx1

% ---- Noise (MDS) stays as you wrote (independent of freq in current model) ----
Tsys_s = (10.^(NF./10) -1).*T0 + gs_s;
Tsys_w = (10.^(NF./10) -1).*T0 + gs_w;

MDS_s = SNR_req + 10.*log10(K_boltz .* Tsys_s .* BW.*1e3); % [dBm]
MDS_w = SNR_req + 10.*log10(K_boltz .* Tsys_w .* BW.*1e3); % [dBm]

% ---- Distance (as you wrote) ----
Elev = linspace(0,90,100);           % [deg], 1xE
Re = 6371; GS_pos = Re + HOSL;
slant_dist = sqrt(GS_pos.^2 .* sind(Elev).^2 + 2 * GS_pos * alt + alt.^2) - GS_pos .* sind(Elev); % [km]
slant_dist_m = slant_dist * 1e3;     % [m], 1xE

% Colors for later plotting
colSummer = [0.9, 0.3, 0.4];
colWinter = [0.3, 0.55, 0.85];

% Helper to compute p_rx matrix for one atmosphere
function p_rx_NE = compute_prx_for_atm(atmName, freqGHz, ElevDeg, slant_m, ...
    sat_tx, geff_satAnt, geff_gsAnt, numElem, distElem, D, sat_ptg_error, gs_ptg_error, gs_pol_type, pol_angle, alt, hstep, HOSL, gs_lat)

    N = numel(freqGHz); E = numel(ElevDeg);
    p_rx_NE = nan(N, E);

    % Get absorption loss for all freqs & elevations
    % absLossSlant must accept vector freq (GHz) and vector Elev. If it returns [N x 1 x E], index accordingly:
    l_abs_all = absLossSlant(alt, freqGHz, ElevDeg, hstep, HOSL, atmName, gs_lat);

    for i = 1:N
        for j = 1:E
            % Robustly extract l_abs (works whether l_abs_all is 3D or 2D)
            lij = squeeze(l_abs_all(min(i, size(l_abs_all,1)), 1, j));
            p_rx_NE(i, j) = linkBudget( ...
                sat_tx, geff_satAnt(i), geff_gsAnt(i), freqGHz(i)*1e9, ...
                numElem, distElem, slant_m(j), lij, D, ...
                sat_ptg_error, gs_ptg_error, gs_pol_type, pol_angle);
        end
    end
end

% ---- Compute Boston (interp) cases only, then pick best frequency by earliest MDS crossing ----
% Interp Summer
p_rx_summer = compute_prx_for_atm("InterpSummer", freq, Elev, slant_dist_m, ...
    sat_tx, geff_satAnt, geff_gsAnt, numElem, distElem, D, sat_ptg_error, gs_ptg_error, gs_pol_type, pol_angle, alt, hstep, HOSL, gs_lat);

% Interp Winter
p_rx_winter = compute_prx_for_atm("InterpWinter", freq, Elev, slant_dist_m, ...
    sat_tx, geff_satAnt, geff_gsAnt, numElem, distElem, D, sat_ptg_error, gs_ptg_error, gs_pol_type, pol_angle, alt, hstep, HOSL, gs_lat);

% For each frequency, find first elevation index where p_rx >= MDS
N = numel(freq);
E = numel(Elev);
firstCross_s = nan(N,1);
firstCross_w = nan(N,1);

for i = 1:N
    js = find(p_rx_summer(i,:) >= MDS_s, 1, 'first');
    if ~isempty(js), firstCross_s(i) = Elev(js); end

    jw = find(p_rx_winter(i,:) >= MDS_w, 1, 'first');
    if ~isempty(jw), firstCross_w(i) = Elev(jw); end
end

% Pick frequency (or frequencies) with minimum crossing elevation
minElev_s = min(firstCross_s, [], 'omitnan');
minElev_w = min(firstCross_w, [], 'omitnan');

bestFreqIdx_s = find(firstCross_s == minElev_s);
bestFreqIdx_w = find(firstCross_w == minElev_w);

bestFreqs_s = freq(bestFreqIdx_s);
bestFreqs_w = freq(bestFreqIdx_w);

% ---- Print a small table ----
fprintf('\n=== Boston Interp Atmospheres: Smallest Elevation Crossings ===\n');
if ~isempty(bestFreqIdx_s)
    fprintf('Summer: min crossing elevation = %.2f deg at freq(s) [GHz]: %s\n', ...
        minElev_s, num2str(bestFreqs_s(:).', '%.3g '));
else
    fprintf('Summer: (no crossing with MDS)\n');
end
if ~isempty(bestFreqIdx_w)
    fprintf('Winter: min crossing elevation = %.2f deg at freq(s) [GHz]: %s\n', ...
        minElev_w, num2str(bestFreqs_w(:).', '%.3g '));
else
    fprintf('Winter: (no crossing with MDS)\n');
end

% ---- Plot only the “best” curves (not every frequency) ----
figure; hold on;
if ~isempty(bestFreqIdx_s)
    for k = bestFreqIdx_s(:).'
        plot(Elev, p_rx_summer(k,:), 'LineWidth', 1.8, 'LineStyle','--', 'Color', colSummer, ...
            'DisplayName', sprintf('Boston Summer @ %.0f GHz', freq(k)));
    end
end
if ~isempty(bestFreqIdx_w)
    for k = bestFreqIdx_w(:).'
        plot(Elev, p_rx_winter(k,:), 'LineWidth', 1.8, 'LineStyle',':', 'Color', colWinter, ...
            'DisplayName', sprintf('Boston Winter @ %.0f GHz', freq(k)));
    end
end

% MDS lines (hidden from legend)
yline(MDS_s, '--', sprintf('MDS Summer = %.2f dBm', MDS_s), 'Color', [0.6 0.6 0.6], 'LabelVerticalAlignment','top');
yline(MDS_w, '--', sprintf('MDS Winter = %.2f dBm', MDS_w), 'Color', [0.3 0.3 0.3], 'LabelVerticalAlignment','top');

xlabel('Elevation Angle [deg]');
ylabel('Received Power [dBm]');
title('Received Power vs Elevation — Best Boston Frequencies Only');
xlim([0 90]); ylim([-200 -40]); grid on;
legend('Location','southoutside'); box on;
