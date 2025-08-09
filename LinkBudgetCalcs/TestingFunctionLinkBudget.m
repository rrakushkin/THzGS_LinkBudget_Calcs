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
distElem = 0.003125;        % [m] CHECK WITH ALBERT!!!
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
% Calculations & Plotting %

% ---- Antenna Gain Computations ----

% Antenna Directivity (100% efficiency Gain) Calculations
lambda = c ./ freq_Hz;
g_satAnt = (4.*pi.*(distElem.^2)./(lambda.^2)).*numElem.*numElem; %assumes square array and negligble thickness of horn so that adj. horn separation distance = horn aperature length and width dimensions
g_gsAnt = (pi.*D./lambda).^2;

% GS Dish surface roughness loss (Ruze formula)
surface_efficiency = exp(-1 * (4 * pi * (surface_rms * 1e-6) ./ lambda).^2);
%l_rms = 10 * log10(surface_efficiency); % in dB

% Effective Antenna Gain Calculations (including inefficiency) in dBi
otherSat_efficiencies = 0.99; %placeholder
otherGS_efficiencies = 0.7;  %placeholder

geff_satAnt = 10.*log10(g_satAnt.*otherSat_efficiencies)
geff_gsAnt = 10.*log10(g_gsAnt * surface_efficiency.*otherGS_efficiencies)


% ---- Noise Calculations ----
Tsys_s = (10.^(NF./10) -1).*T0 + gs_s;
Tsys_w = (10.^(NF./10) -1).*T0 + gs_w;
 
%Minimum Detectable Signal (MDS)
Power_noisefloor = 10.*log10(K_boltz .* Tsys_s .* BW.*1e3); % [dBm]
MDS_s = SNR_req + Power_noisefloor ; % Minimum Detectable Signal (SNR_req = 3dB margin)

Power_noisefloor = 10.*log10(K_boltz .* Tsys_w .* BW.*1e3); % [dBm]
MDS_w = SNR_req + Power_noisefloor ; % Minimum Detectable Signal (SNR_req = 3dB margin)



% ---- Distance computation ----
Re = 6371;                  % [km] volumetric mean radius 
GS_pos = Re + HOSL;         % [km] geometric positional height of gs 
slant_dist = sqrt(GS_pos.^2 .* sind(Elev).^2 + ...
                  2 * GS_pos * alt + alt.^2) - ...
                  GS_pos .* sind(Elev);     % [km] verified formula to be correct by hand via proof (law of cosines then convert to sine via trig property etc.)

slant_dist_m = slant_dist * 1e3;   % [m]


% ---- Plot setup ----
colors = [
    0.4, 0.8, 0.4;   % Summer 1 – light green
    0.6, 0.4, 0.8;   % Winter 1 – lavender
    0.3, 0.55, 0.85;  % Soft, cool blue
    0.9, 0.3, 0.4;   % Boston - deep rose color
    0.9, 0.3, 0.4;   % Boston - deep rose color
];

figure;
hold on;

% ---- Loop over atmospheric conditions ----
for k = 1:numel(atmTypes)

    atm = atmTypes(k);
    
    % Get atmospheric absorption loss [dB]
    l_abs = absLossSlant(alt, freq, Elev, hstep, HOSL, atm, gs_lat);
    
    % Compute received power over all Elev angles
    p_rx = zeros(size(Elev));
    for j = 1:numel(Elev)
        p_rx(j) = linkBudget(sat_tx, geff_satAnt, geff_gsAnt, freq_Hz, ...
                             numElem, distElem, ...
                             slant_dist_m(j), l_abs(1,1,j), D, ...
                             sat_ptg_error, gs_ptg_error, gs_pol_type, ...
                             pol_angle);
    end
    
    % Plot results for this atmosphere type
    if atm == "InterpWinter"
        linetype = ":";
        displayName = "Boston Winter";
    elseif atm == "InterpSummer"
        linetype = "--";
        displayName = "Boston Summer";
    else
        linetype = "-";
        displayName = atm+"°";
    end
    
    plot(Elev, p_rx, 'DisplayName', displayName, 'LineWidth', 1.5, 'Color', colors(k,:),LineStyle=linetype);

end



% ---- Plot formatting ----
xlabel('Elevation Angle [deg]');
ylabel('Received Power [dBm]');
title(sprintf('Received Power vs Elevation Angle at %.0f GHz', freq));
xlim([0, 90]);
ylim([-200, -40]);
labelStr = sprintf('MDS_s = %.2f dBm', ((MDS_s+MDS_w)./2));
ylineHandle = yline(MDS_s, '--', labelStr, 'Color', [0.6, 0.6, 0.6], 'LabelVerticalAlignment','top');
ylineHandle.Annotation.LegendInformation.IconDisplayStyle = 'off';  % Hide from legend
ylineHandle = yline(MDS_w, '--', labelStr, 'Color', [0.3, 0.3, 0.3], 'LabelVerticalAlignment','top');
ylineHandle.Annotation.LegendInformation.IconDisplayStyle = 'off';  % Hide from legend
lgd = legend;
lgd.Position = [0.65, 0.18, 0.0, 0.25];
lgd.Title.String = 'Latitudinal & Seasonal Atmospheres';
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UNCERTANTIES/CHECK THESE FORMULAS:

% - claculation of pointing loss (for each antenna as well as total)
% - noise calculation section of this script...
% - BW used for noise calculations... 
% - uncertain as to if I am supposed to ALSO add the gain from the RF
%   chain? Because I included the NF as equivalent noise temperature (so
%   this might mean the gain was not included in the final link budget 
%   calculations)...I did not explicitly include gain from any LNA? I didn't
%   want to double count if it is included implicitly in the NF
%   value...???


% NOTE! This doesn't account for losses through cables yet! MDS does
% include SNR of 3dB though...


