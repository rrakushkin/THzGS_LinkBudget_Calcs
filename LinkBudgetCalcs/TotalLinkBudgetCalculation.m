clear, clc, close all;
addpath(genpath('./functions'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLES

% ---- Constants ----
K_boltz = physconst('boltzman');
T0      = 290;
c       = physconst('LightSpeed');

% ---- Link budget params ----
sat_tx        = 25;        % [dBm]
geff_satAnt   = 44;        % [dBi]
numElem       = 16;
distElem      = 0.0045;    % [m]
geff_gsAnt    = 65;        % [dBi]
freq          = 225;       % [GHz]
freq_Hz       = freq * 1e9;
gs_pol_type   = 'linear';
pol_angle     = 15;        % [deg]
D             = 1.5;       % [m]
gs_ptg_error  = 0.01;      % [deg]
sat_ptg_error = 0.1;       % [deg]
surface_rms   = 100;       % [um]
targetBER     = 1e-5;
rolloff       = 0.3;

% ---- Geometry ----
HOSL = 37e-3;              % [km]
alt  = 420;                % [km]
Elev = linspace(0,90,100); % [deg]

% ---- Atmospheres ----
hstep    = 0.1;            % [km]
gs_lat   = 42.3378054237531;
atmTypes = ["Summer 45","Winter 45","Annual 15","InterpWinter","InterpSummer"];

% ---- GS Noise Profiles ----
[T1,P1,e1]     = atmProfile(HOSL,"Annual 15");
[T2S,P2S,e2S]  = atmProfile(HOSL,"Summer 45");
[T2W,P2W,e2W]  = atmProfile(HOSL,"Winter 45");
[gs_s,~,~]     = InterpAtm({T1,P1,e1},{T2S,P2S,e2S},gs_lat);
[gs_w,~,~]     = InterpAtm({T1,P1,e1},{T2W,P2W,e2W},gs_lat);

% ---- Single-BW items used only for Figure 1 MDS lines ----
BW       = 2e9;   % [Hz] (Figure 1 reference only)
NF       = 7;       % [dB]
SNR_req  = 3;       % [dB]
LM_THR   = 3;       % [dB] threshold shown on Figure 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRECOMPUTATIONS

% Antenna directivities (for reference; geff_* are held fixed above)
lambda   = c ./ freq_Hz;
g_satAnt = (4.*pi.*(distElem.^2)./(lambda.^2)).*numElem.*numElem;
g_gsAnt  = (pi.*D./lambda).^2;

% Surface efficiency (not applied because geff_gsAnt provided)
surface_efficiency = exp(-(4*pi*(surface_rms*1e-6)./lambda).^2);

% BPSK required Eb/N0 for target BER
M = 2; b = log2(M);
BER_function = @(x) berawgn(x,"psk",M,"nondiff") - targetBER;
EbNo_min_dB  = fzero(BER_function,10);   % good initial guess

% System noise temperatures for seasonal cases
Tsys_s = (10^(NF/10)-1)*T0 + gs_s;
Tsys_w = (10^(NF/10)-1)*T0 + gs_w;

% Slant range (km->m)
Re          = 6371;
GS_pos      = Re + HOSL;
slant_dist  = sqrt(GS_pos.^2 .* sind(Elev).^2 + 2*GS_pos*alt + alt.^2) ...
              - GS_pos .* sind(Elev);                     % [km]
slant_dist_m = slant_dist*1e3;                            % [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURES

% Figure 1 (received power) axes
f1 = figure(1); clf(f1); ax1 = axes(f1); hold(ax1,'on');

% Figure 2 (link margin vs BW & season) axes
f2 = figure(2); clf(f2); ax2 = axes(f2); hold(ax2,'on');

% Color map for atmospheres (Fig 1)
colors_atm = [ 0.40 0.80 0.40;   % Summer 45
               0.60 0.40 0.80;   % Winter 45
               0.30 0.55 0.85;   % Annual 15
               0.90 0.30 0.40;   % InterpWinter
               0.90 0.30 0.40 ]; % InterpSummer

% Bandwidths (Fig 2): color encodes BW, style encodes season
bwList   = [2e9, 100e6, 10e6];
bwLabels = ["2 GHz","100 MHz","10 MHz"];
colorsBW = [0.20 0.44 0.69;   % blue
            0.85 0.33 0.10;   % orange
            0.47 0.67 0.19];  % green

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP (per atmosphere)

for k = 1:numel(atmTypes)
    atm = atmTypes(k);

    % Atmospheric absorption along slant path
    l_abs = absLossSlant(alt, freq, Elev, hstep, HOSL, atm, gs_lat); % [dB]

    % Received power vs elevation (dBm)
    p_rx = zeros(size(Elev));
    for j = 1:numel(Elev)
        p_rx(j) = linkBudget(sat_tx, geff_satAnt, geff_gsAnt, freq_Hz, ...
                             numElem, distElem, ...
                             slant_dist_m(j), l_abs(1,1,j), D, ...
                             sat_ptg_error, gs_ptg_error, gs_pol_type, ...
                             pol_angle);
    end

    % ---------- Figure 1: received power ----------
    if atm == "InterpWinter"
        ls1 = ":";  name1 = "Boston Winter";
    elseif atm == "InterpSummer"
        ls1 = "--"; name1 = "Boston Summer";
    else
        ls1 = "-";  name1 = atm + "°";
    end

    plot(ax1, Elev, p_rx, 'DisplayName', name1, ...
         'LineWidth', 1.5, 'Color', colors_atm(k,:), 'LineStyle', ls1);

    % ---------- Figure 2: link margin (seasonal only) ----------
    if atm=="InterpWinter" || atm=="InterpSummer"

        % season mapping for style and Tsys
        if atm=="InterpWinter"
            Tsys   = Tsys_w;
            lmLS   = '-';     % winter solid
            season = "Boston Winter";
        else
            Tsys   = Tsys_s;
            lmLS   = '--';    % summer dashed
            season = "Boston Summer";
        end

        for m = 1:numel(bwList)
            BWm = bwList(m);

            % Noise floor (dBm), SNR, Eb/N0, Link Margin
            Pnf_dBm    = 10*log10(K_boltz*Tsys*BWm*1e3);        % dBm
            snr_dB     = p_rx - Pnf_dBm;                         % dB
            EbNo_dB    = snr_dB + 10*log10((1+rolloff)/b);       % dB
            linkMargin = EbNo_dB - EbNo_min_dB;                  % dB

            % First crossing with LM_THR (if any)
            y  = linkMargin - LM_THR;
            zc = find(y(1:end-1).*y(2:end) <= 0, 1, 'first');    % first sign change
            x_cross = NaN;
            if ~isempty(zc)
                i1 = zc; i2 = zc+1;
                x_cross = Elev(i1) - y(i1)*(Elev(i2)-Elev(i1))/(y(i2)-y(i1));
            end

            % Legend label: include crossing elevation
            if ~isnan(x_cross)
                lbl = sprintf('%s – %s (E=%.1f°)', season, bwLabels(m), x_cross);
            else
                lbl = sprintf('%s – %s (E=N/A)',  season, bwLabels(m));
            end

            % Plot curve
            plot(ax2, Elev, linkMargin, ...
                 'DisplayName', lbl, ...
                 'LineWidth', 1.5, ...
                 'Color', colorsBW(m,:), ...
                 'LineStyle', lmLS);

            % Optional marker at the crossing (hidden from legend)
            if ~isnan(x_cross)
                hMk = plot(ax2, x_cross, LM_THR, 'o', ...
                           'MarkerSize', 5, ...
                           'MarkerFaceColor', colorsBW(m,:), ...
                           'MarkerEdgeColor', colorsBW(m,:), ...
                           'HandleVisibility','off');
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORMATTING

% ---- Figure 1 ----
figure(f1)
xlabel('Elevation Angle [deg]');
ylabel('Received Power [dBm]');
title(sprintf('Received Power vs Elevation Angle at %.0f GHz', freq));
xlim([0,90]); ylim([-200,-40]); grid on;
legend(ax1,'show','Location','best');

% MDS reference lines for Fig 1 (keep out of legend)
Pnf_s = 10*log10(K_boltz*Tsys_s*BW*1e3);
Pnf_w = 10*log10(K_boltz*Tsys_w*BW*1e3);
MDS_s = SNR_req + Pnf_s;
MDS_w = SNR_req + Pnf_w;
hs = yline(ax1,MDS_s,'--',sprintf('MDS_s = %.2f dBm',MDS_s), ...
           'Color','black','LabelVerticalAlignment','top', 'LabelHorizontalAlignment','left');
hs.Annotation.LegendInformation.IconDisplayStyle='off';
%hw = yline(ax1,MDS_w,'--',sprintf('MDS_w = %.2f dBm',MDS_w), ...
           %'Color', 'black','LabelVerticalAlignment','top','LabelHorizontalAlignment','left');

%hw.Annotation.LegendInformation.IconDisplayStyle='off';

% ---- Figure 2 ----
figure(f2)
xlabel('Elevation Angle [deg]');
ylabel('Link Margin [dB]');
title(sprintf('Link Margin vs Elevation Angle at %.0f GHz', freq));
xlim([0,90]); ylim([-200,50]); grid on;

% 3 dB reference line (hidden from legend)
hthr = yline(ax2, LM_THR, '-', '3 dB', ...
             'Color','black', 'LabelHorizontalAlignment','left');
hthr.Annotation.LegendInformation.IconDisplayStyle = 'off';

lgd = legend(ax2,'show','Location','best');
lgd.Title.String = 'Color = BW, Line = Season';
