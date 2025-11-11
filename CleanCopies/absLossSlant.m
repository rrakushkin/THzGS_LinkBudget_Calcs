function [L_abs] = absLossSlant(r, f, a, dh, r0, atm, Lat)
    %ABSLOSSSLANT Calculates the absorption loss for Earth-space slant paths.
    %
    % Input:
    %   r         [km] Satellite altitude.
    %   f         [GHz] Frequency. VALID RANGE: (1, 1000) GHz.
    %   a         [deg] Elevation angle.
    %   dh        [km] Height step.
    %   r0        [km] Initial height.
    %   atm       "globalAnnual", "lowAnnual", "midSummer", "midWinter", "highSummer", "highWinter" or "InterpolateBoston"
    %   Lat       [deg] latitude of ground station
    % Output:
    %    L_abs    [dB] Atmospheric absorption loss.

    a(a == 0) = 0.1;

    R = 6371 + r0; % [km] Earth (volumetric mean) radius.

    nLength    = numel(r);
    nFrequency = numel(f);
    nAngle     = numel(a);

    L_abs = zeros(nLength, nFrequency, nAngle); % [dB] Atmospheric absorption loss (ITU-R P.835-7 (08/2024))
    T=zeros(nLength);P=zeros(nLength);e=zeros(nLength);
    for iLength = 1 : nLength
        h = r0 : dh : min(r(iLength), 100); % [km] Height.
        nHeight = numel(h);

        % [m] Propagation distance (Rec. ITU-R S.1257).
        r(iLength) = r(iLength) - r0;
        d = sqrt(R ^ 2 * sind(a) .^ 2 + 2 * R * r(iLength) + r(iLength) ^ 2) - R * sind(a);

        if atm == "AndrewsBostonProfile"
            [T,P,e] = AndrewsbostonAtmProfile(h, "summer");
        elseif atm == "InterpSummer"    
            [T1,P1,e1] = atmProfile(h, "Annual 15");            
            [T2S,P2S,e2S] = atmProfile(h, "Summer 45");
            [T,P,e] = InterpAtm({T1,P1,e1}, {T2S,P2S,e2S}, Lat);
        elseif atm == "InterpWinter"
            [T1,P1,e1] = atmProfile(h, "Annual 15");
            [T2W,P2W,e2W] = atmProfile(h, "Winter 45");
            [T,P,e] = InterpAtm({T1,P1,e1}, {T2W,P2W,e2W}, Lat);
        else
            [T,P,e] = atmProfile(h, atm);
        end

            %this linearly interpolates the low global atmosphere
            %conditions with the mid summer and mid winter atmosphere conditions (45
            %deg summer and winter conditions) to determine atm conditions for GS

        % Atmospheric refractive index (Rec. ITU-R P.453).
        N0 = 315;  % Average value of atmospheric refractivity extrapolated to sea level.
        h0 = 7.35; % Scale height [km].
        ref = 1 + N0 * 1e-6 * exp(-h / h0);

        gamma   = zeros(nFrequency, nHeight, nAngle); % [dB/km] Specific attenuation by atmospheric gases.
        gamma_o = zeros(nFrequency, nHeight);         % [dB/km] Oxygen specific attenuation.
        gamma_w = zeros(nFrequency, nHeight);         % [dB/km] Water vapor specific attenuation.
        slant = zeros(nHeight, nAngle);               % Slant.

        for iFrequency = 1 : nFrequency
            for iHeight = 1 : nHeight
                gamma_o(iFrequency, iHeight) = specificAttenuationDryAir(f(iFrequency), T(iHeight), ...
                    P(iHeight) - e(iHeight), e(iHeight));
                gamma_w(iFrequency, iHeight) = specificAttenuationWaterVapor(f(iFrequency), T(iHeight), ...
                    P(iHeight) - e(iHeight), e(iHeight));

                slant(iHeight, :) = sqrt(1 - (R / (R + h(iHeight)) * ref(1) / ref(iHeight) * cosd(a)) .^ 2);
                gamma(iFrequency, iHeight, :) = (gamma_o(iFrequency, iHeight) ...
                    + gamma_w(iFrequency, iHeight)) ./ slant(iHeight, :);
            end
            L_abs(iLength, iFrequency, :) = trapz(gamma(iFrequency, :, :)) * dh;
        end
    end

    

end

