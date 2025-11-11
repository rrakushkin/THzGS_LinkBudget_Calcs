function [T, P, e] = AndrewsbostonAtmProfile(h, season)
    %BOSTONATMPROFILE Atmospheric profiles for Boston using seasonal averages
    %
    % Provides temperature, pressure, and water vapor pressure profiles
    % based on Boston's seasonal climate conditions using ITU-R P.835 methodology.
    %
    % Input:
    %   h       [km] Geometric height. VALID RANGE: (0, 100) km.
    %   season  "winter", "spring", "summer", "fall"
    %
    % Output:
    %   T       [K] Temperature.
    %   P       [hPa] Pressure.
    %   e       [hPa] Water-vapor pressure.
    %
    % Boston Location: 42.36°N, 71.06°W (mid-latitude)
    %
    % Uses NOAA climate data for Boston with ITU-R P.835 atmospheric structure

    % Validate inputs
    if any(h < 0) || any(h > 100)
        error('Height must be in range [0, 100] km');
    end
    
    valid_seasons = ["winter", "spring", "summer", "fall"];
    if ~any(strcmpi(season, valid_seasons))
        error('Season must be "winter", "spring", "summer", or "fall"');
    end

    nHeight = numel(h);
    T = zeros(1, nHeight); % [K] Temperature
    P = zeros(1, nHeight); % [hPa] Pressure

    % Convert geometric to geopotential height (unchanged from ITU-R P.835)
    hh = 6356.766 * h ./ (6356.766 + h); % [km] Geopotential height

    % Define altitude ranges (unchanged from ITU-R P.835)
    range1 = hh >= 0 & hh <= 11;
    range2 = hh > 11 & hh <= 20;
    range3 = hh > 20 & hh <= 32;
    range4 = hh > 32 & hh <= 47;
    range5 = hh > 47 & hh <= 51;
    range6 = hh > 51 & hh <= 71;
    range7 = hh > 71 & hh <= 84.852;
    range8 = hh > 84.852 & hh <= 89.716;
    range9 = hh > 89.716 & hh <= 98.452;

    switch lower(season)
        case "winter" % Boston winter conditions
            
            % Boston winter climate data
            T0 = 32.5+273.15;  % [K] Boston average winter temperature (32.5°F)
            P0 = 1015;  % [hPa] Boston average winter pressure
            
            % Temperature profile (ITU-R structure, Boston surface temperature)
            T(range1) = T0 - 6.5 * hh(range1);
            T(range2) = 216.65;
            T(range3) = 216.65 + (hh(range3) - 20);
            T(range4) = 228.65 + 2.8 * (hh(range4) - 32);
            T(range5) = 270.65;
            T(range6) = 270.65 - 2.8 * (hh(range6) - 51);
            T(range7) = 214.65 - 2.0 * (hh(range7) - 71);
            T(range8) = 186.8673;
            T(range9) = 263.1905 - 76.3232 * sqrt(1 - ((h(range9) - 91) / 19.9429) .^ 2);

            % Pressure profile (ITU-R methodology, Boston surface pressure)
            P(range1) = P0 * (T0 ./ T(range1)) .^ (-34.1632 / 6.5);
            P(range2) = 226.3226 * exp(-34.1632 * (hh(range2) - 11) ./ T(range2));
            P(range3) = 54.74980 * (216.65 ./ T(range3)) .^ 34.1632;
            P(range4) = 8.680422 * (228.65 ./ T(range4)) .^ (34.1632 / 2.8);
            P(range5) = 1.109106 * exp(-34.1632 * (hh(range5) - 47) ./ T(range5));
            P(range6) = 0.6694167 * (270.65 ./ T(range6)) .^ (-34.1632 / 2.8);
            P(range7) = 0.03956649 * (214.65 ./ T(range7)) .^ (-34.1632 / 2.0);
            P(range8 | range9) = exp(polyval([1.340543e-6, -4.789660e-4, 6.424731e-2, -4.011801, 95.571899], h(range8 | range9)));

            % Water vapor density (ITU-R midWinter model - unchanged)
            rho = 3.4742 * exp(polyval([0.0004489, -0.03604, -0.2697, 0], h)); % [g/m^3]
            rho(h > 10) = 0; % ITU-R midWinter cutoff

        case "spring" % Boston spring conditions
            
            % Boston spring climate data  
            T0 = 48.4+273.15;  % [K] Boston average spring temperature (48.4°F)
            P0 = 1014;  % [hPa] Boston average spring pressure

            % Temperature profile (ITU-R structure, Boston surface temperature)
            T(range1) = T0 - 6.5 * hh(range1);
            T(range2) = 216.65;
            T(range3) = 216.65 + (hh(range3) - 20);
            T(range4) = 228.65 + 2.8 * (hh(range4) - 32);
            T(range5) = 270.65;
            T(range6) = 270.65 - 2.8 * (hh(range6) - 51);
            T(range7) = 214.65 - 2.0 * (hh(range7) - 71);
            T(range8) = 186.8673;
            T(range9) = 263.1905 - 76.3232 * sqrt(1 - ((h(range9) - 91) / 19.9429) .^ 2);

            % Pressure profile (ITU-R methodology, Boston surface pressure)
            P(range1) = P0 * (T0 ./ T(range1)) .^ (-34.1632 / 6.5);
            P(range2) = 226.3226 * exp(-34.1632 * (hh(range2) - 11) ./ T(range2));
            P(range3) = 54.74980 * (216.65 ./ T(range3)) .^ 34.1632;
            P(range4) = 8.680422 * (228.65 ./ T(range4)) .^ (34.1632 / 2.8);
            P(range5) = 1.109106 * exp(-34.1632 * (hh(range5) - 47) ./ T(range5));
            P(range6) = 0.6694167 * (270.65 ./ T(range6)) .^ (-34.1632 / 2.8);
            P(range7) = 0.03956649 * (214.65 ./ T(range7)) .^ (-34.1632 / 2.0);
            P(range8 | range9) = exp(polyval([1.340543e-6, -4.789660e-4, 6.424731e-2, -4.011801, 95.571899], h(range8 | range9)));

            % Water vapor density (interpolated between ITU-R winter/summer)
            % Interpolate between midWinter and midSummer scaling factors
            rho_scale = (3.4742 + 14.3542) / 2; % Average of winter/summer scaling: 8.9142
            % Interpolate polynomial coefficients
            winter_coeffs = [0.0004489, -0.03604, -0.2697, 0];
            summer_coeffs = [0.001007, -0.02290, -0.4174, 0];
            spring_coeffs = (winter_coeffs + summer_coeffs) / 2;
            rho = rho_scale * exp(polyval(spring_coeffs, h)); % [g/m^3]
            rho(h > 12) = 0; % Intermediate cutoff height

        case "summer" % Boston summer conditions
            
            % Boston summer climate data
            T0 = 71.6+273.15;  % [K] Boston average summer temperature (71.6°F)
            P0 = 1008;  % [hPa] Boston average summer pressure

            % Temperature profile (ITU-R structure, Boston surface temperature)
            T(range1) = T0 - 6.5 * hh(range1);
            T(range2) = 216.65;
            T(range3) = 216.65 + (hh(range3) - 20);
            T(range4) = 228.65 + 2.8 * (hh(range4) - 32);
            T(range5) = 270.65;
            T(range6) = 270.65 - 2.8 * (hh(range6) - 51);
            T(range7) = 214.65 - 2.0 * (hh(range7) - 71);
            T(range8) = 186.8673;
            T(range9) = 263.1905 - 76.3232 * sqrt(1 - ((h(range9) - 91) / 19.9429) .^ 2);

            % Pressure profile (ITU-R methodology, Boston surface pressure)
            P(range1) = P0 * (T0 ./ T(range1)) .^ (-34.1632 / 6.5);
            P(range2) = 226.3226 * exp(-34.1632 * (hh(range2) - 11) ./ T(range2));
            P(range3) = 54.74980 * (216.65 ./ T(range3)) .^ 34.1632;
            P(range4) = 8.680422 * (228.65 ./ T(range4)) .^ (34.1632 / 2.8);
            P(range5) = 1.109106 * exp(-34.1632 * (hh(range5) - 47) ./ T(range5));
            P(range6) = 0.6694167 * (270.65 ./ T(range6)) .^ (-34.1632 / 2.8);
            P(range7) = 0.03956649 * (214.65 ./ T(range7)) .^ (-34.1632 / 2.0);
            P(range8 | range9) = exp(polyval([1.340543e-6, -4.789660e-4, 6.424731e-2, -4.011801, 95.571899], h(range8 | range9)));

            % Water vapor density (ITU-R midSummer model - unchanged)
            rho = 14.3542 * exp(polyval([0.001007, -0.02290, -0.4174, 0], h)); % [g/m^3]
            rho(h > 15) = 0; % ITU-R midSummer cutoff

        case "fall" % Boston fall conditions
            
            % Boston fall climate data
            T0 = 55.1+273.15;  % [K] Boston average fall temperature (55.1°F)
            P0 = 1018.1;  % [hPa] Boston average fall pressure

            % Temperature profile (ITU-R structure, Boston surface temperature)
            T(range1) = T0 - 6.5 * hh(range1);
            T(range2) = 216.65;
            T(range3) = 216.65 + (hh(range3) - 20);
            T(range4) = 228.65 + 2.8 * (hh(range4) - 32);
            T(range5) = 270.65;
            T(range6) = 270.65 - 2.8 * (hh(range6) - 51);
            T(range7) = 214.65 - 2.0 * (hh(range7) - 71);
            T(range8) = 186.8673;
            T(range9) = 263.1905 - 76.3232 * sqrt(1 - ((h(range9) - 91) / 19.9429) .^ 2);

            % Pressure profile (ITU-R methodology, Boston surface pressure)
            P(range1) = P0 * (T0 ./ T(range1)) .^ (-34.1632 / 6.5);
            P(range2) = 226.3226 * exp(-34.1632 * (hh(range2) - 11) ./ T(range2));
            P(range3) = 54.74980 * (216.65 ./ T(range3)) .^ 34.1632;
            P(range4) = 8.680422 * (228.65 ./ T(range4)) .^ (34.1632 / 2.8);
            P(range5) = 1.109106 * exp(-34.1632 * (hh(range5) - 47) ./ T(range5));
            P(range6) = 0.6694167 * (270.65 ./ T(range6)) .^ (-34.1632 / 2.8);
            P(range7) = 0.03956649 * (214.65 ./ T(range7)) .^ (-34.1632 / 2.0);
            P(range8 | range9) = exp(polyval([1.340543e-6, -4.789660e-4, 6.424731e-2, -4.011801, 95.571899], h(range8 | range9)));

            % Water vapor density (interpolated between ITU-R summer/winter)
            % Closer to winter values since fall is transitioning to winter
            rho_scale = (14.3542 + 3.4742) / 2; % Average of summer/winter: 8.9142
            % Interpolate polynomial coefficients (closer to winter)
            summer_coeffs = [0.001007, -0.02290, -0.4174, 0];
            winter_coeffs = [0.0004489, -0.03604, -0.2697, 0];
            fall_coeffs = (summer_coeffs + winter_coeffs) / 2;
            rho = rho_scale * exp(polyval(fall_coeffs, h)); % [g/m^3]
            rho(h > 11) = 0; % Intermediate cutoff height

        otherwise
            error("Choose a valid season: 'winter', 'spring', 'summer', or 'fall'.");
    end

    % Calculate water vapor pressure from density (unchanged from ITU-R P.835)
    e = rho .* T / 216.7; % [hPa] Water-vapor pressure
    
    % Apply mixing ratio constraint (unchanged from ITU-R P.835)
    r = (e ./ P < 2e-6); % Mixing ratio limit
    e(r) = 2e-6 * P(r);
end