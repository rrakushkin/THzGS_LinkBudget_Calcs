function [T_interp, P_interp, e_interp] = InterpAtm(ref1, ref2, lat)
%this linearly interpolates the low global atmosphere
%conditions with the mid summer and mid winter atmosphere conditions (45
%deg summer and winter conditions) to determine atm conditions for GS

% GS on Egan Roof Coordinates: (42.3378054237531, -71.08874165317037) (lat,long)
% thus intended for: Reference1 = "Annual 15", Reference2 = "Summer 45" and
% output is interpolation giving you Summer at whatever angle your
% interested in

% ref1 and ref2 are cell arrays: {T, P, e}

% Default latitude
    if nargin < 3 || isempty(lat)
        lat = 42.3378054237531;
    end

    % Interpolation weight
    w = (45 - lat) / (45 - 15);

    % Extract components
    T1 = ref1{1}; P1 = ref1{2}; e1 = ref1{3};
    T2 = ref2{1}; P2 = ref2{2}; e2 = ref2{3};

    % Linear interpolation
    T_interp = w * T1 + (1 - w) * T2;
    P_interp = w * P1 + (1 - w) * P2;
    e_interp = w * e1 + (1 - w) * e2;
end