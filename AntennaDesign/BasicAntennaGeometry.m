function [theta0_deg, f_over_d, f_m, d_m] = BasicAntennaGeometry(varargin)
% ReflectorGeometry  Compute parabolic reflector geometry from θ₀, f/d, or f.
%
%   Usage examples:
%       ReflectorGeometry('theta0_deg', 46.1); <<< here θ₀ is given so f/d is derived from it
%       ReflectorGeometry('f_over_D', 0.42); <<< here f/d is given so θ₀ is derived from it
%       ReflectorGeometry('f_m', 0.68);   % <<< here the focal length f_m
%       is given only so default D = 1.0 m is used to determine f/d and
%       corresponding θ₀ 
%
%   Formula (Balanis Eq. 15-25):
%       f = (D/4) * cot(theta0/2)
%       → f/D = (1/4)*cot(theta0/2)
%       → θ₀ = 2*atan(1 / (4*(f/D)))
%---------------------------------------------------

% --- Parse Inputs ---
p = inputParser;
addParameter(p, 'theta0_deg', []);
addParameter(p, 'f_over_d', []);
addParameter(p, 'f_m', []);
parse(p, varargin{:});

theta0_deg_in = p.Results.theta0_deg;
f_over_d_in   = p.Results.f_over_d;
f_m_in        = p.Results.f_m;

% --- Flags ---
have_theta = ~isempty(theta0_deg_in);
have_fD    = ~isempty(f_over_d_in);
have_f     = ~isempty(f_m_in);

% --- Case validation ---
if sum([have_theta, have_fD, have_f]) == 0
    error('Provide one of "theta0_deg", "f_over_D", or "f_m".');
elseif sum([have_theta, have_fD, have_f]) > 1
    warning('Multiple inputs provided — using first valid priority: theta0 > f/D > f.');
end

% -------------------------------------------------------------------------
% CASE 1: θ₀ given
% -------------------------------------------------------------------------
if have_theta
    theta0_deg = theta0_deg_in;
    f_over_d   = 0.25 * cotd(theta0_deg / 2);
    f_m = [];
    d_m = [];

% -------------------------------------------------------------------------
% CASE 2: f/D given
% -------------------------------------------------------------------------
elseif have_fD
    f_over_d   = f_over_d_in;
    theta0_deg = 2 * atand(1 / (4 * f_over_d));
    f_m = [];
    d_m = [];

% -------------------------------------------------------------------------
% CASE 3: f given (default D = 1.0 m)
% -------------------------------------------------------------------------
elseif have_f
    d_m        = 1.0;           % fixed default diameter
    f_m        = f_m_in;
    f_over_d   = f_m / d_m;
    theta0_deg = 2 * atand(1 / (4 * f_over_d));
end

% -------------------------------------------------------------------------
% Display results
% -------------------------------------------------------------------------
fprintf('---------------------------------------------\n');
fprintf('Reflector Geometry (Balanis Eq. 15-25)\n');
fprintf('---------------------------------------------\n');

if have_theta
    fprintf('Input: θ₀ = %.6f°\n', theta0_deg_in);
elseif have_fD
    fprintf('Input: f/D = %.8f\n', f_over_d_in);
elseif have_f
    fprintf('Input: f = %.6f m (default D = 1.0000 m used)\n', f_m_in);
end

fprintf('---------------------------------------------\n');
fprintf('Computed θ₀ = %.6f°\n', theta0_deg);
fprintf('Computed f/D = %.8f\n', f_over_d);
if have_f
    fprintf('Default D = %.6f m\n', d_m);
    fprintf('Focal length f = %.6f m\n', f_m);
end
fprintf('---------------------------------------------\n');

end
