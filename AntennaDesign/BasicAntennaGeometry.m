function [theta0_deg, f_over_d, f_m, d_m] = BasicAntennaGeometry(varargin)
% BasicAntennaGeometry  Compute parabolic reflector geometry from θ₀, f/D, f, and/or D.
%
% Usage:
%   BasicAntennaGeometry('theta0_deg', 46.1, 'f_m', 0.68)
%   BasicAntennaGeometry('theta0_deg', 46.1, 'd_m', 1.7)
%   BasicAntennaGeometry('f_over_d', 0.42, 'f_m', 0.68)
%   BasicAntennaGeometry('f_over_d', 0.42, 'd_m', 1.7)
%   BasicAntennaGeometry('f_m', 0.68, 'd_m', 1.7)
%
% Balanis Eq. 15-25:
%   f = (D/4)*cot(θ₀/2)  →  f/D = (1/4)*cot(θ₀/2)
%   θ₀ = 2*atan(1 / (4*(f/D)))   [degrees via cotd/atand]

% ------------------------------
% Parse name–value pairs
% ------------------------------
theta0_in = [];
fd_in     = [];
f_in      = [];
d_in      = [];

for k = 1:2:numel(varargin)
    name = lower(string(varargin{k}));
    val  = varargin{k+1};
    switch name
        case "theta0_deg", theta0_in = val;
        case "f_over_d",   fd_in     = val;
        case "f_m",        f_in      = val;
        case "d_m",        d_in      = val;
        otherwise, error('Unknown parameter name "%s".', varargin{k});
    end
end

% ------------------------------
% Presence flags
% ------------------------------
have_theta = ~isempty(theta0_in);
have_fd    = ~isempty(fd_in);
have_f     = ~isempty(f_in);
have_d     = ~isempty(d_in);

% Cannot specify both θ₀ and f/D
if have_theta && have_fd
    error('Specify only one of "theta0_deg" or "f_over_d".');
end

% ------------------------------
% Conversion relationships
% ------------------------------
theta_from_fd = @(fd) 2*atand(1./(4*fd));
fd_from_theta = @(th) 0.25 * cotd(th/2);

% ------------------------------
% Compute missing quantities
% ------------------------------
theta0_deg = NaN; f_over_d = NaN; f_m = NaN; d_m = NaN;

if have_theta
    theta0_deg = theta0_in;
    f_over_d   = fd_from_theta(theta0_deg);
    if have_f
        f_m = f_in;
        d_m = f_m / f_over_d;
    elseif have_d
        d_m = d_in;
        f_m = f_over_d * d_m;
    end

elseif have_fd
    f_over_d   = fd_in;
    theta0_deg = theta_from_fd(f_over_d);
    if have_f
        f_m = f_in;
        d_m = f_m / f_over_d;
    elseif have_d
        d_m = d_in;
        f_m = f_over_d * d_m;
    end

elseif have_f && have_d
    f_m = f_in;
    d_m = d_in;
    f_over_d = f_m / d_m;
    theta0_deg = theta_from_fd(f_over_d);

else
    error('Not enough information provided to compute geometry.');
end

% ------------------------------
% Output summary
% ------------------------------
fprintf('---------------------------------------------\n');
fprintf('Reflector Geometry (Balanis Eq. 15-25)\n');
fprintf('---------------------------------------------\n');
if have_theta, fprintf('Input: θ₀ = %.6f°\n', theta0_in); end
if have_fd,    fprintf('Input: f/D = %.8f\n', fd_in); end
if have_f,     fprintf('Input: f = %.6f m\n', f_in); end
if have_d,     fprintf('Input: D = %.6f m\n', d_in); end
fprintf('---------------------------------------------\n');
fprintf('Computed θ₀   = %.6f°\n', theta0_deg);
fprintf('Computed f/D  = %.8f\n', f_over_d);
fprintf('Computed f    = %.6f m\n', f_m);
fprintf('Computed D    = %.6f m\n', d_m);
fprintf('---------------------------------------------\n');
end
