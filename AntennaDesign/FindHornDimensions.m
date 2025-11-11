function [chi, rho_e, rho_h, a1, b1, p_e, p_h] = FindHornDimensions(G0_in, a_in, b_in, freq_Hz)
% PYRAMIDAL_HORN_FROM_GAIN  Design an optimal pyramidal horn from target gain.
% Implements Balanis eqs. (13-56), (rho_e/λ=χ), (rho_h/λ = G0^2/(8π^3·1/χ)),
% a1 = sqrt(3 λ ρ_h), b1 = sqrt(2 λ ρ_e), and (13-49a,b) for p_e, p_h.
%
% Output struct fields:
%   chi, rho_e, rho_h, a1, b1, p_e, p_h
% Notes:
% - All geometry results are in meters.
% - Real solutions require a1>a and b1>b and the square-root arguments >0.

c = physconst('LightSpeed');
lambda = c/freq_Hz
G0_lin = 10.^(G0_in/10);
%G0_lin = G0_in
% -------- solve for chi (other terms depend on this value)
chi = solve_chi_1356(a_in, b_in, freq_Hz, G0_lin);

% -------- subsequent equations
rho_e = chi * lambda;                 % ρ_e
rho_h = (G0_lin^2)/(8*pi^3) * (1/chi) * lambda;  % ρ_h

a1 = sqrt(3*lambda*rho_h);            % aperture H-plane dimension
b1 = sqrt(2*lambda*rho_e);            % aperture E-plane dimension

% -------- pe, ph (13-49a,b)
arg_e = (rho_e/b1)^2 - 1/4;
arg_h = (rho_h/a1)^2 - 1/4;

if any([a1<=a_in, b1<=b_in, arg_e<=0, arg_h<=0])
    warning(['Geometry may be infeasible for this G0/a/b/λ: ' ...
             'check that a1>a, b1>b and sqrt-arguments are positive.']);
end

p_e = (b1 - b_in) * sqrt(arg_e);
p_h = (a1 - a_in) * sqrt(arg_h);

end





%---------------------Archived Method 1---------------------
%Source: https://www.youtube.com/watch?v=ISKvYdAffhU 
%Note: keep Ph_e < 0.25 and Ph_h < 0.4

%Step 1: define (enter as input arguments) desired beamwidth 
%Step 2: calculate required ApertureE and ApertureH
%Step 3: calculate the corresponding Le based on tolerable Ph_e (< 0.25)
%Step 4: calculate the corresponding Lh based on tolerable Ph_h (< 0.4)
%Step 5: verify via FEKO (NOT included in this script!)
%---------------------Constants---------------------
%c  = physconst('LightSpeed');
%freq = 225e9;
%lambda = c./freq;
%Ph_e = 0.24;
%Ph_h = 0.39;
%---------------------Calculations---------------------
%ApertureE = 56*lambda/HPBW_e;
%ApertureH = 67*lambda/HPBW_h;
%Le = ApertureE^2/(8*lambda*Ph_e);
%Lh = ApertureH^2/(8*lambda*Ph_h);
%end






%---------------------Archived Method 2 (Info From Other Online Source)---------------------
%Le =  %is the slant length of the side in the E-field direction
%Lh = %is the slant length of the side in the H-field direction
%L = %is the slant length of the cone from the apex.
%eff = 0.6 %aperture efficiency (value between 0 and 1)

%ApertureE = sqrt(2*lambda*Le) %rectangular
%ApertureH = sqrt(3*lambda*Lh) %rectangular 
%CHorn_d= sqrt(3*lambda*L) %diameter of connical horn

%Gain_Pyrmd = efficiency*(4*pi*A)/(lambda^2)
%Gain_con = (pi*CHorn_d/lambda)^2*eff