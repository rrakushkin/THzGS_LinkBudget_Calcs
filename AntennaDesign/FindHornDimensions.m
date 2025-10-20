function [ApertureE,ApertureH,Le,Lh] = FindHornDimensions(HPBW_e,HPBW_h)
%---------------------Function Description---------------------
%Source: https://www.youtube.com/watch?v=ISKvYdAffhU 
%Note: keep Ph_e < 0.25 and Ph_h < 0.4

%Step 1: define (enter as input arguments) desired beamwidth 
%Step 2: calculate required ApertureE and ApertureH
%Step 3: calculate the corresponding Le based on tolerable Ph_e (< 0.25)
%Step 4: calculate the corresponding Lh based on tolerable Ph_h (< 0.4)
%Step 5: verify via FEKO (NOT included in this script!)
%---------------------Constants---------------------
c  = physconst('LightSpeed');
freq = 225e9;
lambda = c./freq;
Ph_e = 0.24;
Ph_h = 0.39;
%---------------------Calculations---------------------
ApertureE = 56*lambda/HPBW_e;
ApertureH = 67*lambda/HPBW_h;
Le = ApertureE^2/(8*lambda*Ph_e);
Lh = ApertureH^2/(8*lambda*Ph_h);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------Archived Lines (Info From Other Online Source)---------------------
%Le =  %is the slant length of the side in the E-field direction
%Lh = %is the slant length of the side in the H-field direction
%L = %is the slant length of the cone from the apex.
%eff = 0.6 %aperture efficiency (value between 0 and 1)

%ApertureE = sqrt(2*lambda*Le) %rectangular
%ApertureH = sqrt(3*lambda*Lh) %rectangular 
%CHorn_d= sqrt(3*lambda*L) %diameter of connical horn

%Gain_Pyrmd = efficiency*(4*pi*A)/(lambda^2)
%Gain_con = (pi*CHorn_d/lambda)^2*eff