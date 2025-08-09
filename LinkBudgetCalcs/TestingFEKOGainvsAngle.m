clear, clc, close all;
addpath(genpath('./functions'))
addpath('C:\Users\rache\Projects\CubeSat\TeraLink-1\GroundStation\Simulations\MATLAB Sims\SubTHzLinkBudget\DataSourceFiles')

[hpbw, points] = FEKO_Gain_vs_Angle('0.1_000Final_HPBW_WorkingModel.dat');
