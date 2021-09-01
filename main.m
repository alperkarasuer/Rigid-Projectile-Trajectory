clear all
close all
clc 

% Run the script for initial values and load them as strut
run('givenData.m');
initialVals = load('data.mat');


% Run interpolation and get polynomial coefficients from aerodynamic
% coefficients
polys = interpData(initialVals);

% Get the initial speed of the projectile
states.Vm_mpers(1) = initialVals.Vm_mpers;

% Given Beta value is zero, using the approximation Beta ~ v/u we can see that
% v = 0
states.v(1) = 0;

% We also have Alpha ~ w/u approximation and Alpha is given as 10, we use
% this inside V_t equation to find u and w
states.u(1) = fzero(@(u) sqrt(u^2 + (initialVals.alpha*u)^2) - states.Vm_mpers(1),[0 1E4]);
states.w(1) = states.u(1)*initialVals.alpha;






% Density, Temperature and Speed of Sound in given altitude
function [rho, T, speedOfSound] = altitudeProp(h, consts)
    if((h) <= 10000)
        rho = consts.rho0*(1 - 0.00002256*(h))^4.256;
        T = consts.T0*(1 - 0.0002256*(h));
    else
        rho = 0.412*exp(-0.000151*(h-10000));
        T = 0.7744*T0;
    end
    speedOfSound = sqrt(consts.k*consts.R*T) 
end

% Obtain interpolation coefficients to use with ppval
function interped = interpData(consts)
    interped.Cd = interp1(consts.Machpoints, consts.Cd_data,'spline','pp');
    
    interped.Cza = interp1(consts.Machpoints, consts.Cza_data,'spline','pp');
    
    interped.Czq = interp1(consts.Machpoints, consts.Czq_data,'spline','pp');
    
    interped.Cma = interp1(consts.Machpoints, consts.Cma_data,'spline','pp');
    
    interped.Cmq = interp1(consts.Machpoints, consts.Cmq_data,'spline','pp');
    
    interped.Clp = interp1(consts.Machpoints, consts.Clp_data,'spline','pp');
    
    interped.Czd = interp1(consts.Machpoints, consts.Czd_data,'spline','pp');
    
    interped.Cmd = interp1(consts.Machpoints, consts.Cmd_data,'spline','pp');
    
    interped.Cld = interp1(consts.Machpoints, consts.Cld_data,'spline','pp');
end



    