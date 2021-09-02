clear all
close all
clc 

% Run the script for initial values and load them as struct
run('givenData.m');
initialVals = load('data.mat');
endTime = 20; % seconds
elemCount = size((0:initialVals.dt:endTime)',1);

% NED to Body Transformation Matrix
Lbe = @(psi,theta,phi) [cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta);...
    sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi), sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi), sin(phi)*cos(theta);...
    cos(phi)*sin(theta)*cos(phi)+sin(phi)*sin(psi), cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi), cos(phi)*cos(theta)];

%% Initialization of states

% Run interpolation and get polynomial coefficients from aerodynamic
% coefficients
polys = interpData(initialVals);

% Get the initial speed of the projectile
states.Vm_mpers(1) = initialVals.Vm_mpers;
states.mach(1) = initialVals.mach;

% Initial Angle of Attack and Sideslip Angle
states.alpha(1) = initialVals.alpha;
states.beta(1) = initialVals.beta;

% Given Beta value and known V_t we can find the v value
states.u(1) = cos(initialVals.alpha)*cos(initialVals.beta)*states.Vm_mpers(1);
states.v(1) = -cos(initialVals.alpha)*sin(initialVals.beta)*states.Vm_mpers(1);
states.w(1) = sin(initialVals.alpha)*states.Vm_mpers(1);

% Initial angular rates and Euler angles
states.p(1) = initialVals.p;
states.q(1) = initialVals.q;
states.r(1) = initialVals.r;

states.phi(1) = initialVals.phi;
states.theta(1) = initialVals.theta;
states.psi(1) = initialVals.psi;

% Deflection Angles
states.def_de = initialVals.def_de;
states.def_dr = initialVals.def_dr;
states.def_da = initialVals.def_da;

% NED axis positions
states.earthPosX(1) = initialVals.xm;
states.earthPosY(1) = initialVals.ym;
states.earthPosAlt(1) = initialVals.hm;
states.earthPosZ(1) = initialVals.zm;

[rho, T, sos] = altitudeProp(states.earthPosAlt(1), initialVals);
Qd = 0.5*rho*states.Vm_mpers(1);
aeroCoef = aeroCoefs(states, initialVals, polys, 1);


%% Iterative Solution
for i = 1:elemCount   
   forces = transpose(Qd*initialVals.A.*aeroCoef(1:3));
   moments = transpose(Qd*initialVals.A*initialVals.d*aeroCoef(3:6));
   gravOnBody = Lbe(states.psi(i), states.theta(i), states.phi(i))*[0; 0; initialVals.g];
   forces = forces + initialVals.m*gravOnBody;
   
   angRateDots = [moments(1)/initialVals.Ix;...
       moments(2)/initialVals.Iy + states.r(i)*states.p(i)*(initialVals.Iy - initialVals.Ix)/initialVals.Iy;...
       moments(3)/initialVals.Iy + states.p(i)*states.q(i)*(initialVals.Ix - initialVals.Iy)/initialVals.Iy];
   
   velDots = [forces(1)/initialVals.m - states.w(i)*states.q(i) + states.v(i)*states.r(i);...
       forces(2)/initialVals.m - states.u(i)*states.r(i) + states.w(i)*states.p(i);...
       forces(3)/initialVals.m - states.v(i)*states.p(i) + states.u(i)*states.q(i)];
       
   eulerDots = [(states.q(i)*sin(states.phi(i)) + states.r(i)*cos(states.phi(i)))/cos(states.theta(i));...
                states.q(i)*cos(states.phi(i)) - states.r(i)*sin(states.phi(i));...
                states.p(i) + (states.q(i)*sin(states.phi(i)) + states.r(i)*cos(states.phi(i)))*tan(states.theta(i))];
            
   
   states.psi(i+1) = states.psi(i) + eulerDots(1)*initialVals.dt;
   states.theta(i+1) = states.theta(i) + eulerDots(2)*initialVals.dt;
   states.phi(i+1) = states.phi(i) + eulerDots(3)*initialVals.dt;
   
   states.p(i+1) = states.p(i) + angRateDots(1)*initialVals.dt;
   states.q(i+1) = states.q(i) + angRateDots(2)*initialVals.dt;
   states.r(i+1) = states.r(i) + angRateDots(3)*initialVals.dt;
   
   states.u(i+1) = states.u(i) + velDots(1)*initialVals.dt;
   states.v(i+1) = states.v(i) + velDots(2)*initialVals.dt;
   states.w(i+1) = states.w(i) + velDots(3)*initialVals.dt;
                
   velOnNED = transpose(Lbe(states.psi(i+1), states.theta(i+1), states.phi(i+1)))*[states.u(i+1); states.v(i+1); states.w(i+1)];
   
   states.earthPosX(i+1) = states.earthPosX(i) + velOnNED(1)*initialVals.dt;
   states.earthPosY(i+1) = states.earthPosY(i) + velOnNED(2)*initialVals.dt;
   states.earthPosZ(i+1) = states.earthPosZ(i) + velOnNED(3)*initialVals.dt;
   states.earthPosAlt(i+1) = states.earthPosAlt(i) - velOnNED(3)*initialVals.dt;
   
   states.Vm_mpers(i+1) = sqrt((states.u(i+1))^2+(states.v(i+1))^2+(states.w(i+1))^2);
   
   states.alpha(i+1) = atan(states.w(i+1)/states.u(i+1));
   states.beta(i+1) = asin(states.v(i+1)/states.Vm_mpers(i+1));
   
   [rho, T, sos] = altitudeProp(states.earthPosAlt(i+1), initialVals);
   Qd = 0.5*rho*states.Vm_mpers(i+1);
   states.mach(i+1) = states.Vm_mpers(i+1)/sos;
   aeroCoef = aeroCoefs(states, initialVals, polys, i+1);

end

%% Function Definitions

% Density, Temperature and Speed of Sound in given altitude
function [rho, T, speedOfSound] = altitudeProp(h, consts)
    if((h) <= 10000)
        rho = consts.rho0*(1 - 0.00002256*(h))^4.256;
        T = consts.T0*(1 - 0.00002256*(h));
    else
        rho = 0.412*exp(-0.000151*(h-10000));
        T = 0.7744*consts.T0;
    end
    speedOfSound = sqrt(consts.k*consts.R*T); 
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

% Aerodynamic Coefficients
function coefs = aeroCoefs(state, init, polyCoefs, i)
    Cx = ppval(polyCoefs.Cd, state.mach(i));
    
    Cy = ppval(polyCoefs.Cza, state.mach(i))*state.beta(i) +...
        ppval(polyCoefs.Czd, state.mach(i))*state.def_dr - ...
        ppval(polyCoefs.Czq, state.mach(i))*state.r(i)*init.d/(2*state.Vm_mpers(i));
    
    Cz = ppval(polyCoefs.Cza, state.mach(i))*state.alpha(i) + ppval(polyCoefs.Czd, state.mach(i))*state.def_de + ...
        ppval(polyCoefs.Czq, state.mach(i))*state.q(i)*init.d/(2*state.Vm_mpers(i));
    
    Cl = ppval(polyCoefs.Cld, state.mach(i))*state.def_da + ppval(polyCoefs.Clp, state.mach(i))*...
        state.p(i)*init.d/(2*state.Vm_mpers(i));
    
    Cm = ppval(polyCoefs.Cma, state.mach(i))-state.alpha(i) + ppval(polyCoefs.Cmd, state.mach(i))*state.def_de + ...
        ppval(polyCoefs.Cmq, state.mach(i))*state.q(i)*init.d/(2*state.Vm_mpers(i));
    
    Cn = -ppval(polyCoefs.Cma, state.mach(i))*state.beta(i) - ppval(polyCoefs.Cmd, state.mach(i))*state.def_dr + ...
        ppval(polyCoefs.Cmq, state.mach(i))*state.r(i)*init.d/(2*state.Vm_mpers(i));
    
    coefs = [Cx, Cy, Cz, Cl, Cm, Cn];
end


    