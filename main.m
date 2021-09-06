clear
close
clc 

% If the data file doesn't exist then
% run the script, load the resulting data for initial values as a struct
if(isfile('data.mat'))
    inits = load('data.mat');
else
    run('givenData.m')
    inits = load('data.mat');
end

endTime = 100; % seconds, should reach ground before this
elemCount = size((0:inits.dt:endTime)',1); 

% NED to Body Transformation Matrix
Lbe = @(psi,theta,phi) [cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta);...
    sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi), sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi), sin(phi)*cos(theta);...
    cos(phi)*sin(theta)*cos(phi)+sin(phi)*sin(psi), cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi), cos(phi)*cos(theta)];

%% Initialization of states

% Run interpolation and get polynomial coefficients for aerodynamic
% coefficients
polys = interpData(inits);

% Initialize empty arrays with the maximum expected size
[states.Vm_mpers, states.mach, states.alpha, states.beta, states.u,...
    states.v, states.w, states.p, states.q, states.r, states.phi,...
    states.theta, states.psi, states.earthPosX, states.earthPosY,...
    states.earthPosZ, states.earthPosAlt] = deal(zeros(1,elemCount + 1));

% Get the initial speed of the projectile
states.Vm_mpers(1) = inits.Vm_mpers;
states.mach(1) = inits.mach;

% Initial Angle of Attack and Sideslip Angle
states.alpha(1) = inits.alpha;
states.beta(1) = inits.beta;

% Given Beta value and known V_t we can find the v value
states.u(1) = cos(inits.alpha)*cos(inits.beta)*states.Vm_mpers(1);
states.v(1) = -cos(inits.alpha)*sin(inits.beta)*states.Vm_mpers(1);
states.w(1) = sin(inits.alpha)*states.Vm_mpers(1);

% Initial angular rates and Euler angles
states.p(1) = inits.p;
states.q(1) = inits.q;
states.r(1) = inits.r;

states.phi(1) = inits.phi;
states.theta(1) = inits.theta;
states.psi(1) = inits.psi;

% Fin deflection Angles
states.def_de = inits.def_de;
states.def_dr = inits.def_dr;
states.def_da = inits.def_da;

% NED axis positions
states.earthPosX(1) = inits.xm;
states.earthPosY(1) = inits.ym;
states.earthPosAlt(1) = inits.hm;
states.earthPosZ(1) = inits.zm;

% Initial atmospheric properties and aerodynamic coefficients
[rho, T, sos] = altitudeProp(states.earthPosAlt(1), inits);
Qd = 0.5*rho*states.Vm_mpers(1)^2;
aeroCoef = aeroCoefs(states, inits, polys, 1);

%% Iterative Solution
for i = 1:elemCount
    % Aeropropulsive forces and moments --> [Fx; Fy; Fz] & [L; M; N]
    forces = transpose(Qd*inits.A.*aeroCoef(1:3));
    moments = transpose(Qd*inits.A*inits.d*aeroCoef(4:6));

    % Apply NED -> Body matrix to gravitational force vector 
    gravOnBody = Lbe(states.psi(i), states.theta(i), states.phi(i))*[0; 0; inits.m*inits.g];
    forces = forces + gravOnBody;
   
    % Rate of change for angular rates p q r --> [pDot; qDot; rDot]
    % Found via Eqn 2.37
    angRateDots = inits.I\moments -...
       inits.I\(cross([states.p(i); states.q(i); states.r(i)], inits.I*[states.p(i); states.q(i); states.r(i)]));
   
    % Rate of change for velocities in body axis -> [uDot; vDot; wDot]
    % Found using Eqn 2.29-2.30-2.31
    velDots = [forces(1)/inits.m - states.w(i)*states.q(i) + states.v(i)*states.r(i);...
        forces(2)/inits.m - states.u(i)*states.r(i) + states.w(i)*states.p(i);...
        forces(3)/inits.m - states.v(i)*states.p(i) + states.u(i)*states.q(i)];
    
    % Rate of change for Euler Angles -> [psiDot; thetaDot; phiDot]
    % As given in Eqn 2.19-2.20-2.21
    eulerDots = [(states.q(i)*sin(states.phi(i)) + states.r(i)*cos(states.phi(i)))/cos(states.theta(i));...
                 states.q(i)*cos(states.phi(i)) - states.r(i)*sin(states.phi(i));...
                 states.p(i) + (states.q(i)*sin(states.phi(i)) + states.r(i)*cos(states.phi(i)))*tan(states.theta(i))]; % psi ; theta; phi
            
    % Find the velocities in body axis for the next step using Euler's Method   
    states.u(i+1) = states.u(i) + velDots(1)*inits.dt;
    states.v(i+1) = states.v(i) + velDots(2)*inits.dt;
    states.w(i+1) = states.w(i) + velDots(3)*inits.dt;

    % Total speed
    states.Vm_mpers(i+1) = sqrt((states.u(i+1))^2+(states.v(i+1))^2+(states.w(i+1))^2);
    
    % Convert the velocity in body axis to NED axis using the
    % transformation matrix
    velOnNED = transpose(Lbe(states.psi(i), states.theta(i), states.phi(i)))*[states.u(i); states.v(i); states.w(i)];

    % Find the next positions in NED axis using Euler's Method and
    % previously found velocities in NED axis
    states.earthPosX(i+1) = states.earthPosX(i) + velOnNED(1)*inits.dt;
    states.earthPosY(i+1) = states.earthPosY(i) + velOnNED(2)*inits.dt;
    states.earthPosZ(i+1) = states.earthPosZ(i) + velOnNED(3)*inits.dt;
    states.earthPosAlt(i+1) = -states.earthPosZ(i+1);
   
    % Find the Euler Angles in the next iteration
    states.psi(i+1) = states.psi(i) + eulerDots(1)*inits.dt;
    states.theta(i+1) = states.theta(i) + eulerDots(2)*inits.dt;
    states.phi(i+1) = states.phi(i) + eulerDots(3)*inits.dt;
   
    % Find the angular rates in the next step
    states.p(i+1) = states.p(i) + angRateDots(1)*inits.dt;
    states.q(i+1) = states.q(i) + angRateDots(2)*inits.dt;
    states.r(i+1) = states.r(i) + angRateDots(3)*inits.dt;
   
    % Find the Angle of Attack and Sideslip angle in the next step
    states.alpha(i+1) = atan(states.w(i+1)/states.u(i+1));
    states.beta(i+1) = atan(states.v(i+1)/states.u(i+1));
   
    % Calculate the next atmospheric propoerties, convert the velocity into
    % Mach number and get the new aerodynamic coefficients
    [rho, T, sos] = altitudeProp(states.earthPosAlt(i+1), inits);
    Qd = 0.5*rho*states.Vm_mpers(i+1)^2;
    states.mach(i+1) = states.Vm_mpers(i+1)/sos;
    aeroCoef = aeroCoefs(states, inits, polys, i+1);

    % Terminate the iteration if the projectile is below sea level
    if states.earthPosAlt(i+1) <= 0
        break
    end
end

% Array for time values to use in plots
timeArr = 0:inits.dt:(i-1)*inits.dt;

%% Plots

% 3D Trajectory Plot
figure('Name', 'X-Y-H, 3D Plot')
set(gcf, 'WindowState', 'maximized');
plot3(states.earthPosX(1:i), states.earthPosY(1:i), states.earthPosAlt(1:i))
title('Trajectory','FontSize',14)
xlabel('X Axis (m)','FontSize',12)
ylabel('Y Axis (m)','FontSize',12)
zlabel('Altitude (m)','FontSize',12)

% X-H Plot
figure('Name', 'X-H Plot')
set(gcf, 'WindowState', 'maximized');
plot(states.earthPosX(1:i), states.earthPosAlt(1:i));
title('X-H Plot','FontSize',14)
xlabel('X Axis (m)','FontSize',12)
ylabel('Altitude (m)','FontSize',12)

% Y-H Plot
figure('Name', 'Y-H Plot')
set(gcf, 'WindowState', 'maximized');
plot(states.earthPosY(1:i), states.earthPosAlt(1:i));
title('Y-H Plot','FontSize',14)
xlabel('Y Axis (m)','FontSize',12)
ylabel('Altitude (m)','FontSize',12)

% X-Y Plot
figure('Name', 'X-Y Plot')
set(gcf, 'WindowState', 'maximized');
plot(states.earthPosX(1:i), states.earthPosY(1:i));
title('X-Y Plot','FontSize',14)
xlabel('X Axis (m)','FontSize',12)
ylabel('Y Axis (m)','FontSize',12)

% Body Angular Rates Plots
figure('Name', 'p-q-r Body Angular Rates Plots')
set(gcf, 'WindowState', 'maximized');
subplot(3,1,1)
plot(timeArr, states.p(1:i))
title('p vs Time','FontSize',14)
xlim([0 timeArr(end)])
xlabel('Time (s)','FontSize',12)
ylabel('p (rad/s)','FontSize',12)
hold on

subplot(3,1,2)
plot(timeArr, states.q(1:i))
title('q vs Time','FontSize',14)
xlim([0 timeArr(end)])
xlabel('Time (s)','FontSize',12)
ylabel('q (rad/s)','FontSize',12)

subplot(3,1,3)
plot(timeArr, states.r(1:i))
title('r vs Time','FontSize',14)
xlim([0 timeArr(end)])
xlabel('Time (s)','FontSize',12)
ylabel('r (rad/s)','FontSize',12)

% u-v-w Body Velocities Plots
figure('Name', 'u-v-w Velocities in Body Frame')
set(gcf, 'WindowState', 'maximized');
subplot(3,1,1)
plot(timeArr, states.u(1:i))
title('u vs Time','FontSize',14)
xlim([0 timeArr(end)])
xlabel('Time (s)','FontSize',12)
ylabel('u (m/s)','FontSize',12)
hold on

subplot(3,1,2)
plot(timeArr, states.v(1:i))
title('v vs Time','FontSize',14)
xlim([0 timeArr(end)])
xlabel('Time (s)','FontSize',12)
ylabel('v (m/s)','FontSize',12)

subplot(3,1,3)
plot(timeArr, states.w(1:i))
title('w vs Time','FontSize',14)
xlim([0 timeArr(end)])
xlabel('Time (s)','FontSize',12)
ylabel('w (m/s)','FontSize',12)

% Euler Angles Plots
figure('Name', 'Euler Angles Plot')
set(gcf, 'WindowState', 'maximized');
subplot(3,1,1)
plot(timeArr,states.phi(1:i))
title('\phi vs Time','FontSize',14)
xlim([0 timeArr(end)])
xlabel('Time (s)','FontSize',12)
ylabel('\phi (rad)','FontSize',12)
hold on

subplot(3,1,2)
plot(timeArr,states.theta(1:i))
title('\theta vs Time','FontSize',14)
xlim([0 timeArr(end)])
xlabel('Time (s)','FontSize',12)
ylabel('\theta (rad)','FontSize',12)

subplot(3,1,3)
plot(timeArr,states.psi(1:i))
title('\psi vs Time','FontSize',14)
xlim([0 timeArr(end)])
xlabel('Time (s)','FontSize',12)
ylabel('\psi (rad)','FontSize',12)

% Mach Number - Angle of Attack and Sideslip Angle Plots
figure('Name', 'Mach Number - Alpha - Beta Plots')
set(gcf, 'WindowState', 'maximized');
subplot(3,1,1)
plot(timeArr, states.mach(1:i))
title('Mach Number vs Time','FontSize',14)
xlim([0 timeArr(end)])
xlabel('Time (s)','FontSize',12)
ylabel('Mach','FontSize',12)
hold on

subplot(3,1,2)
plot(timeArr, states.alpha(1:i))
title('\alpha vs Time','FontSize',14)
xlim([0 timeArr(end)])
xlabel('Time (s)','FontSize',12)
ylabel('\alpha (rad)','FontSize',12)

subplot(3,1,3)
plot(timeArr, states.beta(1:i))
title('\beta vs Time','FontSize',14)
xlim([0 timeArr(end)])
xlabel('Time (s)','FontSize',12)
ylabel('\beta (rad)','FontSize',12)

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

% Obtain interpolation coefficients
function interped = interpData(consts)
    interped.Cd = griddedInterpolant(consts.Machpoints, consts.Cd_data);
    
    interped.Cza = griddedInterpolant(consts.Machpoints, consts.Cza_data);
    
    interped.Czq = griddedInterpolant(consts.Machpoints, consts.Czq_data);
    
    interped.Cma = griddedInterpolant(consts.Machpoints, consts.Cma_data);
    
    interped.Cmq = griddedInterpolant(consts.Machpoints, consts.Cmq_data);
    
    interped.Clp = griddedInterpolant(consts.Machpoints, consts.Clp_data);
    
    interped.Czd = griddedInterpolant(consts.Machpoints, consts.Czd_data);
    
    interped.Cmd = griddedInterpolant(consts.Machpoints, consts.Cmd_data);
    
    interped.Cld = griddedInterpolant(consts.Machpoints, consts.Cld_data);
end
    
% Aerodynamic Coefficients
function coefs = aeroCoefs(state, init, polyCoefs, i)
    Cx = polyCoefs.Cd(state.mach(i));
    
    Cy = polyCoefs.Cza(state.mach(i))*state.beta(i) + polyCoefs.Czd(state.mach(i))*state.def_dr - ...
    polyCoefs.Czq(state.mach(i))*state.r(i)*init.d/(2*state.Vm_mpers(i));

    Cz = polyCoefs.Cza(state.mach(i))*state.alpha(i) + polyCoefs.Czd(state.mach(i))*state.def_de + ...
        polyCoefs.Czq(state.mach(i))*state.q(i)*init.d/(2*state.Vm_mpers(i));
    
    Cl = polyCoefs.Cld(state.mach(i))*state.def_da +...
        polyCoefs.Clp(state.mach(i))*state.p(i)*init.d/(2*state.Vm_mpers(i));
    
    Cm = polyCoefs.Cma(state.mach(i))*state.alpha(i) + polyCoefs.Cmd(state.mach(i))*state.def_de + ...
        polyCoefs.Cmq(state.mach(i))*state.q(i)*init.d/(2*state.Vm_mpers(i));
    
    Cn = -polyCoefs.Cma(state.mach(i))*state.beta(i) - polyCoefs.Cmd(state.mach(i))*state.def_dr + ...
        polyCoefs.Cmq(state.mach(i))*state.r(i)*init.d/(2*state.Vm_mpers(i));
    
    coefs = [Cx, Cy, Cz, Cl, Cm, Cn];
    
end 