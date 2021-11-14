clear variables; clc;
addpath('/home/omarkahol/MATLAB/ROTOR_AERODYNAMICS/lab2/lib');

dataFile;
helicopter = buildHelicopter(R,R0,c,Nb);
N = 500;

helicopter = makeMesh(helicopter,N);

h = 0; %altitude
theta = 12; %degrees
Vup = 4; %meters per second

conf = flightConfiguration(h,theta,RPM,Vup);
sol = solver(helicopter,conf);
data = postProcessor(helicopter,conf,sol);

%plots
figure(1);
plot(helicopter.Mesh/helicopter.Radius,sol.ThrustCoefficient,'k-','LineWidth',2);
title('dC_T / dR  in the blade');
xlabel('r [-]');
ylabel('dC_T / dR');

figure(2);
plot(helicopter.Mesh/helicopter.Radius,sol.TorqueCoefficient,'k-','LineWidth',2);
title('dQ_T / dR  in the blade');
xlabel('r [-]');
ylabel('dQ_T / dR');


figure(3);
plot(helicopter.Mesh/helicopter.Radius,sol.InflowRatio,'k-','LineWidth',2);
title('Total inflow ratio \Lambda  in the blade');
xlabel('r [-]');
ylabel('\Lambda');

figure(4);
plot(helicopter.Mesh/helicopter.Radius,sol.SwirlEffect,'k-','LineWidth',2);
title('Swirl Effect parameter in the blade');
xlabel('r [-]');
ylabel('a');

figure(5);
plot(helicopter.Mesh/helicopter.Radius,sol.InducedAngle,'k-','LineWidth',2);
title('Induced angle \phi in the blade');
xlabel('r [-]');
ylabel('\phi');


figure(6);
plot(helicopter.Mesh/helicopter.Radius,sol.IterationsRequired,'k-','LineWidth',2);
title('it necessary for convergence');
xlabel('r [-]');
ylabel('it');
