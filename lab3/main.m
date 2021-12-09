clear variables; clc;
addpath('/home/omarkahol/MATLAB/ROTOR_AERODYNAMICS/lab3/lib');

dataFile;
helicopter = buildHelicopter(R,R0,c,Nb);
N = 50;

helicopter = makeMesh(helicopter,N);

h = 0; %altitude
theta = 8.2391; %degrees
theta1s = 2.53; %degrees
theta1c = 0; %degrees
beta = 1.9308;
alfa_sh = -1.0694;
V = 12.9857;



conf = flightConfiguration(h,theta,RPM,V,theta1s,theta1c,beta,alfa_sh);
hsol = hoveringSolver(helicopter,conf);

sol = forewardFlightSolver(hsol,helicopter,conf,50);
data = postProcessor(helicopter,conf,sol);

