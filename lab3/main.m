clear variables; clc;
addpath('/home/omarkahol/MATLAB/ROTOR_AERODYNAMICS/lab3/lib');

dataFile;
helicopter = buildHelicopter(R,R0,c,Nb);
N = 500;

helicopter = makeMesh(helicopter,N);

h = 0; %altitude
theta = 12; %degrees
theta1s = 5; %degrees
theta1c = 1; %degrees
beta = 3;
alfa_sh = 1;
V = 10;



conf = flightConfiguration(h,theta,RPM,V,theta1s,theta1c,beta,alfa_sh);
hsol = hoveringSolver(helicopter,conf);

sol = forewardFlightSolver(hsol,helicopter,conf,100);

