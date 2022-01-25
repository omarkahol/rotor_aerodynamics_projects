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

[X,Y,L,l] = sweepPlot(helicopter.Mesh,linspace(0,2*pi,50),real(sol.CT));
figure(1);
contourf(X,Y,L,linspace(min(min(L)),max(max(L)),1000),'edgecolor','none');
colormap("jet");
colorbar('southoutside');
daspect([1 1 1]);
title('Inflow Ratio');
xlabel('x');
ylabel('y');

