clear variables; clc;
addpath('/home/omarkahol/MATLAB/ROTOR_AERODYNAMICS/lab4/lib');

%-----------------------------------
%    STRUCTURE
%------------------------------------
R = 1.143; %m
C = 0.1905; %m
R0 = 0.1950; %m cutout 
NB = 2; %number of blades
GAMMA = 4.0; %assumed locke number for mass distribution

%------------------------------------
%    DYNAMICS
%------------------------------------
H = 670; %altitude [m]
V = 12.9857; %foreward flight velocity [m/s]
RPM = 1209.67; %rpm
ASHAFT = -1.0694; %angle of attack of the shaft [°]
THETA = 8.2391; %collective pitch angle [°]
T1S = 2.53; %sinusoidal cyclic pitch angle [°]
T1C = 0; %cosinusoidal cyclic pitch angle [°]

%-----------------------------------
%    SOLVER
%------------------------------------
SWEEPS = 6; %number of complete rotations
NT = 1000; %number of time steps per sweep
NP = 500; %number of blade panels
ITMAX = 1000; %maximum number of iterations for BLADE convergence
ITMAX_FZERO = 1000; %maximum number of iterations for nonlinear problem solution
DAMP = 0.1; %damping function for BLADE convergence
DAMP_FZERO = 0.1; %damping for nonlinear problem solution
TOL = 1e-5; %tolerance for BLADE convergence
TOL_FZERO = 1e-5; %tolerance for nonlinear problem solution

%-----------------------------------
%    BUILD THE STRUCTURES
%------------------------------------
h = buildHelicopter(R,R0,C,NB,NP,GAMMA);
fc = flightConfiguration(H,RPM,V,THETA,T1S,T1C,ASHAFT);

data = problemData(h,fc);
solver = solverConfiguration(NT,SWEEPS,ITMAX,ITMAX_FZERO,DAMP,DAMP_FZERO,TOL,TOL_FZERO);

sol = forewardFlightSolver(data,solver);



