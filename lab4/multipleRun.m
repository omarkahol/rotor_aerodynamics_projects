clear variables; clc;
addpath('/home/omarkahol/MATLAB/ROTOR_AERODYNAMICS/lab4/lib');

%-----------------------------------
%    STRUCTURE
%------------------------------------
R = 1.143; %m
C = 0.1905; %m
R0 = 0.1950; %m cutout 
NB = 2; %number of blades
GAMMA = 4; %assumed locke number for mass distribution

%------------------------------------
%    DYNAMICS
%------------------------------------
H = 0; %altitude [m]
V = 12.9; %foreward flight velocity [m/s]
RPM = 1209.67; %rpm
ASHAFT = 1.0694; %angle of attack of the shaft [°]
THETA = 8.2391; %collective pitch angle [°]
T1S = 2.53; %sinusoidal cyclic pitch angle [°]
T1C = 0; %cosinusoidal cyclic pitch angle [°]

%-----------------------------------
%    SOLVER
%------------------------------------
SWEEPS = 5; %number of complete rotations
NT = 500; %number of time steps per sweep
NP = 250; %number of blade panels
ITMAX = 1000; %maximum number of iterations for BLADE convergence
ITMAX_FZERO = 1000; %maximum number of iterations for nonlinear problem solution
DAMP = 0.1; %damping function for BLADE convergence
DAMP_FZERO = 0.1; %damping for nonlinear problem solution
TOL = 1e-4; %tolerance for BLADE convergence
TOL_FZERO = 1e-4; %tolerance for nonlinear problem solution

%-----------------------------------
%    BUILD THE STRUCTURES
%------------------------------------
h=buildHelicopter(R,R0,C,NB,NP,GAMMA);
solver = solverConfiguration(NT,SWEEPS,ITMAX,ITMAX_FZERO,DAMP,DAMP_FZERO,TOL,TOL_FZERO);

sols = [];
for V = 5:10:105
    fprintf('RUNNING FOR V = %d\n',V);
    fc = flightConfiguration(H,RPM,V,THETA,T1S,T1C,ASHAFT);
    data = problemData(h,fc);
    sol = forewardFlightSolver(data,solver);
    sols = [sols,sol];
end
