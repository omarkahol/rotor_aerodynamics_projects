function [solver] = solverConfiguration(NT,SWEEPS,ITMAX,ITMAX_FZERO,DAMP,DAMP_FZERO,TOL,TOL_FZERO)
    solver = struct();
    solver.NT = NT;
    solver.Sweeps= SWEEPS;
    solver.timeMesh = linspace(0,2*pi*SWEEPS,NT*SWEEPS);
    solver.dt = solver.timeMesh(2)-solver.timeMesh(1);
    
    solver.ItMax = ITMAX;
    solver.ItmaxFzero = ITMAX_FZERO; %maximum number of iterations for nonlinear problem solution
    solver.Damp=DAMP; %damping function for BLADE convergence
    solver.DampFzero = DAMP_FZERO; %damping for nonlinear problem solution
    solver.Tol = TOL; %tolerance for BLADE convergence
    solver.TolFzero = TOL_FZERO; %tolerance for nonlinear problem solution
end

