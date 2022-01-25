function [X,Y,L] = sweepPlot(data,solver,x)
    r = data.mesh;
    psi_data = solver.timeMesh(end-solver.NT-1:end);
    X = [];
    Y = [];
    for psi = psi_data
        X = [X,r'*cos(psi)];
        Y = [Y,r'*sin(psi)];
    end
    L = x(end-solver.NT-1:end,:)';
end

