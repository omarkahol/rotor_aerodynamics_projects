function [X,Y,L,l] = sweepPlot(r,psi_mesh,x)
    X = [];
    Y = [];
    for psi = psi_mesh
        X = [X,r'.*cos(psi)];
        Y = [Y,r'.*sin(psi)];
    end
    L = x';
    l = x(1,:);
end

