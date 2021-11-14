function [R_array, r_array, dR] = makeMesh(N)
%makeMesh Return an arry of discretized radial coordinates

dataFile;

dR = (R-R0)/N;
space = 0:1:N-1;

R_array = R0 + dR.*(space + 0.5);
r_array = R_array ./ R;
end

