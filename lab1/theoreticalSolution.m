function [l] = theoreticalSolution(sigma,Cla,theta,r)

l = (sigma*Cla/16)*(-1 + sqrt(1+32*(theta*r/(sigma*Cla))));
end

