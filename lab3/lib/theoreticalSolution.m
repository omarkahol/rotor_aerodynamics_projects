function [l] = theoreticalSolution(sigma,Cla,theta,r,l_c)

f1 = (1/16)*sigma*Cla - l_c/2;
f2 = 0.125*sigma*Cla*theta*r;
l = sqrt(f1^2 + f2)-f1;
end

