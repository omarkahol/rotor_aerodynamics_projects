clear variables; clc;
addpath('/home/omarkahol/MATLAB/ROTOR_AERODYNAMICS/lab6/lib');
alfa = linspace(-pi,pi,1000);
M = 0.4;
f = arrayfun(@(x) f_kirchoff(x,M),alfa);
cl = [];
cd = [];
for a = alfa
    [cl_d,cd_d,~]= cpcrcm(a,M);
    cl = [cl;cl_d];
    cd = [cd;cd_d];
end
plot(alfa,cl,'LineStyle','-','LineWidth',2);
hold on;
plot(alfa,cd,'LineStyle','-','Color','r','LineWidth',2)
hold off;

