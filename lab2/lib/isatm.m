function [p,d,T] = isatm(h)

% International Standard Atmosphere
%
% syntax: [p,d,T] = isatm(h)
%
% given altitude h
% gives pressure p, density d and temperature T
%
% NB - SI units implied
% ---------------------------------
% DIA - Politecnico di Milano, 2002

% constants
g_0 = 9.80665;
lambda = -.0065;
R = 287.05;
h_s = 11000;
h_m = 20000;

if h <= h_s
   
   % reference values for troposphere
   p_0 = 101325;
   d_0 = 1.225;
   T_0 = 288.15;
   % ISA equations for troposphere
   p = p_0 * (1+lambda*h/T_0)^(-g_0/(R*lambda));
   d = d_0 * (1+lambda*h/T_0)^(-1-g_0/(R*lambda));
   T = T_0 + lambda*h;

elseif h <= h_m
   
   % reference values for stratosphere
   p_s = 22632;
   d_s = .3639;
   T_s = 216.65;
   % ISA equations for stratosphere
   p = p_s * exp(-g_0/(R*T_s)*(h-h_s));
   d = d_s * exp(-g_0/(R*T_s)*(h-h_s));
   T = T_s;

else
   
   % diagnostics
   disp(' ')
   disp('   ***** function isatm:')
   disp('   ***** INPUT DATA ERROR !')
   disp('   ***** output data set to limit values')
   disp('   ***** results may be inaccurate')
   % ISA values for stratopause
   p = 22632 * exp(-g_0/(R*216.65)*9000);
   d = .3639 * exp(-g_0/(R*216.65)*9000);
   T = 216.65;
   
end