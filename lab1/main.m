clear variables; clc;

dataFile;

N = 100;
[R_array, r_array, dR] = makeMesh(N);

T_ext = 20 + 273.15; %outside temperature
sound_speed = sqrt(1.4*287*T_ext);

RPM = 1250;
Omega = RPM*2*pi/60;
Vtip = Omega*R;

theta = deg2rad(8);


ITMAX = 1000;
TOL = 1e-12;
DAMP = 0.5;

lambda_array = [];
ct_array = [];
cq_array = [];

for r = r_array
    
    it = 0;
    err = TOL + 1;
    
    l_old = theoreticalSolution(sigma,2*pi,theta,r);
    l_new = l_old;
    
    while (it < ITMAX) && (err > TOL)
        phi = atan(l_old/r);
        alpha = theta - phi;
        M = (Vtip/sound_speed) * sqrt(r^2 + l_old^2);
        [cl,cd,~] = cpcrcm(alpha, M);
        
        ct = 0.5*sigma*(r^2 + l_old^2)*(cos(phi)*cl - sin(phi)*cd);
        l_new = sqrt(ct/(4*r));
        
        l_new = l_old + DAMP*(l_new-l_old);
        
        err = abs(l_new-l_old);
        l_old = l_new;
        it = it+1;
    end
    
    lambda_array = [lambda_array, l_new];
    cq_array = [cq_array, 0.5*sigma*(l_new^2 + r^2)*(sin(phi)*cl + cos(phi)*cd)*r];
    ct_array = [ct_array, ct];   
    
end

CT = sum(ct_array .* (dR/R));
CQ = sum(cq_array .* (dR/R));

FM = CT^(1.5) / (CQ*sqrt(2));

