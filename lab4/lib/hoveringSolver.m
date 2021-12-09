function solution = hoveringSolver(data)

solution = struct();
ITMAX = 1000;
TOL = 1e-12;
DAMP_l = 0.1;

solution.InflowRatio = zeros(1,data.NP);
solution.ThrustCoefficient = zeros(1,data.NP);

for i = 1:data.NP
    
    it = 0;
    err = TOL + 1;
    
    r = data.mesh(i);
    l_old = theoreticalSolution(data.solidity,2*pi,data.theta,r,0);
    
    while (it < ITMAX) && (err > TOL)
        it = it+1;
        
        %mach number computation
        M = data.Mtip * sqrt(l_old*2 + r^2);
        
        %induced angle and angle of attack
        phi = atan(l_old/r);
        alfa = data.theta - phi;
        
        %cl and cd
        [cl,cd,~] = cpcrcm(alfa,M);
        
        %thrust coeff
        ct = 0.5*data.solidity*(cos(phi)*cl -sin(phi)*cd)*(r^2 + l_old^2)*data.dr;
        
        %recompute lambda
        l_new = sqrt(ct/(4*r*l_old^2));
        l_new = l_old + DAMP_l*(l_new-l_old);
        
        %errors and update
        err = abs(l_new-l_old);
        
        l_old=l_new;
    end 
    
    %save results
    solution.InflowRatio(i) = l_new;
    solution.ThrustCoefficient(i) = ct;
end

end

