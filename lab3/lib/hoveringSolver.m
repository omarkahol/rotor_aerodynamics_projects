function solution = hoveringSolver(helicopter,conf)

solution = struct();
ITMAX = 1000;
TOL = 1e-12;
DAMP_l = 0.5;

N = length(helicopter.Mesh());
solution.InflowRatio = zeros(1,N);
solution.ThrustCoefficient = zeros(1,N);
solution.IterationsRequired = zeros(1,N);
solution.InducedAngle = zeros(1,N);
solution.TipVelocity = conf.Omega*helicopter.Radius;

for i = 1:N
    
    it = 0;
    err = TOL + 1;
    
    r = helicopter.Mesh(i)/helicopter.Radius;
    l_old = theoreticalSolution(helicopter.Solidity,2*pi,conf.AngleOfAttack,r,0);
    
    while (it < ITMAX) && (err > TOL)
        it = it+1;
        
        %mach number computation
        M = (solution.TipVelocity/conf.SpeedOfSound) ...
            * sqrt(l_old*2 + r^2);
        
        %induced angle and angle of attack
        phi = atan(l_old/r);
        alfa = conf.AngleOfAttack - phi;
        
        %cl and cd
        [cl,cd,~] = cpcrcm(alfa,M);
        
        %thrust coeff
        ct = 0.5*helicopter.Solidity*(cos(phi)*cl -sin(phi)*cd)*(r^2 + l_old^2);
        
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
    solution.IterationsRequired(i) = it-1;
    solution.InducedAngle(i) = rad2deg(phi);
end

end

