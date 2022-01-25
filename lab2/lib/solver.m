function solution = solver(helicopter,conf)

solution = struct();
ITMAX = 10000;
TOL = 1e-12;
DAMP_l = 0.1;
DAMP_a = 0.1;

N = length(helicopter.Mesh());
solution.InflowRatio = zeros(1,N);
solution.TheoreticalInflowRatio = zeros(1,N);
solution.SwirlEffect = zeros(1,N);
solution.ThrustCoefficient = zeros(1,N);
solution.TorqueCoefficient = zeros(1,N);
solution.IterationsRequired = zeros(1,N);
solution.InducedAngle = zeros(1,N);
solution.TipVelocity = conf.Omega*helicopter.Radius;
solution.AxialInflowRatio = conf.AxialVelocity/(conf.Omega*helicopter.Radius);

for i = 1:N
    
    it = 0;
    err = TOL + 1;
    
    r = helicopter.Mesh(i)/helicopter.Radius;
    l_old = theoreticalSolution(helicopter.Solidity,2*pi,conf.AngleOfAttack,r,solution.AxialInflowRatio);
    solution.TheoreticalInflowRatio(i) = l_old;
    a_old = 0;
    
    while (it < ITMAX) && (err > TOL)
        it = it+1;
        
        %correction of tangential velocity due to swirl effect
        v = r*(1-a_old);
        
        %mach number computation
        M = (solution.TipVelocity/conf.SpeedOfSound) ...
            * sqrt(l_old*2 + v^2);
        
        %induced angle and angle of attack
        phi = atan(l_old/v);
        alfa = conf.AngleOfAttack - phi;
        
        %cl and cd
        [cl,cd,~] = cpcrcm(alfa,M);
        
        %prandtl factors
        f = 0.5*helicopter.NumberOfBlades*(1-r)/(phi*r);
        F = (2/pi)*acos(exp(-f));
        
        %thrust coeff
        ct = real(0.5*helicopter.Solidity*(cos(phi)*cl -sin(phi)*cd)*(v^2 + l_old^2));
        
        %recompute lambda
        l_new = real((solution.AxialInflowRatio +...
            sqrt(solution.AxialInflowRatio^2 + ct/(r*F)))*0.5);
        
        l_new = l_old + DAMP_l*(l_new-l_old);
        
        %compute cq
        cq = 0.5*helicopter.Solidity*(sin(phi)*cl + cos(phi)*cd)*(v^2 + l_new^2)*r;
        
        %recompute a
        a_new = cq/(4*l_new*r^3);
        a_new = a_old + DAMP_a*(a_new-a_old);
        
        %errors and update
        err = sqrt((l_new-l_old)^2 + (a_new-a_old)^2);
        
        a_old=a_new;
        l_old=l_new;
    end 
    
    %save results
    solution.InflowRatio(i) = l_new;
    solution.SwirlEffect(i) = a_new;
    solution.ThrustCoefficient(i) = ct;
    solution.TorqueCoefficient(i) = cq;
    solution.IterationsRequired(i) = it-1;
    solution.InducedAngle(i) = rad2deg(phi);
end

end

