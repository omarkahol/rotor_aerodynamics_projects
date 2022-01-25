function sol = forewardFlightSolver(data,solver)
    
    %----------------------------------------------------------------------
    % SOLUTION
    %----------------------------------------------------------------------
    sol = struct();
    sol.ct = zeros(solver.Sweeps*solver.NT,data.NP);
    sol.cq = zeros(solver.Sweeps*solver.NT,data.NP);
    sol.inflow_ratio = zeros(solver.Sweeps*solver.NT,data.NP);
    sol.it_required = zeros(1,solver.Sweeps*solver.NT);
    sol.CT = zeros(1,solver.Sweeps*solver.NT);
    sol.CQ = zeros(1,solver.Sweeps*solver.NT);
    sol.beta = zeros(1,solver.Sweeps*solver.NT);
    sol.dbeta = zeros(1,solver.Sweeps*solver.NT);
    sol.ddbeta = zeros(1,solver.Sweeps*solver.NT);
    sol.alfa_e = zeros(solver.Sweeps*solver.NT,data.NP);
    sol.err = zeros(1,solver.Sweeps*solver.NT);
    sol.v = zeros(solver.Sweeps*solver.NT,data.NP);
    sol.vp = zeros(solver.Sweeps*solver.NT,data.NP);
    sol.vt = zeros(solver.Sweeps*solver.NT,data.NP);

    %----------------------------------------------------------------------
    %COMMON INPUTS
    %----------------------------------------------------------------------
    %nondimensional time for Wagner
    s = (2*data.mu./data.cr).*solver.timeMesh;
    ds = (2*data.mu./data.cr).*solver.dt;
    %helicopter mesh
    r = data.mesh;
    %Pradtl deficiency model --> sum dct up to index
    index = floor(0.97*data.NP);
    %Wagner coefficients
    A1 = 0.165;
    A2 = 0.335;
    b1 = 0.0455;
    b2 = 0.3;
    %----------------------------------------------------------------------
    %HOVERINGSOLUTION
    %----------------------------------------------------------------------
    hoveringSolution = hoveringSolver(data);   
    %----------------------------------------------------------------------
    %MEMORY ALLOCATION
    %----------------------------------------------------------------------
    Cl=zeros(1,data.NP);
    Cd=zeros(1,data.NP);
    %pitch
    theta = data.theta-data.theta1s*sin(s)+...
        -data.theta1c*cos(s);
    %wagner deficiency functions
    Xw = zeros(solver.Sweeps*solver.NT,data.NP);
    Yw = zeros(solver.Sweeps*solver.NT,data.NP);
    %----------------------------------------------------------------------
    % START THE TIME ITERATION
    %----------------------------------------------------------------------
    beta_zero = 0;
    beta_180 = 0;
    for n = 3:solver.NT*solver.Sweeps

        %current sweep
        sweep = floor(n/solver.NT)+1;
        sweep_perc = (n/solver.NT - sweep+1)*100;

        if abs(sweep_perc) < 1/solver.NT
            beta_zero = sol.beta(n);
            fprintf('Sweep number %d out of %d is %3.2f %% complete\n',sweep,...
            solver.Sweeps,sweep_perc);
        end

        if abs(sweep_perc-50) < 1/solver.NT
            beta_180 = sol.beta(n);
            fprintf('Sweep number %d out of %d is %3.2f %% complete\n',sweep,...
            solver.Sweeps,sweep_perc);
        end 

        psi = solver.timeMesh(n);
        err = solver.Tol+1;
        it = 0;
        l_old = hoveringSolution.InflowRatio; 
        %------------------------------------------------------------------
        % START CONVERGENCE WHILE LOOP
        %------------------------------------------------------------------
        while err > solver.Tol && it < solver.ItMax
            %--------------------------------------------------------------
            % COMPUTE THE VELOCIETIES
            %--------------------------------------------------------------
            sol.vt(n,:) = r*cos(sol.beta(n)) + data.mu*sin(psi)*cos(data.alfashaft);
            sol.vp(n,:) = data.mu*sin(data.alfashaft)*cos(sol.beta(n))  + data.mu*...
                cos(data.alfashaft)*cos(psi)*sin(sol.beta(n)) + ...
                r.*sol.dbeta(n);
            sol.v(n,:) = sqrt(sol.vt(n,:).^2 + (l_old+sol.vp(n,:)).^2);
            
            %--------------------------------------------------------------
            % AERODYNAMIC MODEL
            %--------------------------------------------------------------
            % Compute the induced angle of attack
            phi = atan((l_old+sol.vp(n,:))./sol.vt(n,:));
            %equivalent angle of attack
            dvp = (sol.vp(n,:)-sol.vp(n-1))./ds;
            sol.alfa_e(n,:) = theta(n) - sol.vp(n,:)./sol.v(n,:);
            
            %step in the angle of attack times the velocity
            dalfa_e = (sol.v(n,:).*sol.alfa_e(n,:)-...
                sol.v(n-1,:).*sol.alfa_e(n-1,:));
            %wagner deficiency functions
            Xw(n,:) = exp(-b1*ds).*Xw(n-1,:) + A1.*dalfa_e...
                .*exp(-b1*0.5*ds);
            Yw(n,:) = exp(-b2*ds).*Yw(n-1,:) + A2.*dalfa_e...
                .*exp(-b2*0.5*ds);
            %wagner angle of attack
            alfa_w = sol.alfa_e(n,:) - ...
                (1./sol.v(n,:)).*(Xw(n,:) + Yw(n,:));
            %Compute the mach number
            M = data.Mtip .* sol.v(n,:);
            %Non circulatory angle of attack
            ddalfa_e = (sol.alfa_e(n)-2*sol.alfa_e(n-1)+sol.alfa_e(n-2))...
                ./(ds^2);
            alfa_nc = 0.25*(data.cr./(sol.v(n,:).^2)).*(dalfa_e./ds -...
               dvp - data.cr*0.5.*ddalfa_e);
            %Compute CL and CD 
            for i = 1:length(alfa_w)
                [cl,cd,~]=cpcrcm(alfa_w(i)+alfa_nc(i)-phi(i),M(i));
                Cl(i) = cl;
                Cd(i) = cd;
            end
            %--------------------------------------------------------------
            % BLADE ELEMENT MOMENTUM 
            %--------------------------------------------------------------
            % evaluation of the loads
            sol.ct(n,:) = real(0.5*data.solidity*(sol.v(n,:).^2)...
                .*(cos(phi).*Cl-sin(phi).*Cd));
            sol.cq(n,:) = real(0.5*data.solidity*(sol.v(n,:).^2)...
                .*(sin(phi).*Cl+cos(phi).*Cd));
            CT_old = sol.CT(n);
            sol.CT(n) = sum(sol.ct(n,1:index))*data.dr;
            err = abs(sol.CT(n)-CT_old);
            % DREES MODEL FOR THE VELOCITY DISTRIBUTION
            atpp = data.alfashaft + 0.5*(beta_zero-beta_180);
            %RECOMPUTE THE INFLOW RATIO
            err_in = solver.TolFzero+1;
            it_in = 0;
            l_old_in = l_old;
            while err_in > solver.TolFzero && it_in < solver.ItmaxFzero
                mux = real(data.mu*cos(atpp));
                muy = real(-data.mu*sin(atpp)+l_old_in);
                chi = atan2(mux,muy);
                kx = (4.0/3.0)*(1-1.8*mux^2)*tan(chi./2);
                ky = -2*mux;
                f = 1+r.*cos(psi).*kx + r.*sin(psi).*ky;
                l_new = 0.5*sol.CT(n).*f./sqrt((data.mu.*cos(atpp)).^2 +...
                    (data.mu*sin(atpp)+l_old_in).^2);
                err_in = norm(l_new-l_old_in)/data.NP;
                it_in = it_in+1;
                l_old_in = real(l_old_in + solver.DampFzero.*(l_new -l_old_in));
            end
            %--------------------------------------------------------------
            % PREPARE NEXT ITERATION 
            %--------------------------------------------------------------
            err = err + norm(l_old-l_new)/data.NP;
            l_old = real(l_old + solver.Damp.*(l_new -l_old));
            it = it + 1;
            %--------------------------------------------------------------
            % DYNAMIC EQUATION
            %--------------------------------------------------------------
%             sol.ddbeta(n) = real(data.mass*sum(r.*(sol.v(n,:).^2)...
%                     .*(cos(phi).*Cl-sin(phi).*Cd))*data.dr - ...
%                     sin(sol.beta(n))*cos(sol.beta(n)) -...
%                     1.5*data.gravity*(cos(data.alfashaft)*cos(sol.beta(n)) - ...
%                     cos(psi)*sin(data.alfashaft)*sin(sol.beta(n))));
            sol.ddbeta(n) = real(data.mass*sum((r-data.r0).*(Cl.*cos(phi)-Cd.*sin(phi)).*...
                sol.v(n,:).^2 ).*data.dr - 3*cos(sol.beta(n))*sin(sol.beta(n))*sum((r-data.r0).*r*data.dr));
            
        end 

        if n ~=solver.NT*solver.Sweeps
            sol.beta(n+1) = sol.beta(n) + (1.5*sol.dbeta(n)-...
                 0.5*sol.dbeta(n-1))*solver.dt;
            sol.dbeta(n+1) = sol.dbeta(n) + (1.5*sol.ddbeta(n)-...
                0.5*sol.ddbeta(n-1))*solver.dt;
        end
        
        %------------------------------------------------------------------
        % SAVE RESULTS
        %------------------------------------------------------------------
        
        sol.inflow_ratio(n,:)=l_new;
        sol.CQ(n) = sum(sol.cq(n,:));
        sol.it_required(n)=it;
        sol.err(n) = err;
    end
end





%--------------------------------------------------------------------------
%DERIVATOR FUNCTIONS
% Compute the derivatives at the previous iteration
%--------------------------------------------------------------------------
function [dx] = d(x,n,h,matrix)
    if matrix
        dx = (x(n,:)-x(n-2,:))./(2*h);
    else 
        dx = (x(n)-x(n-2))./(2*h);
    end
end

function [dx] = dd(x,n,h,matrix)
    if matrix
        dx = (x(n,:)-2*x(n-1,:)+x(n-2,:))./(h*h);
    else
        dx = (x(n)-2*x(n-1)+x(n-2))./(h*h);
    end
end
