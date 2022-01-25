function sol = forewardFlightSolver(data,solver)
    
    sol = struct();
    sol.ct = zeros(solver.Sweeps*solver.NT,data.NP);
    sol.cq = zeros(solver.Sweeps*solver.NT,data.NP);
    sol.inflow_ratio = zeros(solver.Sweeps*solver.NT,data.NP);
    sol.it_required = zeros(solver.Sweeps*solver.NT);
    sol.CT = zeros(1,solver.Sweeps*solver.NT);
    sol.CQ = zeros(1,solver.Sweeps*solver.NT);
    sol.beta = zeros(1,solver.Sweeps*solver.NT);
    sol.dbeta = zeros(1,solver.Sweeps*solver.NT);
    r = data.mesh;
    index = floor(0.97*data.NP);
    %----------------------------------------------------------------------
    %HOVERINGSOLUTION
    %----------------------------------------------------------------------
    hoveringSolution = hoveringSolver(data);
    %-----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    %MEMORY ALLOCATION
    %----------------------------------------------------------------------
    Cl=zeros(1,data.NP);
    Cd=zeros(1,data.NP);
    %derivative of pitch
    tv = -data.theta1s*sin(solver.timeMesh)+...
        -data.theta1c*cos(solver.timeMesh);
    dtv = -data.theta1s*sin(solver.timeMesh) +...
        data.theta1c*cos(solver.timeMesh);
    ddtv = data.theta1s*cos(solver.timeMesh) +...
        data.theta1c*sin(solver.timeMesh);
    
    %saved parameters
    lambda = zeros(solver.NT*solver.Sweeps,data.NP);
    beta = zeros(1,solver.NT*solver.Sweeps+1);
    dbeta = zeros(1,solver.NT*solver.Sweeps+1);
    ddbeta = zeros(1,solver.NT*solver.Sweeps+1);
    lambdap = zeros(solver.NT*solver.Sweeps, data.NP);
    
    %----------------------------------------------------------------------


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
            beta_zero = beta(n);
            fprintf('Sweep number %d out of %d is %3.2f %% complete\n',sweep,...
            solver.Sweeps,sweep_perc);
        end

        if abs(sweep_perc-50) < 1/solver.NT
            beta_180 = beta(n);
            fprintf('Sweep number %d out of %d is %3.2f %% complete\n',sweep,...
            solver.Sweeps,sweep_perc);
        end 

        psi = solver.timeMesh(n);
        
        %------------------------------------------------------------------
        % INITIALIZE OLD SOLUTION FROM HOVERING AND BETA, H FROM PREV IT
        l_old = hoveringSolution.InflowRatio;
        %------------------------------------------------------------------

        err = solver.Tol+1;
        it = 0;
        %------------------------------------------------------------------
        % START CONVERGENCE WHILE LOOP
        %------------------------------------------------------------------
        while err > solver.Tol && it < solver.ItMax
            it = it + 1;

            %--------------------------------------------------------------
            % COMPUTE THE DERIVATIVES
            %BETA FLUCTUATING
            betav = beta(n) - mean(beta(5:n));
            dbetav = d(beta-mean(beta(5:n)),n,solver.dt,false);
            ddbetav = dd(beta-mean(beta(5:n)),n,solver.dt,false);
            %LAMBDA
            dlambda = d(lambda,n,solver.dt,true);
            ddlambda = dd(lambda,n,solver.dt,true);
            %LAMBDA FLUCTUATING
            dlambdav = d(lambda-mean(lambda(5:n,:),1),n,solver.dt,true);
            ddlambdav = dd(lambda-mean(lambda(5:n,:),1),n,solver.dt,true);
            %LAMBDA P FLUCTUATING
            dlambdapv = d(lambdap-mean(lambdap(5:n,:)),n,solver.dt,true);
            %LAMBDA P MEAN
            dlambdap0 = (mean(lambdap(5:n,:),1) - ...
                mean(lambdap(5:n-1,:),1))./solver.dt;
            %--------------------------------------------------------------
            % COMPUTE THE VELOCIETIES
            %--------------------------------------------------------------
            vt = cos(beta(n))*r + data.mu*sin(psi)*cos(data.alfashaft);
            vp = l_old + data.mu*sin(data.alfashaft)*cos(beta(n)) + data.mu*...
                cos(data.alfashaft)*cos(psi)*sin(beta(n)) + ...
                data.mesh.*dbeta(n);
            v = sqrt(vt.^2 + vp.^2);
            %--------------------------------------------------------------
            % AERODYNAMIC MODEL
            %--------------------------------------------------------------
            % Compute the angle of attack
            alfa_atg = (-data.F.*dlambdapv -dlambdap0 +data.G.*ddlambdav+...
                data.G.*lambda(n,:).*(dtv(n)+betav))./vt;
            alfa = real(data.theta + data.F.*tv(n) + atan(alfa_atg));
            phi = real(atan(vp./vt));
            %Compute the mach number
            M = data.Mtip .* v;
            % Noncirculatory angle
            alfa_nc = real(0.25*data.cr.*((dtv(n)+beta(n-1))./lambda(n,:) +...
                ddlambdav./lambda(n,:)  - 0.25.*data.cr.*...
                (ddtv(n)+ddbetav)./lambda(n,:)));
            for i = 1:length(alfa)
                [cl_circ,cd_circ,~]=cpcrcm(alfa(i)-phi(i),M(i));
                [cl_nc, cd_nc,~] = cpcrcm(alfa_nc(i),M(i));
                Cl(i) = cl_circ+cl_nc;
                Cd(i) = cd_circ+cd_nc;
            end
            %--------------------------------------------------------------
            % BLADE ELEMENT MOMENTUM 
            %--------------------------------------------------------------
            % evaluation of the loads
            ct = real(0.5*data.solidity*(v.^2).*(cos(phi).*Cl-sin(phi).*Cd)...
                .*data.dr);
            cq = real(0.5*data.solidity*(v.^2).*(sin(phi).*Cl+cos(phi).*Cd)...
                .*data.dr);
            ct(isnan(ct)) = 0.0;
            CT = sum(ct(1:index));
            % DREES MODEL FOR THE VELOCITY DISTRIBUTION
            atpp = data.alfashaft + 0.5*(beta_zero-beta_180);
            kx = 1.2;
            ky = 0;
            f = 1+r.*cos(psi).*kx + r.*sin(psi).*ky;
            %RECOMPUTE THE INFLOW RATIO
            err_in = solver.TolFzero+1;
            it_in = 0;
            l_old_in = l_old;
            while err_in > solver.TolFzero && it_in < solver.ItmaxFzero
                l_new = 0.5*CT.*f./sqrt((data.mu.*cos(atpp)).^2 +...
                    (data.mu*sin(atpp)+l_old_in).^2);
                err_in = norm(l_new-l_old_in);
                it_in = it_in+1;
                l_old_in = real(l_old_in + solver.DampFzero.*(l_new -l_old_in));
            end

            ddbeta(n) = (data.LockeNumber/(4*pi))*sum((r-data.r0).*(v.^2)...
                    .*(cos(phi).*Cl-sin(phi).*Cd))*data.dr*cos(beta(n)) - ...
                    3*sin(beta(n))*sum((r-data.r0).*r*data.dr);

            %--------------------------------------------------------------
            % PREPARE NEXT ITERATION 
            %--------------------------------------------------------------
            err = norm(l_old-l_new);
            l_old = real(l_old + solver.Damp.*(l_new -l_old));
            lambdap(n,:)=vp;
            lambda(n,:)=v;
        end
        %------------------------------------------------------------------
        % DYNAMIC EQUATION
        %------------------------------------------------------------------
        if n ~= solver.NT*solver.Sweeps
            beta(n+1) = beta(n) + (1.5*dbeta(n)-0.5*dbeta(n-1))*solver.dt;
            dbeta(n+1) = dbeta(n) + (1.5*ddbeta(n)-0.5*ddbeta(n-1))*solver.dt;
        end
        %------------------------------------------------------------------
        % SAVE RESULTS
        %------------------------------------------------------------------
        sol.ct(n,:) = ct;
        sol.cq(n,:) = cq;
        sol.inflow_ratio(n,:)=l_new;
        sol.CT(n) = CT;
        sol.CQ(n) = sum(cq);
        sol.beta(n) = beta(n);
        sol.dbeta(n) = dbeta(n);
        sol.it_required(n)=it;
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
