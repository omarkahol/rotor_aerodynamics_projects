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
    for n = 4:solver.NT*solver.Sweeps

        clc;
        fprintf('Computing iteration %d out of %d\n',n,...
            solver.NT*solver.Sweeps);

        psi = solver.timeMesh(n);
        
        %------------------------------------------------------------------
        % INITIALIZE OLD SOLUTION FROM HOVERING AND BETA, H FROM PREV IT
        l_old = hoveringSolution.InflowRatio;
        lambdap(n,:) = lambdap(n-1,:);
        lambda(n,:) = lambdap(n-1,:);
        beta(n) = beta(n-1);
        dbeta(n)=dbeta(n-1);
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
            betav = beta(n) - mean(beta(3:n));
            dbetav = d(beta-mean(beta(3:n)),n,solver.dt,false);
            ddbetav = dd(beta-mean(beta(3:n)),n,solver.dt,false);
            %LAMBDA
            dlambda = d(lambda,n,solver.dt,true);
            ddlambda = dd(lambda,n,solver.dt,true);
            %LAMBDA FLUCTUATING
            dlambdav = d(lambda-mean(lambda(3:n,:),1),n,solver.dt,true);
            ddlambdav = dd(lambda-mean(lambda(3:n,:),1),n,solver.dt,true);
            %LAMBDA P FLUCTUATING
            dlambdapv = d(lambdap-mean(lambdap(3:n,:)),n,solver.dt,true);
            %LAMBDA P MEAN
            dlambdap0 = (mean(lambdap(3:n,:),1) - ...
                mean(lambdap(3:n-1,:),1))./solver.dt;
            %--------------------------------------------------------------
            % COMPUTE THE VELOCIETIES
            %--------------------------------------------------------------
            vt = r + data.mu*sin(psi)*cos(data.alfashaft);
            vp = l_old + data.mu*sin(data.alfashaft) + data.mu*...
                cos(data.alfashaft)*cos(psi)*sin(beta(n)) + ...
                data.mesh.*dbeta(n);
            v = sqrt(vt.^2 + vp.^2);
            %--------------------------------------------------------------
            % AERODYNAMIC MODEL
            %--------------------------------------------------------------
            % Compute the angle of attack
            alfa_atg = (-data.F.*dlambdapv -dlambdap0 +data.G.*ddlambdav+...
                data.G.*lambda(n,:).*(dtv(n)+betav))./vt;
            alfa = data.theta + data.F.*tv(n) + atan(alfa_atg);
            phi = atan(vp./v);
            %Compute the mach number
            M = data.Mtip .* v;
            %Equivalent lift and drag
            i = 1;
            for a = alfa-phi
                [cl_sect,cd_sect,~]=cpcrcm(a,M(i));
                Cl(i) = cl_sect;
                Cd(i) = cd_sect;
                i=i+1;
            end
            %Correction for the non circulatory part
            Cl = Cl + 0.5*pi*data.cr.*((dtv(n)+beta(n-1))./lambda(n,:) +...
                ddlambdav./lambda(n,:)  - 0.25.*data.cr.*...
                (ddtv(n)+ddbetav)./lambda(n,:));
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
            mux = data.mu*cos(data.alfashaft);
            muy = data.mu*sin(data.alfashaft) + l_old;
            chi = atan(mux./muy);
            kx = (4/3).*(1-cos(chi)-1.8*data.mu^2)./(sin(chi));
            ky = -2*data.mu;
            f = 1+r.*cos(psi).*kx + r.*sin(psi).*ky;
            %RECOMPUTE THE INFLOW RATIO
            err_in = solver.TolFzero+1;
            it_in = 0;
            l_old_in = l_old;
            while err_in > solver.TolFzero && it_in < solver.ItmaxFzero
                l_new = 0.5*CT.*f./sqrt((data.mu.*cos(data.alfashaft)).^2 +...
                    (data.mu*sin(data.alfashaft)+l_old_in).^2);
                err_in = norm(l_new-l_old_in);
                it_in = it_in+1;
                l_old_in = real(l_old_in + solver.DampFzero.*(l_new -l_old_in));
            end

            %--------------------------------------------------------------
            % PREPARE NEXT ITERATION 
            %--------------------------------------------------------------
            err = norm(l_old-l_new);
            l_old = real(l_old + solver.Damp.*(l_new -l_old));
            lambdap(n,:)=vp;
            lambda(n,:)=v;
            %------------------------------------------------------------------
            % DYNAMIC EQUATION
            %------------------------------------------------------------------
            ddbeta(n-1) = 3*sum((1./data.mass).*(data.mesh(1:index)-...
                data.r0).*ct(1:index)) - ...
                3*sum((data.mesh(1:index)-data.r0)...
                .*data.mesh(1:index).*sin(beta(n)));
            beta(n) = beta(n-1) + ((23/12)*dbeta(n-1)-(16/12)*dbeta(n-2)+...
                (5/12)*dbeta(n-3))*solver.dt;
            dbeta(n) = dbeta(n-1) + ((23/12)*ddbeta(n-1)-(16/12)*ddbeta(n-2)+...
                (5/12)*ddbeta(n-3))*solver.dt;
            beta(n) =0.0;
            dbeta(n)=0.0;
        end
        %------------------------------------------------------------------
        % SAVE RESULTS
        %------------------------------------------------------------------
        sol.ct(n,:) = ct;
        sol.cq(n,:) = cq;
        sol.inflow_ratio(n,:)=l_new;
        sol.CT(n) = CT;
        sol.CQ(n) = sum(cq);
        sol.beta(n) = rad2deg(beta(n));
        sol.dbeta(n) = rad2deg(dbeta(n));
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
