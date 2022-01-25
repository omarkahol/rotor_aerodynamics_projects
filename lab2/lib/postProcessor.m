function data = postProcessor(helicopter,conf,sol)
    data = struct();
    data.AdimensionalQuantities = struct();
    data.DimensionalQuantities = struct();
    
    data.AdimensionalQuantities.CT = sum(sol.ThrustCoefficient)*helicopter.dR/helicopter.Radius;
    data.AdimensionalQuantities.CQ = sum(sol.TorqueCoefficient)*helicopter.dR/helicopter.Radius;
    data.AdimensionalQuantities.FM = (1/sqrt(2))*...
        (data.AdimensionalQuantities.CT^1.5)/data.AdimensionalQuantities.CQ;
    
    
    data.DimensionalQuantities.Thrust = conf.Density*...
        helicopter.Area*data.AdimensionalQuantities.CT*sol.TipVelocity^2;
    data.DimensionalQuantities.PowerRequired = conf.Density*...
        helicopter.Area*data.AdimensionalQuantities.CQ*sol.TipVelocity^3;
    data.DimensionalQuantities.Torque = conf.Density*...
        helicopter.Area*data.AdimensionalQuantities.CQ*sol.TipVelocity^2*helicopter.Radius;
end

