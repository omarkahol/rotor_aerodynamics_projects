function data = postProcessor(helicopter,conf,sol)
    data = struct();
    data.AdimensionalQuantities = struct();
    data.DimensionalQuantities = struct();
    
    data.AdimensionalQuantities.CT = sum(sol.ThrustCoefficientMean);
    
    data.DimensionalQuantities.Thrust = conf.Density*...
        helicopter.Area*data.AdimensionalQuantities.CT*sol.TipVelocity^2;
end

