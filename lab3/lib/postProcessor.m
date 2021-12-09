function data = postProcessor(helicopter,conf,sol)
    data = struct();
    data.AdimensionalQuantities = struct();
    data.DimensionalQuantities = struct();

    max_index = floor(0.97*length(sol.ThrustCoefficientMean));
    dr = helicopter.dR/helicopter.Radius;
    
    data.AdimensionalQuantities.CT = sum(sol.ThrustCoefficientMean(1:max_index))*dr;
    
    data.DimensionalQuantities.Thrust = conf.Density*...
        helicopter.Area*data.AdimensionalQuantities.CT*sol.TipVelocity^2;
end

