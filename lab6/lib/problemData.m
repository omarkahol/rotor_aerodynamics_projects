function [data] = problemData(helicopter,configuration)

data = struct();
data.NP = length(helicopter.Mesh);
data.dr = helicopter.dR ./ helicopter.Radius;
data.cr = helicopter.Chord/helicopter.Radius;
data.r0 = helicopter.CutOut/helicopter.Radius;
data.solidity = helicopter.Solidity;
data.mesh = helicopter.Mesh ./ helicopter.Radius; 
data.mass = helicopter.LockeNumber/(4*pi);
data.gravity = 9.81/(helicopter.Radius*configuration.Omega^2);
data.Mtip = configuration.Omega*helicopter.Radius/configuration.SpeedOfSound;
data.mu = configuration.ForewardVelocity/(configuration.Omega*helicopter.Radius);

k = configuration.Omega*helicopter.Chord/(2*configuration.ForewardVelocity);
J0 = besselj(0,k);
J1 = besselj(1,k);
Y0 = bessely(0,k);
Y1 = bessely(1,k);

den = (J1+Y1)^2 + (Y1-J0)^2;
numf = J1*(J1+Y0)+Y1*(Y1-J0);
numg = Y1*Y0+J1*J0;

data.F = numf/den;
data.G = numg/den;

data.theta = configuration.CollectivePitch;
data.theta1s = configuration.PitchSin;
data.theta1c= configuration.PitchCosin;
data.alfashaft= configuration.AlfaShaft;
data.LockeNumber = helicopter.LockeNumber;

end

