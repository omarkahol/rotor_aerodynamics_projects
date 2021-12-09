function conf = flightConfiguration(H,RPM,V,THETA,T1S,T1C,ASHAFT)
conf = struct();

conf.Altitude=H;
[p,rho,T] = isatm(H);
conf.Pressure = p;
conf.Density = rho;
conf.Temperature = T;
conf.SpeedOfSound = sqrt(287*1.4*T);
conf.RPM = RPM;
conf.Omega = RPM*2*pi/60;
conf.ForewardVelocity = V;
conf.CollectivePitch = deg2rad(THETA);
conf.PitchSin=deg2rad(T1S);
conf.PitchCosin=deg2rad(T1C);
conf.AlfaShaft = deg2rad(ASHAFT);

end

