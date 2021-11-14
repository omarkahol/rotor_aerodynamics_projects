function conf = flightConfiguration(h,theta,RPM,Vup)

conf = struct();

conf.Altitude=h;
[p,rho,T] = isatm(h);
conf.Pressure = p;
conf.Density = rho;
conf.Temperature = T;
conf.SpeedOfSound = sqrt(287*1.4*T);
conf.RPM = RPM;
conf.Omega = RPM*2*pi/60;
conf.AxialVelocity = Vup;
conf.AngleOfAttack = deg2rad(theta);

end

