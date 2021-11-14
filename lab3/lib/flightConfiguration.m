function conf = flightConfiguration(h,theta,RPM,V,theta1s,theta1c,beta,alfa_sh)

conf = struct();

conf.Altitude=h;
[p,rho,T] = isatm(h);
conf.Pressure = p;
conf.Density = rho;
conf.Temperature = T;
conf.SpeedOfSound = sqrt(287*1.4*T);
conf.RPM = RPM;
conf.Omega = RPM*2*pi/60;
conf.ForewardVelocity = V;
conf.AngleOfAttack = deg2rad(theta);
conf.Theta1s=deg2rad(theta1s);
conf.Theta1c=deg2rad(theta1c);
conf.Flap = deg2rad(beta);
conf.AlfaShaft = deg2rad(alfa_sh);

end

