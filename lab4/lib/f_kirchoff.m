function [f] = f_kirchoff(alfa,mach)
    [cl,~,~] = cpcrcm(alfa,mach);
    f = real((2*sqrt(cl/(2*pi*alfa))-1)^2);
end

