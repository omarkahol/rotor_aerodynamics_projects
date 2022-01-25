CPs = [];
CQs = [];
CTs = [];
V = 5:10:105;

for sol = sols
    CPs = [CPs,mean(sol.CP(end-solver.NT:end))];
    CQs =[CQs,mean(sol.CQ(end-solver.NT:end))];
    CTs = [CTs,mean(sol.CT(end-solver.NT:end))];
end

mu = V ./ (fc.Omega*h.Radius);
figure(1);
hold on;
title('Total Power');
plot(mu,CPs,'b-','LineWidth',2);
xlabel('mu');
ylabel('CP');
 [~,cd0,~] = cpcrcm(0,0);