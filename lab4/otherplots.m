figure(1);
hold on;
title('CT');
plot(solver.timeMesh,sol.CT,'k-','LineWidth',2);
xlabel('t [-]');
ylabel('CT [-]');

figure(2);
hold on;
title('flap')
plot(solver.timeMesh,rad2deg(sol.beta),'k-','LineWidth',2);
xlabel('t [-]');
ylabel('beta [-]');

% figure(4);
% hold on;
% title('CP');
% plot(solver.timeMesh,sol.CP,'k-','LineWidth',2);
% xlabel('t [-]');
% ylabel('CP [-]');

