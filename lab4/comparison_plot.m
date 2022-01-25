figure(1); hold on; grid on;
title('CT')
xlabel('t [-]')
ylabel('CT [-]')
plot(solver.timeMesh,sol_gl.CT,'k-','LineWidth',2);
plot(solver.timeMesh,sol_bw.CT,'r-','LineWidth',2);
plot(solver.timeMesh,sol_drees.CT,'g-','LineWidth',2);
plot(solver.timeMesh,sol_hw.CT,'b-','LineWidth',2);
legend('Glauert','Blake and White', 'Drees', 'Howlett');

figure(2); hold on; grid on;
title('beta')
xlabel('t [-]')
ylabel('beta [°]')
plot(solver.timeMesh,rad2deg(sol_gl.beta),'k-','LineWidth',2);
plot(solver.timeMesh,rad2deg(sol_bw.beta),'r-','LineWidth',2);
plot(solver.timeMesh,rad2deg(sol_drees.beta),'g-','LineWidth',2);
plot(solver.timeMesh,rad2deg(sol_hw.beta),'b-','LineWidth',2);
legend('Glauert','Blake and White', 'Drees', 'Howlett');
