%This script creates two plots.

%First it creates a plot of element displacements after a single
%iteration perturbation by a boundary block. The displacement is plotted
%for several different fractions of t_c, or several values of f_c.

load('ftc');

figure(1);
plot(f);
text(10,f(10,1),'$f_{tc}=0.1$','interpreter','latex','fontsize',16,'verticalalignment','top');
text(9,f(9,2),'$f_{tc}=0.2$','interpreter','latex','fontsize',16,'verticalalignment','top');
text(8,f(8,3),'$f_{tc}=0.3$','interpreter','latex','fontsize',16,'verticalalignment','top');
text(7,f(7,4),'$f_{tc}=0.4$','interpreter','latex','fontsize',16,'verticalalignment','top');
text(8,.55,'$f_{tc}=0.5$','interpreter','latex','fontsize',16,'verticalalignment','middle');
text(7,.6,'$f_{tc}=0.6$','interpreter','latex','fontsize',16,'verticalalignment','middle');
text(5.5,.65,'$f_{tc}=0.7$','interpreter','latex','fontsize',16,'verticalalignment','middle');
text(3,.65,'$f_{tc}=0.8$','interpreter','latex','fontsize',16,'horizontalalignment','center');
ylabel('Element Displacement (m)','interpreter','latex','fontsize',18);
xlabel('Iteration','interpreter','latex','fontsize',18);
title({'Effect of $f_{tc}$ on Element Responsiveness'},'interpreter','latex','fontsize',20);

%Then it plots the relationships between u_n, k_n, and t_c

load('tc_un_kn');

figure(2);

subplot(1,2,1);
plot(un_percent,kn);
title('$K_n$ vs $\delta_n$\%','interpreter','latex','fontsize',18);
xlabel('Percentage Overlap $\big(100\cdot\frac{\delta_n}{r_j}\big)$','interpreter','latex','fontsize',18);
ylabel('Stiffness Coefficient $K_n$','interpreter','latex','fontsize',18);

subplot(1,2,2);
plot(un_percent,tc);
title('$\Delta t_c$ vs $\delta_n$\%','interpreter','latex','fontsize',18);
xlabel('Percentage Overlap $\big(100\cdot\frac{\delta_n}{r_j}\big)$','interpreter','latex','fontsize',18);
ylabel('Critical Time Step $\Delta t_c$','interpreter','latex','fontsize',18);