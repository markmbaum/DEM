%plots demonstration of the effect of different fvc on a single element

figure(1);
clf;

subplot(3,1,1);
load('fvc_p05');
frictionlimit=max(temp);
plot(temp,'b.');
hold on;
plot([0;length(temp)],[frictionlimit;frictionlimit],'r--');
title({'Block Shear Force, $f_{vc}=0.05$'},'interpreter','latex','fontsize',18);
%xlabel('Iterations','interpreter','latex','fontsize',18);
%ylabel('Shear Force (N)','interpreter','latex','fontsize',18);

subplot(3,1,2);
load('fvc_p25');
frictionlimit=max(temp);
plot(temp,'b.');
hold on;
plot([0;length(temp)],[frictionlimit;frictionlimit],'r--');
title({'Block Shear Force, $f_{vc}=0.25$'},'interpreter','latex','fontsize',18);
%xlabel('Iterations','interpreter','latex','fontsize',18);
ylabel('Shear Force (N)','interpreter','latex','fontsize',18);

subplot(3,1,3);
load('fvc_1');
plot(temp,'b.');
hold on;
plot([0;length(temp)],[frictionlimit;frictionlimit],'r--');
title({'Block Shear Force, $f_{vc}=1$'},'interpreter','latex','fontsize',18);
xlabel('Iteration','interpreter','latex','fontsize',18);
%ylabel('Shear Force (N)','interpreter','latex','fontsize',18);