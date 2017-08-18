%This script plots a demonstration of the fft control on tcfactor

load('fft_variables');

xits=10:10:160;
freq=linspace(0,.5,9);

subplot(2,2,1);
hold on;
plot(xits,tstepsbad,'bo');
plot(xits,tstepsbad,'r-');
title('Unacceptable $\Delta t$ Progression','interpreter','latex','fontsize',16);
xlabel('Iteration','interpreter','latex','fontsize',16);
ylabel('$\Delta t$','interpreter','latex','fontsize',16);

subplot(2,2,2);
hold on;
plot(xits,tstepsgood,'bo');
plot(xits,tstepsgood,'r-');
title('Acceptable $\Delta t$ Progression','interpreter','latex','fontsize',16);
xlabel('Iteration','interpreter','latex','fontsize',16);
ylabel('$\Delta t$','interpreter','latex','fontsize',16);
%set(gca,'ylim',[.8,1.8]);

subplot(2,2,3);
bar(freq,fftsigbad);
title('Unacceptable Frequencies','interpreter','latex','fontsize',16);
xlabel({'Frequency (iterations$^{-1}$)'},'interpreter','latex','fontsize',16);
ylabel('Normalized Magnitude','interpreter','latex','fontsize',16);
set(gca,'xlim',[-.05,.55]);

subplot(2,2,4);
bar(freq,fftsiggood);
title('Acceptable Frequencies','interpreter','latex','fontsize',16);
xlabel({'Frequency (iterations$^{-1}$)'},'interpreter','latex','fontsize',16);
ylabel('Normalized Magnitude','interpreter','latex','fontsize',16);
set(gca,'xlim',[-.05,.55]);