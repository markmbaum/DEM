%Approximation validation script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % making contour plot of percentage approximation error

figure(1);clf;

mult=1; %using to make sure it doesn't change for different values

rl=mult*(1:.01:25);
rj=mult;

R=(rl*rj)./(rl+rj); %not used

minrl=min(rl);
U=minrl/1000:minrl/1000:minrl/10;

exact=zeros(length(rl),length(U));
appx=zeros(length(rl),length(U));

for i=1:length(rl)
    for j=1:length(U)
        appx(i,j)=a_appx(rl(i),rj,U(j));
        exact(i,j)=a_exact(rl(i),rj,U(j));
    end
end

percentagediff=100*((appx-exact)./exact);

figure(1);
contour(100*U./rj,rl./rj,percentagediff,'showtext','on');

%set(gca,'YScale','log');
title('Percentage Approximation Error $\Big(\frac{\sqrt{2 R \delta_n}-a}{a}\cdot 100\Big)$',...
    'Interpreter','Latex','fontsize',22);
xlabel('Percentage Overlap \Big($\frac{\delta_n}{r_j}\cdot 100$\Big)',...
    'Interpreter','Latex','fontsize',20);
ylabel('Ratio of Element Radii \Big($\frac{r_l}{r_j}$\Big)',...
    'Interpreter','Latex','fontsize',20);

axis fill

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotting a demonstration of why it works

figure(2);
x1=0:.001:1.25;
y1=sqrt(1-(x1-1).^2);
x2=0:.001:1;
y2=zeros(1,length(x2));
for i=1:length(x2);
    y2(i)=a_appx(1,1,2*x2(i));
end

plot(x1,y1,x2,y2);

legend({'$f_1=\sqrt{1 - \left(x - 1\right)^2} \approx a$','$f_2=\sqrt{2x}\approx\sqrt{\delta_n}$'},...
    'location','northwest','interpreter','latex','fontsize',16,'box','off');
title({'Shape of the Approximation Curve Compared to a Circle'},...
    'interpreter','latex','fontsize',20);
xlabel('$x$','interpreter','latex','fontsize',20);
ylabel('$f(x)$','interpreter','latex','fontsize',20);
axis equal
