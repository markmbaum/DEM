%Plotting displacement and rotation across an array

%File name or names of the trial. For multiple trials, the data will be
%plotted together. For only one trial input the initial and final arrays
%will be plotted with the selected elements highlighted and the evolution
%of the selected element's rotations throughout the trial will be plotted
%in separate figures as well as the cross section data.
files={'friction_1_random_sh1',...
    'friction_p75_random_sh1',...
    'friction_p5_random_sh1',...
    'friction_p2_random_sh1'};

%Basically controls how many elements are included in the cross section.
%It selects with distance from the middle of the array, using a
%threshold of bandlength*mean(radii) above and below the middle.
bandlength=2.5;

%turn on/off plotting for the variables. With both on (1), it's best only
%to use one file at a time.
dispplot=1; %y-axis displacement plotting
rotplot=1;  %rotation plotting

%plot data on the same plot or in individual subplots
subplots=1;

%manually input titles or go with titles=files
%titles=files;
%titles={'Completed Trial'};
titles={'Random Element Diameters, $\mu_s=1$',...
    'Random Element Diameters, $\mu_s=0.75$',...
    'Random Element Diameters, $\mu_s=0.5$',...
    'Random Element Diameters, $\mu_s=0.2$'};

%---------------------

figure(1);clf;hold on;
cmap=colormap('lines');
markers={'o','*','+','x','s','d','^','v','p','h'};

for i=1:length(files)
    
    trial=load(files{i});

    L=length(trial.Estore);
    N=length(trial.Estore(1).x);

    %find approximate middle vertical point of initial array
    idx=FindMiddleElements(trial.Estore(1),trial.r,bandlength);
    
    %Plot initial array, completed trial, and highlight selected elements
    if(length(files)==1)
        figure(2);
        clf;hold on;
        subplot(1,2,1);
        %TraceCompletedTrial(trial.Estore,trial.r,trial.Bstore,'on','off');
        LineElements(trial.Estore(1),trial.r,'off','k','-');
        LineBlocks(trial.Bstore(1));
        HighlightElements(trial.Estore(1),trial.r,idx,0);
        title('Initial Array','interpreter','latex','fontsize',20)
        xlabel('$x$-dimension (m)','interpreter','latex','fontsize',18)
        ylabel('$y$-dimension (m)','interpreter','latex','fontsize',18);
        subplot(1,2,2);
        LineElements(trial.Estore(L),trial.r,'off','k','-');
        LineBlocks(trial.Bstore(L));
        HighlightElements(trial.Estore(L),trial.r,idx,0);
        title(titles,'interpreter','latex','fontsize',20);
        xlabel('$x$-dimension (m)','interpreter','latex','fontsize',18);
        figure(3);
        clf;
        GrabData(trial.Estore,'theta',idx,'on','simul',trial.ithorzoff,...
            trial.itverton,trial.itsstore,trial.its);
    end

    if(~isempty(idx))
        %calculate displacement/rotation and lateral position of selected elements
        l=trial.Estore(1).x(idx);
        displacement=(trial.Estore(L).y(idx)-trial.Estore(1).y(idx));
        theta=trial.Estore(L).theta(idx);
        figure(1);
        if(subplots)
            subplot(length(files),1,i);
        end
        if(dispplot && rotplot)
            fprintf('%f\n',mean(theta));
            [hAx,hLine1,hLine2]=plotyy(l,theta,l,displacement);
            %set(hAx(1),'xlim',[-1000 10000]);
            %set(hAx(2),'xlim',[-1000 10000]);
            hLine1.LineStyle='none';
            hLine1.Color=cmap(2,:);
            hLine1.Marker='*';
            hAx(1).YColor=cmap(2,:);
            %hLine1.MarkerFaceColor=cmap(2*i-1,:);
            hLine1.MarkerSize=10; 
            hLine2.LineStyle='none';
            hLine2.Color=cmap(1,:);
            hLine2.Marker='o';
            hAx(2).YColor=cmap(1,:);
            %hLine2.MarkerFaceColor=cmap(2*i,:);
            hLine2.MarkerSize=10;
            ylabel(hAx(2),'$y$-displacement (m)','interpreter','latex',...
                'fontsize',16);
            ylabel(hAx(1),'Rotation (rad)','interpreter','latex',...
                'fontsize',16);
            title(titles{i},'interpreter','latex','fontsize',18);
        elseif(dispplot && ~rotplot)
            if(subplots)
                plot(l,displacement,':o','markersize',10,'color',cmap(1,:));
                title(titles{i},'interpreter','latex','fontsize',18);
            else
                plot(l,displacement,[':',markers{i}],'markersize',10);
                legend(titles,'interpreter','latex','fontsize',14,...
                    'location','best');
                title('Element Displacement vs. Initial $x$ Position',...
                    'interpreter','latex','fontsize',18);
            end
            ylabel('$y$-displacement (m)','interpreter','latex','fontsize',16);
        elseif(~dispplot && rotplot)
            if(subplots)
                hold on;
                plot(l,theta,':*','markersize',10,'color',cmap(2,:));
                x=get(gca,'XLim');
                plot(x,[0 0],'k--');
                title(titles{i},'interpreter','latex','fontsize',18);
            else
                plot(l,theta,[':',markers{i}],'markersize',10);
                legend(titles,'interpreter','latex','fontsize',14,...
                    'location','best');
                title('Element Rotation vs. Initial $x$ Position',...
                    'interpreter','latex','fontsize',18);
            end
            ylabel('Rotation (rad)','interpreter','latex','fontsize',16);
        end
    end
    
end

xlabel('Initial $x$-position (m)','interpreter','latex','fontsize',18);

% subplot(length(files)+1,1,i+1);
% 
% data=load('LVVSZ Rotations');
% 
% DSond=data.DistFromLVVSZ_Sonder(:,1)*1e3;
% DUncSond=data.DistFromLVVSZ_Sonder(:,2)*1e3;
% RotSond=data.RotLVVSZ_Sonder(:,1)*(2*pi/360);
% RotUncSond=data.RotLVVSZ_Sonder(:,2)*(2*pi/360);
% 
% DNels=data.DistFromLVVSZ_Nelson(:,1)*1e3;
% RotNels=data.RotLVVSZ_Nelson(:,1)*(2*pi/360);
% RotUncNels=data.RotLVVSZ_Nelson(:,2)*(2*pi/360);
% 
% figure(1);hold on;
% errorbar(DSond,RotSond,RotUncSond,'r*');
% errorbar(DNels,RotNels,RotUncNels,'g*');
% 
% title('Rotation vs. Distance From LVVSZ','interpreter','latex','fontsize',18);
% ylabel('Rotation (rad)','interpreter','latex','fontsize',16);
% xlabel('Distance From LVVSZ Axis (m)','interpreter','latex','fontsize',16);
% legend({'Gale Hills Traverse','Las Vegas Range Traverse'},...
%     'interpreter','latex','fontsize',13);
% set(legend,'Location','Northwest')
% set(gca,'XLim',[-1.5e4 1.5e4]);