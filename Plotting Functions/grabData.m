function [varargout] = grabData(varargin)
%grabData extracts the values of a single variable (x,y,theta,disp.) from
%each iteration of a trial and returns it in a matrix where colums are made
%up of values for a single element over all iterations. Varargout is a set
%of arrays for each input in variable, in the same order. If
%variables=='all' three vectors are returned-x,y,theta in that order. The
%function requires at least the first three inputs. If plotting is on and
%multiple variables are called for, theta is plotted.
%   [varargout] = grabData(Estore,variables,elements,plotting,...
                        %plottype,ithorzoff,itverton,storeInt,totalits)
%   [varargout] = grabData(filename,variables,elements,plotting,plottype)
%
%   E - element array of structures
%   filename - .mat file with saved variables from a completedtrial
%   variable - cell of strings with combo of 'x','y', and 'theta', or 'all'
%   elements - vector specifying which elements to grab data from or 'all'
%   plotting - toggle whether the grabbed data is automatically plotted
%               with 'on' or 'off'
%   plottype - 'subplot' or 'simul'
%   ithorzoff - iteration where the horizontal block velocities zero
%   itverton - iteration where the vertical block velocities turn on
%   itsstore - vector of iterations marking which iterations are saved in
%               Estore and Bstore
%   totalits - total number of iterations of the trial

switch nargin
    case 3
        if(ischar(varargin{1}))
            load(varargin{1});
            variables=varargin{2};
            elements=varargin{3};
            plotting='on';
            plottype='simul';
        else
            Estore=varargin{1};
            variables=varargin{2};
            elements=varargin{3};
            plotting='on';
            plottype='simul';
            ithorzoff=0;
            itverton=0;
            itsstore=0;
        end
    case 5
        if(ischar(varargin{1}))
            load(varargin{1});
            variables=varargin{2};
            elements=varargin{3};
            plotting=varargin{4};
            plottype=varargin{5};
        else
            Estore=varargin{1};
            variables=varargin{2};
            elements=varargin{3};
            plotting=varargin{4};
            plottype=varargin{5};
            ithorzoff=0;
            itverton=0;
            itsstore=0;
        end
    case 8
        Estore=varargin{1};
        variables=varargin{2};
        elements=varargin{3};
        plotting=varargin{4};
        plottype=varargin{5};
        itverton=varargin{6};
        ithorzoff=itverton;
        itsstore=(0:length(Estore)-1)*varargin{7};
        its=varargin{8};
    otherwise
        error('Number of inputs must be 3, 5, or 9');
end

hold on;

%managing certain variables
if(ischar(variables))
    variables={variables};
end

L=length(Estore);

if(strcmp(elements,'all'))
    elements=1:length(Estore(1).x);
end

N=length(elements);

%compiling data from Estore
plotnum=1;
varargout=cell(1,length(variables));
for q=1:length(variables)
    if(strcmp(variables{q},'x'))
        A=zeros(L,N);
        for i=1:N
            for j=1:L
                A(j,i)=Estore(j).x(elements(i));
            end
        end
    elseif(strcmp(variables{q},'y'))
        A=zeros(L,N);
        for i=1:N
            for j=1:L
                A(j,i)=Estore(j).y(elements(i));
            end
        end
    elseif(strcmp(variables{q},'theta'));
        plotnum=q;
        A=zeros(L,N);
        for i=1:N
            for j=1:L
                A(j,i)=Estore(j).theta(elements(i));
            end
        end
    else
        error('Invalid variable input. Must be ''x'', ''y'', or ''theta''');
    end
    varargout{q}=A;
end

%plotting compiled data
if(strcmp(plotting,'on'))
    if(strcmp(plottype,'simul'))
        clf;plot(repmat(itsstore',1,size(varargout{plotnum},2)),...
            varargout{plotnum});
        for i=1:N
            text(its-mod(its,itsstore(end)),...
                varargout{plotnum}(end,i),num2str(elements(i)),...
                'verticalalignment','bottom','horizontalalignment',...
                'right');
        end
        if(strcmp(variables{plotnum},'theta'))
            str='Rotation Angles';
            if(ithorzoff~=0)
                ithorzoffvirtual=ceil((ithorzoff/its)*length(Estore));
                L=[max(Estore(ithorzoffvirtual).theta(elements));...
                    min(Estore(ithorzoffvirtual).theta(elements))];
                M=mean(abs(L))*.25;
                L=L+[M;-M];
                if(L(1)<max(Estore(end).theta(elements)))
                    L(1)=mean([L(1),max(Estore(end).theta(elements))]);
                end
                if(L(2)>min(Estore(end).theta(elements)))
                    L(2)=mean([L(2),min(Estore(end).theta(elements))]);
                end
                line(ones(2,1)*ithorzoff,L,'color','k');
                text(ithorzoff,L(1),'v_{comp} off','fontsize',12);
            end
            if(itverton~=0)
                itvertonvirtual=ceil((itverton/its)*length(Estore));
                if(itverton==ithorzoff)
                    text(itverton,L(2),'v_{shear} on','fontsize',12,...
                        'verticalalignment','top');
                elseif(ithorzoff~=0 && itverton~=ithorzoff)
                    line(ones(2,1)*itverton,L,'color','k')
                    text(itverton,L(2),'v_{shear} on','fontsize',12,...
                        'verticalalignment','top');
                elseif(ithorzoff==0)
                    L=[max(Estore(itvertonvirtual).theta(elements));...
                        min(Estore(itvertonvirtual).theta(elements))];
                    M=mean(abs(L))*.25;
                    L=L+[M;-M];
                    if(L(1)<max(Estore(end).theta(elements)))
                    L(1)=mean([L(1),max(Estore(end).theta(elements))]);
                    end
                    if(L(2)>min(Estore(end).theta(elements)))
                        L(2)=mean([L(2),min(Estore(end).theta(elements))]);
                    end
                    line(ones(2,1)*itverton,L,'color','k')
                    text(itverton,L(2),'v_{shear} on','fontsize',12,...
                        'verticalalignment','top');
                end
            end
        elseif(strcmp(variables{plotnum},'x'))
            str='X-position';
        elseif(strcmp(variables{plotnum},'y'))
            str='Y-position';
        end
        str=['Element ',str,' vs. Iteration Number'];
        title(str,'fontsize',20,'interpreter','latex');
        if(strcmp(variables{plotnum},'theta'))
            ylabel('$\theta$ (rad)','fontsize',18,'interpreter','latex');
        else
            ylabel(variables{plotnum},'fontsize',18,'interpreter','latex');
        end
        xlabel('Iteration','fontsize',18,'interpreter','latex');
    elseif(strcmp(plottype,'subplot'))
        clf;
        for i=1:N
            subplot(N,1,i);
            plot(itsstore,varargout{plotnum}(:,i));
            title(['Element ',num2str(elements(i))],...
                'fontsize',20,'interpreter','latex');
            if(strcmp(variables{plotnum},'theta'))
                ylabel('$\theta$ (rad)','fontsize',18,'interpreter','latex');
            else
                ylabel(variables{plotnum},'fontsize',18,'interpreter','latex');
            end
        end
    end
    drawnow;
end 

end