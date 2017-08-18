function replayTrial(varargin)
%replayTrial plots the movement of elements and blocks over a
%trial.
%   replayTrial(E,r,B,frames,framerate,highlight)
%   replayTrial(filename,frames,framerate,highlight)
%
%   filename - name of trial file
%   E - element info array of structures
%   r - radius array
%   B - block info array of structures
%   frames - the number of different snapshots plotted. Set frames to 'all'
                %to plot all iterations.
%   framerate - set the framerate for previously saved movie replay
%   highlight - toggle whether middle elements are highlighted or not

%SPECIFY FOLDER THAT MOVIE FILES ARE SAVED TO
path='/Volumes/BAUMPRIMARY/Thesis/DEM-v2-Parallel/Trial Movies';

switch nargin
    case 1
        framerate=5;
        frames=25;
        loaded=1;
        blocks=false;
        highlight=false;
    case 2
        if(ischar(varargin{1}))
            frames=varargin{2};
            framerate=5;
            loaded=1;
        else
            Estore=varargin{1};
            r=varargin{2};
            Bstore=[];
            frames=25;
            framerate=5;
            loaded=0;
        end
        blocks=true;
        highlight=false;
    case 3
        if(ischar(varargin{1}))
            frames=varargin{2};
            framerate=varargin{3};
            loaded=1;
        else
            Estore=varargin{1};
            r=varargin{2};
            Bstore=varargin{3};
            frames=25;
            framerate=5;
            loaded=0;
        end
        blocks=true;
        highlight=false;
    case 4
        if(ischar(varargin{1}))
            frames=varargin{2};
            framerate=varargin{3};
            highlight=varargin{4};
            loaded=1;
        else
            Estore=varargin{1};
            r=varargin{2};
            Bstore=varargin{3};
            frames=varargin{4};
            framerate=5;
            highlight=false;
            loaded=0;
        end
        blocks=true;
    case 5
        Estore=varargin{1};
        r=varargin{2};
        Bstore=varargin{3};
        blocks=true;
        frames=varargin{4};
        framerate=varargin{5};
        highlight=false;
        loaded=0;
    case 6
        Estore=varargin{1};
        r=varargin{2};
        Bstore=varargin{3};
        blocks=true;
        frames=varargin{4};
        framerate=varargin{5};
        highlight=varargin{6};
        loaded=0;
    otherwise
        error('Inputs incorrect or unsupported');
end

if(loaded)
    load(varargin{1});
end

L=length(Estore);
if(strcmp(frames,'all') || frames>L)
    frames=L;
end

madealready=false;
if(loaded && exist(['Movie_',varargin{1},'.mat'],'file'))
    moviefile=load(['Movie_',varargin{1}]);
    if((moviefile.frames>=frames) && (moviefile.highlight==highlight))
        madealready=true;
    end
end

fdim=get(gcf,'OuterPosition');
scrsz=get(groot,'ScreenSize');
if(madealready && scrsz(3)>=moviefile.FigureDim(3) &&...
        scrsz(4)>=moviefile.FigureDim(4))
    set(gcf,'OuterPosition',moviefile.FigureDim);
    title('Trial Replay','interpreter','latex','fontsize',20);
    movie(moviefile.Mov,1,framerate);
    set(gcf,'OuterPosition',fdim);
    idx=FindMiddleElements(Estore(1),r,1.5);
else
    N=length(Estore(1).x);
    
    frameinterval=round(L/frames);
    frameits=[1,frameinterval:frameinterval:L];

    %generate circles and anglebars
    circleres=100; %number of steps taken in drawing the circles
    sincircle=sin(0:(2*pi/circleres):2*pi)';
    coscircle=cos(0:(2*pi/circleres):2*pi)';
    elementx=zeros(length(sincircle),N);
    elementy=elementx;
    for i=1:N
        elementx(:,i)=r(i)*sincircle;
        elementy(:,i)=r(i)*coscircle;
    end
    movie_Ex=zeros(length(sincircle),N,length(frameits));
    movie_Ey=movie_Ex;
    movie_barx=zeros(2,N,length(frameits));
    movie_bary=movie_barx;
    movie_shrblkx=zeros(4,2,length(frameits));
    movie_shrblky=movie_shrblkx;
    movie_otrblkx=zeros(2,2,length(frameits));
    movie_otrblky=movie_otrblkx;
    for i=1:length(frameits)
        %elements and anglebars
        for j=1:N
            movie_Ex(:,j,i)=elementx(:,j)+Estore(frameits(i)).x(j);
            movie_Ey(:,j,i)=elementy(:,j)+Estore(frameits(i)).y(j);
            movie_barx(:,j,i)=[Estore(frameits(i)).x(j);...
                Estore(frameits(i)).x(j)+r(j)*sin(Estore(frameits(i)).theta(j))];
            movie_bary(:,j,i)=[Estore(frameits(i)).y(j);...
                Estore(frameits(i)).y(j)+r(j)*cos(Estore(frameits(i)).theta(j))];
        end
        %modified shear block coordinates
        horz=(Bstore(frameits(i)).wb(1,2)-Bstore(frameits(i)).wb(2,2))/10;
        movie_shrblkx(:,1,i)=[Bstore(frameits(i)).wb(1,1)-horz;...
            Bstore(frameits(i)).wb(:,1); Bstore(frameits(i)).wb(2,1)-horz];
        movie_shrblkx(:,2,i)=[Bstore(frameits(i)).eb(1,1)+horz;...
            Bstore(frameits(i)).eb(:,1); Bstore(frameits(i)).eb(2,1)+horz];
        movie_shrblky(:,1,i)=[Bstore(frameits(i)).wb(1,2);...
            Bstore(frameits(i)).wb(:,2); Bstore(frameits(i)).wb(2,2)];
        movie_shrblky(:,2,i)=[Bstore(frameits(i)).eb(1,2);...
            Bstore(frameits(i)).eb(:,2); Bstore(frameits(i)).eb(2,2)];
        %top and bottom block coordinates
        movie_otrblkx(:,1,i)=[Bstore(frameits(i)).nb(1,1);...
            Bstore(frameits(i)).nb(2,1)];
        movie_otrblkx(:,2,i)=[Bstore(frameits(i)).sb(1,1);...
            Bstore(frameits(i)).sb(2,1)];
        movie_otrblky(:,1,i)=[Bstore(frameits(i)).nb(1,2);...
            Bstore(frameits(i)).nb(2,2)];
        movie_otrblky(:,2,i)=[Bstore(frameits(i)).sb(1,2);...
            Bstore(frameits(i)).sb(2,2)];
    end

    clf;
    axis equal;
    title('Trial Replay','interpreter','latex','fontsize',20);
    xlabel('$x$-dimension (m)','interpreter','latex','fontsize',20);
    ylabel('$y$-dimension (m)','interpreter','latex','fontsize',20);
    
    idx=FindMiddleElements(Estore(1),r,1.5);
    
    he=line(movie_Ex(:,:,1),movie_Ey(:,:,1),'color','k');
    hbar=line(movie_barx(:,:,1),movie_bary(:,:,1),'color','b');
    if(blocks)
        hb1=line(movie_shrblkx(:,:,1),movie_shrblky(:,:,1),'color','b','linewidth',2);
        hb2=line(movie_otrblkx(:,:,1),movie_otrblky(:,:,1),'color','b');
    end
    if(highlight)
        hhl=line(movie_Ex(:,idx,1),movie_Ey(:,idx,1),'color','r','linewidth',1.5);
    end
    drawnow;
    Mov(length(frameits)+1)=struct('cdata',[],'colormap',[]);
    Mov(1)=getframe;
    
    for i=2:length(frameits)
        for j=1:N
            set(he(j),'XData',movie_Ex(:,j,i));
            set(he(j),'YData',movie_Ey(:,j,i));
            set(hbar(j),'XData',movie_barx(:,j,i));
            set(hbar(j),'YData',movie_bary(:,j,i));
        end
        if(blocks)
            for j=1:2
                set(hb1(j),'XData',movie_shrblkx(:,j,i));
                set(hb1(j),'YData',movie_shrblky(:,j,i));
                set(hb2(j),'XData',movie_otrblkx(:,j,i));
                set(hb2(j),'YData',movie_otrblky(:,j,i));
            end
        end
        if(highlight)
            for j=1:length(idx)
                set(hhl(j),'XData',movie_Ex(:,idx(j),i));
                set(hhl(j),'YData',movie_Ey(:,idx(j),i));
            end
        end
        drawnow;
        Mov(i)=getframe;
    end

    fprintf('\n');

    Mov(i+1)=getframe;
    FigureDim=get(gcf,'OuterPosition');
    %save([path,'/Movie_',varargin{1}],...
    %    'Mov','frames','FigureDim','highlight');
end

% temp=zeros(1,length(idx));
% Etemp=struct('x',temp,'y',temp,'theta',temp);
% Etemp=repmat(Etemp,length(Estore),1);
% for i=1:length(Estore)
%     Etemp(i).x=Estore(i).x(idx);
%     Etemp(i).y=Estore(i).y(idx);
%     Etemp(i).theta=Estore(i).theta(idx);
% end

clf;
% if(isempty(Bstore) || ~blocks)
%     TraceCompletedTrial(Estore,r,Bstore,'off','on')
% elseif(length(Estore(1).x)>100)
%     TraceCompletedTrial(Estore,r,Bstore,'on','off');
% else
%     TraceCompletedTrial(Estore,r,Bstore,'on','on');
% end
lineElements(Estore(L),r,'off','k','-');
lineBlocks(Bstore(L));
if(highlight)
    highlightElements(Estore(L),r,idx,0);
end
title('Completed Trial','fontsize',20,'interpreter','latex');
xlabel('$x$-dimension (m)','interpreter','latex','fontsize',20);
ylabel('$y$-dimension (m)','interpreter','latex','fontsize',20);

end