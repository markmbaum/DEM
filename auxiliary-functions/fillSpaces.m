function [x,y,r,N] = fillSpaces(x,y,r,nrows,ncols,minRadius,demo)
%fillSpaces places elements in the spaces between other elements
%automatically adjacent to them without overlaps.
%   [x,y,r,N] = fillSpaces(x,y,r,nrows,minradius,demo)
%
%
%   x,y,r - full element location and radius arrays
%   nrows - number of rows
%   ncols - number of columns
%   minradius - minimum radius allowed for filling elements
%   demo - toggle for illustrative plotting commands

%DEMO
    if(demo~=0)
        plotElements(x,y,r,'off','k-');
        pause(demo)
    end
%DEMO

initialN=length(x);

%generating a grid of points over the elements
spacing=minRadius*.9;
row=min(x):spacing:max(x);
maxy=max(y);
[maxyfirstcol,idx]=max(y(1:nrows));
if(maxy>(maxyfirstcol+r(idx)))
    maxy=maxyfirstcol;
end
col=min(y):spacing:maxy;
[xx,yy]=meshgrid(row,col);
xx=reshape(xx,[1,length(row)*length(col)]);
yy=reshape(yy,[1,length(row)*length(col)]);
fprintf('Fill grid size = %d points\n',length(xx));

%DEMO plotting full grid
    if(demo~=0)
        plot(xx,yy,'k.')
        pause(demo);
    end
%DEMO

%removing grid points contained by elements
[xx,yy]=removeContained(xx,yy,x,y,r);

%DEMO plotting grid without contained points
    if(demo~=0)
        plotElements(x,y,r,'off','k-');
        plot(xx,yy,'k.')
        pause(demo);
    end
%DEMO

%splitting the grid into sections
S=Sectioning(xx,yy,spacing,0);

%DEMO plotting sections in color
    if(demo~=0)
        for j=1:length(S)
            plot(S{j}.x,S{j}.y,'.','markersize',12);
        end
        pause(demo);
        %getting rid of grid
        plotElements(x,y,r,'off','k-');
    end
%DEMO
    
%placing elements
miny=min(y);
maxx=max(x);
minx=min(x);
i=length(x);
xlocprev=-1;
ylocprev=-1;
[xsize(1),xloc(1),ysize(1),yloc(1)]=maxSectionDimensions(S{1},spacing,0);
while(~isempty(S))
    %checking repeat (xloc,yloc)
    if((xloc(1)==xlocprev) && (yloc(1)==ylocprev))
        if(length(S)>1)
            S=S(2:end);
            [xsize(1),xloc(1),ysize(1),yloc(1)]=...
                maxSectionDimensions(S{1},spacing,0);
        else
            S=[];
        end
    else
        %make new element in the space
        [D,I]=closestElements(xloc(1),yloc(1),x,y,r);
        i=i+1;
        [x(i),y(i),r(i)]=...
            autoNewElement(x,y,r,xloc(1),yloc(1),I,D,[mean(D),mean(D)],0);
        [x(i),y(i),r(i)]=AutoAdjust(i,x,y,r,I,0);
        %check if the new element is adequately sized and if it is within the
        %confines of the original array. Deleting it and (xloc,yloc).
        if((r(i)<minRadius)||(y(i)>maxy)||(y(i)<miny)||(x(i)>maxx)||(x(i)<minx))
            %deleting element
            r(i)=[];
            x(i)=[];
            y(i)=[];
            i=i-1;
            %deleting section
            if(length(S)>1)
                S=S(2:end);
                [xsize(1),xloc(1),ysize(1),yloc(1)]=...
                    maxSectionDimensions(S{1},spacing,0);
            else
                S=[];
            end
        else
            xlocprev=xloc(1);
            ylocprev=yloc(1);
            %remove contained points
            [S{1}.x,S{1}.y]=removeContained(S{1}.x,S{1}.y,x(i),y(i),r(i));

            %DEMO
                if(demo)
                    lineElements(x(i),y(i),r(i),'off','k','-');
                    pause(demo);
                end
            %DEMO

            %resection
            resection=Sectioning(S{1}.x,S{1}.y,spacing,0);
            %add new sections to the list
            if(~isempty(resection) && length(S)>1)
                S=[resection,S(2:end)];
            elseif(~isempty(resection) && length(S)==1)
                S=resection;
            elseif(isempty(resection) && length(S)>1)
                S=S(2:end);
            elseif(isempty(resection) && length(S)==1)
                S=[];
            end
            %get section info for next round
            if(~isempty(S))
                [xsize(1),xloc(1),ysize(1),yloc(1)]=...
                    maxSectionDimensions(S{1},spacing,0);
            end
        end
    end
end

N=length(x);

%DEMO
    if(demo~=0)
        clf;
        lineElements(x,y,r,'off','k','-');
        if((N-initialN)>0)
            title('Initial array, fill elements highlighted',...
                'interpreter','latex','fontsize',18);
            highlightElements(x,y,r,initialN+1:N,3);
        end
    end
%DEMO

%potentially interesting stats, exit message
if((N-initialN)>0)
fprintf('N=%d, %d fill elements added, N-nrows*ncols=%+d, avgr=%g\n',...
    N,N-initialN,N-nrows*ncols,sum(r)/N);
end

end