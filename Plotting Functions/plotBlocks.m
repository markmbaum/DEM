function [varargout]=plotBlocks(varargin)
%plotBlocks draws the boundary blocks
%   [blockhandles] = plotBlocks(nb,eb,sb,wb,bv,blockVert)
%   [blockhandles] = plotBlocks(Bstore)
%
%   nb,eb,sb,wb - north block, east block, etc.
%   bv - block velocity
%   blockVert - 1 for blocks moving vertically, 0 for horizontally
%   Bstore - array of structures containing block positions   

if(length(varargin) == 1)
    nb = varargin{1}.nb;
    eb = varargin{1}.eb;
    sb = varargin{1}.sb;
    wb = varargin{1}.wb;
    arrows = 0;
else
    nb = varargin{1};
    eb = varargin{2};
    sb = varargin{3};
    wb = varargin{4};
    bv = varargin{5};
    arrows = 1;
    blockVert = varargin{6};
end

blockcolor = 'b';

horz=abs(wb(1,2)-wb(2,2))/10;

Lx=[wb(1,1)-horz; wb(:,1); wb(2,1)-horz];
Ly=[wb(1,2); wb(:,2); wb(2,2)];
Rx=[eb(1,1)+horz; eb(:,1); eb(2,1)+horz];
Ry=[eb(1,2); eb(:,2); eb(2,2)];

hold on;

%plot the blocks
h=plot(Lx,Ly,blockcolor,'LineWidth',2);
h=[h;plot(Rx,Ry,blockcolor,'LineWidth',2)];
h=[h;plot([nb(1,1);nb(2,1)],[nb(1,2);nb(2,2)],...
    blockcolor,'LineWidth',1)];
h=[h;plot([sb(1,1);sb(2,1)],[sb(1,2);sb(2,2)],...
    blockcolor,'LineWidth',1)];

%plot the velocity arrows
x=[wb(1,1)-1.5*horz;eb(1,1)+1.5*horz];
y=[mean([wb(1,2),wb(2,2)]);mean([eb(1,2),eb(2,2)])];
if(blockVert)
    u = [0;0];
    v = [1;-1];
else
    u = [1;-1];
    v = [0;0];
end
h=[h;quiver(x,y,u,v,.075,'color','blue')];

if(blockVert)
    h = [h; text(x(1),y(1),['(0,',num2str(bv),')'],...
        'verticalalignment','top','horizontalalignment','center')];
    h = [h; text(x(2),y(2),['(0,',num2str(bv),')'],...
        'verticalalignment','bottom','horizontalalignment','center')];
else
    h = [h; text(x(1),y(1),['(',num2str(bv),',0)'],...
        'verticalalignment','bottom','horizontalalignment','center')];
    h = [h; text(x(2),y(2),['(',num2str(bv),',0)'],...
        'verticalalignment','bottom','horizontalalignment','center')];
end

varargout{1}=h;

end