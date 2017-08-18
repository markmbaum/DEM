function [varargout]=lineBlocks(varargin)
%lineBlocks draws the boundary blocks
%   [blockhandles] = lineBlocks(nb,eb,sb,wb,bv,blockVert)
%   [blockhandles] = lineBlocks(Bstore)
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

horz=(wb(1,2)-wb(2,2))/10;

Lx=[wb(1,1)-horz; wb(:,1); wb(2,1)-horz];
Ly=[wb(1,2); wb(:,2); wb(2,2)];
Rx=[eb(1,1)+horz; eb(:,1); eb(2,1)+horz];
Ry=[eb(1,2); eb(:,2); eb(2,2)];

h(1)=line(Lx,Ly,'linewidth',2,'color',blockcolor);
h(2)=line(Rx,Ry,'linewidth',2,'color',blockcolor);
h(3)=line([nb(1,1);nb(2,1)],[nb(1,2);nb(2,2)],...
    'linewidth',1,'color',blockcolor);
h(4)=line([sb(1,1);sb(2,1)],[sb(1,2);sb(2,2)],...
    'linewidth',1,'color',blockcolor);

x=[wb(1,1)-1.5*horz;eb(1,1)+1.5*horz];
y=[mean([wb(1,2),wb(2,2)]);mean([eb(1,2),eb(2,2)])];

if(arrows)
    if(blockVert)
        u = [0;0];
        v = [1;-1];
    else
        u = [1;-1];
        v = [0;0];
    end
    h(5)=quiver(x,y,u,v,.075,'color','blue');
    if(blockVert)
        h(6) = text(x(1),y(1),['(0,',num2str(bv),')'],...
            'verticalalignment','top','horizontalalignment','center');
        h(7) = text(x(2),y(2),['(0,',num2str(bv),')'],...
            'verticalalignment','bottom','horizontalalignment','center');
    else
        h(6) = text(x(1),y(1),['(',num2str(bv),',0)'],...
            'verticalalignment','bottom','horizontalalignment','center');
        h(7) = text(x(2),y(2),['(',num2str(bv),',0)'],...
            'verticalalignment','bottom','horizontalalignment','center');
    end
end

varargout{1}=h;

end