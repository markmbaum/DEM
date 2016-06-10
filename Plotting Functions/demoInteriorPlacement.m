function demoInteriorPlacement(xloc,yloc,x,y,r,i,c,w,Iclose,demo)
%demoInteriorPlacement draws lines to the anchor elements used to generate
%a new interior element.
%   [h] = demoInteriorPlacement(xloc,yloc,x,y,r,i,c,w,Iclose,demo)

h = cell(1,6 + length(c) + length(Iclose));

%crosshair
h{1} = line(xloc,yloc,'color','k','marker','x','markerSize',14);
h{2} = line(xloc,yloc,'color','k','marker','+','markerSize',16);
h{3} = line(xloc,yloc,'color','k','marker','o','markerSize',16);

%lines to position averaged elements
h{4} = line([xloc,x(i)],[yloc,y(i)],...
    'color','k','linewidth',4*w(end));
for j = 1:length(c)
    h{4+j} = line([xloc,x(c(j))],[yloc,y(c(j))],...
        'color','k','linewidth',4*w(j));
end
%the position averaged elements themselves
h{4+length(c)+1} = lineElements(x(c),y(c),r(c),'off','k','-');

if(~isempty(Iclose))
    pause(demo);
    %lines to "close" elements
    for j = 1:length(Iclose)
        h{5+length(c)+j} = line([xloc,x(Iclose(j))],[yloc,y(Iclose(j))],...
            'color',[0.63 0.07 0.18],'linestyle','-','linewidth',max(w)*4);
    end
    %the close elements themselves
    hElements = lineElements(x(Iclose),y(Iclose),r(Iclose),'off',...
        [0.63 0.078 0.18],'-');
end

pause(demo);

%draw element
h{end} = lineElements(x(i+1),y(i+1),r(i+1),'off','k','-');

pause(demo);

%delete the graphics
if(~isempty(Iclose))
    delete(hElements);
end
for j = 1:length(h)-1
    delete(h{j});
end

lineElements(x(i),y(i),r(i),'off',[.75 .75 .75],'-');

end