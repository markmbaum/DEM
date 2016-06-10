function [xsize,xloc,ysize,yloc] = maxSectionDimensions(S,spacing,demo)
%maxSectionDimensions finds the length of the longest row and column in a
%section that intersect
%   [xsize,xloc,ysize,yloc] = maxSectionDimensions(S,spacing,demo)
%
%   S - a section, or a structure array with grid points separated into x
%       and y fields
%   spacing - spacing between the grid points
%   demo - toggle for illustrative plotting commands

%DEMO
    if(demo~=0)
        demofrac=10;
        h=[];
    end
%DEMO

%sorting by x-coordinate
[xx,idx]=sort(S.x);
yy=S.y(idx);
ysize=0;
place=[0,0];
xloc=0;
%finding column lengths by where the x-coordinate changes
for i=1:length(xx)-1
    if(xx(i+1)~=xx(i))
        place(1)=place(2);
        place(2)=i;
        if((place(2)-place(1))>ysize)
            ysize=place(2)-place(1);
            xloc=xx(i);
        end
    end
end
            
%DEMO plotting the longest column
    if(demo~=0)
        if(xloc~=0)
            h{1}=plot(xloc,yy(xx==xloc),'b.','MarkerSize',12);
            pause(demo/demofrac);
        end
    end
%DEMO

%sorting by y-coordinate
[yy,idx]=sort(S.y);
xx=S.x(idx);
xsize=0;
place=[0,0];
yloc=0;
%finding row lengths by where the y-coordinate changes
for i=1:length(yy)-1
    if(yy(i+1)~=yy(i))
        place(1)=place(2);
        place(2)=i;
        %added condition that the selected row must intersect xloc
        if(((place(2)-place(1))>xsize) &&...
                any(abs(xx(place(1)+1:place(2))-xloc)<(.1*spacing)))
            xsize=place(2)-place(1);  
            yloc=yy(i);
        end
    end
end

            
%DEMO plotting current longest row
    if(demo~=0)
        if(yloc~=0)
            h{2}=plot(xx(yy==yloc),yloc,'r.','MarkerSize',12);
            pause(demo/demofrac);
            %plotting (xloc,yloc) with a noticeable combo of markers
            if(xloc~=0)
                h{3}=plot(xloc,yloc,'rx','MarkerSize',14);
                h{4}=plot(xloc,yloc,'b+','MarkerSize',16);
                h{5}=plot(xloc,yloc,'ko','MarkerSize',16);
                pause(demo/demofrac)
            end
        end
        pause(demo/demofrac);
%         for i=1:length(h)
%             if(~isempty(h{i}))
%                 delete(h{i});
%             end
%         end
    end
%DEMO

xsize=xsize*spacing;
ysize=ysize*spacing;
 
end