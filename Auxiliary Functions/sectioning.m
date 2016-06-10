function [S] = sectioning(xx,yy,spacing,demo)
%   [S]= sectionGrid(xx,yy,spacing,demo)
%
%   xx - array of x-coordinates for the grid points to be sectioned
%   yy - array of y-coordinates for the grid points to be sectioned
%   spacing - spacing between the grid points
%   demo - toggle for illustrative plotting commands
%   S - a cell array of structures. Each cell entry is a structure
%       containing 'x' and 'y' fields with correspoding coordinates for the
%       points in that section.
% 
% sectioning takes grid points in the form of two vectors of x and y
% coordinates and separates them into sections. A section is a subgroup of
% points where each point is adjacent to at least one other point in the
% subgroup, so that no two points in different sections are adjacent. 
% First the function separates the entire grid into rows by sorting the y 
% vector and arranging the x vector in the same order. The rows contain only 
% the indices of the points from the initial x and y vectors, now sorted 
% accordingly. An initial point is selected (the first point in the first 
% row, but it could be any point). The initial point is added to the first 
% section, S{1}, which now only contains that point. Then the point is used 
% as a ?search point? to find adjacent points, following these steps:
% -Find rows adjacent to the row containing the search point by comparing 
% the y coordinates of the first points all of the rows
% -Compare the coordinates of each point in each adjacent row to the search 
% point
% -If any of those points are adjacent to the search point, add their 
% indices to the section and delete their index from the row containing 
% them.
% Once all adjacent points are found, each of the points that were just 
% added to the section are used as the next round of search points to find 
% other adjacent points that belong in the same section, following the above
% process. When the search for adjacent points returns nothing, a new 
% section is created and a point from one of the rows that isn?t in any 
% section yet is used as the initial point, and the process repeats. The 
% sectioning stops when each point in the entire grid has been added to a 
% section. Finally, the sections are converted back from indices to x and y
% coordinates corresponding to the oringinal grid.
% By deleting the index of each point added to a section from its respective
% row, none of the points are found more than once during the search 
% process. Each set of search points contains points that are adjacent to 
% each other. Naturally, when each search point is being used to find 
% adjacent points independently, they would find the same points at least 
% twice and maybe more. Deleting the index of a found point from its row 
% means that once it's been found, it can't be found again by other search 
% points in the same group, which speeds up the process significantly. 
% I think this process could be adapted to work on more irregular grids of 
% points, or even grids of other kinds of information instead of just 
% spatial coordinates. This version assumes adjacent points are within 
% sqrt(2) times the spacing of the grid, essentially defining adjacent 
% points as those that are within the diagonal distance between points in 
% the grid. Adjacency could be defined by another value, which would be an 
% easy change. It could also be defined by something like spectral 
% qualities, which would be more difficult.

if(~isempty(xx))

%DEMO
    if(demo~=0)
        demofrac=15;
    end
%DEMO

%separating the grid into partial rows
[yy,idx]=sort(yy,'ascend');
xx=xx(idx);
place=[0,0];
tol1=.1*spacing;
count=1;
rows=cell(0);
for j=1:length(xx)-1
    if((yy(j+1)-yy(j))>tol1)
        place(1)=place(2);
        place(2)=j;
        rows{count}=place(1)+1:place(2);
        count=count+1;
    end
end
rows{count}=place(2)+1:length(xx);

%DEMO plotting the separated rows
    if(demo~=0)
        for j=1:length(rows)
            plot(xx(rows{j}),yy(rows{j}),'.','MarkerSize',11);
        end
        pause(demo/demofrac);
    end
%DEMO

%moving points into sections
tol2=1.01*spacing;
S=cell(1,5); %allocation, even though the size of S varies a lot
section=1;
S{1}=1;
searchpts{1}=1;
rows{1}(1)=[];
if(isempty(rows{1}))
    rows=rows(2:end);
end
searchpts{2}=[];
matched=zeros(1,length(xx));
matched(1)=1;
previous=1;
place=1;
while(matched(end)==0)
    for k=1:length(searchpts{1})
        %finding nearby rows
        nbrrows=zeros(1,5);
        count=1;
        l=length(rows);
        m=1;
        while((length(nbrrows(nbrrows~=0))<3) && m<=l)
            if(abs(yy(rows{m}(1))-yy(searchpts{1}(k)))<tol2)
                nbrrows(count)=m;
                count=count+1;
            end
            m=m+1;
        end
        nbrrows=nbrrows(nbrrows~=0);
        %finding neighbor pts in those rows
        nbrpts=zeros(1,8);
        count=1;
        l=length(nbrrows);
        m=1;
        while(m<=l)%checking each row
            n=length(rows{nbrrows(m)});
            o=1;
            while(o<=n)%checking each pt in the row
                if(abs(xx(rows{nbrrows(m)}(o))-xx(searchpts{1}(k)))<tol2)
                    nbrpts(count)=rows{nbrrows(m)}(o);
                    count=count+1;
                    %removing found points so they're not refound later
                    rows{nbrrows(m)}(o)=[];
                    if(isempty(rows{nbrrows(m)}))
                        nbrrows(m)=[];
                        l=l-1;
                        m=m-1;
                    end
                    n=n-1;
                    o=o-1;
                end
                o=o+1;
            end
            m=m+1;
        end
        nbrpts=nbrpts(nbrpts~=0);
        %removing rows completely emptied when finding nbrpts
        rows=rows(~cellfun('isempty',rows));
        %adding the neighboring points to the list
        searchpts{2}=[searchpts{2},nbrpts];
    end

        %DEMO plotting new searchpoints in r and base searchpoints in g
            if(demo~=0)
                plot(xx(searchpts{1}),yy(searchpts{1}),'g.',...
                    'MarkerSize',12);
                if(~isempty(searchpts{2}))
                    plot(xx(searchpts{2}),yy(searchpts{2}),'r.',...
                    'MarkerSize',12);
                end
                pause(demo/(10*demofrac));
            end
        %DEMO

    %Adding the new pts to the current section or creating a new one if 
    %there are none
    if(~isempty(searchpts{2}))
        S{section}=[S{section},searchpts{2}];
        L=length(searchpts{2});
        matched(place+1:place+L)=searchpts{2};
        place=place+L;
        searchpts{1}=searchpts{2};
        searchpts{2}=[];
    else

        %DEMO plotting a completed section when a new one begins
            if(demo~=0)
                plot(xx(S{section}),yy(S{section}),'.','MarkerSize',12);
                pause(demo/(demofrac*10));
            end
        %DEMO

        if(rows{1}(1)==previous)
            S{section}=previous;
            rows{1}(1)=[];
            if(isempty(rows{1}))
                rows=rows(2:end);
            end
            searchpts{1}=rows{1}(1);
        else
            searchpts{1}=rows{1}(1);
            rows{1}(1)=[];
            if(isempty(rows{1}))
                rows=rows(2:end);
            end
            previous=searchpts{1};
        end
        section=section+1;
        S(section)=searchpts(1);
        matched(place+1)=searchpts{1};
        place=place+1;
    end
end

%DEMO plotting the last section after sectioning is completed
    if(demo~=0)
        plot(xx(S{section}),yy(S{section}),'.','MarkerSize',12);
        pause(demo/demofrac);
    end
%DEMO

%removing empty cells if S was overallocated
S=S(~cellfun('isempty',S));

%Separating the indices of points in each section into the x and y
%coordinates to be returned as a structure
T=S;
L=length(T);
S=cell(1,L);
for j=1:L
    S{j}=struct('x',xx(T{j}),'y',yy(T{j}));
end

else
    %Rarely, the last point in a section will be removed after an element
    %that contains it is placed. In this case xx and yy are passed in empty
    %and the sectioning process results in error. So, S is simply passed
    %back empty and will get deleted by BuildRandomRadiusArray.
    S=cell(0);
end

end