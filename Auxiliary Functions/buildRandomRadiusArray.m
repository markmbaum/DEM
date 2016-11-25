function[x,y,r]=buildRandomRadiusArray(nrows,ncols,randBounds,demo)
%buildRandomRadiusArray generates the locations and radii of element
%centers when the radii are randomly generated within bounds.
%   [x,y,r,N] = buildRandomRadiusArray(nrows,ncols,randBounds,demo)
%
%   nrows - number of rows
%   ncols - number of columns
%   randBounds - new element size constraints
%   demo - toggle for illustrative plotting commands

%------------------------------------------------------nrw-------------------
%Tuning parameters for the interior element generation process

minRadius = min(randBounds); maxRadius = max(randBounds);
%delete new elements if r(i) < minRadius or r(i) > maxRadius. When the
%radius bounds are close to each other, like RandBounds=[1,2], this should
%be lowered a little to prevent the array from being partially filled. If
%it's equal to min(RandBounds), spaces in the interior usually develop that
%can't be filled. If minRadius were closer to min(RandBounds)/2 the smaller
%elements would be allowed in those spaces.

maxReach = max(randBounds);%compared against distance to closest elements
%to determine which elements the new element is anchored to. Should be at
%least max(RandBounds) or the likelihood of elements failing to adjust
%increases a lot. New elements might have a radius large enough to
%overlap with elements that aren't selected for I by the comparison
%with maxReach.

%-------------------------------------------------------------------------

%DEMO turning pause and hold on, clearing figure
    if(demo ~= 0)
        fprintf('Assembling array with randombly distributed radii...\n');
        pause on;
        hold on;
        clf;
        mov = repmat(struct('cdata',[],'colormap',[]), 1000, 1);
        mov_count = 1;
    end
%DEMO

if(ncols == 1 && nrows == 1)
    x = 0;
    y = 0;
    r = boundedRandom(randBounds);
    N = 1;
elseif(xor(ncols==1,nrows==1))
    error('2x1 or 1x2');
else

%constructing the leftmost column
x = zeros(3*nrows*ncols,1);
y = zeros(3*nrows*ncols,1);
r = zeros(3*nrows*ncols,1);
r(1) = boundedRandom(randBounds);
%DEMO
    if(demo~=0)
        lineElements(x(1),y(1),r(1),'off','k','-');
        mov(mov_count) = getframe;
        mov_count = mov_count + 1;
        pause(demo);
    end
%DEMO
for i = 2:nrows
    r(i) = boundedRandom(randBounds);
    x(i) = r(i) - r(1);
    y(i) = y(i-1) - sqrt((r(i-1) + r(i))^2 - (r(i) - r(i-1))^2);
    [x(i),y(i),r(i)] = autoEdgeAdjust(i,x,y,r,'l',demo);
    %DEMO
        if(demo~=0)
            lineElements(x(i),y(i),r(i),'off','k','-');
            mov(mov_count) = getframe;
            mov_count = mov_count + 1;
            pause(demo);
        end
    %DEMO
    %first element cannot be trapped by the second
    if(i==2)
        if((y(2) + r(2))>y(1))
            r(2) = r(2)/2;
            x(2) = r(2) - r(1);
            y(2) = y(1) - sqrt((r(1) + r(2))^2 - (r(2) - r(1))^2);
        end
    end
end

%checking last two elements in the first column. If the center of the last
%element is above the lowest extent of the second to last column, the
%bottom row can't be properly executed.
while(y(i)>(y(i-1) - r(i-1)))
    r(i) = boundedRandom(randBounds);
    x(i) = r(i) - r(1);
    y(i) = y(i-1) - sqrt((r(i-1) + r(i))^2 - (r(i) - r(i-1))^2);
    [x(i),y(i),r(i)] = autoEdgeAdjust(i,x,y,r,'l',demo);
    %DEMO
        if(demo~=0)
            lineElements(x(i),y(i),r(i),'off','k','-');
            mov(mov_count) = getframe;
            mov_count = mov_count + 1;
            pause(demo);
        end
    %DEMO
end

y = y - y(i);
x = x - x(i);

%DEMO
    if(demo~=0)
        plotElements(x(1:i),y(1:i),r(1:i),'off','k-');
        mov(mov_count) = getframe;
        mov_count = mov_count + 1;
    end
%DEMO

corner1 = i;

countError = 4*ncols*nrows;

%constructing the bottom row
for i = i+1:i+ncols-1
    r(i) = boundedRandom(randBounds);
    y(i) = y(i-1) - r(i-1) + r(i);
    x(i) = x(i-1) + sqrt((r(i-1) + r(i))^2 - (r(i) - r(i-1))^2);
    [x(i),y(i),r(i)] = autoEdgeAdjust(i,x,y,r,'b',demo);
    %DEMO
        if(demo~=0)
            lineElements(x(i),y(i),r(i),'off','k','-');
            mov(mov_count) = getframe;
            mov_count = mov_count + 1;
            pause(demo);
        end
    %DEMO
    if(i>countError)
        plotElements(x(1:i),y(1:i),r(1:i),'on','k-');
        error('Number of edge elements grew unusually large, bottom row');
    end
end

%ensuring the bottom row is wider than the width of the first column
%otherwise building the next column will result in error
[maxx,idx] = max(x(1:nrows));
arraywidth = maxx + r(idx);
additions = 0;
while((x(i) + r(i))<arraywidth)
    i = i + 1;
    r(i) = boundedRandom(randBounds);
    y(i) = y(i-1) - r(i-1) + r(i);
    x(i) = x(i-1) + sqrt((r(i-1) + r(i))^2 - (r(i) - r(i-1))^2);
    [x(i),y(i),r(i)] = autoEdgeAdjust(i,x,y,r,'b',demo);
    %DEMO
        if(demo~=0)
            lineElements(x(i),y(i),r(i),'off','k','-');
            mov(mov_count) = getframe;
            mov_count = mov_count + 1;
            pause(demo);
        end
    %DEMO
    additions = additions + 1;
    if(i>countError)
        plotElements(x(1:i),y(1:i),r(1:i),'on','k-');
        error('Number of edge elements grew unusually large, bottom row');
    end
end
if(additions>0)
    fprintf('%d additions were made to the bottom row\n',additions);
end

%also ensuring the last element in the bottom row can meet the right column
%the previous element can't extend farther in the x direction than it
while(x(i-1) + r(i-1) > x(i) + r(i))
    r(i) = r(i)*2;
    y(i) = y(i-1) - r(i-1) + r(i);
    x(i) = x(i-1) + sqrt((r(i-1) + r(i))^2 - (r(i) - r(i-1))^2);
end

%DEMO
    if(demo~=0)
        plotElements(x(1:i),y(1:i),r(1:i),'off','k-');
        mov(mov_count) = getframe;
        mov_count = mov_count + 1;
    end
%DEMO

corner2 = i;

arraywidth = x(i) - x(nrows) + r(i) + r(nrows);

%constructing the rightmost column
%placing the top element first so that it's even with the top element in
%the left column
i = i + 1;
r(i) = boundedRandom(randBounds);
x(i) = x(i-1) + r(i-1) - r(i);
y(i) = y(1) + r(1) - r(i);
[x(i),y(i),r(i)] = autoEdgeAdjust(i,x,y,r,'c3',demo);
%DEMO
    if(demo~=0)
        lineElements(x(i),y(i),r(i),'off','k','-');
        mov(mov_count) = getframe;
        mov_count = mov_count + 1;
        pause(demo);
    end
%DEMO
corner3 = i;
%building the column top down
I = [];
D = sqrt((x(corner2) - x(i))^2 + (y(corner2) - y(i))^2) - r(corner2) - r(i);
while(isempty(I) && (D > (.1*mean(r))))
    i = i + 1;
    r(i) = boundedRandom(randBounds);
    x(i) = x(i-1) + r(i-1) - r(i);
    y(i) = y(i-1) - sqrt((r(i-1) + r(i))^2 - (r(i) - r(i-1))^2);
    if(i == corner3 + 1)
        while((y(i) + r(i)) > y(i-1))
            r(i) = r(i)/2;
            x(i) = x(i-1) + r(i-1) - r(i);
            y(i) = y(i-1) - sqrt((r(i-1) + r(i))^2 - (r(i) - r(i-1))^2);
        end
    end
    [x(i),y(i),r(i)] = autoEdgeAdjust(i,x,y,r,'r',demo);
    %DEMO
        if(demo~=0)
            lineElements(x(i),y(i),r(i),'off','k','-');
            mov(mov_count) = getframe;
            mov_count = mov_count + 1;
            pause(demo);
        end
    %DEMO
    if(r(i)>arraywidth)
        plotElements(x(1:i),y(1:i),r(1:i),'on','k-');
        error('Unstable element position forced in right column');
    end
    if(i>countError)
        plotElements(x(1:i),y(1:i),r(1:i),'on','k-');
        error('Number of edge elements grew unusually large, right col');
    end
    idx = [corner2,i];
    I = adjacency(2,x(idx),y(idx),r(idx));
    D = sqrt((x(corner2) - x(i))^2 + (y(corner2) - y(i))^2) - r(corner2) - r(i);
end
%rearranging right column to correctly sequence the edge elements
x(corner3:i) = myFlip(x(corner3:i));%Have to use custom function 'myFlip'
y(corner3:i) = myFlip(y(corner3:i));%instead of built-in 'flip' because the
r(corner3:i) = myFlip(r(corner3:i));%latter is not supported in earlier
corner3 = i;                        %version of Matlab.

%DEMO
    if(demo~=0)
        plotElements(x(1:i),y(1:i),r(1:i),'off','k-');
        mov(mov_count) = getframe;
        mov_count = mov_count + 1;
    end
%DEMO

%creating the top row
done = 0;
while(~done)
    i = i + 1;
    r(i) = boundedRandom(randBounds);
    y(i) = y(i-1) + r(i-1) - r(i);
    x(i) = x(i-1) - sqrt((r(i-1) + r(i))^2 - (r(i) - r(i-1))^2);
    [x(i),y(i),r(i)] = autoEdgeAdjust(i,x,y,r,'t',demo);
    %DEMO
        if(demo~=0)
            lineElements(x(i),y(i),r(i),'off','k','-');
            mov(mov_count) = getframe;
            mov_count = mov_count + 1;
            pause(demo);
        end
    %DEMO
    if(r(i)>arraywidth)
        plotElements(x(1:i),y(1:i),r(1:i),'on','k-');
        error('Unstable element position was forced in top row');
    end
    if(i>countError)
        plotElements(x(1:i),y(1:i),r(1:i),'on','k-');
        error('Number of edge elements grew unusually large, top row');
    end
    %complete when new elements are running up against element 1, forcing
    %them to be small
    temp = distance(x(i),y(i),x(1),y(1)) - r(i) - r(1);
    if((r(i) < minRadius) && (temp < minRadius*.1))
        done = 1;
        x(i) = 0; y(i) = 0; r(i) = 0;
        i = i - 1;
    end
end

last = i;

%indices of elements around the edges and combinations of them
if(i > 15)
    temp = 1:round(i/15):i;
else
    temp = 1:round(i/5):i;
end
%go heavy on edge elements near the corners
temp = [4,3,2,1,last,last-1,last-2,...
        corner1+(-3:3),corner2+(-3:3),corner3+(-3:3),temp];
%can't include out of bounds indices
temp = temp(temp <= last);
temp = temp(temp > 0);
%only need unique elements selected into 'temp'
c1 = unique(temp);

%%%%%
% c1 = [1,round(corner1/2.),...
%     corner1,round(mean(corner1,corner2)),...
%     corner2,round(mean(corner2,corner3)),...
%     corner3,round(mean(corner3,i))];
% c1 = unique(c1);
% c1 = 1:i;
%%%%%

%DEMO
    if(demo~=0)
        plotElements(x(1:i),y(1:i),r(1:i),'off','k-');
        highlightElements(x,y,r,c1,2);
        mov(mov_count) = getframe;
        mov_count = mov_count + 1;
        clf;
        lineElements(x(1:i),y(1:i),r(1:i),'off',[.75 .75 .75],'-');
        mov(mov_count) = getframe;
        mov_count = mov_count + 1;
        pause(demo);
    end
%DEMO

%get indices of elements in the edge rows and colums, without overlaps
lcol = c1(c1(c1 < corner1) > 1);
lcolL = length(lcol);
brow = c1(c1(c1 < corner2) > corner1);
browL = length(brow);
rcol = c1(c1(c1 < corner3) > corner2);
rcolL = length(rcol);
trow = c1(c1 > corner3);
trowL = length(trow);

%assemble all pairs of those elements
c2 = zeros(lcolL*(browL+rcolL+trowL) + browL*(rcolL+trowL) + rcolL*trowL,2);
count = 1;
for j = lcol
    for k = [brow,rcol,trow];
        c2(count,1) = j;
        c2(count,2) = k;
        count = count + 1;
    end
end
for j = brow
    for k = [rcol,trow]
        c2(count,1) = j;
        c2(count,2) = k;
        count = count + 1;
    end
end
for j = rcol
    for k = trow
        c2(count,1) = j;
        c2(count,2) = k;
        count = count + 1;
    end
end

%shuffle c1 and c2
c1 = c1(randperm(length(c1)));
c2 = c2(randperm(length(c2)),:);

%weights for averaging the positions of elements in c1
c1w = zeros(9,2);
count = 1;
for j = [.01,.05,.1,.25,.5,.75,.9,.95,.99]
    c1w(count,1) = j;
    c1w(count,2) = 1 - j;
    count = count + 1;
end

%weights for averaging the positions of elements in c2
c2w = zeros(6,3);
count = 1;
for j = 0.2:0.2:0.8
    for k = 0.2:0.2:(0.8 - j + 0.01)
        c2w(count,1) = j;
        c2w(count,2) = k;
        c2w(count,3) = 1 - j - k;
        count = count + 1;
    end
end

c1L = length(c1);
c2L = length(c2);

placed = true;
while(placed)
    placed = false;
    %crawl through weighted sums of 2 edge elements with the ith element
    for j = 1:c2L
        for k = 1:6
            %take weighted average of various element positions
            xloc = c2w(k,1)*x(c2(j,1)) + c2w(k,2)*x(c2(j,2)) + c2w(k,3)*x(i);
            yloc = c2w(k,1)*y(c2(j,1)) + c2w(k,2)*y(c2(j,2)) + c2w(k,3)*y(i);
            %check if (xloc,yloc) is inside any existing elements
            if(~isContained(xloc,yloc,x,y,r))
                [D,Iclose] = closestElements(xloc,yloc,x(1:i),y(1:i),r(1:i));
                if(any(D > minRadius))%avoid repeatedly placing tiny elements
                    %pare down D and I to elements within an acceptable reach
                    temp = D < maxReach;
                    D = D(temp);
                    Iclose = Iclose(temp);
                    %place a new element
                    i = i + 1;
                    x(i) = xloc;
                    y(i) = yloc;
                    r(i) = boundedRandom(randBounds);
                    [x(i),y(i),r(i),fail] = autoAdjust(...
                        i,x(1:i),y(1:i),r(1:i),D,Iclose);
                    if(~fail && (r(i) <= maxRadius) && (r(i) >= minRadius))
                        %successful placement
                        %DEMO
                            if(demo ~= 0)
                                [mov, mov_count] = demoInteriorPlacement(xloc,yloc,x,y,r,i-1,...
                                    c2(j,:),c2w(k,:),Iclose,demo,mov,mov_count);
                                demo = 0.975*demo;
                            end
                        %DEMO
                        placed = true;
                        break;
                    else %failed element placement, delete it
                        x(i) = 0;
                        y(i) = 0;
                        r(i) = 0;
                        i = i - 1;
                    end
                end
            end
        end
    end
    if(~placed)
        %crawl through weighted sums of 1 edge element and the ith element
        for j = 1:c1L
            for k = 1:9
                %take weighted average of various element positions
                xloc = c1w(k,1)*x(c1(j)) + c1w(k,2)*x(i);
                yloc = c1w(k,1)*y(c1(j)) + c1w(k,2)*y(i);
                %check if (xloc,yloc) is inside any existing elements
                if(~isContained(xloc,yloc,x,y,r))
                    [D,Iclose] = closestElements(xloc,yloc,x(1:i),y(1:i),r(1:i));
                    if(any(D > minRadius))%avoid repeatedly placing tiny elements
                        %pare down D and I to elements within an acceptable reach
                        temp = D < maxReach;
                        D = D(temp);
                        Iclose = Iclose(temp);
                        %place a new element
                        i = i + 1;
                        x(i) = xloc;
                        y(i) = yloc;
                        r(i) = boundedRandom(randBounds);
                        [x(i),y(i),r(i),fail] = autoAdjust(...
                            i,x(1:i),y(1:i),r(1:i),D,Iclose);
                        if(~fail && (r(i)<=maxRadius) && (r(i)>=minRadius))
                            %successful placement
                            %DEMO
                                if(demo ~= 0)
                                    [mov, mov_count] = demoInteriorPlacement(xloc,yloc,x,y,r,...
                                        i-1,c1(j),c1w(k,:),Iclose,demo,mov,mov_count);
                                    demo = 0.975*demo;
                                end
                            %DEMO
                            placed = true;
                            break;
                        else %failed element placement, delete it
                            x(i) = 0;
                            y(i) = 0;
                            r(i) = 0;
                            i = i - 1;
                        end
                    end
                end
            end
        end
    end
end

%complex numbers indicate a failure
if(~(isreal(x) && isreal(y) && isreal(r)))
    error('Imaginary numbers generated');
end

%getting rid of empty spots
for i=1:length(r)
    if(r(i)==0)
        r(i:end)=[];
        x(i:end)=[];
        y(i:end)=[];
        break
    end
end

%DEMO
    if(demo ~= 0)
        plotElements(x,y,r,'off','k-');
        mov(mov_count) = getframe;
        mov = mov(1:mov_count);
        save('demo_movie', 'mov');
    end
%DEMO

%final count of elements
N=length(x);

%a few quick stats
fprintf('Initialized array complete\n');
fprintf('N=%d, nrows*ncols=%d, N-nrows*ncols=%+d, Bounds=[%d,%d], avgr=%g\n',...
    N,nrows*ncols,N - nrows*ncols,randBounds(1),randBounds(2),sum(r)/N);

end

end
