function [x,y,r,N] = autoNewElement(xloc,yloc,x,y,r,demo)
%autoNewElement Summary of this function goes here
%   Detailed explanation goes here

N = length(x);

%check if (xloc,yloc) is inside any existing elements
contained = isContained(xloc,yloc,x,y,r);
if(~contained)
    [D,Iclose] = closestElements(xloc,yloc,x(1:N),y(1:N),r(1:N));
    %place a new element
    i = N+1;
    x(i) = xloc;
    y(i) = yloc;
    r(i) = mean(r);
    [x(i),y(i),r(i),fail] = autoAdjust(...
        i,x(1:i),y(1:i),r(1:i),D,Iclose);                  
    if(fail)
        %failed element placement, delete it
        x(i) = [];
        y(i) = [];
        r(i) = [];
    else
        %DEMO
            if(demo ~= 0)
                plotElements(x,y,r,'off','k-');
            end
        %DEMO
        N = N + 1;
    end
else
    fprintf('(xloc,yloc) is in an element, refine the initial location\n');
end

end