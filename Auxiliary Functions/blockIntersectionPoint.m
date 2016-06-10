function [xint,yint] = blockIntersectionPoint(x1,y1,x2,y2,line)
%blockIntersectionPoint finds the intersection point of a line between two
%input points and an input line.
%
%   x1,y1,x2,y2 - point coordinates
%   line - coefficients of line [slope,y-intercept]

if(x1==x2)
    %if the secondary line is vertical
    xint=x1;
    yint=line(1)*x1+line(2);
else
    %make line between input points
    P(1)=(y1-y2)/(x1-x2);
    P(2)=y1-P(1)*x1;
    %find x-coordinate of intersection point
    xint=(P(2)-line(2))/(line(1)-P(1));
    %find y-coordinate of intersection point
    yint=line(1)*xint+line(2);
end

end