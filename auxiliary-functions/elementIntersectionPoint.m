function [xint,yint] = elementIntersectionPoint(x1,y1,r1,x2,y2,r2)
%elementIntersectionPoint geometrically calculates the point that lies at the
%middle of the imaginary line connecting the two intersection points of two
%circles.
%   [xint,yint] = elementIntersectionPoint(x1,y1,r1,x2,y2,r2)

%solution adapted from:
%http://2000clicks.com/mathhelp/GeometryConicSectionCircleIntersection.aspx

C = (r1^2 - r2^2)/((x1 - x2)^2 + (y1 - y2)^2);

xint = .5*((x1+x2) + (x2-x1)*C);

yint = .5*((y1+y2) + (y2-y1)*C);

end
