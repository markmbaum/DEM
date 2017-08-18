function [a] = perpendicularOverlap(x1,y1,r1,x2,y2,r2,threshold)
%perpendicularOverlap geometrically computes HALF the maximum width of the 
%overlapping section of two circles with the radii of the circles and the 
%distance of the overlap/intersection. If a is imaginary or less than the
%input threshold percentage of the average radius, it is set to zero.
%   %perpendicularOverlap(x1,y1,r1,x2,y2,r2)
%
%   x1,y1,r1 - position and radius of the first element
%   x1,y1,r1 - position and radius of the second element

%formula from http://mathworld.wolfram.com/Circle-CircleIntersection.html
d=(x1-x2)^2+(y1-y2)^2;
a=(1/2)*sqrt((4*d*r1^2-(d-r2^2+r1^2)^2)/d);
%a=(1/(2*d))*sqrt((-d+r1-r2)*(-d-r1+r2)*(-d+r1+r2)*(d+r1+r2));

if((~isreal(a))||(a<(threshold*(r1+r2)/2)))
    a=0;
end

end