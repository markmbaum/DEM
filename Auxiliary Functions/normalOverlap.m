function [u] = normalOverlap(x1,y1,r1,x2,y2,r2,threshold)
%mormalOverlap calculates the normal overlap of two circles, setting the
%return to zero if it's below a threshold percentage of the average radius
%or if it's imaginary (the circles don't overlap).
%   [u] = mormalOverlap(x1,y1,r1,x2,y2,r2)
%
%   x1,y1,r1 - position and radius of the first element
%   x2,y2,r2 - position and radius of the second element
%   threshold - 'high' or 'low'

u=r1+r2-sqrt((x1-x2)^2+(y1-y2)^2);

if(u<(threshold*r1))
    u = 0;
elseif(~isreal(u))
    u = 0;
end

end