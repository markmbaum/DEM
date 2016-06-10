function [found] = isContained(xp,yp,x,y,r)
%isContained checks if a given point is contained within any
%elements, and returns true if it is, false if it isn't.
%   [found] = isContained(xp,yp,x,y,r)
%
%   xp - x-coordinate of the point under investigation
%   yp - y-coordinate of the point under investigation
%   x,y,r - full element location and radius arrays

found = false;

for i = 1:length(x)
    if(xp >= (x(i) - r(i)))
        if(yp <= (y(i) + r(i)))
            if((sqrt((xp - x(i))^2 + (yp - y(i))^2)) <= r(i))
                found = true;
                break;
            end
        end
    end
end

end