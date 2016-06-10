function [xn,yn,rn] = bisectionLeftEdge(~,~,rn,xp,yp,rp,xu,yu,ru,tol)
%bisectionLeftEdge uses a bisection root finding method to improve the
%position of an element placed on the left edge that overlaps another
%element.
%   bisectionLeftEdge(~,~,rn,xp,yp,rp,xu,yu,ru,tol)
%
%   rn - radius of the element being adjusted
%   xp,yp,rp - information of previously placed element
%   xu,yu,ru - information of overlapping element
%   tol - overlap tolerance using overlap relative to the element's radius

rhigh = rn;
rn = rhigh/2;
rlow = 0;
xn = rn + xp - rp;
yn = yp - sqrt((rp + rn)^2 - (rn - rp)^2);
u = rn + ru - sqrt((xn - xu)^2 + (yn - yu)^2);

while(abs(u/rn) > tol)
    if(u > 0)
        rhigh = rn;
        rn = (rhigh + rlow)/2;
    elseif(u < 0)
        rlow = rn;
        rn = (rhigh + rlow)/2;
    end
    xn = rn + xp - rp;
    yn = yp - sqrt((rp + rn)^2 - (rn - rp)^2);
    u = rn + ru - sqrt((xn - xu)^2 + (yn - yu)^2);
end

end