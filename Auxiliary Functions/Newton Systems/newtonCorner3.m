function [xn,yn,rn] = newtonCorner3(xn,yn,rn,xp,~,rp,xu,yu,ru,tol)
%newtonBottomEdge uses generalized Newton's method to adjust the position and
%size of an element being placed along the left edge. Corner 3 refers to
%the upper right corner of the array because it's sort of the third corner
%completed.
%   [xn,yn,rn] = newtonCorner3(xn,yn,rn,xp,yp,rp,xu,yu,ru,tol)
%
%   xn,yn,rn - information for the element being adjusted
%   xp,yp,rp - information for the previously placed element
%   xu,yu,ru - information for the overlapping element
%   tol - termination tolerance using relative radius improvement

J = ones(2,2);
z = rn*ones(2,1);
while(abs(z(2))/rn > tol)
    %evaluate each functions
    F = [rn + ru - sqrt((xn - xu)^2 + (yn - yu)^2);
        xn + rn - xp - rp];
    %evaluate the Jacobian
    J(1,1) = -(xn - xu)*((xn - xu)^2 + (yn - yu)^2)^(-.5);
    %solve J*z = -F
    z = gaussianEliminationSolve(J,-F);
    %set improved estimate
    xn = xn + z(1);
    rn = rn + z(2);
end

end