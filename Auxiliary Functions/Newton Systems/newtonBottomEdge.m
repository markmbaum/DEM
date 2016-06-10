function [xn,yn,rn] = newtonBottomEdge(xn,yn,rn,xp,yp,rp,xu,yu,ru,tol)
%newtonBottomEdge uses generalized Newton's method to adjust the position and
%size of an element being placed along the left edge.
%   [xn,yn,rn] = newtonBottomEdge(xn,yn,rn,xp,yp,rp,xu,yu,ru,tol)
%
%   xn,yn,rn - information for the element being adjusted
%   xp,yp,rp - information for the previously placed element
%   xu,yu,ru - information for the overlapping element
%   tol - termination tolerance using relative radius improvement

J = ones(3,3);
J(3,3) = -1;
J(3,1) = 0;
z = rn*ones(3,1);
while(max(abs(z))/rn > tol)
    %evaluate each functions
    F = [rn + ru - sqrt((xn - xu)^2 + (yn - yu)^2);
        rn + rp - sqrt((xn - xp)^2 + (yn - yp)^2);
        yn - rn - yp + rp];
    %evaluate the Jacobian
    temp = ((xn - xu)^2 + (yn - yu)^2)^(-.5);
    J(1,1) = -(xn - xu)*temp;
    J(1,2) = -(yn - yu)*temp;
    temp = ((xn - xp)^2 + (yn - yp)^2)^(-.5);
    J(2,1) = -(xn - xp)*temp;
    J(2,2) = -(yn - yp)*temp;
    %solve J*z = -F
    z = gaussianEliminationSolve(J,-F);
    %set improved estimate
    xn = xn + z(1);
    yn = yn + z(2);
    rn = rn + z(3);
end

end