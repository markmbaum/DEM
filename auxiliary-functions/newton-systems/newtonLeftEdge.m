function [xn,yn,rn] = newtonLeftEdge(xn,yn,rn,xp,yp,rp,xu,yu,ru,tol)
%newtonLeftEdge uses generalized Newton's method to adjust the position and
%size of an element being placed along the left edge.
%   [xn,yn,rn] = newtonLeftEdge(xn,yn,rn,xp,yp,rp,xu,yu,ru,tol)
%
%   xn,yn,rn - information for the element being adjusted
%   xp,yp,rp - information for the previously placed element
%   xu,yu,ru - information for the overlapping element
%   tol - termination tolerance using relative radius improvement

J = ones(3,3);
J(3,3) = -1;
J(3,2) = 0;
z = rn*ones(3,1);
while(abs(z(3))/rn > tol)
    %evaluate each function
    F = [rn + ru - sqrt((xn - xu)^2 + (yn - yu)^2);
        rn + rp - sqrt((xn - xp)^2 + (yn - yp)^2);
        xn - rn - xp + rp];
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