function [xn,yn,rn,fail] = newtonInterior2xr(xn,yn,rn,x,y,r,tol)
%newtonInterior2 uses generalized Newton's method to adjust the position and
%size of an element being placed in the interior. The new element is being
%set to achor 2 other elements, with its x-coordinate and radius free.
%   [xn,yn,rn] = newtonInterior2(xn,yn,rn,x,y,r,tol)
%
%   xn,yn,rn - information for the element being adjusted
%   x,y,r - information of anchor elements
%   tol - termination tolerance using relative radius improvement

J = ones(2,2);
z = rn*ones(2,1);
fail = false;
i = 0;
while(max(abs(z))/rn > tol)
    %evaluate each functions
    F = [rn + r(1) - sqrt((xn - x(1))^2 + (yn - y(1))^2);
        rn + r(2) - sqrt((xn - x(2))^2 + (yn - y(2))^2)];
    %evaluate the Jacobian
    J(1,1) = -(xn - x(1))*((xn - x(1))^2 + (yn - y(1))^2)^(-.5);
    J(2,1) = -(xn - x(2))*((xn - x(2))^2 + (yn - y(2))^2)^(-.5);
    %solve J*z = -F
    z = gaussianEliminationSolve(J,-F);
    %set improved estimate
    xn = xn + z(1);
    rn = rn + z(2);
    
    i = i + 1;
    if(i == 10)
        fail = true;
        break;
    end
end

end