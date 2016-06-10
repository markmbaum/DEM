function [xn,yn,rn,fail] = newtonInterior3(xn,yn,rn,x,y,r,tol)
%newtonInterior3 uses generalized Newton's method to adjust the position and
%size of an element being placed in the interior. The new element is being
%set to achor 3 other elements.
%   [xn,yn,rn,fail] = newtonInterior3(xn,yn,rn,x,y,r,tol)
%
%   xn,yn,rn - information for the element being adjusted
%   x,y,r - information of anchor elements
%   tol - termination tolerance using relative radius improvement

J = ones(3,3);
z = rn*ones(3,1);
fail = false;
i = 0;
while(max(abs(z))/rn > tol)
    %evaluate each functions
    F = [rn + r(1) - sqrt((xn - x(1))^2 + (yn - y(1))^2);
        rn + r(2) - sqrt((xn - x(2))^2 + (yn - y(2))^2);
        rn + r(3) - sqrt((xn - x(3))^2 + (yn - y(3))^2)];
    %evaluate the Jacobian
    temp = ((xn - x(1))^2 + (yn - y(1))^2)^(-.5);
    J(1,1) = -(xn - x(1))*temp;
    J(1,2) = -(yn - y(1))*temp;
    temp = ((xn - x(2))^2 + (yn - y(2))^2)^(-.5);
    J(2,1) = -(xn - x(2))*temp;
    J(2,2) = -(yn - y(2))*temp;
    temp = ((xn - x(3))^2 + (yn - y(3))^2)^(-.5);
    J(3,1) = -(xn - x(3))*temp;
    J(3,2) = -(yn - y(3))*temp;
    %solve J*z = -F
    z = gaussianEliminationSolve(J,-F);
    %set improved estimate
    xn = xn + z(1);
    yn = yn + z(2);
    rn = rn + z(3);
    
    i = i + 1;
    if(i == 10)
        fail = true;
        break;
    end
end



end