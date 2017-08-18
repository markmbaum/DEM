function [xn,yn,rn,fail] = newtonInterior3_cubic(xn,yn,rn,x,y,r,tol)
%newtonInterior3 uses a modification of generalized Newton's method to
%adjust the position and size of an element being placed in the interior.
%The modification is intended to generate a cubically convergent iteration,
%and is described here:
%
%http://www.researchgate.net/publication/220564348_A_third-order_Newton-type_method_to_solve_systems_of_nonlinear_equations
%
%However, it was slower in practice to use this scheme instead of regular
%old quadratically convergent Newton's method.
%The new element is being set to achor 3 other elements.
%   [xn,yn,rn,fail] = newtonInterior3_cubic(xn,yn,rn,x,y,r,tol)
%
%   xn,yn,rn - information for the element being adjusted
%   x,y,r - information of anchor elements
%   tol - termination tolerance using relative radius improvement

J = ones(3,3);
z = rn*ones(3,1);
fail = false;
i = 0;
while(max(abs(z))/rn > tol)
    %evaluate each function
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
    %solve J*z = F for the intermediate estimate
    z = gaussianEliminationSolve(J,F);
    z(1) = xn - z(1);
    z(2) = yn - z(2);
    z(3) = rn - z(3);
    %evaluate F(z)
    Fstar = [z(3) + r(1) - sqrt((z(1) - x(1))^2 + (z(2) - y(1))^2);
        z(3) + r(2) - sqrt((z(1) - x(2))^2 + (z(2) - y(2))^2);
        z(3) + r(3) - sqrt((z(1) - x(3))^2 + (z(2) - y(3))^2)];
    %solve J*z = F + Fstar for the final improved guess
    z = gaussianEliminationSolve(J,F+Fstar);
    xn = xn - z(1);
    yn = yn - z(2);
    rn = rn - z(3);
    
    i = i + 1;
    if(i == 10)
        fail = true;
        break;
    end
end

end