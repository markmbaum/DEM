function [a] = a_appx(r1,r2,U)
%test function for the approximate value of a

R=(r1*r2)/(r1+r2);

a=sqrt(2*U*R);

end

