function [a] = a_exact(r1,r2,U)
%test function for the exact value of a

d=(r1+r2-U)^2;

a=(1/2)*sqrt((4*d*r1^2-(d-r2^2+r1^2)^2)/d);

end

