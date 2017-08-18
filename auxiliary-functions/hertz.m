function [kn] = hertz(r1,r2,C,un)
%hertz calculates the normal force stiffness constant for two overlapping
%elements. It utilizes the approximation a~=sqrt(2*R*un). The overlap value
%is assumed to be that which is calculated in one of the FindOverlaps
%functions, where it is the difference between the distance between the
%element centers and the sum of the radii. The real strain would be half
%this value.
%   [kn] = hertz(r1,r2,C,un)
%
%   r1,r2 - element radii
%   C - combined elastic constant
%   un - normal overlap

R = (r1*r2)/(r1 + r2);

kn = C*((log(8*R/un)-1)^-1);

end 