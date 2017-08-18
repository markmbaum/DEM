function [ks] = mindlin(r1,r2,C,un)
%mindlin computes the shear stiffness constant for two elastic cylinders
%moving past eachother
%   [ks] = Mindlin(r1,r2,C,un)
%
%   r1,r2 - element radii
%   C - combined elastic constant
%   un - normal overlap

R = (r1*r2)/(r1 + r2);

ks=C*(un/R);
    
end