function [idx] = map6(i,j,k,l,m,n)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

s = sort([i,j,k,l,m,n],'descend');
a=s(1);
b=s(2);
c=s(3);
d=s(4);
e=s(5);
f=s(6);

idx = (1/720)*(a-6)*(a-5)*(a-4)*(a-3)*(a-2)*(a-1) +...
    (1/120)*(b-5)*(b-4)*(b-3)*(b-2)*(b-1) +...
    (1/24)*(c-4)*(c-3)*(c-2)*(c-1) +...
    (1/6)*(d-3)*(d-2)*(d-1) +...
    (1/2)*(e-2)*(e-1) +...
    f;
end

