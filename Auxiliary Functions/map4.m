function [idx] = map4(i,j,k,l)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[s] = sort([i,j,k,l],'descend');
a=s(1);
b=s(2);
c=s(3);
d=s(4);

idx = (1/24)*(a-4)*(a-3)*(a-2)*(a-1) + (1/6)*(b-3)*(b-2)*(b-1) + ...
            (1/2)*(c-2)*(c-1) + d;

end

