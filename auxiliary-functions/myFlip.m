function [f] = myFlip(v)
%myFlip rearranges a vector so that the elements are in the reverse order.
%Writing this function was necessary because Matlab's flip function is only
%supported in 2014, not earlier, and I figured the whole program should be
%portable to 2013 at least.
%   [f]=Flip(v)
%   
%   v - vector to flip
%   f - flipped vector

L=length(v);
f = zeros(1,L);

for i=1:L
    f(i) = v(L+1-i);
end

end