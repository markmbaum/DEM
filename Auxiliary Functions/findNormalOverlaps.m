function [I,U] = findNormalOverlaps(i,x,y,r,threshold)
%findNormalOverlaps searches the element array with respect to a specific 
%element for other elements that overlap with it, returning the indices of 
%the overlapping elements in I and the normal overlaps in U.
%   [I,U] = findNormalOverlaps(i,x,y,r,threshold);
%
%   i - index of element being investigated
%   x,y,r - element location and radius arrays

I = zeros(1,5);
U = zeros(1,5);

count = 1;
for n = [1:i-1,i+1:length(x)]
    temp = r(i)+r(n);
    if(abs(x(i) - x(n)) < temp)
        if(abs(y(i) - y(n)) < temp)
            u = r(n) + r(i) - sqrt((x(n) - x(i))^2+(y(n) - y(i))^2);
            if(u > (threshold*r(i)))
                U(count) = u;
                I(count) = n;
                count = count + 1;
            end
        end
    end
end

U(count:5) = [];
I(count:5) = [];

end