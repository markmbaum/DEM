function [D,I] = closestElements(xp,yp,x,y,r)
%closestElements takes a grid point and finds the three element's
%edges that are closest to the point. Inputs are the point's coordinates,
%separately, and the element info arrays. D will contain negative values if
%the point lies inside of any elements.
%   [D,I] = closestElements(xp,yp,x,y,r)
%   
%   xp - x-coordinate of the point under investigation
%   yp - y-coordinate of the point under investigation
%   x,y,r - full element location and radius arrays

l = length(x);

edgedist=zeros(1,l);
for i = 1:l
    edgedist(i) = sqrt((xp - x(i))^2+(yp - y(i))^2) - r(i);
end

%get the indices of the sorted first three edge distances with optimal
%sorting of three elements
if(edgedist(1) < edgedist(2)) %first pair in correct order
    if(edgedist(2) > edgedist(3)) %second pair in incorrect order
        if(edgedist(3) < edgedist(1)) %third should be in first
            I = [3,1,2];
        else %third should be in second
            I = [1,3,2];
        end
    else %do nothing, they're already in order
        I = [1,2,3];
    end
else %first pair in incorrect order
    if(edgedist(1) > edgedist(3)) %first should be last
        if(edgedist(2) < edgedist(3)) %second should be first
            I = [2,3,1];
        else %second should be second
            I = [3,2,1];
        end
    else %first should be second
        I = [2,1,3];
    end
end

%It would suffice to simply sort the entire array of edge distances and
%return the first three distances and indices. However, for very large
%arrays the sorting algorithm is actually slower than the one used below.
%This is because time is wasted on sorting all of the entries beyond those
%first three desired ones. The process below will only require a single
%comparison for almost all of the elements by running through a logic tree,
%making it significantly faster.

for i=4:l
    if(edgedist(i) < edgedist(I(3)))
        if(edgedist(i) < edgedist(I(2)))
            if(edgedist(i) < edgedist(I(1)))
                I(3) = I(2);
                I(2) = I(1);
                I(1) = i;
            else
                I(3) = I(2);
                I(2) = i;
            end
        else
            I(3) = i;
        end
    end
end

D = [edgedist(I(1)),edgedist(I(2)),edgedist(I(3))];

end