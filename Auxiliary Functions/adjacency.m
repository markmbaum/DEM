function [I] = adjacency(n,x,y,r)
%adjacency finds elements that are adjacent to a particular element and
%returns the indices of those elements. Because there is often a tiny
%amount of roundoff error or other error associated with the element sizes
%and positions, adjacency is determined by whether the distance between the
%elements is within 0.01% of the sum of their radii.
%   I = adjacency(n,x,y,r)
%
%   n - index of element being adjusted
%   x,y,r - full element location and radius arrays

I=zeros(1,5);
found=1;
for i=1:length(x)
    if(i~=n)
        a=(r(i)+r(n)-sqrt((x(i)-x(n))^2+(y(i)-y(n))^2));
        if(a>=0 || (abs(a)/((r(i)+r(n))/2))<.0001)
            I(found)=i;
            found=found+1;
        end
    end
end
I(I==0)=[];

end