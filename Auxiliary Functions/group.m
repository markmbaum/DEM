function [G] = group(x,y,r,maxr,radiusfactor)
%group finds the group of elements nearest an element, for each element,
%and returns the groups in a cell array G. Elements are in a group if their
%coordinates are within some number of multiples of the maximum radius of
%the array away from the element being considered. The first element in
%each group is the element being used as the center. It has to be included
%for the functions that find overlaps in the main program.
%   [G] = group(x,y,r,maxr,maxrdistance);
%
%   x,y - element locations
%   r - element radii
%   maxr - maximum radius in the array
%   radiusfactor - controls how big the groups are. Elements within
%                   radiusfactor*max(r) of the edge of the element are
%                   grouped together.

N=length(x);
G=cell(N,1);

d=maxr*radiusfactor;

for i=1:N
    temp=d+r(i);
    G(i) = {zeros(20,1)};
    G{i}(1)=i;
    count=2;
    for j=1:N
        if(i~=j)
           if(abs(x(j)-x(i))<temp)
               if(abs(y(j)-y(i))<temp)
                   G{i}(count)=j;
                   count=count+1;
               end
           end
        end
    end
    G{i}(G{i} == 0) = [];
end

end