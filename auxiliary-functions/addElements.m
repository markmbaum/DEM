function [x,y,r,N] = addElements(x,y,r,I)
%addElements allows the user to put in extra elements after the array has
%been generated. 3 element indices are given in I and the new element is
%automatically placed in contact with them and adjusted if it overlaps with
%any other elements. For putting in multiple elements, I can be an Nx3
%matrix.
%   [x,y,r,N] = addElements(x,y,r,I)
%
%   x,y,r - full element location and radius arrays
%   I - indices of elements that the adjusted element must be adjacent to

if(size(I,2)~=3)
    error('length(I) must be 3');
end

for j=1:size(I,1)
    xloc=mean(x(I(j,:)));
    yloc=mean(y(I(j,:)));
    [D,~]=closestElements(xloc,yloc,x(I(j,:)),y(I(j,:)),r(I(j,:)));
    L=length(x)+1;
    [x(L),y(L),r(L)]=autoNewElement(x,y,r,xloc,yloc,I(j,:),D,[],0);
    [x(L),y(L),r(L)]=autoAdjust(L,x,y,r,I(j,:),1);
end

plotElements(x,y,r,'on','k-');

N=L;

end