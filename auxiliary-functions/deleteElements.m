function [x,y,r,N] = deleteElements(x,y,r,I)
%deleteElements removes a set of elements from the array
%   [x,y,r,N] = deleteElements(x,y,r,I)
%
%   x,y,r - full element location and radius arrays
%   I - indices of elements to be deleted

for i=1:length(I)
    x(I(i))=-1;
    y(I(i))=-1;
    r(I(i))=-1;
end
x=x(x~=-1);
y=y(y~=-1);
r=r(r~=-1);

plotElements(x,y,r,'on','k-');

N=length(x);

end