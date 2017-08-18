function [I,U,A] = findAllOverlaps(varargin)
%findAllOverlaps searches the element array with respect to a specific element
%for other elements that overlap with it, returning the indices of the
%overlapping elements in I, the normal overlapping distances in U, and the
%perpendicular overlap distances in A.
%   [I,U,A] = findAllOverlaps(i,x,y,r,threshold);
%   [I,U,A] = findAllOverlaps(i,E,threshold);
%
%   i - index of element being investigated
%   x,y,r - full element location and radius arrays
%   E - structure containing position and size info of the elements

I=zeros(1,10);
U=zeros(1,10);
A=zeros(1,10);

i=varargin{1};
switch nargin
    case 5
        x=varargin{2};
        y=varargin{3};
        r=varargin{4};
        threshold=varargin{5};
    case 3
        E=varargin{2};
        x=E.x;
        y=E.y;
        r=E.r;
        threshold=varargin{3};
    otherwise
        error('FindOverlaps doesn''t support that number of inputs');
end

found=1;
for n=1:length(x)
    if(n~=i)
        u=normalOverlap(x(i),y(i),r(i),x(n),y(n),r(n),threshold);
        a=perpendicularOverlap(x(i),y(i),r(i),x(n),y(n),r(n),threshold);
        if(u>0)
            I(found)=n;
            U(found)=u;
            A(found)=a;
            found=found+1;
        end
    end
end

U(I==0)=[];
A(I==0)=[];
I(I==0)=[];

end