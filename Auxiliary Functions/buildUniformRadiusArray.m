function [x,y,r,N] = buildUniformRadiusArray(nrows,ncols,R,packing)
%buildUniformRadiusArray generates the locations of elements in an array
%where all the elements are the same size.
%   [x,y,r,N] = buildUniformRadiusArray(nrows,ncols,R,packing)
%
%   nrows - number of rows
%   ncols - number of columns
%   R - element radii
%   packing - toggle for cubic or hexagonal packing
%             'h' is hexagonal packing, anything else is cubic

%Assigning the x and y coordinates for cubic packing
N=nrows*ncols;
r=R*ones(N,1);
x=zeros(N,1);
y=x;
for i=1:nrows
    for j=1:ncols
        x((i-1)*ncols+j)=(j-1)*2*R;
        y((i-1)*ncols+j)=(i-1)*2*R;
    end
end

%Hexagonal packing
if(strcmp(packing,'h'))
    x=(sqrt(3)/2)*x;
    for i=1:2:ncols
        y(i:ncols:N)=y(i:ncols:N)+r(1);
    end
end

fprintf('N=%d, R=%g\n',N,R);
end