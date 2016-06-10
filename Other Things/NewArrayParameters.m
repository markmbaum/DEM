function [r2,nrows2] = NewArrayParameters(r1,ncols1,nrows1,ncols2)
%NewArrayParameters calculates the size and number of rows necessary for a
%new array of uniformly sized elements in closest packing to be as close to
%the same size as another array as possible, using the element radii and
%dimensions of the other array.
%   [r,nrows] = NewArrayParameters(r1,ncols1,ncols2)
%
%   r1 - the radius of the elements in the previous array
%   ncols1 - number of columns in the previous array
%   ncols2 - number of columns that will be used in the new array

r2=r1*(2+sqrt(3)*(ncols1-1))/(2+sqrt(3)*(ncols2-1));

nrows2=round((r1*nrows1+r1-r2)/r2);

end

