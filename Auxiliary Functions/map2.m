function [idx] = map2(i,j)
%map2 accepts the indices of two elements and maps them to one index.
%It's used to reduce the size of the array storing previous intersection
%points between two elements from N^2 - N to (N^2 - N)/2, which is
%reduction by at least 50%, by making the order of the indices non-unique.
%Identical incoming indices are not allowed because an element can't exert
%a force on itself.
%   [idx] = map2(i,j)
%
%   i,j - element indices
%   idx - return index

if(i>j)
    a = i;
    b = j;
else
    a = j;
    b = i;
end

%Horner's method for efficient polynomial evaluation
idx = a*(.5*a - 1.5) + 1 + b;

end

