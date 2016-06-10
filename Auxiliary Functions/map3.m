function [idx] = map3(i,j,k)
%map3 is a non-reversable hash from three integers to a single
%integer, I think. It can be used to store information for combinations of
%indices where the order of the indices doesn't matter and no indices can
%be repeated. If all combinations of the indices were stored it would take
%an N^3 array to store information for all of the combos. Mapping the
%indices reduces the array size to N^3/6 - N^2/2 + N/3, more 
%than an 80% reduction. It basically works by sorting the input indices, 
%so this tecnique might not be ideal for many more than 3 inputs.
%   [idx] = map3(i,j,k)
%
%   i,j,k - integers to be mapped
%   idx - mapped index

%optimal (I think) sorting of the three inputs, descending order
if(i > j) %i and j in right order
    if(j > k) %j and k in tight order
        a = i; b = j; c = k;
    else %j and k in wrong order
        if(k > i) %k goes first
            a = k; b = i; b = j;
        else %k goes second
            a = i; b = k; c = j;
        end
    end
else %i and j in wrong order
    if(k > j) %k goes first
        a = k; b = j; c = i;
    else
        if(k > i) %k goes second
            a = j; b = k; c = i; 
        else %k goes last
            a = j; b = i; c = k;
        end
    end
end

%Horner's method for efficient polynomial evaluation?
idx = ((a/6 - 1)*a + 11/6)*a + (.5*b - 1.5)*b + c;
            
end