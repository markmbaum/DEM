function idx = mapGroup(g)
%Evaluate a map_k polynomial for an arbitrary number of input indices,
%where g is a vector of the indices.

g = sort(g);

idx = g(1);

for i = 2:length(g)
    p = 1/factorial(i);
    for j = 1:i
        p = p*(g(i) - j);
    end
    idx = idx + p;
end

end

