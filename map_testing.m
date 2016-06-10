%test functions that map combinations of integers into indices, minimizing
%the size of the array needed to store information for all possible
%combinations of those integers.

max_n = 20;
max_k = 6;

%iterate over n, the number of elements in the target list
for n = 3:max_n
    %iterate over k, the number of elements in combinations
    for k = 2:max_k
        %generate all combinations of k elements in a list of n elements
        I = nchoosek(1:n, k);
        a = zeros(1, size(I,1));
        %map each combination
        for i = 1:size(I, 1)
            a(i) = mapGroup(I(i,:));
        end
        %test for uniqueness
        if(length(a) ~= length(unique(a)))
            fprintf('Non-unique index returned for n = %d, k = %d', n, k);
        end
    end
end

fprintf('Everything was checked.\n');

%make a plot?
n = 16;
k = 4;
I = nchoosek(1:n, k);
a = zeros(1, size(I,1));
%map each combination
for i = 1:size(I, 1)
    a(i) = mapGroup(I(i,:));
end
plot(a, 'k.');
axis tight;