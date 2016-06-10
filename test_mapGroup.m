%test

k = 70; n = 72;

c = nchoosek(1:n,k);
l = size(c,1);

a = zeros(1,l);

for i = 1:l
    a(i) = mapGroup(c(i,:));
end

plot(a,'k.');
axis tight;
title('Unsorted Mapped Values');

title('Sorted Mapped Values');
fprintf('%d indices generated\n',nchoosek(n,k));
fprintf('%d non-unique indices generated\n\n',l - length(unique(a)));
pause(1);