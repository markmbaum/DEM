%%%Time testing BuildRandomRadiusArray

randBounds = [1,10];
demo = 0;

%test ncols and nrows up to 30 elements each
n = 8:15;
l = length(n);
t = zeros(1,l);

for i = n
    for temp = 1:10
        tic;
        buildRandomRadiusArray(i,i,randBounds,demo);
        t(i-l+1) = t(i-l+1) + toc;
    end
    t(i-l+1) = t(i-l+1)/10;
end