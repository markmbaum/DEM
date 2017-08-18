function [xmin, xmax, ymin, ymax] = arrayEdges(x, y, r)
%Find the edges of an array of elements

[m, idx] = min(x);
xmin = m - r(idx);

[m, idx] = max(x);
xmax = m + r(idx);

[m, idx] = min(y);
ymin = m - r(idx);

[m, idx] = max(y);
ymax = m + r(idx);

end