function [d] = distance(x1,y1,x2,y2)
%distance computes the euclidian distance between two points
%   distance(x1,y1,x2,y2)

d = sqrt((x1 - x2)^2 + (y1 - y2)^2);

end

