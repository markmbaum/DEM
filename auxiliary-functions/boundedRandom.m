function [r] = boundedRandom(in)
%boundedRandom returns a random number scaled between two input bounds.
%   [r] = boundedRandom(in)
%
%   in - vector with two numbers, the bounds of the desired random number

r = in(1) + (in(2) - in(1))*rand;

end