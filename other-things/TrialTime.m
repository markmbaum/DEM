function [hours] = TrialTime(iterationspeed,blockvertspeed,distvertend,tstep)
%TrialTime estimates the number of hours required for a trial with the
%given parameters.
%   [hours] = TrialTime(iterationspeed,blockvertspeed,distvertend,tstep)

hours=(distvertend/blockvertspeed)*(1/tstep)*(1/iterationspeed)*(1/60^2);

end