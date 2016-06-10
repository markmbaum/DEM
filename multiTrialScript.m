%This script can be used to run several different DEM trials consecutively
%without supervision by passing in different parts of the cell array
%ModelVariables. The array has to be edited so that the first column
%contains all of the necessary variables and the rest of the columns
%contain the values of those variables for separate trials. The loop will
%use the first column of values for a trial (the second column), then the
%third, etc. until there are no more columns.

starttime=clock;

%loading trial parameters, the model input variables
ModelVariables=load('ModelVariables.mat');
ModelVariables=ModelVariables.ModelVariables;
L=size(ModelVariables,2);

% %calculating time estimates for each trial
% maxsize=1:max([ModelVariables{1,2:end}].*[ModelVariables{2,2:end}]);
% speed=5704./(maxsize+8.477); %using curve fit to model speed measurements
% EstimatedTime=zeros(1,L-1);
% for i=2:L
%     if(ischar(ModelVariables{33,i}))
%         temp=str2double(ModelVariables{33,i}(2:end-1));
%     else
%         temp=ModelVariables{33,i};
%     end
%     EstimatedTime(i-1)=TrialTime(speed(ModelVariables{1,i}*...
%         ModelVariables{2,i}),abs(ModelVariables{21,i}-...
%         ModelVariables{23,i}),temp,ModelVariables{17,i});
% end
% fprintf('\nTOTAL TIME ESTIMATE = %g hours\n\n',sum(EstimatedTime));

%running each trial
for i=2:L
%     msg=['-Trial ',num2str(i-1),'/',num2str(L-1),', Estimated Time: ',...
%         num2str(EstimatedTime(i-1)),' hours\n'];
%     fprintf(msg);
    DistinctElementModel_Trial(ModelVariables(:,[1,i]));
end

endtime=clock;
Hours=etime(endtime,starttime)/3600;

fprintf('MultiTrialScript Complete, Total Time: %g\n',Hours);