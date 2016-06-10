function [idx] = findMiddleElements(E,r,bandlength)
%findMiddleElements assesses the dimensions of an array and finds elements
%that are within some distance away from a line across the middle of the
%aray.
%   findMiddleElements(E,bandlength)
%
%   E - single structure, if an array of structures of passed the
%           first one is ised
%   r - element radius array
%   bandlength - multiplied by the average radius and used as a threshold
%                   for which elements are included

if(length(E)>1)
    E=E(1);
end

N=length(E.x);

%find approximate middle vertical point of initial array
middle=(max(E.y)+min(E.y))/2;

%define the size of the band of included elements
verticalband=[middle+mean(r)*bandlength,...
    middle-mean(r)*bandlength];

%find elements initially within the band
idx=zeros(1,length(E.x));
count=1;
for j=1:N
    if(E.y(j)<verticalband(1) && E.y(j)>verticalband(2))
        idx(count)=j;
        count=count+1;
    end
end
idx(idx==0)=[];

%bubble sort element indices by horizontal locations
l=length(idx);
done=0;
while(~done)
    done=1;
    for j=1:l-1
        if(E(1).x(idx(j))>E.x(idx(j+1)))
            temp=idx(j);
            idx(j)=idx(j+1);
            idx(j+1)=temp;
            done=0;
        end
    end
end

end