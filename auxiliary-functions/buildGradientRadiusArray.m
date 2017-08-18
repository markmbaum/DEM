function[x,y,r]=buildGradientRadiusArray(nrows,ncols,meanRadius,inversion)
%buildGradientRadiusArray generates the locations and radii of element
%centers when the size of the elements is constant for each column but
%varies horizontally.
%   [x,y,r,N]=buildGradientRadiusArray(nrows,ncols,meanRadius,inversion)
%
%   nrows - number of rows
%   ncols - number of columns
%   meanRadius - the average radius of the array
%   inversion - toggle for the direction of the "gradient"
%               '+' means larger elements in the center
%               '-' means smaller elements in the center

%Counting the number of elements so that adjustments to radii can be made
%to conform to Mean Radius
    % N - total number of elements
N=0;
if(inversion=='-')
    for i=1:fix(ncols/2)
        N=N+nrows*2^(i-1);
    end
    N=N*2;
    if(mod(ncols,2)~=0)
        N=N+nrows*(2^i);
    end
elseif(inversion=='+')
    for i=1:fix(ncols/2)
        N=N+ceil(nrows*.5^(i-1));
    end
    N=N*2;
    if(mod(ncols,2)~=0)
        N=N+ceil(nrows*.5^i);
    end
end

%allocating space now that elements are counted
x=zeros(N,1);
y=zeros(N,1);
r=zeros(N,1);

%Finding the appropriate radii for MeanRadius
    %r1 - radius of the elements in the outside columns
    %temp - temporary variable
if(mod(ncols,2)==0 && inversion=='-')
    r1=(meanRadius*N)/(nrows*ncols);
elseif(mod(ncols,2)~=0 && inversion=='-')
    r1=(meanRadius*N)/(2*nrows*(ceil(ncols/2)-1/2));
elseif(mod(ncols,2)==0 && inversion=='+')
    temp=zeros(1,ncols/2);
    for i=1:ncols/2
        temp(i)=ceil(nrows/(2^(i-1)))*2^(i-1);
    end
    temp=sum(temp);
    r1=(meanRadius*N)/(2*temp);
elseif(mod(ncols,2)~=0 && inversion=='+');
    temp=zeros(1,ceil(ncols/2));
    for i=1:ceil(ncols/2)
        temp(i)=ceil(nrows/(2^(i-1)))*2^(i-1);
    end
    temp=sum(temp)-ceil(nrows/(2^(i-1)))*2^(i-2);
    r1=(meanRadius*N)/(2*temp);
else
    error('Calculation of initial radius r1 failed')
end

%Assigning the x,y coordinates of the element centers and their radii

syms q; %symbolic variable for use in the sum function symsum()
z=r1*sqrt(8); %geometric factor separating the x-coordinates of columns
place=0; %index placeholder variable
vec=[]; %index vector

if(inversion=='-')
    for i=1:ceil(ncols/2)
        vec=place+1:place+nrows*2^(i-1);  %assigning whole vectors instead
            x(vec)=z*double(symsum(2^(1-q),q,2,i));   % of using loops
            r(vec)=r1*2^(1-i);
            y(vec)=(vec-min(vec))*2*r(vec(1))+r(vec(1));
        place=vec(end);
    end
    if(mod(ncols,2)==0)
        vec=place+1:place+nrows*2^(i-1);
            x(vec)=x(place)+2*r(place);
            r(vec)=r(place);
            y(vec)=y(vec-length(vec));
        place=vec(end);
    end
    for j=1:fix(ncols/2)-abs(mod(ncols,2)-1) %The abs(mod(ncols,2)-1) term
        vec=place+1:place+nrows*2^(i-1-j);    %prevents the even middle
            x(vec)=x(place)+r(place)*sqrt(8); %column from being counted
            r(vec)=2*r(place);                %twice
            y(vec)=(vec-min(vec))*2*r(vec(1))+r(vec(1));
        place=vec(end);
    end
else
    for i=1:ceil(ncols/2)
        vec=place+1:place+ceil(nrows*2^(1-i));
            x(vec)=z*double(symsum(2^(q-2),q,2,i));
            r(vec)=r1*2^(i-1);
            y(vec)=(vec-min(vec))*2*r(vec(1))+r(vec(1));
        place=vec(end);
    end
    if(mod(ncols,2)==0)
        vec=place+1:place+ceil(nrows*2^(1-i));
            x(vec)=x(place)+2*r(place);
            r(vec)=r(place);
            y(vec)=y(vec-length(vec));
        place=vec(end);
    end
    for j=1:ncols/2-abs(mod(ncols,2)-1)
        vec=place+1:place+ceil(nrows*2^(1-i+j));
            x(vec)=x(place)+.5*r(place)*sqrt(8);
            r(vec)=.5*r(place);
            y(vec)=(vec-min(vec))*2*r(vec(1))+r(vec(1));
        place=vec(end);
    end
end

fprintf('N=%d, avgr=%g\n',N,mean(r));

end
