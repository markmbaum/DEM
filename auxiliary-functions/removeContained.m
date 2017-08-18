function [xx,yy] = removeContained(xx,yy,x,y,r)
%removeContained deletes grid points that are inside of elements from x
%and y coordinate vectors
%   [xx,yy] = removeContained(xx,yy,x,y,r)
%
%   xx - vector of x-coordinates for the points
%   yy - vector of y-coordinates for the points
%   x,y,r - element information

if(length(x)==1)
    cand=zeros(1,length(xx));
    count=1;
    for j=1:length(xx)
        if((xx(j)>(x-r))&&(xx(j)<(x+r))&&(yy(j)>(y-r))&&(yy(j)<(y+r)))
            cand(count)=(j);
            count=count+1;
        end
    end
    cand=cand(cand~=0);
else
    cand=1:length(xx);
end

L=length(cand);
j=1;
M=length(x);
rem=zeros(1,length(cand));
remcount=1;
while(j<=L)
    found=0;
    k=1;
    while(~found && k<=M)
        d=sqrt((xx(cand(j))-x(k))^2+(yy(cand(j))-y(k))^2);
        if(d<r(k))
            found=true;
            rem(remcount)=cand(j);
            remcount=remcount+1;
        end
        k=k+1;
    end
    j=j+1;
end
rem=rem(rem~=0);
xx(rem)=[];
yy(rem)=[];
    
end