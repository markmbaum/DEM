function highlightElements(varargin)
%highlightElements noticeably plots a set of elements and removes the plot
%soon after, so that specific elements can be found in large arrays.
%   highlightElements(E(i),r,n,pauselength)
%   highlightElements(x,y,r,n,pauselength)
%   highlightElements('filename',n,pauselength)
%
%   n - set of element indices to highlight
%   x,y,r,theta - full element location, radius, and angle arrays
%   E - structure containing position and size info of the elements
%   filename - saved .mat file with trial data, uses E(end)
%   pauselength - time that highlighting will last, use 0 for permanent

switch nargin
    case 3
        load(varargin{1});
        n=varargin{2};
        l=length(Estore);
        x=Estore(l).x(n);
        y=Estore(l).y(n);
        r=r(n);
        theta=Estore(l).theta(n);
        pauselength=varargin{3};
    case 4
        n=varargin{3};
        x=varargin{1}.x(n);
        y=varargin{1}.y(n);
        theta=varargin{1}.theta(n);
        r=varargin{2}(n);
        pauselength=varargin{4};
    case 5
        n=varargin{4};
        x=varargin{1}(n);
        y=varargin{2}(n);
        r=varargin{3}(n);
        theta=zeros(length(n));
        pauselength=varargin{5};
    otherwise
        error('Incompatible number of inputs');
end

steps=250;

sincircle=sin(0:(2*pi/steps):2*pi)';
coscircle=cos(0:(2*pi/steps):2*pi)';
a=zeros(length(sincircle),length(x));
b=a;
c=zeros(2,1);
d=c;
for i=1:length(x)
    a(:,i)=x(i)+r(i)*sincircle;
    b(:,i)=y(i)+r(i)*coscircle;
    c(:,i)=[x(i), x(i)+r(i)*sin(theta(i))];
    d(:,i)=[y(i), y(i)+r(i)*cos(theta(i))];
end

ElementHandles=line(a,b,'color',[0.85 0.32 0.09],'linewidth',3,'linestyle','-');
%ElementHandles=[ElementHandles;line(c,d,'color','r','linewidth',2)];


if(pauselength~=0)
    pause on;
    pause(pauselength);
    delete(ElementHandles);
end

end