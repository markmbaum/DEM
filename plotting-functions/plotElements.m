function [varargout] = plotElements(varargin)
%PlotElements draws circular line elements and angle bars to represent the
%discrete elements at a given iteration. PlotElements will clear the figure
%and replot.
%   [ElementHandles,LabelHandles]=plotElements(x,y,r,theta,labeling,attributes)
%   [ElementHandles,LabelHandles]=plotElements(x,y,r,labeling,attributes)
%   [ElementHandles,LabelHandles]=plotElements(E{i},r,labeling,attributes)
%
%   x,y,r,theta - full element location, radius, and angle arrays
%   E - structure containing position info of the elements
%   labeling - toggle whether the element numbers are shows next to their
%              centers and the text handles are returned
%   attributes - circle circumference color and style (dashed, dotted, etc)

switch nargin
     case 4
        x=varargin{1}.x;
        y=varargin{1}.y;
        theta=varargin{1}.theta;
        r=varargin{2};
    case 5
        x=varargin{1};
        y=varargin{2};
        r=varargin{3};
        theta=zeros(1,length(x));
    case 6
        x=varargin{1};
        y=varargin{2};
        r=varargin{3};
        theta=varargin{4};
    otherwise
        error('# of inputs incompatible. Remember labeling and attributes.');
end
labeling=varargin{end-1};
attributes=varargin{end};

%Removing elements with zero radius that are passed in if the array is
%still being assembled.
remove=(r==0);
if(~isempty(remove))
    x(remove)=[];
    y(remove)=[];
    r(remove)=[];
    theta(remove)=[];
end 

clf;

steps=100; %number of steps taken in drawing the circles

hold on;
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

ElementHandles=plot(a,b,attributes);
ElementHandles=[ElementHandles;plot(c,d,'color',[0 0.44 0.74])];

varargout{1}=ElementHandles;
    
axis equal;

if(strcmp(labeling,'on'))
    LabelHandles=zeros(1,length(x));
    for i=1:length(x)
        label=num2str(i);
        if(sin(theta)>=0)
            LabelHandles(i)=text(x(i),y(i),label,'color','r');
        else
            LabelHandles(i)=text(x(i),y(i),label,'color','r',...
                'HorizontalAlignment','right');
        end
    end
    varargout{2}=LabelHandles;
end

end