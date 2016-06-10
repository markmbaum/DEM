function [varargout] = lineElements(varargin)
%lineElements draws circular line elements and angle bars to represent the
%discrete elements at a given iteration. LineElements will not clear the
%figure.
%   [ElementHandles,LabelHandles]=lineElements(x,y,r,theta,labeling,color,
%        style)
%   [ElementHandles,LabelHandles]=lineElements(x,y,r,labeling,color,style)
%   [ElementHandles,LabelHandles]=lineElements(E(i),r,labeling,color,style)
%
%   x,y,r,theta - full element location, radius, and angle arrays
%   E - structure containing position and size info of the elements
%   labeling - toggle whether the element numbers are shows next to their
%              centers and the text handles are returned
%   color - circle circumference line colors
%   style - circle circumference line styles (dashed, dotted, etc)

switch nargin
     case 5
        x=varargin{1}.x;
        y=varargin{1}.y;
        theta=varargin{1}.theta;
        r=varargin{2};
    case 6
        x=varargin{1};
        y=varargin{2};
        r=varargin{3};
        theta=zeros(1,length(x));
    case 7
        x=varargin{1};
        y=varargin{2};
        r=varargin{3};
        theta=varargin{4};
    otherwise
        error('# of inputs incompatible. Remember labeling and style.');
end
labeling=varargin{end-2};
color=varargin{end-1};
style=varargin{end};

%Removing elements with zero radius that are passed in if the array is
%still being assembled.
remove=(r==0);
if(~isempty(remove))
    x(remove)=[];
    y(remove)=[];
    r(remove)=[];
    theta(remove)=[];
end

steps=100; %number of steps taken in drawing the circles

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

ElementHandles=[line(a,b,'color',color,'linestyle',style);
                line(c,d,'color',[0 0.44 0.74])];

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