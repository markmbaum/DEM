function traceCompletedTrial(varargin)
%traceCompletedTrial plots the final positions of elements and draws tracer
%lines along the element center's paths and the outside block corners
%   traceCompletedTrial(E,r,B,blocks,labeling)
%   traceCompletedTrial(filename,blocks,labeling)
%
%   E - element array of structs
%   r - radius array
%   B - block array of structs
%   filename - name of saved .mat file to load completed trial from
%   blocks - toggle whether the boundary blocks are plotted in their final
%               position.
%   labeling - toggle whether elements are labeled with their index numbers

switch nargin
    case 1
        load(varargin{1});
        blocks='off';
        labeling='off';
    case 2
        if(ischar(varargin{1}))
            load(varargin{1});
            blocks=varargin{2};
        else
            Estore=varargin{1};
            r=varargin{2};
            Bstore=[];
            blocks='on';
        end
        labeling='off';
    case 3
        if(ischar(varargin{1}))
            load(varargin{1});
            blocks=varargin{2};
            labeling=varargin{3};
        else
            Estore=varargin{1};
            r=varargin{2};
            Bstore=varargin{3};
            blocks='off';
            labeling='off';
        end
    case 5
        Estore=varargin{1};
        r=varargin{2};
        Bstore=varargin{3};
        blocks=varargin{4};
        labeling=varargin{5};
    otherwise
        error('Incompatible number of inputs');
end

clf;
L=length(Estore);
N=length(Estore(1).x);
xx=zeros(length(Estore),N);
yy=xx;
hold on;
for j=1:N
    for k=1:L
        xx(k,j)=Estore(k).x(j);
        yy(k,j)=Estore(k).y(j);
    end
end
line(xx,yy);
lineElements(Estore(end),r,labeling,'k','-');
if(strcmp(blocks,'on'))
    xx=zeros(L,4);
    yy=xx;
    for j=1:L
        xx(j,1)=Bstore(j).wb(1,1);
        xx(j,2)=Bstore(j).wb(2,1);
        xx(j,3)=Bstore(j).eb(1,1);
        xx(j,4)=Bstore(j).eb(2,1);
        yy(j,1)=Bstore(j).wb(1,2);
        yy(j,2)=Bstore(j).wb(2,2);
        yy(j,3)=Bstore(j).eb(1,2);
        yy(j,4)=Bstore(j).eb(2,2);
    end
     lineBlocks(Bstore(end));
    line(xx,yy,'color','k');
end

drawnow;

end