function[S,spacing,sectioning]=makeGrid(x,y,r,corner2,corner3,nrows,ncols,i,demo)
%makeGrid returns a grid of points filling in the edge elements of the
%array being generated in BuildRandomRadiusArray. It generates line
%equations between adjacent edge elements, excluding ones that are trapped
%on the outside by larger elements, and used those lines as boundaries to
%fill in the array. Points that lie within elements are filtered out, then
%the grid is divided into sections if any of the elements in the left
%column touch any elements in the right column.
%   [S,spacing,sect]=makeGrid(x,y,r,corner2,corner3,nrows,ncols,i,demo)
%
%   x,y,r - full element location and radius array
%   corner2 - index of the element at the bottom of the leftmost column
%   corner3 - index of the element at the right end of the bottom row
%   nrows - number of rows
%   ncols - number of columns
%   i - index of element being adjusted
%   demo - toggle for illustrative plotting commands

%DEMO
    if(demo~=0)
        demofrac=1;
    end
%DEMO

%finding elements that are trapped on the outside of the array by larger
%elements and discounting them from subsequent steps
% trapped=zeros(1,i);
% place=1;
% for j=1:i-2 %using j for loops because i is the element counter
%     a=Adjacency(j,x,y,r);
%     a=a(a>j+1);
%     if(any(a>j+1))
%         if(any(j==1:(nrows-1))) 
%             for k=corner2+1:i %removing adjacencies in opposite columns
%                 if(any(a==k))
%                     a(a==k)=[];
%                 end
%             end
%         end
%         if(~isempty(a))
%             a=max(a);
%             insert=j+1:a-1;
%             trapped(place+1:place+length(insert))=insert;
%             place=place+length(insert);
%         end
%     end
% end
% trapped(trapped==0)=[];

%defining convenient variables
leftcolumn=1:nrows;
bottomrow=nrows:corner2;
rightcolumn=corner2:corner3;
toprow=[corner3:i,1];
edgeelements=[leftcolumn,bottomrow(2:end-1),rightcolumn,toprow(2:end-1)];

%generating linear equations between each pair of consecutive edge elements
%The polynomials are functions of y, and not x, because it is most useful
%to find the x-value of the lines between edge elements for a given y-value
%for filling in the boundary created by all the lines with a grid.
edgeelements=[edgeelements,1];
P=zeros(length(edgeelements)-1,2);
for j=1:length(edgeelements)-1
    xx=[x(edgeelements(j)), x(edgeelements(j+1))];
    yy=[y(edgeelements(j)), y(edgeelements(j+1))];
    if(diff(yy)==0)
        P(j,:)=[0,0];
    else
        P(j,:)=polyfit(yy,xx,1);
    end
end
edgeelements(end)=[];

%DEMO plotting edge element connection lines
    if(demo~=0)
        PlotElements(x,y,r,'off','k-');
        for k=1:length(P)-1
            yy=[y(edgeelements(k)),y(edgeelements(k+1))];
            xx=polyval(P(k,:),yy);
            line(xx,yy,'color','r');
        end
        line([x(edgeelements(end)),x(1)],[y(edgeelements(end)),y(1)],...
            'color','r');
        pause(demo/demofrac);
    end
%DEMO

%storing boundaries
leftin=max(x(leftcolumn));
rightin=min(x(rightcolumn));
topin=min(y(toprow));
bottomin=max(y(bottomrow));
leftout=min(x(leftcolumn));
rightout=max(x(rightcolumn));
topout=max(y(toprow));
bottomout=min(y(bottomrow));

%filling in the space defined by the polyfit lines with gridpoints
spacing=(mean(r)/10)*log(exp(1)*ncols)*.75;
xrow=leftin:spacing:rightin;
ycolumn=bottomin:spacing:topin;
fgrid=zeros(length(xrow)*length(ycolumn),2);
count=1;
for i=1:length(ycolumn)
    for j=1:length(xrow)
        fgrid(count,:)=[xrow(j),ycolumn(i)];
        count=count+1;
    end
end

%DEMO plotting grid that fits inside the lines
    if(demo~=0)
        fprintf('Uncontained interior grid size = %d points\n',length(fgrid));
        plot(fgrid(:,1),fgrid(:,2),'k.');
        pause(demo/demofrac);
    end
%DEMO

    %left flap
xflaprow=xrow(1)-spacing:-spacing:leftout;
yflapcolumn=[ycolumn(1):-spacing:y(leftcolumn(end)),...
    ycolumn(1)+spacing:spacing:y(leftcolumn(1))];
if(~isempty(xflaprow) && ~isempty(yflapcolumn))
flap=zeros(length(xflaprow)*length(yflapcolumn),2);
count=1;
for i=1:length(yflapcolumn)
    for j=1:length(xflaprow)
        
        %DEMO
            if(demo~=0)
                h=plot(xflaprow(j),yflapcolumn(i),'r.');
                %pause(demo);
            end
        %DEMO
        
        %find appropriate polynomial
        done=0;
        k=0;
        while(~done && k<(length(leftcolumn)-1))
            k=k+1;
            if(yflapcolumn(i)<=y(leftcolumn(k))&&...
                    yflapcolumn(i)>=y(leftcolumn(k+1)))
                done=1;
            end
        end
        %test if the point is inside the polynomial
        if(done)
            if(xflaprow(j)>polyval(P(leftcolumn(k),:),...
                    yflapcolumn(i)))
                flap(count,:)=[xflaprow(j),yflapcolumn(i)];
                
                %DEMO
                    if(demo)
                        delete(h);
                        plot(flap(count,1),flap(count,2),'g.');
                    end
                %DEMO
                
                count=count+1;
            end
        end
    end
end
i=0;
done=0;
l=size(flap,1);
while(~done && i<l)
    i=i+1;
    if((flap(i,1)==0) && (flap(i,2)==0))
        done=1;
    end
end
fgrid=[fgrid;flap(1:i-1,:)];
end

%DEMO
    if(demo~=0)
        pause(demo/demofrac);
    end
%DEMO

    %bottom flap
xflaprow=xrow;
yflapcolumn=ycolumn(1)-spacing:-spacing:bottomout;
if(~isempty(xflaprow) && ~isempty(yflapcolumn))
flap=zeros(length(xflaprow)*length(yflapcolumn),2);
count=1;
for i=1:length(yflapcolumn)
    for j=1:length(xflaprow)
        
        %DEMO
            if(demo~=0)
                h=plot(xflaprow(j),yflapcolumn(i),'r.');
                %pause(demo);
            end
        %DEMO
        
        %find appropriate polynomial
        done=0;
        k=0;
        while(~done && k<(length(bottomrow)-1))
            k=k+1;
            if(xflaprow(j)>=x(bottomrow(k))&&...
                    xflaprow(j)<=x(bottomrow(k+1)))
                done=1;
            end
        end
        %test if the point is inside the polynomial
        if(done)
        if(((xflaprow(j)>polyval(P(bottomrow(k),:),...
                yflapcolumn(i))) && P(bottomrow(k),1)<0) ||...
                ((xflaprow(j)<polyval(P(bottomrow(k),:),...
                yflapcolumn(i))) && P(bottomrow(k),1)>0))
            flap(count,:)=[xflaprow(j),yflapcolumn(i)];
            
            %DEMO
                if(demo~=0)
                    delete(h);
                    plot(flap(count,1),flap(count,2),'g.');
                end
            %DEMO
            
            count=count+1;
        end
        end
    end
end
i=0;
done=0;
l=size(flap,1);
while(~done && i<l)
    i=i+1;
    if((flap(i,1)==0) && (flap(i,2)==0))
        done=1;
    end
end
fgrid=[fgrid;flap(1:i-1,:)];
end

%DEMO
    if(demo~=0)
        pause(demo/demofrac);
    end
%DEMO

    %right flap
xflaprow=xrow(end)+spacing:spacing:rightout;
yflapcolumn=[ycolumn(1):-spacing:y(rightcolumn(1)),...
    ycolumn(1)+spacing:spacing:y(rightcolumn(end))];
if(~isempty(xflaprow) && ~isempty(yflapcolumn))
flap=zeros(length(xflaprow)*length(yflapcolumn),2);
count=1;
for i=1:length(yflapcolumn)
    for j=1:length(xflaprow)
        
        %DEMO
            if(demo~=0)
                h=plot(xflaprow(j),yflapcolumn(i),'r.');
                %pause(demo);
            end
        %DEMO
        
        %find appropriate polynomial
        done=0;
        k=0;
        while(~done && k<(length(rightcolumn)-1))
            k=k+1;
            if(yflapcolumn(i)>=y(rightcolumn(k))&&...
                    yflapcolumn(i)<=y(rightcolumn(k+1)))
                done=1;
            end
        end
        %test if the point is inside the polynomial
        if(done)
            if(xflaprow(j)<polyval(P(rightcolumn(k),:),...
                    yflapcolumn(i)))
                flap(count,:)=[xflaprow(j),yflapcolumn(i)];
               
                %DEMO
                    if(demo~=0)
                        delete(h);
                        plot(flap(count,1),flap(count,2),'g.');
                    end
                %DEMO
                
                count=count+1;
            end
        end
    end
end
i=0;
done=0;
l=size(flap,1);
while(~done && i<l)
    i=i+1;
    if((flap(i,1)==0) && (flap(i,2)==0))
        done=1;
    end
end
fgrid=[fgrid;flap(1:i-1,:)];
end

%DEMO
    if(demo~=0)
        pause(demo/demofrac);
    end
%DEMO

    %top flap
xflaprow=xrow;
yflapcolumn=ycolumn(end)+spacing:spacing:topout;
if(~isempty(xflaprow) && ~isempty(yflapcolumn))
flap=zeros(length(xflaprow)*length(yflapcolumn),2);
count=1;
for i=1:length(yflapcolumn)
    for j=1:length(xflaprow)
        
        %DEMO
            if(demo~=0)
                h=plot(xflaprow(j),yflapcolumn(i),'r.');
                %pause(demo);
            end
        %DEMO
        
        %find appropriate polynomial
        done=0;
        k=0;
        while(~done && k<(length(toprow)-1))
            k=k+1;
            if(xflaprow(j)<=x(toprow(k))&&...
                    xflaprow(j)>=x(toprow(k+1)))
                done=1;
            end
        end
        %test if the point is inside the polynomial
        if(done)
        if(((xflaprow(j)>polyval(P(toprow(k),:),...
                yflapcolumn(i))) && P(toprow(k),1)>0) ||...
                ((xflaprow(j)<polyval(P(toprow(k),:),...
                yflapcolumn(i))) && P(toprow(k),1)<0))
            flap(count,:)=[xflaprow(j),yflapcolumn(i)];
            
            %DEMO
                if(demo~=0)
                    delete(h);
                    plot(flap(count,1),flap(count,2),'g.');
                end
            %DEMO
            
            count=count+1;
        end
        end
    end
end
i=0;
done=0;
l=size(flap,1);
while(~done && i<l)
    i=i+1;
    if((flap(i,1)==0) && (flap(i,2)==0))
        done=1;
    end
end
fgrid=[fgrid;flap(1:i-1,:)];
end

%DEMO
    if(demo~=0)
        pause(demo/demofrac);
    end
%DEMO

%removing grid points that lie within elements and
[xx,yy]=removeContained(fgrid(:,1),fgrid(:,2),x,y,r);
fgrid=[xx,yy];

%DEMO plotting grid after containment check
    if(demo~=0)
    fprintf('Contained interior grid size = %d points\n',length(fgrid));
        plotElements(x,y,r,'off','k-');
        plot(fgrid(:,1),fgrid(:,2),'k.');
        pause(demo/demofrac);
    end
%DEMO

%checking adjacencies in the vertical columns to see if the grid has to be
%sectioned
sectioning=false;
j=1;
while(~sectioning && (j<=length(leftcolumn)))
    I=Adjacency(leftcolumn(j),x,y,r);
    for k=1:length(I)
        if(any(I(k)==rightcolumn(2:end)))
            sectioning=true;
        end
    end
    j=j+1;
end
j=1;
while(~sectioning && (j<=length(rightcolumn)))
    I=Adjacency(rightcolumn(j),x,y,r);
    for k=1:length(I)
        if(any(I(k)==leftcolumn(1:end-1)))
            sectioning=true;
        end
    end
    j=j+1;
end

%sorting the gridpoints into sections
if(sectioning)
    S=sectioning(fgrid(:,1),fgrid(:,2),spacing,demo);
    
    %DEMO
        if(demo~=0)
        fprintf('Interior grid separated into %d sections\n',length(S));
        end
    %DEMO
    
else
    S{1}.x=fgrid(:,1);
    S{1}.y=fgrid(:,2);
    
    %DEMO plotting the entire grid as one section if it isn't sectioned
        if(demo~=0)
            plot(S{1}.x,S{1}.y,'.','MarkerSize',12);
            pause(demo/demofrac);
        end
    %DEMO
    
end

end