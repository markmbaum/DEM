function [xn,yn,rn] = autoEdgeAdjust(n,x,y,r,spec,demo)
%autoEdgeAdjust re-adjusts the size and location of an element based 
%on the overlaps it makes with other elements. The adjustments are made
%with solutions to nonlinear systems of equations that describe zero
%overlap with other elements and an adherance to the the edge in question.
%   [xn,y(n),rn] = autoEdgeAdjust(n,x,y,r,spec,demo)
%
%   x,y,r - information for the elements
%   spec - string denoting which edge the element is in
%          ('l', 'b', 'r', 'c3', or 't')
%   demo - toggle for illustrative plotting commands

[I,U] = findNormalOverlaps(n,x,y,r,1e-2);
if(~isempty(U))
    u = U(1)/r(n);
else
    u = 0;
end

while(u > 0.01)
    
    %DEMO plotting element before adjustments. Defining demofrac.
        if(demo~=0)
            lineElements(x(n),y(n),r(n),'off','k','-');
            pause(demo);
        end
    %DEMO
    
    if(spec == 'l')%left edge

        if(u > 0.5)
            [x(n),y(n),r(n)] = bisectionLeftEdge(...
                x(n),y(n),r(n),x(n-1),y(n-1),r(n-1),x(I(1)),y(I(1)),r(I(1)),1e-1);
        end
        [x(n),y(n),r(n)] = newtonLeftEdge(x(n),y(n),r(n),x(n-1),y(n-1),r(n-1),x(I(1)),y(I(1)),r(I(1)),1e-3);

    elseif(spec == 'b')%bottom edge

        if(u > 0.5)
            [x(n),y(n),r(n)] = bisectionBottomEdge(...
                x(n),y(n),r(n),x(n-1),y(n-1),r(n-1),x(I(1)),y(I(1)),r(I(1)),1e-1);
        end
        [x(n),y(n),r(n)] = newtonBottomEdge(x(n),y(n),r(n),x(n-1),y(n-1),r(n-1),x(I(1)),y(I(1)),r(I(1)),1e-3);

    elseif(strcmp(spec,'c3'))%"corner 3" the upper right or northeastern corner

        [x(n),y(n),r(n)] = newtonCorner3(x(n),y(n),r(n),x(n-1),y(n-1),r(n-1),x(I(1)),y(I(1)),r(I(1)),1e-3);

    elseif(spec == 'r')%right edge

        if(u > 0.5)
            [x(n),y(n),r(n)] = bisectionRightEdge(...
                x(n),y(n),r(n),x(n-1),y(n-1),r(n-1),x(I(1)),y(I(1)),r(I(1)),1e-1);
        end
        [x(n),y(n),r(n)] = newtonRightEdge(x(n),y(n),r(n),x(n-1),y(n-1),r(n-1),x(I(1)),y(I(1)),r(I(1)),1e-3);

    elseif(spec == 't')%top edge

        if(u > 0.5)
            [x(n),y(n),r(n)] = bisectionTopEdge(...
                x(n),y(n),r(n),x(n-1),y(n-1),r(n-1),x(I(1)),y(I(1)),r(I(1)),1e-1);
        end
        [x(n),y(n),r(n)] = newtonTopEdge(x(n),y(n),r(n),x(n-1),y(n-1),r(n-1),x(I(1)),y(I(1)),r(I(1)),1e-3);

    end

    [I,U] = findNormalOverlaps(n,x,y,r,1e-2);
    if(~isempty(U))
        u = U(1)/r(n);
    else
        u = 0;
    end
end

xn = x(n); yn = y(n); rn = r(n);

end