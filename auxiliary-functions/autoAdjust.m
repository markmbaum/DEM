function [xn,yn,rn,fail] = autoAdjust(n,x,y,r,D,I)
%autoAdjust re-adjusts the location and size of an element based
%on the overlaps it makes with other elements. Solutions to nonlinear
%equations that describe zero overlap with other elements are solved using
%generalized Newton's method.
%   [xn,yn,rn,fail] = autoAdjust(n,x,y,r,D,I)
%
%   n - index of element being adjusted
%   x,y,r - information for all the elements placed so far, including the
%           one being adjusted
%   D - distances to three closest elements within reach, could be empty
%   I - indices of elements within reach of the new element, for 
%            anchoring, could be empty

xn = x(n); yn = y(n); rn = r(n);

fail = false;
count = 0;
go = true;
while(go)
    switch length(I)
        case 0
            %do nothing, the new element has no anchors
            
        case 1 %only a single anchor element, easy
            rn = D(1);
            
        case 2 %anchor the new element to 2 others, with xn and yn free
            [xn,yn,rn,fail] = newtonInterior2xy(xn,yn,rn,x(I),y(I),...
                r(I),1e-3);
            if(fail)
                [xn,yn,rn,fail] = newtonInterior2xr(xn,yn,rn,x(I),...
                    y(I),r(I),1e-3);
                if(fail)
                    [xn,yn,rn,fail] = newtonInterior2yr(xn,yn,rn,x(I),...
                        y(I),r(I),1e-3);
                    if(fail)
                        rn = D(1);
                    end
                end
            end
            
        otherwise %anchor to 3 other elements with xn, yn, and rn free
            [xn,yn,rn,fail] = newtonInterior3(xn,yn,rn,x(I(1:3)),...
                y(I(1:3)),r(I(1:3)),1e-3);
            if(fail)
                [xn,yn,rn,fail] = newtonInterior2xy(xn,yn,rn,x(I(1:2)),...
                    y(I(1:2)),r(I(1:2)),1e-3);
                if(fail)
                    [xn,yn,rn,fail] = newtonInterior2xr(xn,yn,rn,x(I),...
                        y(I),r(I),1e-3);
                    if(fail)
                        [xn,yn,rn,fail] = newtonInterior2yr(xn,yn,rn,...
                            x(I),y(I),r(I),1e-3);
                        if(fail)
                            rn = D(1);
                        end
                    end
                end
            end
    end
    
    %check for overlaps
    [D,I] = closestElements(xn,yn,x(1:n-1),y(1:n-1),r(1:n-1));
    if(D(1) > (rn - 0.03*rn))
        go = false;
    end
    
    %loop limit
    count = count + 1;
    if(count > 5)
        fail = true;
        go = false;
    end
end

end