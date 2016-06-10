 function [varargout] = DEM_v2(varargin)
%
%A Distinct Element Model to Investigate Crustal Rotation in Shear Zones
%
%Written by Mark Baum, 6/5/15 - present
%Continued from a Senior Thesis for the Dartmouth Department of Earth
%   Sciences, advised by Professor Leslie Sonder
%
%This distinct element code is meant to simulate crustal deformation in
%strike-slip deformation zones. It arranges an array of cylindrical
%elements that are subject to shear stresses by boundary blocks moving
%parallel to each other and in opposite directions. The elements and blocks
%behave like elastic solids, conforming to basic elastic theory that
%predicts the normal and shear forces exerted on the elements by other
%elements and by the boundary blocks.
%Code for the program PALLA by Marco Antonellini, written 1993-1994 in
%Fortran, was referenced in building the model.

if(nargin == 0)   %Without an input variables are defined in the menus below
%-------------------------------------------------------------------------
%Number, configuration, and size distribution of the elements

nrows = 2;            %Number of initial rows
ncols = 2;            %Number of initial columns
                        %For irregular arrangements, the number of rows and
                        %columns are defined along the outer edges. If the
                        %elements are randomly sized (arrangement==5),
                        %these inputs define the approximate number of
                        %elements in the outside edges only.
arrangement = 2;      %Different modes for the distribution of element sizes:
                        %1 - uniform radii, cubic packing
                        %2 - uniform radii, hexagonal (closest) packing
                        %3 - column sizes DEcreasing by factors of 2 into
                            %the center of the array
                        %4 - column sizes INcreasing by factors of 2 into
                            %the center of the array
                        %5 - randomly sized elements
meanRadius = 1.5e3;   %Average radius size (doesn't apply for arrangement 5)
%------------------------------------------
%For arrangement 5, randomly sized elements

randBounds = 1e3*[1,15];      %Size limits for the generation of random
                            %   element radii for arrangement 5
demo = 0;                     %Toggle for plotting commands placed throughout
                            %   the functions that generate the randomly
                            %   sized arrays. These functions are:
                            %   auto_adjust, auto_edge_adjust,
                            %   auto_new_element,
                            %   build_random_radius_array, fill_spaces,
                            %   make_grid, max_section_dimensions,
                            %   sectioning
%-------------------------------------------------------------------------
%Physical constants and characteristics

rho = 2800;           %Rock density
youngs_mod = 3.4e8;       %Young's modulus
shear_mod = 1;
poissonsq = .25;      %Poisson's ratio squared
staticFriction = .5;  %Static friction coefficient
kineticFriction = .4; %Kinetic friction coefficient
%-------------------------------------------------------------------------
%Boundary conditions

ftc = .4;             %factor scaling the time step from the critical time
                        %step at each iteration
                        %For arrangement 1-4, recommend .3-.4
                        %For arrangement 5, recommend .05-.2
fvc = .01;            %percentage of the critical shear velocity that adaptive
                        %blocks will automatically assume. Recommend .001-.1
                        %with the smaller values for larger arrays.
percHorzOff = .01;    %horizontal distance traveled by shearing blocks before
                        %their horizontal velocity is set to zero, as a
                        %percentage of the total original width of the
                        %array.
percVertEnd = 1;      %percentage of the total length of the array that, once
                        %traveled by one shearing block (relative to the
                        %other block and not necessarily relative to the
                        %array), the iterations are complete.
%-------------------------------------------------------------------------
%Program controls, or non-physical parameters

storeInterval = 1e1;  %Control how often element and block data are stored
                        %in a non-temporary array that will be saved
plotInterval = 1e2;   %Number of iterations between each replotting. Use
                        %zero to turn in-program plotting off.
progInterval = 1e2;   %Iterations between progress updates
regroupInterval = 100;%Iterations between nearest neighbor regroupings
regroupRadii = 1.5;   %Control the size of search groups. It's the factor
                        %multiplied by max(r) to use as the distance
                        %limit for finding nearest neighbors/groups
matSaving = 1;        %Control saving variables to .mat file
matOverwrite = 1;     %Control whether the .mat file is overwritten, 1=yes
ofile = 'completedtrial'; %Output filename, NO EXTENSION
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%Automatic input definitions with varargin.
%If a varaible is passed in, it must be a Nx2 cell array with the input
%variable names in the first column and the input variable values in the
%second column. This enables running trials with different parameters,
%unsupervised. However, to use eval() string variables need to be contained
%with triple apostrophes ''' in the cell array, and variable values that
%depend on previously evaluated variable values need only single
%apostrophes. With single apostrophes the string is evaluated as a variable
%name and with triple it remains a string. Numeric values are recognized
%and can just be doubles, but strings are the safest and will work too, so
%use strings for each variable that isn't a single number, like RandBounds.

elseif(nargin == 1)
    modelVariables = varargin{1};
    for i = 1:size(modelVariables,1)
        if(ischar(modelVariables{i,2}))
            evalc([modelVariables{i,1},'=',modelVariables{i,2}]);
        else
            evalc([modelVariables{i,1},'=',num2str(modelVariables{i,2})]);
        end
    end
else
    error('Multiple input variables not supported');
end
%-------------------------------------------------------------------------
%varargin must be cleared so that it isn't saved and later loaded into
    %other functions
clear varargin;
figure(1); %default plotting destination
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%MAIN PROGRAM

%Initializing the element array
switch arrangement
    case 1
        [xi,yi,r,N] = buildUniformRadiusArray(nrows,ncols,meanRadius,'c');
    case 2
        [xi,yi,r,N] = buildUniformRadiusArray(nrows,ncols,meanRadius,'h');
    case 3
        [xi,yi,r,N] = buildGradientRadiusArray(nrows,ncols,meanRadius,'-');
    case 4
        [xi,yi,r,N] = buildGradientRadiusArray(nrows,ncols,meanRadius,'+');
    case 5
        %put the try-catch code back in from version 1, but it messes with
        %the error output and is not helpful with debugging
                [xi,yi,r,N] = buildRandomRadiusArray(...
                    nrows,ncols,randBounds,demo);

end

%Allocate memory for Estore, an array of structures that stores data after
%an interval of iterations specified by storeinterval. Also declare xf and
%yf, which are used to advance positions to the next iteration, with theta
%variables to store element rotation angles.
temp = zeros(1,N);
Estruct = struct('x',temp,'y',temp,'theta',temp);
Estore = repmat(Estruct,1000,1);
Estore(1).x = xi;
Estore(1).y = yi;
Estore(1).theta = temp;
xf = temp;
yf = temp;
thetai = temp;
thetaf = temp;

%Allocating memory for and generating boundary blocks with another array of
%structures Bstore to record block positions over the trial. Blocks are
%named wb for west block, eb for east block, etc.
temp = zeros(2,2);
nb = temp; eb = temp; sb = temp; wb = temp;
[~,temp] = max(yi);
temp = yi(temp) + r(temp);
nb(:,2) = temp; wb(1,2) = temp; eb(1,2) = temp;
[~,temp] = min(yi);
temp = yi(temp) - r(temp);
sb(:,2) = temp; wb(2,2) = temp; eb(2,2) = temp;
[~,temp] = min(xi);
temp = xi(temp) - r(temp);
wb(:,1) = temp; nb(1,1) = temp; sb(1,1) = temp;
[~,temp] = max(xi);
temp = xi(temp) + r(temp);
eb(:,1) = temp; nb(2,1) = temp; sb(2,1) = temp;
wbTop0 = wb(1,2);
ebTop0 = eb(1,2);
width0 = eb(1,1) - wb(1,1);
temp = zeros(2,2);
Bstruct = struct('nb',temp,'eb',temp,'sb',temp,'wb',temp);
Bstore = repmat(Bstruct,1000,1);
Bstore(1).nb = nb;
Bstore(1).eb = eb;
Bstore(1).sb = sb;
Bstore(1).wb = wb;

%pre-calculate physical quantities
inv2m = zeros(N,1);
inv2I = zeros(N,1);
for i = 1:N
    inv2m(i) = 1/(2*rho*pi*r(i)^2);     %1/2m for displacement calculations
    inv2I(i) = 1/(rho*pi*r(i)^4);       %1/2I for rotation calculations
end
Chertz = ((pi*youngs_mod)/(1 - poissonsq))/2;
Cmindlin = (pi*youngs_mod)/(4*(1-poissonsq));
ravg = mean(r);

%calculate total horizontal and vertical distance the shearing blocks move
distVertEnd = percVertEnd*(nb(1,2) - sb(1,2));
distHorzOff = percHorzOff*(eb(1,1) - wb(1,1));

%Initialize time step and horizonatal block velocity. The horizontal block
%velocity is set to .01 times the critical shear block velocity because,
%even though the blocks will move horizontally first, the critical shear
%block velocity is still a good reference point.
temp = abs(hertz(mean(r),mean(r),Chertz,distHorzOff));
t = 10;%.25*sqrt((1/mean(inv2m))/temp);
bv = .01;%(staticFriction*temp*distHorzOff)/...
    %(100*mindlin(mean(r),mean(r),Cmindlin,distHorzOff)*t);
blockVert = 0;  %0 for blocks moving horizontally, 1 for vertically

%finding groups of neighbors for each element
maxr = max(r);
G = group(xi,yi,r,maxr,regroupRadii);
regroupCount = 0;

%Plotting the initialized array, storing handles for in-loop replotting
elementHandles = plotElements(xi,yi,r,'off','k-');
blockHandles = plotBlocks(nb,eb,sb,wb,bv,blockVert);
title('Initial array','fontsize',18,'interpreter','latex');
xlabel('x-dimension (m)','fontsize',16,'interpreter','latex');
ylabel('y-dimension (m)','fontsize',16,'interpreter','latex');
drawnow;
if(plotInterval ~= 0)
    plotInterval = round(plotInterval);
    fprintf('Plotting every %d iterations.\n',plotInterval);
end
plotCount = 0;
title('Trial In Progress','fontsize',18,'interpreter','latex');

%progress updates
if(progInterval ~= 0)
    progInterval = round(progInterval);
    fprintf('Running main program iterations...\n');
    progstr = ['Waiting for progress update at iteration ',...
        num2str(progInterval)];
    fprintf(progstr);
end
progCount = 0;

%Allocate space for the previous "intersection points" between elements,
%which must be stored for calculation of the shear forces and torques.
%Elements cannot exert forces on themselves and T(i,j) = T(j,i), so the
%number of (x,y) pairs stored needs to be sum(1:N-1), or N*(N-1)/2. It
%could be done easily with an NxN array and lots of empty space, but the
%amount of information stored can be minimized.
prevInt = zeros((N^2 - N)/2,2);

%storage counters
storeCount = 0;
storedIts = 1;

%initializing loop variables
distVert = 0;
i = 0;

%timing
startClock = clock;
clock1 = startClock;

%MAIN PROGRAM LOOP
while(distVert < distVertEnd)

    %move blocks
    temp = t*bv;
    if(blockVert)
        temp = t*bv;
        nb(1,2) = nb(1,2) + temp;
        nb(2,2) = nb(2,2) - temp;
        eb(:,2) = eb(:,2) - temp;
        sb(1,2) = sb(1,2) + temp;
        sb(2,2) = sb(2,2) - temp;
        wb(:,2) = wb(:,2) + temp;
    else
        nb(1,1) = nb(1,1) + temp;
        nb(2,1) = nb(2,1) - temp;
        eb(:,1) = eb(:,1) - temp;
        sb(1,1) = sb(1,1) + temp;
        sb(2,1) = sb(2,1) - temp;
        wb(:,1) = wb(:,1) + temp;
        if((width0 - (eb(1,1) - wb(1,1))) >= distHorzOff)
            blockVert = 1;
            itVertOn = i + 1;
        end
    end
    %record angle and slope of the north and south blocks for subsequent
    %force calculations
    bSlope = (nb(2,2) - nb(1,2))/(nb(2,1) - nb(1,1));
    bAngle = atan(bSlope);

    %move elements
    for j = 1:N

        %find overlapping elements, 'I' contains contacting element indices,
        %and U contains the respective overlap magnitudes.
        [I,U] = findNormalOverlaps(1,xi(G{j}),yi(G{j}),r(G{j}),0);
        I=G{j}(I);

        fx = 0;
        fy = 0;

        %forces between elements
        for k = 1:length(I)
            %NORMAL FORCE
            kn = hertz(r(j),r(I(k)),Chertz,U(k));
            fn = kn*U(k);
            %alpha is the angle between element centers
            alpha = abs(atan((yi(I(k)) - yi(j))/(xi(I(k)) - xi(j))));
            %x component
            if(xi(j) < xi(I(k)))
                fx = fx - fn*cos(alpha);
            else
                fx = fx + fn*cos(alpha);
            end
            %y component
            if(yi(j) < yi(I(k)))
                fy = fy - fn*sin(alpha);
            else
                fy = fy + fn*sin(alpha);
            end

%             %SHEAR FORCE
%             [xint,yint] = elementIntersectionPoint(xi(j),yi(j),r(j),...
%                 xi(I(k)),yi(I(k)),r(I(k)));
%             %find the index of the 1D storage array for the pair of
%             %elements under consideration, elements 'j' and 'I(k)'.
%             idx = map2(j,I(k));
%             %Interserction points initialize at zero, but without a real
%             %previous value, the shear force can't be computed.
%             if((prevInt(idx,1) ~= 0) && (prevInt(idx,2) ~= 0))
%                 
%                 %shear overlap
%                 s = 0;
%                 
%                 %first the component of S from intersection point migration
%                 a = xint - prevInt(idx,1);
%                 b = yint - prevInt(idx,2);
%                 
%                 %alpha is now angle of the intersection point's migration
%                 if(a == 0)
%                     alpha = pi/2;
%                 else
%                     alpha = abs(atan(b/a));
%                 end
%                 
%                 %phi is the angle from element center to intersection point
%                 phi = atan2(yint - yi(j),xint - xi(j));
%                 
%                 %gamma is the angle from element center to the previous 
%                 %intersection point
%                 gamma = atan2(prevInt(idx,2) - yi(j),prevInt(idx,1) - xi(j));
%                 
%                 %figure out which direction the intersecting element is 
%                 %moving around element j, clockwise (positive) or
%                 %counterclockwise (negative), and assign it to the
%                 %calculation for shear overlap 's'.
%                 temp = sqrt(a^2 + b^2)*abs(sin(alpha + phi));
%                 if(phi > (pi/2) && gamma < (-pi/2))
%                     s = s + 
%             end
        end

        %forces between elements and blocks

        %NORTH BLOCK
        if(yi(j) + r(j) > nb(2,2))
            %NORMAL FORCE
            U = r(j) - cos(bAngle)*(nb(1,2)+bSlope*(xi(j)-nb(1,1))-yi(j));
            if(U > 0)
                kn = hertz(r(j),ravg,Chertz,U);
                fn = kn*U;
                fx = fx + fn*sin(bAngle);
                fy = fy - fn*cos(bAngle);
            end

        end

        %EAST BLOCK
        U = (xi(j) + r(j)) - eb(1,1);
        if(U > 0)
            %NORMAL FORCE
            kn = hertz(r(j),ravg,Chertz,U);
            fn = kn*U;
            fx = fx - fn;

            %SHEAR FORCE

        end

        %SOUTH BLOCK
        if((yi(j) - r(j)) < sb(1,2))
            %NORMAL FORCE
            U = r(j) - cos(bAngle)*(yi(j)-sb(1,2)-bSlope*(xi(j)-sb(1,1)));
            if(U > 0)
                kn = hertz(r(j),ravg,Chertz,U);
                fn = kn*U;
                fx = fx - fn*sin(bAngle);
                fy = fy + fn*cos(bAngle);
            end
        end

        %WEST BLOCK
        U = wb(1,1) - (xi(j) - r(j));
        if(U > 0)
            %NORMAL FORCE
            kn = hertz(r(j),ravg,Chertz,U);
            fn = kn*U;
            fx = fx + fn;

            %SHEAR FORCE

        end

        xf(j) = xi(j) + fx*(t^2)*inv2m(j);
        yf(j) = yi(j) + fy*(t^2)*inv2m(j);
        thetaf(j) = thetaf(j) + .001*j*(2*rand(1) - 1);

    end
    %Tasks that must be done outside the main parfor loop to maintain
    %iteration independence
    for j = 2:N
        for k = 1:j-1
            [xint, yint] = elementIntersectionPoint(...
                xi(j),yi(j),r(j),xi(k),yi(k),r(k));
            prevInt(map2(j,k),:) = [xint, yint];
        end
    end

    %advance the positions of the elements
    xi = xf;
    yi = yf;
    thetai = thetaf;

    %regrouping (re-finding nearest neighbors to use in the search for
    %overlaps)
    regroupCount=regroupCount+1;
    if(regroupCount==regroupInterval)
        regroupCount=0;
        G=group(xi,yi,r,maxr,regroupRadii);
    end

    %update loop criteria
    i = i + 1;
    distVert = (wb(1,2) - wbTop0) + (ebTop0 - eb(1,2));

    %storing data, expanding storage variables when necesarry

    %JUST WRITE THIS SHIT TO A TEXT FILE, PREVENTING THE ARRAY FROM NEEDING TO
    %BE EXPANDED EVERY 1000 ITERATIONS, WHICH INVOLVES COPYING HUGE AMOUNTS OF
    %DATA FOR LARGE TRIALS
    storeCount = storeCount + 1;
    if(storeCount == storeInterval)
        storeCount = 0;
        storedIts = storedIts + 1;
        if(storedIts > length(Estore))
            temp = length(Estore);
            Estore(temp + 1 : temp + 1000) = repmat(Estruct,1000,1);
            Bstore(temp + 1 : temp + 1000) = repmat(Bstruct,1000,1);
        end
        Estore(storedIts).x = xi;
        Estore(storedIts).y = yi;
        Estore(storedIts).theta = thetai;
        Bstore(storedIts).nb = nb;
        Bstore(storedIts).eb = eb;
        Bstore(storedIts).sb = sb;
        Bstore(storedIts).wb = wb;
    end

    %progress display
    progCount = progCount + 1;
    if(progCount == progInterval)
        progCount = 0;
        fprintf(repmat('\b',1,length(progstr)));
        clock2 = clock;
        temp = progInterval/etime(clock2,clock1);
        clock1 = clock2;
        progstr = ['Iterations=',num2str(i),', ',...
            num2str(round(100*distVert/distVertEnd)),' percent shear, ',...
            'Speed=',num2str(temp),' its/s'];
        fprintf(progstr);
    end

    %in-loop plotting, can be disabled with plotInterval = 0
    plotCount = plotCount+1;
    if(plotCount == plotInterval);
        plotCount = 0;
        delete(elementHandles);
        delete(blockHandles);
        elementHandles = lineElements(xi,yi,r,thetai,'off','k','-');
        blockHandles = lineBlocks(nb,eb,sb,wb,bv,blockVert);
        drawnow;
    end
end

%calculate average speed and print progress update
averageSpeed=i/etime(clock,startClock);
if(progInterval~=0)
    fprintf(repmat('\b',1,length(progstr)));
end
fprintf('Main program iterations complete\n');
fprintf('%d Iterations, Average Speed=%g its/second\n',i,averageSpeed);

%getting rid of empty space in the storage arrays
Estore(storedIts+1:end)=[];
Bstore(storedIts+1:end)=[];

%tracing movement in a plot
traceCompletedTrial(Estore,r,Bstore,'on','on');
title('Completed Trial','fontsize',18,'interpreter','latex');
drawnow;

%set output variables, if any
varargout=cell(0);

%saving the current variables to a .mat file, if desired
if(matSaving)
    fprintf('Saving variables');
    if(matOverwrite)
        save(ofile);
    else
        ofile = matSaveName(ofile);
        save(ofile);
    end
    fprintf(' to: %s.mat\n',ofile);
end

%figure(2);
%grabData(Estore,'theta','all','on','simul',itVertOn,storeInterval,i);

totalTime = etime(clock,startClock);
hours = floor(totalTime/3600);
totalTime = totalTime - 3600*hours;
minutes = floor(totalTime/60);
totalTime = totalTime - 60*minutes;
seconds = totalTime;
fprintf('Total Time = %d hours, %d minutes, %f seconds\n',...
    hours,minutes,seconds);
fprintf('---Program complete\n\n');
end