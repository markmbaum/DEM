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

nrows = 25;            %Number of initial rows
ncols = 15;            %Number of initial columns
                        %For irregular arrangements, the number of rows and
                        %columns are defined along the outer edges. If the
                        %elements are randomly sized (arrangement == 5),
                        %these inputs define the approximate number of
                        %elements in the outside edges only.
arrangement = 5;      %Different modes for the distribution of element sizes:
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

randBounds = 1e3*[1,10];      %Size limits for the generation of random
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
shear_mod = 1;          %Shear modulus
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

plotInterval = 1e4;   %Number of iterations between each replotting. Use
                        %zero to turn in-program plotting off.
progInterval = 1e3;   %Iterations between progress updates
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
        [xi, yi, r] = buildUniformRadiusArray(nrows, ncols, meanRadius, 'c');
    case 2
        [xi, yi, r] = buildUniformRadiusArray(nrows, ncols, meanRadius, 'h');
    case 3
        [xi, yi, r] = buildGradientRadiusArray(nrows, ncols, meanRadius, '-');
    case 4
        [xi, yi, r] = buildGradientRadiusArray(nrows, ncols, meanRadius, '+');
    case 5
        %put the try-catch code back in from version 1, but it messes with
        %the error output and is not helpful with debugging
        [xi, yi, r] = buildRandomRadiusArray(nrows, ncols, randBounds, demo);
end

%number of elements
N = length(xi);
%number of possible unique contacts between elements
N_contacts = (N^2 - N)/2;

%Allocate arrays used to advance element positions and rotations
xf = zeros(1,N);
yf = zeros(1,N);
thetai = zeros(1,N);
thetaf = zeros(1,N);

%Allocating memory for and generate boundary blocks, automatically at the exact
%edges of the element array. Blocks are named wb for west, eb for east, etc.
%Each block has four numbers stored in a 2x2 array. The rows of those arrays
%are x,y pairs stored from top to bottom and from left to right.
%So for the east block the array stores:
%  [top x coordinate, top y coordinate;
%   bottom x coordinate, bottom y coordinate]
%For the north block the array stores:
%  [left x coordinate, left y coordinate;
%   right x coorinate, right y coordinate]
%and so on...
%Extra variables like wbTop0 are used to track block progress over the duration
%of the trial
%allocate arrays
nb = zeros(2,2); eb = zeros(2,2); sb = zeros(2,2); wb = zeros(2,2);
%populate northern y coordinates of blocks
[~,idx] = max(yi);
y = yi(idx) + r(idx);
nb(:,2) = y; wb(1,2) = y; eb(1,2) = y;
%populate southern y coordinates of blocks
[~,idx] = min(yi);
y = yi(idx) - r(idx);
sb(:,2) = y; wb(2,2) = y; eb(2,2) = y;
%populate western x coordinates of blocks
[~,idx] = min(xi);
x = xi(idx) - r(idx);
wb(:,1) = x; nb(1,1) = x; sb(1,1) = x;
%populate eastern x coordinates of blocks
[~,idx] = max(xi);
x = xi(idx) + r(idx);
eb(:,1) = x; nb(2,1) = x; sb(2,1) = x;
%set tracking variables
wbTop0 = wb(1,2);
ebTop0 = eb(1,2);
width0 = eb(1,1) - wb(1,1);

%pre-calculate physical quantities, masses and moments of inertia
inv2m = zeros(N,1);
inv2I = zeros(N,1);
for i = 1:N
    inv2m(i) = 1/(2*rho*pi*r(i)^2);     %1/2m for displacement calculations
    inv2I(i) = 1/(rho*pi*r(i)^4);       %1/2I for rotation calculations
end
Chertz = ((pi*youngs_mod)/(1 - poissonsq))/2;
Cmindlin = (pi*youngs_mod)/(4*(1-poissonsq));
ravg = mean(r);

%calculate total horizontal and vertical distance the shearing blocks will move
distVertEnd = percVertEnd*(nb(1,2) - sb(1,2));
distHorzOff = percHorzOff*(eb(1,1) - wb(1,1));

%Initialize time step and horizonatal block velocity. The horizontal block
%velocity is set to .01 times the critical shear block velocity because,
%even though the blocks will move horizontally first, the critical shear
%block velocity is still a good reference point.
t = 5;
bv = .1;
blockVert = 0;  %0 for blocks moving horizontally, 1 for vertically

%find groups of neighbors for each element and set the regrouping counter
maxr = max(r);
G = group(xi, yi, r, maxr, regroupRadii);
regroupCount = 0;

%Plotting the initialized array, storing handles for in-loop replotting
elementHandles = plotElements(xi, yi, r, 'off', 'k-');
blockHandles = plotBlocks(nb, eb, sb, wb, bv, blockVert);
title('Initial array', 'fontsize',18, 'interpreter', 'latex');
xlabel('x-dimension (m)', 'fontsize',16, 'interpreter', 'latex');
ylabel('y-dimension (m)', 'fontsize',16, 'interpreter', 'latex');
drawnow;
if(plotInterval ~= 0)
    plotInterval = round(plotInterval);
    fprintf('Plotting every %d iterations.\n', plotInterval);
end
plotCount = 0;
title('Trial In Progress', 'fontsize', 18, 'interpreter', 'latex');

%progress updates
if(progInterval ~= 0)
    progInterval = round(progInterval);
    fprintf('Running main program iterations...\n');
    progstr = ['Waiting for progress update at iteration ',...
        num2str(progInterval)];
    fprintf(progstr);
end
progCount = 0;

%Allocate x,y pairs for the previous "intersection points" between elements,
%which must be stored for calculation of the shear forces and torques.
%Elements cannot exert forces on themselves and T(i,j) = T(j,i), so the
%number of (x,y) pairs stored needs to be sum(1:N-1), or N*(N-1)/2. It
%could be done easily with an NxN array and lots of empty space, but the
%amount of information stored can be minimized.
prevInt = zeros(N_contacts, 2);

%initializing loop variables
distVert = 0;
i = 0;

%timing
startClock = clock;
clock1 = startClock;

%MAIN PROGRAM LOOP
%advance until the blocks have moved vertically to the desired limit
while(distVert < distVertEnd)

    %move blocks
    d = t*bv;
    if(blockVert)
        nb(1,2) = nb(1,2) + d; %north block, left y coordinate
        nb(2,2) = nb(2,2) - d; %north block, right y coordinate
        eb(:,2) = eb(:,2) - d; %east block, both y coordinates
        sb(1,2) = sb(1,2) + d; %south block, left y coordinate
        sb(2,2) = sb(2,2) - d; %south block, right y coordinate
        wb(:,2) = wb(:,2) + d; %west block, both y coordinates
    else
        nb(1,1) = nb(1,1) + d; %north block, left x coordinate
        nb(2,1) = nb(2,1) - d; %north block, right x coordinate
        eb(:,1) = eb(:,1) - d; %east block, both x coordinates
        sb(1,1) = sb(1,1) + d; %south block, left x coordinate
        sb(2,1) = sb(2,1) - d; %south block, right x coordinate
        wb(:,1) = wb(:,1) + d; %west block, both x coordinates
        %when the blocks move in enough, switch to vertical motion by toggling
        %blockVert
        if((width0 - (eb(1,1) - wb(1,1))) >= distHorzOff)
            blockVert = 1;
        end
    end
    %record angle and slope of the north and south blocks for subsequent
    %force and torque calculations
    bSlope = (nb(2,2) - nb(1,2))/(nb(2,1) - nb(1,1));
    bAngle = atan(bSlope);

    %iterate over each element
    for j = 1:N

        %find overlapping elements, 'I' contains contacting element indices,
        %and U contains the respective overlap magnitudes.
        [I,U] = findNormalOverlaps(1, xi(G{j}), yi(G{j}), r(G{j}), 0);
        I=G{j}(I);

        %use running sums for the force components due to several contacts
        fx = 0;
        fy = 0;

        %forces between contacting elements
        for k = 1:length(I)
            %NORMAL FORCE
            fn = youngs_mod*U(k);
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

            %SHEAR FORCE
            [xint,yint] = elementIntersectionPoint(xi(j), yi(j), r(j),...
                xi(I(k)), yi(I(k)), r(I(k)));
            %find the index of the 1D storage array for the pair of
            %elements under consideration, elements 'j' and 'I(k)'.
            idx = map2(j,I(k));
            %Interserction points initialize at zero, but without a real
            %previous value, the shear force can't be computed. At the first
            %iteration the blocks contact, no shear force is computed. In
            %subsequent iterations there will be previous intersection points
            %and the shear force is computed.
            if((prevInt(idx,1) ~= 0) && (prevInt(idx,2) ~= 0))

                %shear overlap
                s = 0;

                %first the component of S from intersection point migration
                a = xint - prevInt(idx,1);
                b = yint - prevInt(idx,2);

                %alpha is now angle of the intersection point's migration
                if(a == 0)
                    alpha = pi/2;
                else
                    alpha = abs(atan(b/a));
                end

                %phi is the angle from element center to intersection point
                phi = atan2(yint - yi(j),xint - xi(j));

                %gamma is the angle from element center to the previous
                %intersection point
                gam = atan2(prevInt(idx,2) - yi(j),prevInt(idx,1) - xi(j));

                %figure out which direction the intersecting element is
                %moving around element j, clockwise (positive) or
                %counterclockwise (negative), and assign it to the
                %calculation for shear overlap 's'.
                temp = sqrt(a^2 + b^2)*abs(sin(alpha + phi));
                if(phi > (pi/2) && gam < (-pi/2))
                    %s = s +
                end
            else
                %store the intersection
                prevInt(idx,1) = xint;
                prevInt(idx,2) = yint;
            end
        end

        %forces between elements and blocks

        %NORTH BLOCK
        if(yi(j) + r(j) > nb(2,2))
            %NORMAL FORCE
            U = r(j) - cos(bAngle)*(nb(1,2) + bSlope*(xi(j) - nb(1,1)) - yi(j));
            if(U > 0)
                fn = youngs_mod*U;
                fx = fx + fn*sin(bAngle);
                fy = fy - fn*cos(bAngle);
            end

            %SHEAR FORCE
        end

        %EAST BLOCK
        U = (xi(j) + r(j)) - eb(1,1);
        if(U > 0)
            %NORMAL FORCE
            fn = youngs_mod*U;
            fx = fx - fn;

            %SHEAR FORCE
        end

        %SOUTH BLOCK
        if((yi(j) - r(j)) < sb(1,2))
            %NORMAL FORCE
            U = r(j) - cos(bAngle)*(yi(j) - sb(1,2) - bSlope*(xi(j) - sb(1,1)));
            if(U > 0)
                fn = youngs_mod*U;
                fx = fx - fn*sin(bAngle);
                fy = fy + fn*cos(bAngle);
            end

            %SHEAR FORCE
        end

        %WEST BLOCK
        U = wb(1,1) - (xi(j) - r(j));
        if(U > 0)
            %NORMAL FORCE
            fn = youngs_mod*U;
            fx = fx + fn;

            %SHEAR FORCE
        end

        %calculate element positions and rotations at next time step
        xf(j) = xi(j) + fx*(t^2)*inv2m(j);
        yf(j) = yi(j) + fy*(t^2)*inv2m(j);

    end

    %advance the positions of the elements
    xi = xf;
    yi = yf;
    thetai = thetaf;

    %regrouping (re-finding nearest neighbors to use in the search for overlaps)
    regroupCount = regroupCount + 1;
    if(regroupCount == regroupInterval)
        regroupCount = 0;
        G = group(xi, yi, r, maxr, regroupRadii);
    end

    %update loop end criteria
    i = i + 1;
    distVert = (wb(1,2) - wbTop0) + (ebTop0 - eb(1,2));

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
        elementHandles = lineElements(xi, yi, r, thetai, 'off', 'k',' -');
        blockHandles = lineBlocks(nb, eb, sb, wb, bv, blockVert);
        drawnow;
    end
end

%calculate average speed and print progress update
averageSpeed = i/etime(clock,startClock);
if(progInterval ~= 0)
    fprintf(repmat('\b', 1, length(progstr)));
end
fprintf('Main program iterations complete\n');
fprintf('%d Iterations, Average Speed=%g its/second\n', i, averageSpeed);

%set output variables, if any
varargout = cell(0);


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
