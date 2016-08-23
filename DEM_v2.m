function [varargout] = DEM_v2(varargin)
%
%A Distinct Element Model to Investigate Crustal Rotation in Shear Zones
%
%Written by Mark Baum, 12/2/14 - present
%Continued from a Senior Thesis in the Dartmouth Department of Earth
%   Sciences, advised by Professor Leslie Sonder
%
%This distinct element code is meant to simulate crustal deformation in
%strike-slip deformation zones. It arranges an array of cylindrical
%elements that are subject to shear stresses by boundary blocks moving
%parallel to each other and in opposite directions. The elements and blocks
%behave like elastic solids, conforming to basic elastic theory that
%predicts the normal and shear forces exerted on the elements by other
%elements and by the boundary blocks.

if(nargin == 0)   %Without an input, variables are defined in the menus below
%-------------------------------------------------------------------------
%Number, configuration, and size distribution of the elements
nrows = 10;            %Number of initial rows
ncols = 6;            %Number of initial columns
                        %For irregular arrangements, the number of rows and
                        %columns are defined along the outer edges. If the
                        %elements are randomly sized (arrangement == 5),
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
mean_radius = 5e3;   %Average radius size (doesn't apply for arrangement 5) (m)
%------------------------------------------
%For arrangement 5, randomly sized elements

rand_bounds = 5e3*[1,10];      %Size limits for the generation of random
                            %   element radii for arrangement 5 (m)
demo = 0.1;                   %Toggle for plotting commands placed throughout
                            %   the functions that generate the randomly
                            %   sized arrays. These functions are:
                            %   auto_adjust, auto_edge_adjust,
                            %   auto_new_element,
                            %   build_random_radius_array, fill_spaces,
                            %   make_grid, max_section_dimensions,
                            %   sectioning
%-------------------------------------------------------------------------
%Physical constants and characteristics

depth_frac = 5;         %Brittle crust depth as fraction of average radius (-)
slip_rate = 3e-10
strain_rate = 1e-15     %Shear strain rate of underlying flow (1/s)
viscosity = 1e21;       %Viscosity of underlying ductile flow (kg/m*s)
rho = 2800;             %Rock density (kg/m^3)
youngs_mod = 3.4e8;     %Young's modulus (kg/m*s^2)
poissons = .5;          %Poisson's ratio (-)
shear_mod = youngs_mod/(2*(1 + poissons)); %Shear modulus (kg/m*s^2)
static_friction = .9;   %Static friction coefficient (-)
kinetic_friction = .4;  %Kinetic friction coefficient (-)
%-------------------------------------------------------------------------
%Boundary

perc_horz_off = .005;   %horizontal distance traveled by shearing blocks before
                        %their horizontal velocity is set to zero, as a
                        %percentage of the total original width of the
                        %array.
perc_vert_end = 1;      %percentage of the total length of the array that, once
                        %traveled by one shearing block (relative to the
                        %other block and not necessarily relative to the
                        %array), the iterations are complete.
%-------------------------------------------------------------------------
%Program controls
plot_interval = 1e2;   %Number of iterations between each replotting. Use
                       %zero to turn in-program plotting off.
prog_interval = 1e2;   %Iterations between progress updates
regroup_interval = 100;%Iterations between nearest neighbor regroupings
regroup_radii = 1.5;   %Control the size of search groups. It's the factor
                        %multiplied by max(r) to use as the distance
                        %limit for finding nearest neighbors/groups
saving = 0;             %toggle text file saving
save_interval = 1e3;    %Iterations between output file writing
file_name = 'completedtrial'; %Output filename, NO EXTENSION
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
%use strings for each variable that isn't a single number, like rand_bounds.

elseif(nargin == 1)
    modelVariables = varargin{1};
    for i = 1:size(modelVariables,1)
        if(ischar(modelVariables{i,2}))
            evalc([modelVariables{i,1},'=',modelVariables{i,2}]);
        else
            evalc([modelVariables{i,1},'=',num2str(modelVariables{i,2})]);
        end
    end
    %varargin must be cleared so that it isn't saved and later loaded into
    %other functions
    clear varargin;
else
    error('Multiple input variables not supported');
end
%-------------------------------------------------------------------------
%Initialize the element array
switch arrangement
    case 1
        [xi, yi, r] = buildUniformRadiusArray(nrows, ncols, mean_radius, 'c');
    case 2
        [xi, yi, r] = buildUniformRadiusArray(nrows, ncols, mean_radius, 'h');
    case 3
        [xi, yi, r] = buildGradientRadiusArray(nrows, ncols, mean_radius, '-');
    case 4
        [xi, yi, r] = buildGradientRadiusArray(nrows, ncols, mean_radius, '+');
    case 5
        %put the try-catch code back in from version 1, but it messes with
        %the error output and is not helpful with debugging
        [xi, yi, r] = buildRandomRadiusArray(nrows, ncols, rand_bounds, demo);
end
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%MAIN

N = length(xi); %number of elements
R_0 = max(r); %max radius of original array used as scaling length
[~, ~, ymin, ymax] = arrayEdges(xi, yi, r);
T = (ymax - ymin)/slip_rate %total time for 100 percent shear

%nondimensionalize
T_hat = T*(youngs_mod/viscosity)
slip_rate_hat = slip_rate*(viscosity/(youngs_mod*R_0))
r = r/R_0;
xi = xi/R_0;
yi = yi/R_0;

%Allocate arrays used to advance element positions and rotations, including
%arrays to store the previous iterations values
xl = xi(:,:);
xf = zeros(1,N);
yl = yi(:,:);
yf = zeros(1,N);
thetal = zeros(1,N);
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
%Extra variables like wbTop0 are used to track block progress over the
%duration of the trial
%allocate block arrays
nb = zeros(2,2); eb = zeros(2,2); sb = zeros(2,2); wb = zeros(2,2);
%get element array edges
[xmin,xmax,ymin,ymax] = arrayEdges(xi, yi, r);
%populate northern y coordinates of blocks
nb(:,2) = ymax; wb(1,2) = ymax; eb(1,2) = ymax;
%populate southern y coordinates of blocks
sb(:,2) = ymin; wb(2,2) = ymin; eb(2,2) = ymin;
%populate western x coordinates of blocks
wb(:,1) = xmin; nb(1,1) = xmin; sb(1,1) = xmin;
%populate eastern x coordinates of blocks
eb(:,1) = xmax; nb(2,1) = xmax; sb(2,1) = xmax;
%set tracking variables
wbTop0 = wb(1,2);
ebTop0 = eb(1,2);
width0 = eb(1,1) - wb(1,1);

%calculate physical quantities
max_r = max(r);
depth_hat = max_r*depth_frac;
M = zeros(N,1);
inv2m = zeros(N,1);
inv2I = zeros(N,1);
for i = 1:N
    M(i) = rho*depth_hat*pi*r(i)^2;
    inv2m(i) = 1/(2*depth_hat*rho*pi*r(i)^2);   %1/2m for displacement calculations
    inv2I(i) = 1/(depth_hat*rho*pi*r(i)^4);     %1/2I for rotation calculations
end

%calculate total horizontal and vertical distance the shearing blocks will move
dist_vert_end = perc_vert_end*(nb(1,2) - sb(1,2));
dist_horz_off = perc_horz_off*(eb(1,1) - wb(1,1));

%find groups of neighbors for each element and set the regrouping counter
G = group(xi, yi, r, max_r, regroup_radii);
regroup_count = 0;

%Initialize time step and horizonatal block velocity. The horizontal block
%velocity is set to .01 times the critical shear block velocity because,
%even though the blocks will move horizontally first, the critical shear
%block velocity is still a good reference point.
t = 0.5*min(r)*sqrt((2*pi*rho)/(youngs_mod))
block_vert = 0;  %0 for blocks moving horizontally, 1 for vertically
array_width = -1; %declaration

%Allocate variables for tracking contacts
%number of possible unique contacts between elements
N_contacts = (N^2 - N)/2;
%storing previous intersection points
prev_int = zeros(N_contacts,2); %element-element intersection points
%storing previous shear overlaps/displacements
S_e = zeros(N_contacts,1); %total shear overlap of element-element contacts
S_nb = zeros(N,1); %total shear overlap of north block contacts
S_eb = zeros(N,1); %total shear overlap of east block contacts
S_sb = zeros(N,1); %total shear overlap of south block contacts
S_wb = zeros(N,1); %total shear overlap of west block contacts

%S_store = zeros(1e5, N_contacts);
%fx_store = zeros(1e5, N);
%fy_store = zeros(1e5, N);
%tau_store = zeros(1e5, N);

%Plotting the initialized array, storing handles for in-loop replotting
figure(1); %default plotting destination
[e_handles, l_handles] = plotElements(xi, yi, r, 'on', 'k-');
blockHandles = plotBlocks(nb, eb, sb, wb, slip_rate_hat, block_vert);
title('Initial array', 'fontsize', 18, 'interpreter', 'latex');
xlabel('x-dimension (m)', 'fontsize', 16, 'interpreter', 'latex');
ylabel('y-dimension (m)', 'fontsize', 16, 'interpreter', 'latex');
drawnow;
if(plot_interval ~= 0)
    plot_interval = round(plot_interval);
    fprintf('Replotting every %d iterations.\n', plot_interval);
end
plot_count = 0;
title('Trial In Progress', 'fontsize', 18, 'interpreter', 'latex');

%progress updates
if(prog_interval ~= 0)
    prog_interval = round(prog_interval);
    fprintf('Running main program iterations...\n');
    progstr = ['Waiting for progress update at iteration ',...
        num2str(prog_interval)];
    fprintf(progstr);
end
prog_count = 0;

%output text file
save_count = 0;
if(saving)
    %open file object
    file_name = [file_name, '.txt'];
    ofile = fopen(file_name, 'w');
    %write the radii and headers
    fprintf(ofile, 'RADII\n');
    for n = 1:N-1
        fprintf(ofile, '%f,', r(n));
    end
    fprintf(ofile, '%f\n\nX,Y,THETA\n', r(N));
end

%initializing loop variables
dist_vert = 0;
i = 0;

%timing
start_clock = clock;
clock1 = start_clock;

%MAIN PROGRAM LOOP
%advance until the blocks have moved vertically to the desired limit
while(dist_vert < dist_vert_end)

    %move blocks
    d = t*slip_rate_hat;
    if(block_vert)
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
        %block_vert
        if((width0 - (eb(1,1) - wb(1,1))) >= dist_horz_off)
            block_vert = 1;
            array_width = eb(1,1) - wb(1,1);
        end
    end
    %record angle and slope of the north and south blocks for subsequent
    %force and torque calculations
    b_slope = (nb(2,2) - nb(1,2))/(nb(2,1) - nb(1,1));
    b_angle = atan(b_slope);

    %set flags to to avoid duplicating calculations for contacts
    %when a contact shared by elements is treated, the flag goes to false
    contact_flags = true(N_contacts,1);
    %arrays to store forces and torques between elements
    fx_e = zeros(N_contacts,1);
    fy_e = zeros(N_contacts,1);
    tau_e = zeros(N_contacts,1);

    %iterate over each element
    for j = 1:N

        %find overlapping elements, 'I' contains contacting element indices,
        %and U contains the respective overlap magnitudes.
        [I,U] = findNormalOverlaps(1, xi(G{j}), yi(G{j}), r(G{j}), 0);
        I = G{j}(I);

        %use running sums for the forces/torques over mulitple contacts
        fx_sum = 0;
        fy_sum = 0;
        tau_sum = 0;

        %forces between contacting elements
        for k = 1:length(I)

            %find the index of the 1D storage array for the pair of
            %elements under consideration, elements 'j' and 'I(k)'.
            idx = map2(j,I(k));
            %only perform calculations for one side of the contact, applying
            %equal and opposite forces for the contacting element
            if(contact_flags(idx))

                %NORMAL FORCE
                fn = youngs_mod*U(k)*depth_hat*R_0^2;
                %alpha is the angle between element centers
                alpha = abs(atan((yi(I(k)) - yi(j))/(xi(I(k)) - xi(j))));
                %x component
                if(xi(j) < xi(I(k)))
                    fx_e(idx) = -fn*cos(alpha);
                else
                    fx_e(idx) = fn*cos(alpha);
                end
                %y component
                if(yi(j) < yi(I(k)))
                    fy_e(idx) = -fn*sin(alpha);
                else
                    fy_e(idx) = fn*sin(alpha);
                end
                %add force components to running sums
                fx_sum = fx_sum + fx_e(idx);
                fy_sum = fy_sum + fy_e(idx);

                %SHEAR FORCE
                %get the "intersection point"
                [x_int,y_int] = elementIntersectionPoint(xi(j), yi(j), r(j),...
                  xi(I(k)), yi(I(k)), r(I(k)));
                %Interserction points initialize at zero, but without a real
                %previous value, the shear force can't be computed. At the first
                %iteration the blocks contact, no shear force is computed. In
                %subsequent iterations there will be previous intersection
                %points and the shear force is computed.
                if(all(prev_int(idx,:)))

                    %components of shear displacement due to intersection point
                    %migration, relative to the element
                    dx = x_int - prev_int(idx,1) - (xi(j) - xl(j));
                    dy = y_int - prev_int(idx,2) - (yi(j) - yl(j));
                    %angle of intersection point migration
                    alpha = atan2(dy, dx);
                    %components of vector from element center to intersection pt
                    a = x_int - xi(j);
                    b = y_int - yi(j);
                    %distance from element center to intersection point
                    dist = sqrt(a^2 + b^2);
                    %angle from element center to intersection point
                    phi = atan2(b, a);
                    %shear displacement parallel to the element surface
                    S_e(idx) = S_e(idx) - sqrt(dx^2 + dy^2)*sin(alpha + phi);

                    %component of shear disp due to the element's rotation
                    dtheta = thetai(j) - thetal(j);
                    S_e(idx) = S_e(idx) - dtheta*dist;
                    %and the contacing element's rotation
                    dtheta = thetai(I(k)) - thetal(I(k));
                    S_e(idx) = S_e(idx) - dtheta*sqrt((x_int - xi(I(k)))^2 + ...
                        (y_int - yi(I(k)))^2);

                    %check that the friction shear limit isn't exceeded
                    f_parallel = S_e(idx)*shear_mod*depth_hat*R_0^2;
                    if(abs(f_parallel) > fn*static_friction)
                        f_parallel = f_parallel*kinetic_friction;
                        S_e(idx) = 0;
                    end

                    %store the torque
                    tau_e(idx) = dist*f_parallel;
                    %add to the cumulative torque for the element
                    tau_sum = tau_sum + tau_e(idx);
                    %add to the cumulative fx and fy

                end
                %store the intersection point
                prev_int(idx,1) = x_int;
                prev_int(idx,2) = y_int;
                %set the contact flag
                contact_flags(idx) = false;
            else
                fx_sum = fx_sum - fx_e(idx);
                fy_sum = fy_sum - fy_e(idx);
                tau_sum = tau_sum + tau_e(idx);
            end
        end

        %forces between elements and blocks

        %NORTH BLOCK
        if(yi(j) + r(j) > nb(2,2))

            U = r(j) - cos(b_angle)*(nb(1,2) + b_slope*(xi(j) - nb(1,1)) - yi(j));
            if(U > 0)

                %NORMAL FORCE
                fn = youngs_mod*U*depth_hat*R_0^2;
                fx_sum = fx_sum + fn*sin(b_angle);
                fy_sum = fy_sum - fn*cos(b_angle);

                %SHEAR FORCE
            end

        end

        %EAST BLOCK
        U = (xi(j) + r(j)) - eb(1,1);
        if(U > 0)
            %NORMAL FORCE
            fn = youngs_mod*U*depth_hat*R_0^2;
            fx_sum = fx_sum - fn;

            %SHEAR FORCE
            %displacement
            S_eb(j) = S_eb(j) + d + (yi(j) - yl(j));

            %rotation
            dtheta = thetai(j) - thetal(j);
            S_eb(j) = S_eb(j) - dtheta*(r(j) - U);

            %compute shear force and check friction limit
            f_parallel = S_eb(j)*shear_mod*depth_hat*R_0^2;
            if(abs(f_parallel) > fn*static_friction)
              f_parallel = f_parallel*kinetic_friction;
              S_eb(j) = 0;
            end

            tau_sum = tau_sum + (r(j) - U)*f_parallel;
        end

        %SOUTH BLOCK
        if((yi(j) - r(j)) < sb(1,2))
            %NORMAL FORCE
            U = r(j) - cos(b_angle)*(yi(j) - sb(1,2) - b_slope*(xi(j) - sb(1,1)));
            if(U > 0)
                fn = youngs_mod*U*depth_hat*R_0^2;
                fx_sum = fx_sum - fn*sin(b_angle);
                fy_sum = fy_sum + fn*cos(b_angle);
            end

            %SHEAR FORCE
        end

        %WEST BLOCK
        U = wb(1,1) - (xi(j) - r(j));
        if(U > 0)
            %NORMAL FORCE
            fn = youngs_mod*U*depth_hat*R_0^2;
            fx_sum = fx_sum + fn;

            %SHEAR FORCE

            %intersection point migration
            S_wb(j) = S_wb(j) + d - (yi(j) - yl(j));

            %rotation
            dtheta = thetai(j) - thetal(j);
            S_wb(j) = S_wb(j) - dtheta*(r(j) - U);

            %compute shear force and check friction limit
            f_parallel = S_wb(j)*shear_mod*depth_hat*R_0^2;
            if(abs(f_parallel) > fn*static_friction)
              f_parallel = f_parallel*kinetic_friction;
              S_wb(j) = 0;
            end

            tau_sum = tau_sum + (r(j) - U)*f_parallel;
        end

        %UNDERLYING DUCTILE FLOW
        if(block_vert)

            %relative velocity
            %rel_vel = slip_rate_hat - (yi(j) - yl(j))/t;

            %y-axis traction
            %coef = 
            %fy = 
            %fy_sum = fy_sum + fy;

            %torque
            %coef = 
            %tau = 
            %tau_sum = tau_sum + tau;

            %if(dist_vert > dist_vert_end*0.25)
            %    temp = 0;
            %end
        end

        %store forces for varargout and post-trial inspection
        %if(i <= length(fx_store))
            %fx_store(i,idx) = fx_e(idx);
            %fy_store(i,idx) = fy_e(idx);
            %tau_store(i,idx) = tau_e(idx);
            %S_store(i,idx) = S_e(idx);
        %end

        %calculate element positions and rotations for next time step
        xf(j) = xi(j) + fx_sum*(t^2)*inv2m(j);
        yf(j) = yi(j) + fy_sum*(t^2)*inv2m(j);
        thetaf(j) = thetai(j) + tau_sum*(t^2)*inv2I(j);

    end

    %advance the positions of the elements, swapping arrays
    xl = xi;
    xi = xf;
    yl = yi;
    yi = yf;
    thetal = thetai;
    thetai = thetaf;

    %regrouping (re-finding nearest neighbors to use in the search for overlaps)
    regroup_count = regroup_count + 1;
    if(regroup_count == regroup_interval)
        regroup_count = 0;
        G = group(xi, yi, r, max_r, regroup_radii);
    end

    %update loop end criteria
    i = i + 1;
    dist_vert = (wb(1,2) - wbTop0) + (ebTop0 - eb(1,2));

    %output file writing
    if(saving)
        save_count = save_count + 1;
        if(save_count == save_interval)
            save_count = 0;
            writeVars(ofile, xf, yf, thetaf);
        end
    end

    %progress display
    prog_count = prog_count + 1;
    if(prog_count == prog_interval)
        prog_count = 0;
        fprintf(repmat('\b',1,length(progstr)));
        clock2 = clock;
        temp = prog_interval/etime(clock2,clock1);
        clock1 = clock2;
        progstr = ['Iterations = ',num2str(i),', ',...
            num2str(round(1000*dist_vert/dist_vert_end)/10),...
            ' percent shear, ', 'Speed = ',num2str(temp),' its/s'];
        fprintf(progstr);
    end

    %in-loop plotting, can be disabled with plot_interval = 0
    plot_count = plot_count + 1;
    if(plot_count == plot_interval);
        plot_count = 0;
        delete(e_handles);
        delete(l_handles)
        delete(blockHandles);
        [e_handles,l_handles] = lineElements(xi, yi, r, thetai,...
            'on', 'k', '-');
        blockHandles = lineBlocks(nb, eb, sb, wb, slip_rate_hat, block_vert);
        drawnow;
    end
end

%final plot
delete(e_handles);
delete(l_handles);
delete(blockHandles);
lineElements(xi, yi, r, thetai, 'on', 'k', '-');
lineBlocks(nb, eb, sb, wb, slip_rate_hat, block_vert);
drawnow;

%calculate average speed and print progress update
average_speed = i/etime(clock,start_clock);
if(prog_interval ~= 0)
    fprintf(repmat('\b', 1, length(progstr)));
end
fprintf('Main program iterations complete\n');
fprintf('%d Iterations, Average Speed = %g its/second\n', i-1, average_speed);

%set output variables, if any
%varargout = {fx_store(1:i,:), fy_store(1:i,:), tau_store(1:i,:)};

%print the final progress update
total_time = etime(clock,start_clock);
hours = floor(total_time/3600);
total_time = total_time - 3600*hours;
minutes = floor(total_time/60);
total_time = total_time - 60*minutes;
seconds = total_time;
fprintf('Total Time = %d hours, %d minutes, %f seconds\n',...
    hours,minutes,seconds);

%close output file
if(saving)
    fclose(ofile);
    fprintf('Results written to %s\n', file_name);
end

fprintf('---Program complete\n\n');

end
