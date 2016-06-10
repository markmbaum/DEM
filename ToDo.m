%TO DO

%sonder, england paper, comparison with width-rotation ratio, 
%compare with thin sheet prediction 
%investigation of geometry

%friction coefficient, block sizes in relation to vel of zone and
%dimensions
%simple shear vs pure shear
%suggest possible field studies and other ways the model could be
%used, maybe with some added features

%are the element-block interactions stronger than the element-element
%interactions, by design?

%did I model half of a shear zone or the whole width? Where would the major
%central fault be?

%OTHER things the model could potentially do

%TO TAKE CARE OF THE NONVERTICAL FAULTING PROBLEM, RANDOMLY GENERATE ANGLES 
%BETWEEN 0-2PI FOR EACH ELEMENT THAT SCALE THE NORMAL FORCES. SO THAT 
%CONTACTS AT EXACTLY THE ANGLE INCREASE IN STRENGTH, THOSE EXACTLY OPPOSITE 
%THE ANGLE DECREASE, WITH A SMOOTH TRANSITION IN BETWEEN. THERE COULD BE 
%SEVERAL OF THESE ANGLES (MAYBE EVEN A RANDOM NUMBER OF ANGLES) AND THEY 
%COULD ROTATE WITH THE ELEMENT, SO THAT THE SIMULATED NONVERTICAL FAULT 
%PLANES MOVE WITH THE ELEMENTS. THIS IS KINDA COOL. A RANDOM NUMBER OF FAULT
%PLANES, EACH WITH A RANDOM DIP ANGLE (BOUNDED), AND EACH LOCATED AT A 
%RANDOM LOCATION (ANGLE) AROUND THE ELEMENTS (WITH SOME CONTROL).

%Be scaled to look at seismic time scales
%Implement a constant velocity field somehow, so that the elements have
%   forces and torques on them that correspond to a vector field and its
%   vorticity
%Modify blocks to look at transpression/transtension
%Build in gravity and look at something
%Build in a normal force slip criterion, like thrust faulting?, but how to
%   do normal faulting in reverse?
%arrangements 3&4, different zones of element sizes

%THIS COULD BE EASY, SORT OF LIKE ^
%use the arrays to look at earthquakes triggering other earthquakes? Move
%one of the boundary blocks until contacts between elements and block slip,
%see if it causes slip between elements in the array?