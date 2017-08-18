%This script plots rotations observed in the Las Vegas Valley Shear Zone
%against distance from the shear zone axis. Data comes from two sources:
%1) Sonder, L. J., C. H. Jones, S. L. Salyards, and K. M. Murphy, Vertical-
    %axis rotations in the Las Vegas Valley Shear Zone, southern Nevada: 
    %Paleomagnetic constraints on kinematics and dynamics of block 
    %rotations, Tectonics, 13 (4), 769-788, 1994.
%2) Nelson, M. R., and C. H. Jones, Paleomagnetism and crustal rotations 
    %along a shear zone, Las Vegas Range, southern Nevada, Tectonics, 6, 
    %13-33, 1987.
%The first set was taken from a traverse of the eastern part of the shear
%zone in the Gale Hills. The second set comes from farther west, in the Las
%Vegas Range.

data=load('LVVSZ Rotations');

DSond=data.DistFromLVVSZ_Sonder(:,1)*1e3;
DUncSond=data.DistFromLVVSZ_Sonder(:,2)*1e3;
RotSond=data.RotLVVSZ_Sonder(:,1)*(2*pi/360);
RotUncSond=data.RotLVVSZ_Sonder(:,2)*(2*pi/360);

DNels=data.DistFromLVVSZ_Nelson(:,1)*1e3;
RotNels=data.RotLVVSZ_Nelson(:,1)*(2*pi/360);
RotUncNels=data.RotLVVSZ_Nelson(:,2)*(2*pi/360);

figure(1);clf;hold on;
errorbar(DSond,RotSond,RotUncSond,'r*');
errorbar(DNels,RotNels,RotUncNels,'g*');
x=get(gca,'XLim');
line(x,[0 0],'color','k','linestyle','--');

title('Rotation vs. Distance From LVVSZ','interpreter','latex','fontsize',20);
ylabel('Rotation (radians)','interpreter','latex','fontsize',18);
xlabel('Distance From LVVSZ Axis (m)','interpreter','latex','fontsize',18);
legend({'Gale Hills Traverse (Eastern)','Las Vegas Range Traverse (Western)'},...
    'interpreter','latex','fontsize',13);