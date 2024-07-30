clc;
clear;
close all;
densityWater = 998;%kg/m^3
dynVisWater = 0.0010016	;%at room temp!
densityAir = 1.204;%kg/m^3 at room temp!!!
dynVisAir = 1.81*10^-5;%at room temp!
gravity = 9.81;%m/s^2
Lmintom3 = 1/60000;
mmtom = 1/1000;
bartoPa = 100000;
inTo_m = 0.0254;
miPerHr_to_m_s = 5280*12*0.0254/3600;
lb2N = 9.81/2.2;
resolution = 100;

sphereDiam = 2*inTo_m;%in m
areaCrossSphere = pi/4 * sphereDiam^2;

windSpeedRange = 'A4:E12';
dragShapesRange = 'A4:E8';
PAirfoilRange = 'A5:K16';
PCylinderRange = 'A5:F23';

closedFlapsClosedSlatsRange = 'A4:N6';
openFlapsClosedSlatsRange = 'A12:N14';
closedFlapsOpenedslatsRange = 'A20:N22';
liftDragAirfoilData = {readtable("EGME306B_Ex04_FinalData",'Sheet',"LiftDragAirfoil",'Range',closedFlapsClosedSlatsRange),readtable("EGME306B_Ex04_FinalData",'Sheet',"LiftDragAirfoil",'Range',openFlapsClosedSlatsRange),readtable("EGME306B_Ex04_FinalData",'Sheet',"LiftDragAirfoil",'Range',closedFlapsOpenedslatsRange)};

labData = {readtable("EGME306B_Ex04_FinalData",'Sheet',"WindSpeed",'Range',windSpeedRange),readtable("EGME306B_Ex04_FinalData",'Sheet',"DragShapes",'Range',dragShapesRange),liftDragAirfoilData,readtable("EGME306B_Ex04_FinalData",'Sheet',"PAirfoil",'Range',PAirfoilRange),readtable("EGME306B_Ex04_FinalData",'Sheet',"PCylinder",'Range',PCylinderRange)};
%manometer heights are in inches
windTunnelSpecificGravities = [labData{1,1}{8,2},labData{1,1}{8,4}];
windTunnelAngles = [labData{1,1}{9,2},labData{1,1}{9,4}];%in degrees
labData{1,1}(8:9,:) = []; %already have SG and angle, so can delete them from table to make things easier

%%converting to SI
labData{1,1}{:,1} = labData{1,1}{:,1} * miPerHr_to_m_s;
labData{1,1}{:,2:5} = labData{1,1}{:,2:5} * inTo_m;
labData{1,2}{:,1} = labData{1,2}{:,1} * miPerHr_to_m_s;
labData{1,2}{:,2:5} = labData{1,2}{:,2:5} * lb2N;
labData{1,3}{1,1}{:,2:end} = labData{1,3}{1,1}{:,2:end} * lb2N;
labData{1,3}{1,2}{:,2:end} = labData{1,3}{1,2}{:,2:end} * lb2N;
labData{1,3}{1,3}{:,2:end} = labData{1,3}{1,3}{:,2:end} * lb2N;
labData{1,4}{:,2:end} = labData{1,4}{:,2:end}*inTo_m;
labData{1,5}{:,3:end} = labData{1,5}{:,3:end}*inTo_m;

%%Finding Actual Wind Tunnel Velocities
windTunnel1Velocities = [sqrt(2*gravity*windTunnelSpecificGravities(1).*abs(labData{1,1}{:,2}-labData{1,1}{:,3}).*sind(windTunnelAngles(1)) * (densityWater/densityAir) )]';
windTunnel2Velocities = [sqrt(2*gravity*windTunnelSpecificGravities(2).*abs(labData{1,1}{:,4}-labData{1,1}{:,5}).*sind(windTunnelAngles(2)) * (densityWater/densityAir) )]';
%Tunnel 1 is pressure
%Tunnel 2 is LD, which stands for  "Lift Drag" wind tunnel
%plotting actual air speed vs meter read air speed
dummyMeterReadAirSpeeds = linspace(min(labData{1,1}{:,1}),max(labData{1,1}{:,1}),resolution);
bflVelocites1 = polyfit(labData{1,1}{:,1},windTunnel1Velocities,1);
bflVelocites2 = polyfit(labData{1,1}{:,1},windTunnel2Velocities,1);
figure;
pointPlot = plot(labData{1,1}{:,1},windTunnel1Velocities,"o");
hold on;
plot(dummyMeterReadAirSpeeds,polyval(bflVelocites1,dummyMeterReadAirSpeeds),'color',get(pointPlot,'color'),'HandleVisibility','off');
hold on;
pointPlot = plot(labData{1,1}{:,1},windTunnel2Velocities,"o");
plot(dummyMeterReadAirSpeeds,polyval(bflVelocites2,dummyMeterReadAirSpeeds),'color',get(pointPlot,'color'),'HandleVisibility','off');
hold off;
legendStuff = [strcat("$Actual \ Velocity \ Wind \ Tunnel \ 1, \ V_{actual} = ",num2str(bflVelocites1(1)), " V_{meter} \ + \ " ,num2str(bflVelocites1(2)) ," $"),strcat("$Actual \ Velocity \ Wind \ Tunnel \ 2, \ V_{actual} = ",num2str(bflVelocites2(1)), " V_{meter} \ + \ " ,num2str(bflVelocites2(2)) ," $")];
lgnd= legend(legendStuff,'Location','best');
set(lgnd, 'Interpreter','latex')
%lgnd.Location = 'best';
title("Actual \ Wind \ Tunnel \ Velocities \ vs \ Meter \ Read \ Velocities",'Interpreter','latex')
xlabel("$ Meter \ Read \ Velocity \ (m/s) $",'Interpreter','latex') 
ylabel("$Actual \ Velocity \ (m/s)$",'Interpreter','latex') 

plot1 = "Actual Velocity vs Meter Read Velocity";
print('-r600','-dpng',plot1);


%%Calculating Drag on Shapes
%We just needed something to attach the wings to, but the attachment itself caused some drag, 
% which isn't what we were interested in. This is what F_D_offset is.
F_D_Sphere = [labData{1,2}{:,3} - labData{1,2}{:,2}]' ;
F_D_Hemisphere = [labData{1,2}{:,5} - labData{1,2}{:,4}]' ;
%the air speeds are 20-40mph, and need actual. So just using last 5 velocities of actual velocities 2 (since wind tunnel 2 is LD)
%LD stands for  "Lift Drag" wind tunnel, so using LD velocities for drag calculations
C_D_Sphere = F_D_Sphere./(densityAir*windTunnel2Velocities(3:7).^2*areaCrossSphere); 
C_D_Hemisphere = F_D_Hemisphere./(densityAir*windTunnel2Velocities(3:7).^2*areaCrossSphere);
%Calculating Re. Because both have the same diam, Re of sphere is same as Re of hemisphere 
Re = densityAir*sphereDiam/dynVisAir * windTunnel2Velocities(3:7);
figure;
plot(Re,C_D_Sphere,"o-");
hold on;
plot(Re,C_D_Hemisphere,"o-");
hold off;
legendStuff = ["$C_{D,Sphere}$","$C_{D,Hemiphere}$"];
lgnd= legend(legendStuff,'Location','best');
set(lgnd, 'Interpreter','latex')
%lgnd.Location = 'best';
title("Sphere \ \& \ Hemisphere \ Drag \ Coefficients",'Interpreter','latex')
xlabel("$ \Re $",'Interpreter','latex') 
ylabel("$C_{D}$",'Interpreter','latex') 

plot1 = "Sphere and Hemisphere Drag vs Re";
print('-r600','-dpng',plot1);

%%Calculating Lift and Drag on Airfoil
velocitiesLD = [windTunnel2Velocities(3),windTunnel2Velocities(5),windTunnel2Velocities(7)];
chordAirFoil = 3.976 * inTo_m;
spanAirfoil = 10.5*inTo_m;
areaAirfoil = chordAirFoil*spanAirfoil;


F_D_airfoil_closedFlaps_closedSlats = [];
F_D_airfoil_openFlaps_closedSlats = [];
F_D_airfoil_closedFlaps_openSlats = [];
C_D_airfoil_closedFlaps_closedSlats = [];
C_D_airfoil_openFlaps_closedSlats = [];
C_D_airfoil_closedFlaps_openSlats = [];
C_L_airfoil_closedFlaps_closedSlats = [];
C_L_airfoil_openFlaps_closedSlats = [];
C_L_airfoil_closedFlaps_openSlats = [];

x = 1:6;
x = 2*x+1; %now, get entries 3,5,7,9,11,13, which are F_D_total values. Don't need a loop :)

%We just needed something to attach the wings to, but the attachment itself caused some drag, 
% which isn't what we were interested in. This is what F_D_offset is.
F_D_airfoil_closedFlaps_closedSlats = labData{1,3}{1,1}{:,x} - labData{1,3}{1,1}{:,2}; %column 2 is offset drag
F_D_airfoil_openFlaps_closedSlats = labData{1,3}{1,2}{:,x} - labData{1,3}{1,1}{:,2}; %column 2 is offset drag
F_D_airfoil_closedFlaps_openSlats = labData{1,3}{1,3}{:,x} - labData{1,3}{1,1}{:,2}; %column 2 is offset drag

%F_L was measured directly. No offset.
F_L_airfoil_closedFlaps_closedSlats = labData{1,3}{1,1}{:,x+1};
F_L_airfoil_openFlaps_closedSlats = labData{1,3}{1,2}{:,x+1}; 
F_L_airfoil_closedFlaps_openSlats = labData{1,3}{1,3}{:,x+1};

for i = 1:6
    C_D_airfoil_closedFlaps_closedSlats(:,i) = (2/(densityAir*areaAirfoil)) * F_D_airfoil_closedFlaps_closedSlats(i)./(velocitiesLD.^2);
    C_D_airfoil_openFlaps_closedSlats(:,i) = (2/(densityAir*areaAirfoil)) * F_D_airfoil_openFlaps_closedSlats(i)./(velocitiesLD.^2);
    C_D_airfoil_closedFlaps_openSlats(:,i) = (2/(densityAir*areaAirfoil)) * F_D_airfoil_closedFlaps_openSlats(i)./(velocitiesLD.^2);
    C_L_airfoil_closedFlaps_closedSlats(:,i) = (2/(densityAir*areaAirfoil)) * F_L_airfoil_closedFlaps_closedSlats(i)./(velocitiesLD.^2);
    C_L_airfoil_openFlaps_closedSlats(:,i) = (2/(densityAir*areaAirfoil)) * F_L_airfoil_openFlaps_closedSlats(i)./(velocitiesLD.^2);
    C_L_airfoil_closedFlaps_openSlats(:,i) = (2/(densityAir*areaAirfoil)) * F_L_airfoil_closedFlaps_openSlats(i)./(velocitiesLD.^2);
end

anglesOfAttack = (0:5)*4;
%Just plotting C_D and C_L against angleOfAttack at 20mph (so just first
%row for all), for all air foil configurations

for i = 1:3
    figure;
    plot(anglesOfAttack,C_L_airfoil_closedFlaps_closedSlats(i,:))
    hold on
    plot(anglesOfAttack,C_L_airfoil_openFlaps_closedSlats(i,:))
    hold on
    plot(anglesOfAttack,C_L_airfoil_closedFlaps_openSlats(i,:))
    hold off
    legendStuff = ["$C_{L,airfoil_closedFlaps_openSlats}$","$C_{L,airfoil_openFlaps_closedSlats}$","$C_{L,airfoil_closedFlaps_openSlats}$"];
    lgnd= legend(legendStuff,'Location','best');
    set(lgnd, 'Interpreter','latex')
    %lgnd.Location = 'best';
    title( strcat("Lift Coefficients \ at \ ", num2str(10*i+10), "mph \ vs \ Angle \ of \ Attack"),'Interpreter','latex')
    xlabel("$ Angle \ Of \ Attack \ (deg^{\circ}) $",'Interpreter','latex') 
    ylabel("$C_{L}$",'Interpreter','latex')

    plot1 = strcat("Lift on Airfoil vs Angles of Attack ",num2str(10*i+10), "mph" );
    print('-r600','-dpng',plot1);

    figure;
    plot(anglesOfAttack,C_D_airfoil_closedFlaps_closedSlats(i,:))
    hold on
    plot(anglesOfAttack,C_D_airfoil_openFlaps_closedSlats(i,:))
    hold on
    plot(anglesOfAttack,C_D_airfoil_closedFlaps_openSlats(i,:))
    hold off
    legendStuff = ["$C_{D,airfoil_closedFlaps_openSlats}$","$C_{D,airfoil_openFlaps_closedSlats}$","$C_{D,airfoil_closedFlaps_openSlats}$"];
    lgnd= legend(legendStuff,'Location','best');
    set(lgnd, 'Interpreter','latex')
    %lgnd.Location = 'best';
    title(strcat("Drag Coefficients \ at \ ", num2str(10*i+10), "mph \ vs \ Angle \ of \ Attack"),'Interpreter','latex')
    xlabel("$ Angle \ Of \ Attack \ (deg^{\circ}) $",'Interpreter','latex') 
    ylabel("$C_{D}$",'Interpreter','latex') 

    plot1 = strcat("Drag on Airfoil vs Angles of Attack",num2str(10*i+10), "mph");
    print('-r600','-dpng',plot1);

end

plot1 = "Drag and Lift on Airfoil vs Angles of Attack";
print('-r600','-dpng',plot1);

%%Finding Pressure Coefficient on Airfoil
velocitiesP = [windTunnel1Velocities(3),windTunnel1Velocities(5),windTunnel1Velocities(7)];
anglesOfAttackPressure = [0,10,-10]; %in deg
cpAirfoils = {};
for i = 1:3
    for ii = 1:3
        refHeights = ones(height(labData{1,4})-1,1);
        refHeights = refHeights*labData{1,4}{end,3*ii+i-1};
        cpAirfoils{1,i}(:,ii) = ( 2*gravity*windTunnelSpecificGravities(1) * densityWater * (labData{1,4}{1:end-1,3*ii+i-1}-refHeights)...
            *sind(anglesOfAttackPressure(i)) )/(densityAir*velocitiesP(i)^2);             %starts on column 3,6,9 so i-1=0,  3*ii = 3,6,9
    end
     %for 0 degrees, C_P is all 0, since sind(0) = 0
end

%plotting at 20mph for each angle of attack, so first column of each

for i = 1:3
    figure;
    for ii = 1:3
        plot(labData{1,4}{1:end-1,2},cpAirfoils{1,ii}(:,i))
        hold on;
    end
    hold off;
    legendStuff = ["$C_{P,0^{/circ}}$","$C_{P,10^{/circ}}$","$C_{P,-10^{/circ}}$"];
    lgnd= legend(legendStuff,'Location','best');
    set(lgnd, 'Interpreter','latex')
    %lgnd.Location = 'best';
    title(strcat("Pressure Coefficients \ at \ ", num2str(10*i+10), "mph \ vs \ Axial \ Position"),'Interpreter','latex')
    xlabel("$ Axial \ Position (m) $",'Interpreter','latex') 
    ylabel("$C_{P}$",'Interpreter','latex') 

    plot1 = strcat("Airfoil Pressure Coefficient vs Axial Position ", num2str(10*i+10), "mph");
    print('-r600','-dpng',plot1);
end



%%Plotting Pressure Distribution Over Circular Cylinder
anglesOfAttackCylinder = 90;%deg
velocitiesPCylinder = [20,26]*miPerHr_to_m_s;

cpCylinder = {};
for i = 1:2
    for ii = 1:2
        refHeights = ones(height(labData{1,5})-1,1);
        refHeights = refHeights*labData{1,5}{end,1 + 2*ii+i-1};
        cpCylinder{1,i}(:,ii) = ( 2*gravity*windTunnelSpecificGravities(1) * densityWater * (labData{1,5}{1:end-1,1 + 2*ii+i-1}-refHeights)...
            *sind(anglesOfAttackCylinder) )/(densityAir*velocitiesPCylinder(i)^2);             
    end
     %for 0 degrees, C_P is all 0, since sind(0) = 0
end

%Plotting C_p vs angular position
figure;
for i =1:2
    for ii = 1:2
        plot(labData{1,5}{1:end-1,2},cpCylinder{1,i}(:,ii))
        hold on;
    end
end

cpCylinder_theo = 1 - 4*sind(labData{1,5}{1:end-1,2}).^2;
plot(labData{1,5}{1:end-1,2},cpCylinder_theo)
hold off;

legendStuff = ["$C_{P,Right,20mph,90^{\circ}angleofattack}$","$C_{P,Right,26mph,90^{\circ}angleofattack}$","$C_{P,Left,20mph,90^{\circ}angleofattack}$","$C_{P,Left,26mph,90^{\circ}angleofattack}$","$C_{P,Theo,90^{\circ}angleofattack}$"];
lgnd= legend(legendStuff,'Location','best');
set(lgnd, 'Interpreter','latex')
%lgnd.Location = 'best';
title("Pressure \ Coefficients \ On \ Circular \ Cylinder \ vs \ Angular \ Position",'Interpreter','latex')
xlabel("$ Angular \ Position, \ \theta (deg^{\circ}) $",'Interpreter','latex') 
ylabel("$C_{P}$",'Interpreter','latex') 

plot1 = "Cylinder Pressure Coefficient vs Angular Position";
print('-r600','-dpng',plot1);