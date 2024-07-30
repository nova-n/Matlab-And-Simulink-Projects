clc;
clear;
clear all;
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

%test subject specifications
widthSubject = 0.125;%in m
heightSubject = 0.125;%in m
diam = 0.0125;%in m
lengthSubject = 0.0951;%in m
areaSurface = pi * diam * lengthSubject;
mass = 0.11295;%in kg
c_specHeat = 380;%in J/(kg*k)
test_Conductivity = 398;%in W/(m*k)

%air properties
air_Conductivity = 0.0259;%in W/(m*k)
densityAir = 1.204;%kg/m^3 at room temp!!!
dynVisAir = 1.81*10^-5;%at room temp!


upstreamRange = 'A6:M10';
%will just have to delete bottom row after
fileNameRange = 'B11:M11';
labData(1,1) = {readtable("Ex05_Data\EGME306B_Ex05_DataSheet_02.xlsx",'Sheet',"Upstream",'Range',upstreamRange)};
[~ , fileNamesList]  = xlsread("Ex05_Data\EGME306B_Ex05_DataSheet_02.xlsx",fileNameRange);

for i = 1:width(fileNamesList)
    excellFileName = strcat("Ex05_Data\",fileNamesList{1,i}, ".xlsx");
    labData(1,1+i) = { readtable(excellFileName) };
end


for i = 2:length(labData) %lab data for the tests starts on index 2

    %%Cooling Curves
    %ALL TEMPS ARE IN CELCIUS
    T_SR_test = labData{1,i}{:,3};%is thermocouple #2, on test subject
    T_SR_airinlet = labData{1,i}{:,4};%is thermocouple #3, on air inlet
    T_SR_airinlet_avg = mean(T_SR_airinlet);%is average of inlet air temp
    T_SR_i = T_SR_test(1);%initial test subject temp
    
    t_SR = labData{1,i}{:,1};%are timestamps of test
    
    dimensionlessTemp = log( (T_SR_test - T_SR_airinlet) / (T_SR_i - T_SR_airinlet_avg) ); %in matlab, log is actually ln
    
    fs = 1e3;
    filteredTestTemp = lowpass(dimensionlessTemp,10,fs);
    
    figure;
    plot(t_SR,dimensionlessTemp);
    hold on
    plot(t_SR,filteredTestTemp);
    hold on
    
    %to get best fit line, finding middle third of filtered line
    %rounds to the lowest number divisible by 3 by subtracting remainder
    roundedIndeciesCount = length(T_SR_test) - mod(length(T_SR_test),3);
    boundingIndecies = [roundedIndeciesCount/3 , 2*roundedIndeciesCount/3];
    
    approxLineEqn = polyfit( t_SR( boundingIndecies(1):boundingIndecies(2) ) , ...
        dimensionlessTemp( boundingIndecies(1):boundingIndecies(2) ) , 1 );
    linEqns(i-1,:) = approxLineEqn;
    t_FinerMesh = linspace(0,t_SR(end),resolution); 
    bestFitTemp = approxLineEqn(1)*t_FinerMesh+approxLineEqn(2);
    plot(t_FinerMesh,bestFitTemp,'--');
    
    legendStuff = ["$Actual \ Test \ Data, \theta (t) $", ...
        "$ Noise \ Filtered \ Test \ Data$","$ Linear \ Approximation \ of \ Test \ Data $"];
    lgnd= legend(legendStuff,'Location','best');
    set(lgnd, 'Interpreter','latex')
    %lgnd.Location = 'best';

    %Biot Number Verification for Lumped Capacitance Assumption
    %in the h*A_s, this is actually equal to the slope, so approxLineEqn(1)
    %will be solving for h to get the biot number
    h(i-1) = -1*approxLineEqn(1)*mass*c_specHeat/areaSurface;
    Bi(i-1) = h(i-1)*diam/(4*test_Conductivity);
    if Bi(i-1) < 0.1
        validity = "Valid";
        signChar = " \leq 0.1 ";
    else
        validity = "Invalid";
        signChar = " > 0.1 ";
    end

    %title(strcat("Dimensionless \ Temperature \ of \ Test \ ", fileNamesList{1,i-1} ," \ vs \ time"),'Interpreter','latex')
    title(strcat("$Test \ ", fileNamesList{1,i-1} ,", \ Bi \ = \ ", num2str(Bi(i-1)), signChar, ...
        ", \ Lumped \ Capacitance \ is \ ", validity , "$"),'Interpreter','latex')
    xlabel("$ Time, \ t \ (s) $",'Interpreter','latex') 
    ylabel("$Dimensionless \ Temperature: \theta (t) =  \ ln( \frac{ T(t) - T_{ \infty } }{ T_i - T_{ \infty ,avg } } ) $",...
        'Interpreter','latex') 
    
    plot1 = strcat("Dimensionless Temperature of Test ", fileNamesList{1,i-1} ,  " vs time");
    print('-r600','-dpng',plot1);
    hold off
end

Re = [];
figure;
for i = 1:3 %lab data for the tests starts on index 2
    %%Calculating Nusselts Number and Reynolds Number
    Nu(i,:) = h( 4*(i-1)+1:4*i ) *diam/air_Conductivity;
    
    airVelocities = [];
    if fileNamesList{1,i*4}(end-3:end-2) == "SR"
        airVelocities = sqrt(2*gravity*(densityWater/densityAir)*( abs(labData{1,1}{2:end,1:4}-labData{1,1}{1,1:4})*mmtom ));
    elseif fileNamesList{1,i*4}(end-3:end-2) == "FB"
        airVelocities = sqrt(2*gravity*(densityWater/densityAir)*( abs(labData{1,1}{2:end,5:8}-labData{1,1}{1,5:8})*mmtom ));
    else
        airVelocities = sqrt(2*gravity*(densityWater/densityAir)*( abs(labData{1,1}{2:end,9:12}-labData{1,1}{1,9:12})*mmtom ));
    end
    avgAirVelocity = mean(airVelocities,2);%want the average of each column (average each percentage opening, vs distance)

    %airVelocities = sqrt(2*gravity*(densityWater/densityAir)*( abs(labData{1,1}{2:end,i}-labData{1,1}{1,i})*mmtom ));
    %avgAirVelocity = mean(airVelocities);
    if fileNamesList{1,i*4}(end-3:end-2) == "SR"
        Re(i,:) = (densityAir*avgAirVelocity*diam)/(dynVisAir*(1 - diam/heightSubject) );
    else
        Re(i,:) = (densityAir*avgAirVelocity*diam)/(dynVisAir*(1 - 5*diam/heightSubject) );
    end
    plot(Re(i,:),Nu(i,:),"o-")
    hold on;
end
hold off;
title("$ Nu vs \Re $",'Interpreter','latex')
xlabel("$ \Re $",'Interpreter','latex') 
ylabel("$ Nu $",'Interpreter','latex') 
legendStuff = ["$SR$","$FB$","$FR$"];
lgnd= legend(legendStuff,'Location','best');
set(lgnd, 'Interpreter','latex')
 plot1 = 'Nu vs Re';
 print('-r600','-dpng',plot1);

