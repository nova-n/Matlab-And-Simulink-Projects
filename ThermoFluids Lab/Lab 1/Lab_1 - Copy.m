clc;
clear;
close all;
densityWater = 1.94;%slugs/ft^3s
gravity = 32.2;%ft/s^2
stationB_Data = {readtable("EGME306B_Ex01_DataSheet_08B.xlsx",'Sheet',"Heights"),readtable("EGME306B_Ex01_DataSheet_08B.xlsx",'Sheet',"Weights")};
stationB_Data{1}{:,[2:end]} = stationB_Data{1}{:,[2:end]}/12; %converting in to ft \
actualThroatDiamB = stationB_Data{1}{14,3}/12;%in ft
waterTempB = stationB_Data{1}{15,3};%in F
stationB_Data{1}([12:end],:) = []; %deleted unused spaces
stationB_Data{1}(:,10) = []; %deleting comments column

stationB_Data{2}([1,7],:) = []; %deleted last and first rows, since didn;t use that time
stationB_Data{2}(:,9) = []; %deleting comments column
%replacing comments column with volume filled
%vol = weight/(density * gravity)
bucketVolumes = [stationB_Data{2}{:,2}/(densityWater*gravity)]';
%stationB_Data{2}{:,3}; %time

%will need entire array of these, for each run
% [row of Initial Pressures, row of flowrates]
systemVarsB = [;] ;
width(stationB_Data{2}{:,[3:end]}); %number of tests, since tests start on column 3

for i = 1:width(stationB_Data{2}{:,[3:end]})
    systemVarsB(1,i) = densityWater*gravity*stationB_Data{1}{1,i+3}; %heights start on col 4, and i starts at 1. Am using initial height, since is same diam as inlet
    fillFlowEqn = polyfit(stationB_Data{2}{:,i+2},bucketVolumes,1);%fill times start on col 3, and i starts at 1
    Q_volflow = fillFlowEqn(1); %in ft^3/s
    systemVarsB(2,i) = Q_volflow;
end
areasB = pi* stationB_Data{1}{:,2}.^2 /4;

%h_theoreticalB: columns are each run, and height is each row in each column
%h_theoreticalB = (systemVarsB(1,:) + systemVarsB(2,:).^2 *(densityWater/2).*(1./stationB_Data{1}{1,2} - 1./stationB_Data{1}{:,2}))./(densityWater*gravity);
h_theoreticalB = [;] ;
for i = 1:length(systemVarsB)
    h_theoreticalB(:,i) = (systemVarsB(1,i) + systemVarsB(2,i).^2 *(densityWater/2).*(1./areasB(1)^2 - 1./areasB(:).^2))./(densityWater*gravity);
end

figure;
for i = 1:width(stationB_Data{2}{:,[3:end]})
    theoreticlPlot = plot(stationB_Data{1}{:,3},h_theoreticalB(:,i),"--",'HandleVisibility','off'); %only set equal so I can get the color
    hold on;
    plot(stationB_Data{1}{:,3},stationB_Data{1}{:,i+3},"o",'color',get(theoreticlPlot,'color'))
end
legendStuff = {'$ 10in \ H_{2}0 $','$ 8in \ H_{2}0 $','$ 6in \ H_{2}0 $','$ 4in \ H_{2}0 $','$ 2in \ H_{2}0 $','$ 0.6in \ H_{2}0 $'};
lgnd= legend(legendStuff);
set(lgnd, 'Interpreter','latex')
title("$ Problem \ 1B \ Section \ 08B $",'Interpreter','latex')
xlabel("$ Axial \ Distance \ (ft) $",'Interpreter','latex') 
ylabel("$Manometer \ Height \ (ft)$",'Interpreter','latex') 
hold off;