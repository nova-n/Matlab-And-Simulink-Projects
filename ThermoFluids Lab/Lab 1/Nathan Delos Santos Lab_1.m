% --------------------------------------------------------
% Experiment No. 01 - Data reduction
%
% Input files :
% Experiment No. 01 rawdata : 
%
% Output files :
% Processed results of Exp. No 01: exp01_processed.xlsx
% Plots of pressure distribution
% Plot of discharge coefficient vs Reynolds number
% Plot of head lossflow rate
%
% Code written by : Nathan Delos Santos
% Last update : 26-FEB -2023 by Nathan Delos Santos
% --------------------------------------------------------
% INITIALIZE CODE
% --------------------------------------------------------
% Clear memory
clc;
clear;
close all;
graphResolution = 100;
densityWater = 998;%kg/m^3
gravity = 9.81;%m/s^2
dynViscosity = 0.0009545; % Pa*s
allData = { {readtable("EGME306B_Ex01_DataSheet_08A.xlsx",'Sheet',"Heights"),readtable("EGME306B_Ex01_DataSheet_08A.xlsx",'Sheet',"Weights")} , {readtable("EGME306B_Ex01_DataSheet_08B.xlsx",'Sheet',"Heights"),readtable("EGME306B_Ex01_DataSheet_08B.xlsx",'Sheet',"Weights")}};

for I = 1:length(allData) 
    station_Data =allData{I};
    [station_Data, testsCount, waterTemp, actualThroatDiam, diams, throatDiam, axialDistances,area,actualHeights] = dataProcessing(station_Data);
    runs = [[1:3];[4:testsCount]];

    heightDiffs = [10,8,6,4,2,0.6]*0.0254;
    Q_theo = pi/4 * diams(1)^2 * throatDiam^2 * sqrt( (2*gravity .* (heightDiffs)) / (diams(1)^4 - throatDiam^4) );
    %h_theo = [h1_Q1,h2_Q1,h3_Q1 ; h1_Q2,h2_Q2,h3_Q2 ; ...]
    h_theo = ones(testsCount,height(station_Data{1})); %preallocating space this time
    % 1./( (stationB_Data{1}{4,2} * (ones(height(stationB_Data{1}),1)) ).^4) - 1./( ( stationB_Data{1}{:,2} ).^4) 
    for i = 1:testsCount
        %ones(height(stationB_Data{1}),1) * 1/stationB_Data{1}{4,2}^2 just makes an array of throat diam,
        %so I can subtract it from another array.
        h_theo(i,:) = 8/(gravity * pi^2) * Q_theo(i)^2 * ( 1./( (throatDiam * ones( height(station_Data{1}),1) ).^4) - 1./(diams.^4) );
    end
    
    %vol = weight/(density * gravity)
    bucketVolumes = [station_Data{2}{:,2}/(densityWater*gravity)]';
    for i = 1:testsCount
        fillFlowEqn = polyfit(station_Data{2}{:,i+2},bucketVolumes,1);%fill times start on col 3, and i starts at 1
        Q_actual(i) = fillFlowEqn(1); %in m^3/s
    end%
    % stationB_Data{2}{:,3}; %time
    
    runsToPlot = runs(I,:);
    fig1 = figure;
    for i = runsToPlot(1):runsToPlot(end)
        theoreticalPlot = plot(axialDistances,h_theo(i,:));
        hold on;
        %have to normalize the plot of actual heights
        plot(axialDistances, actualHeights(:,i) - min(actualHeights(:,i)) ,"o",'color',get(theoreticalPlot,'color'),'HandleVisibility','off')
    end
    legendStuff = [{'$ 10in \ H_{2}0 $','$ 8in \ H_{2}0 $','$ 6in \ H_{2}0 $'};{'$ 4in \ H_{2}0 $','$ 2in \ H_{2}0 $','$ 0.6in \ H_{2}0 $'}];
    lgnd= legend(legendStuff(I,:));
    set(lgnd, 'Interpreter','latex')
    titles = ["$ Problem \ 1A \ Section \ 08A $","$ Problem \ 1B \ Section \ 08B $"];
    title(titles(I),'Interpreter','latex')
    xlabel("$ Axial \ Distance \ (m) $",'Interpreter','latex') 
    ylabel("$Manometer \ Height \ (m)$",'Interpreter','latex') 
    hold off;
    manometerPlotNames = ["ManometersGroupA","ManometersGroupB"];
    plot1 = manometerPlotNames(I)
    print('-r600','-dpng',plot1);
    
    dischargeCoefficients = Q_actual./Q_theo;
    reynoldsNumbers = (4*densityWater)/(pi * dynViscosity * throatDiam) * Q_actual;
    reynoldsAxis = linspace(0,1.1*max(reynoldsNumbers),graphResolution);%want the graph to extend a bit, so thats why is 1.1*max
    figure;
    theoreticalPlot = plot(reynoldsNumbers,dischargeCoefficients,"o");
    hold on;
    bestFit = polyfit(reynoldsNumbers,dischargeCoefficients,3);
    plot(reynoldsAxis,polyval(bestFit,reynoldsAxis),'color',get(theoreticalPlot,'color'))
    hold off;
    xlim([0,max(reynoldsAxis)])
    legendStuff = {'$Experimental \ Data$','$Theoretical \ Curve$'};
    lgnd= legend(legendStuff);
    set(lgnd, 'Interpreter','latex')
    titles = ["$ Problem \ 2A \ Section \ 08A $","$ Problem \ 2B \ Section \ 08B $"];
    title(titles(I),'Interpreter','latex')
    xlabel("$ Reynolds \ Number, \Re $",'Interpreter','latex') 
    ylabel("$Discharge \ Coefficient, \ C_{d}$",'Interpreter','latex') 
    hold off;
    manometerPlotNames = ["Reynolds_VS_Discharge_GroupA","Reynolds_VS_Discharge_GroupB"];
    plot1 = manometerPlotNames(I)
    print('-r600','-dpng',plot1);
    
    %beginning heights - end heights
    headLossAlongPipe = actualHeights(1,:) - actualHeights(end,:);
    figure;
    theoreticalPlot = plot(Q_actual,headLossAlongPipe,"o");
    hold on;
    bestFit = polyfit(Q_actual,headLossAlongPipe,2);
    Q_axis = linspace(0,1.1*max(Q_actual),100); %want the graph to extend a bit, so thats why is 1.1*max
    plot(Q_axis,polyval(bestFit,Q_axis),'color',get(theoreticalPlot,'color'))
    hold off;
    xlim([0,max(Q_axis)])
    legendStuff = {'$Experimental \ Data$','$Theoretical \ Curve$'};
    lgnd= legend(legendStuff);
    set(lgnd, 'Interpreter','latex')
    titles = ["$ Problem \ 3A \ Section \ 08A $","$ Problem \ 3B \ Section \ 08B $"];
    title(titles(I),'Interpreter','latex')
    xlabel("$ Actual \ Volumetric \ Flowrate, \ Q_{actual} \ ( \frac{m^3}{s} ) $",'Interpreter','latex') 
    ylabel("$ Head \ Loss, \ h_{loss} \ (m) $",'Interpreter','latex') 
    hold off;
    manometerPlotNames = ["Head_Loss_GroupA","Head_Loss_GroupB"];
    plot1 = manometerPlotNames(I)
    print('-r600','-dpng',plot1);

end 

function [cleanedStationData, runsCount, waterTemperature, actualThroatDiameter,diamArray, listedThroatDiam,axialDists,crossSectionAreas,experimentalHeights] = dataProcessing(inputTable)
    cleanedData = inputTable;
    waterTemperature = (cleanedData{1}{15,3} - 32 ) * 5/9;%in C
    cleanedData{1}{:,[2:end]} = cleanedData{1}{:,[2:end]}*0.0254; %converting in to m
    actualThroatDiameter = cleanedData{1}{14,3}*0.0254;%in m
    cleanedData{1}([12:end],:) = []; %deleted unused spaces
    cleanedData{1}(:,10) = []; %deleting comments column
    cleanedData{2}([1,7],:) = []; %deleted last and first rows, since didn;t use that time
    cleanedData{2}(:,9) = []; %deleting comments column
    cleanedData{2}{:,2} = cleanedData{2}{:,2} * 4.4482189159; %converting to N 
    %replacing comments column with volume filled\
    diamArray = cleanedData{1}{:,2};
    listedThroatDiam = diamArray(4);
    cleanedData{1}{:,3} = cleanedData{1}{:,3} - cleanedData{1}{1,3}; %normalizing axial distances
    axialDists = cleanedData{1}{:,3};
    crossSectionAreas = pi* cleanedData{1}{:,2}.^2 /4;
    experimentalHeights = cleanedData{1}{:,4:end};
    runsCount = width(cleanedData{2}{:,[3:end]}); %number of tests, since tests start on column 3
    cleanedStationData = cleanedData;
end
