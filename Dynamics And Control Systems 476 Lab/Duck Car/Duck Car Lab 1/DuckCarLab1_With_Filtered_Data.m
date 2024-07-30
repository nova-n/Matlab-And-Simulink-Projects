% Duck Car Lab 1
% Ian Bautista, James Varghese, Nathan Delos Santos

clear;
clc;
close all;

ThreeV_PosData = readtable("3V_Ducky_Position.txt");
FourV_PosData = readtable("4v_Ducky_position.txt");

ThreeVData = table2array([readtable("3V_Ducky_Velocity.txt"),ThreeV_PosData(:,2)]); %both pos and vel have same timestamps
FourVData = table2array([readtable("4v_Ducky_velocity.txt"), FourV_PosData(:,2) ]);
FiveVData = table2array(readtable("Duck Car5V.txt"));

clearvars ThreeV_PosData FourV_PosData; %clears them from memory

%movmean(ThreeVData{:,:},3,1)


filteredThreeVData = dataNoiseRemover(ThreeVData,10);
filteredThreeVData = dataSmoother(filteredThreeVData,50,10);
filteredFourVData = dataNoiseRemover(FourVData,10);
filteredFourVData = dataSmoother(filteredFourVData,20,3);
filteredFiveVData = dataNoiseRemover(FiveVData,10);
filteredFiveVData = dataSmoother(filteredFiveVData,150,3);


rawAndFilteredPlotter(ThreeVData,filteredThreeVData,"$3 V \ Velocity$","$Velocity \ ( \frac{m}{s} )$")
rawAndFilteredPlotter(FourVData,filteredFourVData,"$4 V \ Velocity$","$Velocity \ ( \frac{m}{s} )$")
rawAndFilteredPlotter(FiveVData,filteredFiveVData,"$5 V \ Velocity$","$Velocity \ ( \frac{m}{s} )$")

% %must ALWAYS skip first and last velocity entry, since are NaN, and so if want to
% %plot it, also omit those entries for position

function [filteredData] = dataNoiseRemover(dataIn,cutoffFreq)
    dataIn = NaNRemover(dataIn);
    samplesCount = height(dataIn);
    filteredData = zeros(height(dataIn),width(dataIn));
    filteredData(:,1) = dataIn(:,1);%all time entries for position and velocity are the same
        for i = 2:width(dataIn)
        filteredData(:,i) = lowpass(dataIn(:,i),cutoffFreq,samplesCount);
    end
end

function [smoothedOutData] = dataSmoother(dataIn,sandingPasses,movingAvgWindow) 
    %think of the movingAvgWindow as the grit of sandpaper
    smoothedOutData = dataIn; %No way to pass by reference :(
    for i = 1:sandingPasses
        smoothedOutData = [smoothedOutData(:,1),movmean(smoothedOutData(:,2)...
            ,movingAvgWindow) , movmean(smoothedOutData(:,3),movingAvgWindow)];
    end
end

function [dataWithoutNaNs] = NaNRemover(dataInput)
    dataWithoutNaNs = dataInput;
    [a,~] = find(isnan(dataWithoutNaNs));
    a = unique(a(:).'); %removes any duplicate values, incase there is a NaN twice in one row
    for i = 1:length(a)
        dataWithoutNaNs(a(i),:) = [];
        a = a - 1;
    end
end

function rawAndFilteredPlotter(originalData,filteredData,givenTitle,measuredQuantity)
    figure;
    plot(originalData(:,1),originalData(:,2))
    hold on;
    plot(filteredData(:,1),filteredData(:,2))
    hold off;
    title(givenTitle,'Interpreter','latex')
    xlabel("$ Time, \ t \ (s) $",'Interpreter','latex') 
    ylabel(measuredQuantity,'Interpreter','latex')
    legendStuff = [givenTitle, strcat("$Filtered \ And \ Smoothed \ " , givenTitle, " \ $")];
    lgnd= legend(legendStuff,'Location','best');
    set(lgnd, 'Interpreter','latex')
end