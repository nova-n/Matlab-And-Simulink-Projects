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

filteredThreeVData = dataNoiseRemover(ThreeVData,10);
filteredThreeVData = dataSmoother(filteredThreeVData,50,10);
filteredFourVData = dataNoiseRemover(FourVData,10);
filteredFourVData = dataSmoother(filteredFourVData,20,3);
filteredFiveVData = dataNoiseRemover(FiveVData,10);
filteredFiveVData = dataSmoother(filteredFiveVData,150,3);

ThreeV_fv = finalValueFinder(filteredThreeVData,1*10^(-2.5));
FourV_fv = finalValueFinder(filteredFourVData,1*10^(-2));
FiveV_fv = finalValueFinder(filteredFiveVData,1*10^(-2.5));

rawAndFilteredPlotter(ThreeVData,filteredThreeVData,"$3 V \ Velocity$","$Velocity \ ( \frac{m}{s} )$",ThreeV_fv)
rawAndFilteredPlotter(FourVData,filteredFourVData,"$4 V \ Velocity$","$Velocity \ ( \frac{m}{s} )$",FourV_fv)
rawAndFilteredPlotter(FiveVData,filteredFiveVData,"$5 V \ Velocity$","$Velocity \ ( \frac{m}{s} )$",FiveV_fv)

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
    dataWithoutNaNs = dataInput; % can't pass by reference :(
    [a,~] = find(isnan(dataWithoutNaNs));
    a = unique(a(:).'); %removes any duplicate values, incase there is a NaN twice in one row
    for i = 1:length(a)
        dataWithoutNaNs(a(i),:) = [];
        a = a - 1;
    end
end

function [finalValue] = finalValueFinder(dataInputted,steadyStateTolerance)
    %will only look at last third of time values for final value
    lastThirdOfTime = floor(2*height(dataInputted)/3);
    a = diff(dataInputted(lastThirdOfTime:end,2))./diff(dataInputted(lastThirdOfTime:end,1));
    steadyestLastIndecies = find(abs(a)<=steadyStateTolerance) + lastThirdOfTime - 1; 
    %checks for any changes smaller than the tolerance, so that change is close to 0
    consecutiveSteadyIndecies = [;]; %  [ start1,start2... ; end1,end2... ]
    itt = 1;
    %will now check for the longest streak of consecutive indecies that fall within tolerance
    while itt <= length(steadyestLastIndecies)
        a = steadyestLastIndecies(itt);
        iitt = itt+1; %check for consecutive
        while iitt <= length(steadyestLastIndecies) && steadyestLastIndecies(iitt) - steadyestLastIndecies(itt) < 2
            iitt = iitt+1;
            itt = itt + 1;
        end
        b = steadyestLastIndecies(itt);
        consecutiveSteadyIndecies = [consecutiveSteadyIndecies ; [a,b]];
        itt = itt+1;
    end
    consecutiveSteadyIndecies;
    longestStreakPairs = find(max( consecutiveSteadyIndecies(:,2) - consecutiveSteadyIndecies(:,1) )); %might return multiple indecies if there is a tie
    longestStreakPairs = longestStreakPairs(end); %use the last longest steady streak
    finalValue = mean(dataInputted(consecutiveSteadyIndecies(longestStreakPairs,1):consecutiveSteadyIndecies(longestStreakPairs,1),2));
end

function rawAndFilteredPlotter(originalData,filteredData,givenTitle,measuredQuantity,finalVal)
    figure;
    plot(originalData(:,1),originalData(:,2))
    hold on;
    plot(filteredData(:,1),filteredData(:,2))
    hold on;
    yline(finalVal,"--b");
    hold off;
    title(givenTitle,'Interpreter','latex')
    xlabel("$ Time, \ t \ (s) $",'Interpreter','latex') 
    ylabel(measuredQuantity,'Interpreter','latex')
    legendStuff = [givenTitle, strcat("$Filtered \ And \ Smoothed \ " , givenTitle, " \ $"),strcat("$Final \ Value \ " , givenTitle, " \ $")];
    lgnd= legend(legendStuff,'Location','best');
    set(lgnd, 'Interpreter','latex')
end