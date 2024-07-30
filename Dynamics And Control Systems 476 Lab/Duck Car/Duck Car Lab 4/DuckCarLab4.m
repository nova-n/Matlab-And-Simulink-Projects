clear;
clc;
close all;

%%Experiment Constants
standOffDistInch = 7.00;
standoffV = 1.335;
finalDistInch = 13.00;

K_a = 1.2;
K_m = mean([0.3297,0.3789,0.44694]);
K_m_inch = K_m/0.0254;
K_s = -0.180720588235294;
%K_s = -0.0568;

nonLinearSensor = [0.0864192841833102 , 0.140817267512297];
nonLinearSensor_metric = [3.4023,5.5440];

%%Importing Data
labGainData = table2array(readtable("Lab4_slomo_vel_pos.txt"));

K_p = 7.8;
K_d = (78*10^3)*(13.3*10^-6)

labGainData = NaNRemover(labGainData);

%dampGainData(:,1:2) = dataNoiseRemover(dampGainData(:,1:2),50);
%just removes spikes, thats all. Thats why have a low number for the second argument
labGainData(:,1:2) = medfilt1(labGainData(:,1:2),10); 

posFinalValue = rawDataFinalValueFinder(labGainData(:,3),[0.75,1]);
velFinalValue = round(rawDataFinalValueFinder(labGainData(:,2),[0.75,0.8]) ,2);%is practically 0 in set, and should be 0 in real life

%%finding time constants, settling time, peak time, etc...
[tau,Ts] = TwoPercentSettlingTimeTau(labGainData,3,posFinalValue,10^-2,0.03)

[posMax,posPeakTime,posMaxOvershoot,posMaxPO] = overshootFinder(labGainData,3,posFinalValue)

dampingValue = -log(posMaxPO/100) / sqrt( pi^2 + (log(posMaxPO/100))^2 )
natFreq = pi / ( posPeakTime * sqrt( 1 - dampingValue^2 ) )

%%plotting
sameSetDataPlotter(labGainData,[velFinalValue,posFinalValue],strcat("$ \ K_{p} = " ,...
    num2str(K_p)," , " , " K_{d} = " ,num2str(K_d), ", \ Raw  \ Data \ $"), "$ Time \ (s) $" ,...
    ["$Velocity \ v(t) $", strcat("$Final \ Velocity \ v_{ss} = " ,num2str(velFinalValue) ," \frac{m}{s} $") ,...
    "$Positon \ x(t)$" ,strcat("$Final \ Position \ x_{ss} = " ,num2str(posFinalValue) ," m $") ])


%%functions
function [dataWithoutNaNs] = NaNRemover(dataInput)
    dataWithoutNaNs = dataInput; % can't pass by reference :(
    [a,~] = find(isnan(dataWithoutNaNs));
    a = unique(a(:).'); %removes any duplicate values, incase there is a NaN twice in one row
    for i = 1:length(a)
        dataWithoutNaNs(a(i),:) = [];
        a = a - 1;
    end
end

function [filteredData] = dataNoiseRemover(dataIn,cutoffFreq)
    %Runs a low pass filter for each non-time column of the data, to filter
    %out nise
    dataIn = NaNRemover(dataIn);
    samplesCount = height(dataIn);
    filteredData = zeros(height(dataIn),width(dataIn));
    filteredData(:,1) = dataIn(:,1);%all time entries for position and velocity are the same
    for i = 2:width(dataIn)
        filteredData(:,i) = lowpass(dataIn(:,i),cutoffFreq,samplesCount);
    end
end

function [rawDataFinalValue] = rawDataFinalValueFinder(dataInputted,timeSection)
    %will only look at last set of time values for final value
    lastSectionOfTime = [floor( height(dataInputted)*timeSection(1)) , ceil( height(dataInputted)*timeSection(2))];
    rawDataFinalValue = mean(dataInputted(lastSectionOfTime(1):lastSectionOfTime(2)));
    %runs these indecies through the filtered function, and averages the values at these indecies
end

function sameSetDataPlotter(datas,finalValues,titleStuff,xAxisLabel,legendStuff)
    figure;
    for i = 2:width(datas) %final values length is the same as datas width
        pointPlot = plot(datas(:,1),datas(:,i));
        hold on;
        %just gets it the same color
        yline(finalValues(i-1),"--",'color',get(pointPlot,'color'))
        hold on;
    end
    
    hold off
    title(titleStuff,'Interpreter','latex')
    xlabel(xAxisLabel,'Interpreter','latex') 
    lgnd= legend(legendStuff,'Location','best');
    set(lgnd, 'Interpreter','latex')
end

function [settlingTimeVersion_Tau,settlingTime] = TwoPercentSettlingTimeTau(dataIn,dataColumn,finalVal,tol,streakFraction)
    minContinuousLength = floor( streakFraction * length(dataIn) ); %has to be at least spanning 3% of the graph
    %minContinuousLength = floor( 0.03 * length(dataIn) ); %has to be at least spanning 3% of the graph
    
    % will find the first continous set of values within 2% that is at
    % least "minContinuousLength" indecies long 

    chosenIndecies = find(abs(finalVal*0.982 - dataIn(:,dataColumn)) <= tol); %if the data is 98% of the final val within a tolerance

    consecutiveToleranceIndecies = [;]; %  [ start1,start2 ; end1,end2 ; ... ]
    itt = 1;
    %will now check for the longest streak of consecutive indecies that fall within tolerance
    while itt <= length(chosenIndecies)
        a = chosenIndecies(itt);
        iitt = itt+1; %check for consecutive
        while iitt <= length(chosenIndecies) && chosenIndecies(iitt) - chosenIndecies(itt) < 2
            iitt = iitt+1;
            itt = itt + 1;
        end
        b = chosenIndecies(itt);
        consecutiveToleranceIndecies = [consecutiveToleranceIndecies ; [a,b]];
        itt = itt+1;
    end
    consecutiveToleranceIndecies;

    % will find the first continous set of values within 2% that is at
    % least "minContinuousLength" indecies long 
    for i = 1:height(consecutiveToleranceIndecies)
        if consecutiveToleranceIndecies(i,2) - consecutiveToleranceIndecies(i,1) >= minContinuousLength
            chosenIndecies = consecutiveToleranceIndecies(i,1):1:consecutiveToleranceIndecies(i,2);
            break;
        end
    end
    chosenIndecies;
    %Uses the first index of that streak to get the time
    settlingTime = dataIn(chosenIndecies(1),1); 
    settlingTimeVersion_Tau = settlingTime/4;
end

function [maximumValue,peakTime,overshoot,percentOvershoot] = overshootFinder(dataIn,dataColumn,finalVal)
    [maximumValue , indOfMax] = max(dataIn(:,dataColumn) );
    maximumValue;
    peakTime = dataIn(indOfMax,1);
    overshoot = (maximumValue - finalVal);
    percentOvershoot = overshoot/finalVal *100;
end