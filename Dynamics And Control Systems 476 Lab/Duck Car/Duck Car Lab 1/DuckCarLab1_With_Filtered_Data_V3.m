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

filteredThreeVData = dataNoiseRemover(ThreeVData,10); %Arguments: Raw Data, Cutoff Freq
filteredThreeVData = dataSmoother(filteredThreeVData,50,10); %Arguments: LowPassed Data, How Many Times Moving Avg Taken, Window
filteredFourVData = dataNoiseRemover(FourVData,5);
filteredFourVData = dataSmoother(filteredFourVData,20,3);
filteredFiveVData = dataNoiseRemover(FiveVData,15);
filteredFiveVData = dataSmoother(filteredFiveVData,150,3);

ThreeV_fvVel = finalValueFinder(filteredThreeVData,1*10^(-2.5),[2/3,1]); % Smoothed Out Data, Tolerance, What sections of the graph you want
FourV_fvVel = finalValueFinder(filteredFourVData,1*10^(-2.25),[0.5,1]);
FiveV_fvVel = finalValueFinder(filteredFiveVData,1*10^(-2.5),[0.5,1]);

ThreeV_TimeConst = timeConstantFinder(filteredThreeVData,ThreeV_fvVel,0.06,0.325,"$3 V \ Dimnensionless \ Velocity$","$Dimnensionless \ Velocity$","3V")
FourV_TimeConst = timeConstantFinder(filteredFourVData,FourV_fvVel,0.075,0.485,"$4 V \ Dimnensionless \ Velocity$","$Dimnensionless \ Velocity$","4V")
FiveV_TimeConst = timeConstantFinder(filteredFiveVData,FiveV_fvVel,0.075,0.45,"$5 V \ Dimnensionless \ Velocity$","$Dimnensionless \ Velocity$","5V")

[ThreeV_FirstOrderApprox,ThreeV_Shift] = firstOrderApproximation(filteredThreeVData,ThreeVData,ThreeV_fvVel,ThreeV_TimeConst);
[FourV_FirstOrderApprox,FourV_Shift] = firstOrderApproximation(filteredFourVData,FourVData,FourV_fvVel,FourV_TimeConst);
[FiveV_FirstOrderApprox,FiveV_Shift] = firstOrderApproximation(filteredFiveVData,FiveVData,FiveV_fvVel,FiveV_TimeConst);

rawAndFilteredPlotter(ThreeVData,filteredThreeVData,"$3 V \ Velocity$","$Velocity \ ( \frac{m}{s} )$",ThreeV_fvVel,ThreeV_TimeConst,ThreeV_FirstOrderApprox,ThreeV_Shift,"3V")
rawAndFilteredPlotter(FourVData,filteredFourVData,"$4 V \ Velocity$","$Velocity \ ( \frac{m}{s} )$",FourV_fvVel,FourV_TimeConst,FourV_FirstOrderApprox,FourV_Shift,"4V")
rawAndFilteredPlotter(FiveVData,filteredFiveVData,"$5 V \ Velocity$","$Velocity \ ( \frac{m}{s} )$",FiveV_fvVel,FiveV_TimeConst,FiveV_FirstOrderApprox,FiveV_Shift,"5V")

multiPlotter({ThreeVData(:,1:2),FourVData(:,1:2),FiveVData(:,1:2)},"$Raw \ Velocity \ Data$" , ["$3 V $","$4 V $","$5 V $"] , "$Velocity \ ( \frac{m}{s} )$","Raw Velocity Data")
multiPlotter({filteredThreeVData(:,1:2),filteredFourVData(:,1:2),filteredFiveVData(:,1:2)},"$Filtered \ Velocity \ Data$" ,...
    ["$3 V \ (Filtered)$","$4 V \ (Filtered) $","$5 V \ (Filtered) $"] , "$Velocity \ ( \frac{m}{s} )$","Filtered Velocity Data")
multiPlotter({ThreeV_FirstOrderApprox(:,1:2),FourV_FirstOrderApprox(:,1:2),FiveV_FirstOrderApprox(:,1:2)},...
    "$First \ Order \ Approximations Of \ Velocity$" ,...
    [strcat("$3 V ,", "V(t) = " ,  num2str(ThreeV_fvVel)," (1- e^{ " , num2str(-1/ThreeV_TimeConst),...
        "(t -" , num2str(-1*ThreeV_Shift)," ) } ) \ , \" , " \tau = \ ", num2str(ThreeV_TimeConst), "s, \ $")...
        ,strcat("$4 V ,", "V(t) = " ,  num2str(FourV_fvVel)," (1- e^{ " , num2str(-1/FourV_TimeConst),...
        "(t -" , num2str(-1*FourV_Shift)," ) } ) \ , \" , " \tau = \ ", num2str(FourV_TimeConst), "s, \ $")...
        ,strcat("$5 V ,", "V(t) = " ,  num2str(FiveV_fvVel)," (1- e^{ " , num2str(-1/FiveV_TimeConst),...
        "(t -" , num2str(-1*FiveV_Shift)," ) } ) \ , \" , " \tau = \ ", num2str(FiveV_TimeConst), "s, \ $")] ...
        , "$Velocity \ ( \frac{m}{s} )$","First Order Velocity Approximations")

ThreeV_TwoPercent_Tau =TwoPercentSettlingTimeTau(filteredThreeVData(:,1:2),ThreeV_fvVel,10^-2.5,0.03)
FourV_TwoPercent_Tau =TwoPercentSettlingTimeTau(filteredFourVData(:,1:2),FourV_fvVel,10^-2.5,0.05)
FiveV_TwoPercent_Tau =TwoPercentSettlingTimeTau(filteredFiveVData(:,1:2),FiveV_fvVel,10^-2.5,0.03)

% %must ALWAYS skip first and last velocity entry, since are NaN, and so if want to
% %plot it, also omit those entries for position

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

function [smoothedOutData] = dataSmoother(dataIn,sandingPasses,movingAvgWindow)
    %Runs several moving averages on the noise-removed data in one call, to
    %further smooth it out.
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

function [finalValue] = finalValueFinder(dataInputted,steadyStateTolerance,timeSection)
    %will only look at last set of time values for final value
    lastSectionOfTime = [floor( timeSection *height(dataInputted)),ceil( timeSection *height(dataInputted))];
    if lastSectionOfTime(1) == 0
        lastSectionOfTime(1) = 1;
    end
    a = diff(dataInputted(lastSectionOfTime(1):lastSectionOfTime(2),2))./diff(dataInputted(lastSectionOfTime(1):lastSectionOfTime(2),1));
    steadyestLastIndecies = find(abs(a)<=steadyStateTolerance) + lastSectionOfTime(1) - 1; 
    %checks for any changes smaller than the tolerance, so that change is close to 0
    consecutiveSteadyIndecies = [;]; %  [ start1,start2 ; end1,end2 ; ... ]
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
    %runs these indecies through the filtered function, and averages the values at these indecies
end

function [timeConstant] = timeConstantFinder(dataIn,finalVelocity,firstSectionPercentage,endSectionPercentage,givenTitle,measuredQuantity,voltageString)
    dimensionlessVelocity = log((finalVelocity - dataIn(:,2))/finalVelocity);%log without subscript in matlab is ln
    graphSection = [floor(firstSectionPercentage * length(dataIn(:,1)) ) , ceil(endSectionPercentage * length(dataIn(:,1)) )];
    if graphSection(1) == 0
        graphSection(1) = 1;
    end
    slope = polyfit( dataIn(graphSection(1):graphSection(2),1), ...
        dimensionlessVelocity(graphSection(1):graphSection(2)) ,1 ); 
    %using polyfit to get the slope of the linear part of the plot
    timeConstant = -1/slope(1); 
    figure;
    plot(dataIn(:,1),dimensionlessVelocity); 
    hold on;
    t = linspace(dataIn(1,1),dataIn(end,1),100);
    plot(t,t*slope(1));
    hold off;
    title(givenTitle,'Interpreter','latex')
    xlabel("$ Time, \ t \ (s) $",'Interpreter','latex') 
    ylabel(measuredQuantity,'Interpreter','latex')
    legendStuff = [givenTitle, strcat("$Estimated \ " , givenTitle, ", \ \tau = \ ", num2str(timeConstant), "s \ $"),strcat("$Final \ Value \ " , givenTitle, " \ $")];
    lgnd= legend(legendStuff,'Location','best');
    set(lgnd, 'Interpreter','latex')
    plot1 = strcat("Dimensionless Velocity And Time Constant For ", voltageString ,  " vs time");
    print('-r600','-dpng',plot1);
end

function [firstOrder,shiftPoint] = firstOrderApproximation(filteredData,originalData,finalVal,timeConst)
    pointAtInitTime = originalData(2,2); %is NaN velocity at t=0,
    %y0 = fv * (1 - e^(-(1/tau) * (t - x0) ) ) ----> 1 - (y0/fv) = e^(-(1/tau) * (t - x0) )
    %ln( 1 - (y0/fv) ) / (-1/tau) = (t - x0) ----> x0 = t + tau * ln( 1 - (y0/fv) )
    shiftPoint = originalData(2,1) + timeConst*log( 1 - pointAtInitTime/finalVal);
    %Accounts for the fact that the recorded data from tracker may not have
    %started at 0 velocity, so just shifts the graph left or right to try
    %getting the fit better.
    t = linspace(filteredData(1,1),filteredData(end,1),1000);
    firstOrderApprox = finalVal * (1 - exp( (-1/timeConst) * (t - shiftPoint) ) );
    firstOrder = [t',firstOrderApprox'];
end

function rawAndFilteredPlotter(originalData,filteredData,givenTitle,measuredQuantity,finalVal,timeConst,firstOrder,timeShift,voltageString)
    figure;
    plot(originalData(:,1),originalData(:,2))
    hold on;
    plot(filteredData(:,1),filteredData(:,2))
    hold on;
    yline(finalVal,"--b");
    hold on;
    plot(firstOrder(:,1),firstOrder(:,2));
    hold off;
    
    title(givenTitle,'Interpreter','latex')
    xlabel("$ Time, \ t \ (s) $",'Interpreter','latex') 
    ylabel(measuredQuantity,'Interpreter','latex')
    legendStuff = [givenTitle, strcat("$Filtered \ And \ Smoothed \ " , givenTitle, " \ $")...
        ,strcat("$Final \ Value \ " , givenTitle, " \ = \ ", num2str(finalVal) ," \ ( \frac{m}{s} ) $" ),...
        strcat("$First \ Order \ Approximation \ Of \ " ,...
        givenTitle , " \ $" , newline + " $ V(t) = " ,  num2str(finalVal)," (1 - e^{ " , num2str(-1/timeConst),...
        " (t + " , num2str(timeShift)," ) } ) \ , \" , " \tau = \ ", num2str(timeConst), "s, \ $")];
    lgnd= legend(legendStuff,'Location','best');
    set(lgnd, 'Interpreter','latex')
    plot1 = strcat("All Velocity Plots For", voltageString ,  " vs time");
    print('-r600','-dpng',plot1);
end

function multiPlotter(datas,givenTitle,legendItems,measuredQuantity,fileNameString)
    figure;
    for i = 1:length(datas)
        plot(datas{i}(:,1),datas{i}(:,2))
        hold on;
    end
    hold off;
    title(givenTitle,'Interpreter','latex')
    xlabel("$ Time, \ t \ (s) $",'Interpreter','latex') 
    ylabel(measuredQuantity,'Interpreter','latex')
    lgnd= legend(legendItems,'Location','best');
    set(lgnd, 'Interpreter','latex')
    plot1 = strcat("Multiplot Of", fileNameString ,  " vs time");
    print('-r600','-dpng',plot1);
end

function [settlingTimeVersion_Tau] = TwoPercentSettlingTimeTau(dataIn,finalVal,tol,streakFraction)
    minContinuousLength = floor( streakFraction * length(dataIn) ); %has to be at least spanning 3% of the graph
    %minContinuousLength = floor( 0.03 * length(dataIn) ); %has to be at least spanning 3% of the graph
    
    % will find the first continous set of values within 2% that is at
    % least "minContinuousLength" indecies long 

    chosenIndecies = find(abs(finalVal*0.982 - dataIn(:,2)) <= tol); %if the data is 98% of the final val within a tolerance

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

    %Uses the first index of that streak to get the time
    settlingTime = dataIn(chosenIndecies(1),1); 
    settlingTimeVersion_Tau = settlingTime/4;
end