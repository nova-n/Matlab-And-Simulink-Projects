clc;
clear;
close all;

BoilingWater2Air = readtable("DataQuench_BoilingWater2Air_2023.xlsx");
BoilingWater2IceWater = readtable("DataQuench_BoilingWater2IceWater_2023.xlsx");
IceWater2BoilingWater = readtable("DataQuench_IceWater2BoilingWater_2023.xlsx");

experimentNames = ["$ Quench: \ Boiling \ Water \ To \ Air $","$ Quench: \ Boiling \ Water \ To \ Ice \ Water $" "$ Quench: \ Ice \ Water \ To \ Boiling \ Water $"];

%looks like the 101 columns are the actual data
data = {BoilingWater2Air,BoilingWater2IceWater,IceWater2BoilingWater};
j=1;
t = data{j}{:,2}; %using the time step of the experiment
T_t = data{j}{:,3}; %temperature at each time
figure;
plot(t,T_t)
hold on;
[TF,S1,S2] = ischange(T_t,'linear'); 
S1
% this will find the "abrupt" changes in the graph from curved to linear
plot(t,S1.*t+S2)
hold off;
breakpt=t(TF==1) % gives the values of time where the abrupt changes happen
%in my case, it seems to happen between the first and second abrupt changes
breakIndecies = find(TF==1);

for i = 1:length(data)
    t = data{i}{:,2}; %using the time step of the experiment
    T_t = data{i}{:,3}; %temperature at each time
    
    [T_0,T_final,t,T_t] = trueStartGraph(t,T_t,data{i}); %Resets the graph so that it starts when the quench does
    
    figure; %plotting the logarithmic line
    [timeConst,t,T_t] = timeConstantSolver(T_t,T_final,T_0,t,experimentNames(i)); %solves the time constant, but loses some increments of t and T_t
    hold off;
    
    figure;
    experimentalAndTheoreticalPlotter(T_final,T_0,t,timeConst,T_t,experimentNames(i))
    hold off;
end

function [initTemp,finalTemp,newTime,newTemps] = trueStartGraph(time,Temp,dataTable)
    %Seems like it wasn't quenched initially at the first time entry
    %Will be detecting when the temp suddenly changes to exponential decay,
    %Because working with large numbers, can just round up to see when the
    %exponential part begins

    a = ceil(Temp(1)); %its initial starting temp
    ii = 2;
    while a == ceil(Temp(ii))
        ii = ii+1;
    end
    
    %% note, this part is only good for small changes, such as straight-ish sections
    ceil(Temp(ii)); %this is when it detected something different
    startingInd = ii-1; %so really starts at ii-1, so need to omit all indecies before startingInd, meaning startingInd-1
    startingTime = time(startingInd);
    
    initTemp = dataTable{startingInd,3};
    finalTemp = dataTable{end,3};
    
    time(1:startingInd - 1) = []; %removing all indecies before the quench even starts
    Temp(1:startingInd - 1) = [];
    newTemps = Temp;
    newTime = time - startingTime; %setting the starting index to t = 0 
end


function [timeConstant,filtered_t,filtered_T_t] = timeConstantSolver(Temps,finalTemp,initTemp,time,graphTitle)
    logarithmicPlot = log(abs((Temps - finalTemp)./(initTemp - finalTemp)));
    logarithmicPlot = smooth(logarithmicPlot); %smoothing out the plot with a moving average
    badIndecies = find(abs(logarithmicPlot) == Inf);%removes temp indecies before is exponential
    %must use abs, since some temps before the final in steady state dip below
    %the final temp, and it causes a negative inside the ln, which isn't possible
    %around the steady state point, can ignore values of inf. Is a
    %good approximation because is around steady state
    usefulIndecies = setdiff(1:height(Temps) , badIndecies); %set diff returns all 
    %elements of input 1 that aren't in input 2
    time = time(usefulIndecies); % how to use only the good indecies
    filtered_t = time;
    Temps = Temps(usefulIndecies);
    filtered_T_t = Temps;
    logarithmicPlot = logarithmicPlot(usefulIndecies);
    
    %%solving fot the time constant
    %the test data's ln plot isn't very linear at the ends or beginning, so I
    %will have to filter out the section that is most linear
    [TF,S1,S2] = ischange(logarithmicPlot,'linear'); 
    % this will find the "abrupt" changes in the graph from curved to linear
    brkpt=time(TF==1); % gives the values of time where the abrupt changes happen
    %in my case, it seems to happen between the first and second abrupt changes
    breakIndecies = [find(time == brkpt(1)), find(time == brkpt(2)) ];
    
    bestFitLine = polyfit(time(breakIndecies(1):breakIndecies(2)),logarithmicPlot(breakIndecies(1):breakIndecies(2)),1); %slope is in bestFitLine(1)
    %bestFitLine = polyfit(t,logarithmicPlot,1); %slope is in bestFitLine(1)
    timeConstant = -1/bestFitLine(1);

    %%dealing with the plot
    plot(time,logarithmicPlot)
    hold on;
    plot(time,bestFitLine(1)*time+bestFitLine(2))
    lgnd= legend({"$ ln() \ plot $"," $ best \ fit \ line  =  \frac{-1}{\tau}$"});
    set(lgnd, 'Interpreter','latex')
    graphTitle = erase(graphTitle,"$"); %removes both the $ so I can add the strings together
    title("$ Natural \ Log \ plot \ for \ -- \" +graphTitle +"$",'Interpreter','latex')
    xlabel("$ Time \ In \ Seconds $",'Interpreter','latex') 
end

function [graph] = experimentalAndTheoreticalPlotter(finalTemp,initTemp,time,timeConstant,Temps,graphTitle)
    modelT_t = finalTemp + (initTemp - finalTemp)*exp(1).^(-time/timeConstant);
    y_plots = [Temps , modelT_t];
    [rows,cols]=size(y_plots);
    
    for i = 1:cols
        plot(time,y_plots(:,i),LineWidth=2)
        hold on;
    end

    graphTitles = [graphTitle,"$ Time \ In \ Seconds $","$ Temp \ in \ C^{\circ} $"];
    legendStuff = {'$ Experiment \ Temperature \ Data $',sprintf('$ T(t)_{theoretical} = %f + (%f - %f)e^{-t/%f} $',finalTemp,initTemp,finalTemp,timeConstant)};
    lgnd= legend(legendStuff);
    set(lgnd, 'Interpreter','latex')
    lgnd.FontSize = 12;
    lgnd.Location = 'best';
    title(graphTitles(1),'Interpreter','latex')
    xlabel(graphTitles(2),'Interpreter','latex') 
    ylabel(graphTitles(3),'Interpreter','latex') 
    yline(0,'LineWidth',3,'HandleVisibility','off')
end