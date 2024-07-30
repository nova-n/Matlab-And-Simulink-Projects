clc;
clear;
close all;

BoilingWater2Air = readtable("DataQuench_BoilingWater2Air_2023.xlsx");
BoilingWater2IceWater = readtable("DataQuench_BoilingWater2IceWater_2023.xlsx");
IceWater2BoilingWater = readtable("DataQuench_IceWater2BoilingWater_2023.xlsx");

%looks like the 101 columns are the actual data
data = {BoilingWater2Air,BoilingWater2IceWater,IceWater2BoilingWater};

%Seems like it wasn't quenched initially at the first time entry
%Will be detecting when the temp suddenly changes to exponential decay,
%Because working with large numbers, can just round up to see when the
%exponential part begins
t = data{1}{:,2}; %using the time step of the experiment
T_t = data{1}{:,3}; %temperature at each time

a = ceil(T_t(1)); %its initial starting temp
ii = 2;
while a == ceil(T_t(ii))
    ii = ii+1;
end
ii
round(T_t(ii)); %this is when it detected something different
startingInd = ii-1; %so really starts at ii-1, so need to omit all indecies before startingInd, meaning startingInd-1
startingTime = t(startingInd)

T_0 = data{1}{startingInd,3}
T_final = data{1}{end,3}

t(1:startingInd - 1) = []; %removing all indecies before the quench even starts
T_t(1:startingInd - 1) = [];
t = t - startingTime; %setting the starting index to t = 0 

logarithmicPlot = log(abs((T_t - T_final)./(T_0 - T_final)));
logarithmicPlot = smooth(logarithmicPlot); %smoothing out the plot with a moving average
badIndecies = find(abs(logarithmicPlot) == Inf);%removes temp indecies before is exponential
%must use abs, since some temps before the final in steady state dip below
%the final temp, and it causes a negative inside the ln, which isn't possible
%around the steady state point, can ignore values of inf. Is a
%good approximation because is around steady state
%T_t(height(t)-3) - T_final
usefulIndecies = setdiff(1:height(T_t) , badIndecies); %set diff returns all 
%elements of input 1 that aren't in input 2
t = t(usefulIndecies); % how to use only the good indecies
T_t = T_t(usefulIndecies);
logarithmicPlot = logarithmicPlot(usefulIndecies);

%the test data's ln plot isn't very linear at the ends or beginning, so I
%will have to filter out the section that is most linear
[TF,S1,S2] = ischange(logarithmicPlot,'linear'); 
% this will find the "abrupt" changes in the graph from curved to linear
brkpt=t(TF==1) % gives the values of time where the abrupt changes happen
%in my case, it seems to happen between the first and second abrupt changes
breakIndecies = [find(t == brkpt(1)), find(t == brkpt(2)) ]

bestFitLine = polyfit(t(breakIndecies(1):breakIndecies(2)),logarithmicPlot(breakIndecies(1):breakIndecies(2)),1); %slope is in bestFitLine(1)
%bestFitLine = polyfit(t,logarithmicPlot,1); %slope is in bestFitLine(1)
timeConst = -1/bestFitLine(1)
%timeConst = 200;

figure;
plot(t,logarithmicPlot)
hold on;
plot(t,bestFitLine(1)*t+bestFitLine(2))
lgnd= legend("ln plot","best fit line");
title("Natural Log plot")
xlabel("time") 
hold off;

figure;
modelT_t = T_final + (T_0 - T_final)*exp(1).^(-t/timeConst);
y_plots = [T_t , modelT_t];
[rows,cols]=size(y_plots);

for i = 1:cols
    plot(t,y_plots(:,i),LineWidth=2)
    hold on;
end
graphTitles = ["$ Experimental \ Temperature \ vs \ Theoretical \ Temperature $","$ Time \ In \ Seconds $","$ Temp \ in \ C^{\circ} $"];
legendStuff = {'$ Experiment \ Temperature \ Data $',sprintf('$ T(t)_{theoretical} = %f + (%f - %f)e^{-t/%f} $',T_final,T_0,T_final,timeConst)};
lgnd= legend(legendStuff);
set(lgnd, 'Interpreter','latex')
lgnd.FontSize = 12;
lgnd.Location = 'best';
title(graphTitles(1),'Interpreter','latex')
xlabel(graphTitles(2),'Interpreter','latex') 
ylabel(graphTitles(3),'Interpreter','latex') 
yline(0,'LineWidth',3,'HandleVisibility','off')
hold off;

function [timeConstant,filtered_t,filtered_T_t] = timeConstantSolver(time, Temp)
end