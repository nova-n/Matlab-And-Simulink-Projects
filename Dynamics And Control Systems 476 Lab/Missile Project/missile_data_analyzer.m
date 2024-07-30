%clear;
clc;
close all;

resol = 1000;
t = linspace(0,10,resol);
t = t';

%Original Comparison Inputs
%alphaDesired = (5*pi/180)*ones(1,resol);
%betaDesired = (0*pi/180)*ones(1,resol);

%Step Inputs
%alphaDesired = (5*pi/180)*ones(1,resol);
%betaDesired = (0*pi/180)*ones(1,resol);

%Stronger Step Inputs
alphaDesired = (45*pi/180)*ones(1,resol);
betaDesired = (0*pi/180)*ones(1,resol);

%Linear Inputs
%alphaDesired = (5*pi/180)*t;
%betaDesired = (5*pi/180)*t;

%Stronger Linear Inputs
%alphaDesired = (45*pi/180)*t;
%betaDesired = (45*pi/180)*t;

%Sinusoidal Inputs
%alphaDesired = (45*pi/180)*sin(2*pi*t);
%betaDesired = (45*pi/180)*cos(2*pi*t);

%Faster Sinusoidal Inputs
%alphaDesired = (45*pi/180)*sin(8*pi*t);
%betaDesired = (45*pi/180)*cos(8*pi*t);


alphaDesiredSimIn = timeseries(alphaDesired,t);
betaDesiredSimIn = timeseries(betaDesired,t);

alpha_vs_time = out.alphaData;
alphaDesired = out.desiredAlpha; %gives out an array of all the same value, but only need one of them
alphaDesired = alphaDesired(1);
%stepinfo(alpha_vs_time)

beta_vs_time = out.betaData;
betaDesired = out.desiredBeta; %gives out an array of all the same value, but only need one of them
betaDesired = betaDesired(1);
%stepinfo(beta_vs_time)

[alpha_Tau,alpha_settlingTime] = TwoPercentSettlingTimeTau(alpha_vs_time,2,alphaDesired,10^-0.75,0.10)
[alphaMax,alphaPeakTime,alphaMaxOvershoot,alphaMaxPO] = overshootFinder(alpha_vs_time,2,alphaDesired)

%problematic, since fv=0
%but can get data if just offset it. Ts, tau, peak time, and PO will be accurate, just
%not the max value
%just offset fv, and data by +10^-3
[beta_vs_time(:,1),beta_vs_time(:,2)]
[beta_Tau,beta_settlingTime] = TwoPercentSettlingTimeTau([beta_vs_time(:,1),beta_vs_time(:,2)],2,betaDesired,10^-4,0.15)
[betaMax,betaPeakTime,betaMaxOvershoot,betaMaxPO] = overshootFinder([beta_vs_time(:,1),beta_vs_time(:,2)],2,betaDesired)
%betaMax = betaMax - 10^-3;

figure;
pointPlot1 = plot(alpha_vs_time(:,1),alpha_vs_time(:,2),'LineWidth',2);
hold on;
yline(alphaDesired,"--",'LineWidth',2)
grid on;
plot(alphaPeakTime,alphaMax,"hexagram",'color',get(pointPlot1,'color'),'MarkerSize',12,'linewidth',2);
plot(alphaPeakTime,alphaMax,".",'color',get(pointPlot1,'color'),'HandleVisibility','off','MarkerSize',16);
xline(alphaPeakTime,"-.",'LineWidth',1)
xline(alpha_settlingTime,"-.",'LineWidth',1)
hold off;
xAxisLabel = "$Time \ (s)$";
xlabel(xAxisLabel,'Interpreter','latex') 
yAxisLabel = "$Angle \ of \ Attack, \ \alpha  \ ( deg^{\circ} ) $";
ylabel(yAxisLabel,'Interpreter','latex') 
legendStuff = ["$\alpha$",strcat("$\alpha_{desired} =" ,num2str(alphaDesired) , " ^{\circ} $"),strcat("$P.O = ",num2str(alphaMaxPO) ," \% $"),...
    strcat("$T_{p} = ", num2str(alphaPeakTime) ," s $"),strcat("$T_{s} = " ,num2str(alpha_settlingTime) , " s $")];
lgnd= legend(legendStuff,'Location','best');
set(lgnd, 'Interpreter','latex')
plot1 = "alpha";
print('-r600','-dpng',plot1);

figure;
pointPlot1 = plot(beta_vs_time(:,1),beta_vs_time(:,2),'LineWidth',2);
hold on;
yline(betaDesired,"--",'LineWidth',2)
grid on;
plot(betaPeakTime,betaMax,"hexagram",'color',get(pointPlot1,'color'),'MarkerSize',12,'linewidth',2);
plot(betaPeakTime,betaMax,".",'color',get(pointPlot1,'color'),'HandleVisibility','off','MarkerSize',16);
xline(betaPeakTime,"-.",'LineWidth',1)
hold off;
xAxisLabel = "$Time \ (s)$";
xlabel(xAxisLabel,'Interpreter','latex') 
yAxisLabel = "$Angle \ of \ Side Step, \ \beta  \ ( deg^{\circ} ) $";
ylabel(yAxisLabel,'Interpreter','latex') 
legendStuff = ["$\beta$",strcat("$\beta_{desired} =" ,num2str(betaDesired) , " ^{\circ} $"),strcat("$beta_{max} = ",num2str(betaMax) ," $"),...
    strcat("$T_{p} = ", num2str(betaPeakTime) ," s $")];
lgnd= legend(legendStuff,'Location','best');
set(lgnd, 'Interpreter','latex')
plot1 = "beta";
print('-r600','-dpng',plot1);

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