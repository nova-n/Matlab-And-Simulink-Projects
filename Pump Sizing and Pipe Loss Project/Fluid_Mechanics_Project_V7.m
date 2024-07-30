clear;
clc;
yearsInSeconds = 6.307e+8; %20 years in seconds
hoursADay = 8;
maxTime = yearsInSeconds / ( 24 / hoursADay);
powerCostRate = 0.10;%10 cents an hour
pumpInitialCost = 3500;
impellerSizeCostRate = 1500*12; % $1500 per inch, in feet
valveInitialCost = 300;
valveSizeCostRate = 200 * 12; % $200 per inch, in feet
elbowInitialCost = 50;
elbowSizeCostRate = 50 * 12; % $50 per inch, in feet 
pipeSizeCostRate = 1*12; % $1 per inch of diam per foot of length, in feet

testCurvesPerVar = 50;
meshPoints = testCurvesPerVar*5;
flowRateMin = 6.96252894;
flowRateMax = flowRateMin*3;
Q =linspace(0.001,flowRateMax,meshPoints);
Dguess=linspace(0.5,1.5,testCurvesPerVar);%make sure the minumum is greater than the resolution
density = 1.94;
dynamicViscosity = 0.000018588;
deltaP = 0;
b = 4224;
H = 150;
L = sqrt(H^2 + b^2);
E=4.92e-4;
g=32.2;

%https://www.engineeringtoolbox.com/minor-loss-coefficients-pipes-d_626.html
%minorLossComponents = ["90 deg elbow","butterfly valves"];
minorLossCoefficients = [0.3, 0.2];
minorLossComponentsQuantities = [10,4];
minorLossSumCoefficients = minorLossCoefficients.*minorLossComponentsQuantities;
minorLossMatrix = minorLossMatrixFunc(minorLossSumCoefficients,g, Dguess,Q);

resolution = 0.00001;
Re = [];
Re = reynoldsNumberMatrix(density, dynamicViscosity, Dguess, Q);
fFrictionMatrix = coleBrookMatrix( E, Dguess, Re, resolution);
majorLossMatrix = majorLossMatrixFunc(L,g,fFrictionMatrix,Dguess,Q); 

totalHeadLossMatrix = majorLossMatrix + minorLossMatrix + H;

%%Pump Data
pump_A = readmatrix('Pump_A_Data');
pump_A(:,1) = pump_A(:,1) * (1/7.48)*(1/60);
impellerDiam_A = 0.454;
wGiven_A = 1760*pi/30;
Qfine_A = linspace(0, pump_A(end,1) ,100);
headCurveCoefficeints_A = polyfit(pump_A(:,1),pump_A(:,2), 2);
headCurve_A = polyval(headCurveCoefficeints_A,Qfine_A);
efficiencyCurveCoefficeints_A = polyfit(pump_A(:,1),pump_A(:,3), 2);
efficiencyCurve_A = polyval(efficiencyCurveCoefficeints_A,Qfine_A);

wInputs = linspace(900*pi/30 , 1800*pi/30, testCurvesPerVar);
newImpellerDiams = linspace(impellerDiam_A*2,2,testCurvesPerVar);

figure;
plot(pump_A(:,1) , pump_A(:,2),"b");
hold on;
plot(pump_A(:,1) , pump_A(:,3));
hold on;
plot(Qfine_A,headCurve_A,"b--")
hold on;
plot(Qfine_A,efficiencyCurve_A,"r--")
hold on;
headTrendLineEqn = sprintf('$%fQ^{2} + %fQ + %f $', headCurveCoefficeints_A(1), headCurveCoefficeints_A(2),headCurveCoefficeints_A(3));
text( mean([pump_A(1,1) , pump_A(end,1)]), max( pump_A(:,2) ) + 5  ,"$h_{a}$ "+"$\approx$ " + headTrendLineEqn, 'FontSize', 16, 'Color', 'b','HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Interpreter','latex');
legend('Pump Head in ft','Pump Efficiency in %');
xlabel('Q , Volumetric Flow Rate (ft^{3}/s)')
title("Pump A Curves")
hold off;

%max efficiency coefficients
maxEffIdx = find(pump_A(:,3)==max(pump_A(:,3)));
flowCoeffAtMaxEff = pump_A(maxEffIdx,1)/(wGiven_A * impellerDiam_A^3);
headRiseCoeffAtMaxEff = g*pump_A(maxEffIdx,2)/(wGiven_A^2 * impellerDiam_A^2);
powerMostEfficient_A = density * g *pump_A(maxEffIdx,1)*pump_A(maxEffIdx,2)/(pump_A(maxEffIdx,3)/100);
powerCoeffAtMaxEff = powerMostEfficient_A/(density * wGiven_A^3 * impellerDiam_A^5);


%Actually, each combo of diameter and speed is its OWN PUMP CURVE!!!
%will store data for each pump: its w, D
count = 1;
newPumpHeadMatrix = [];
pumpCombos = []; %[w,D]
for i = 1:length(newImpellerDiams)
    for ii = 1:length(wInputs)
        newPumpHeadMatrix(count,:) =  ( ( newImpellerDiams(i) * wInputs(ii) ) / (wGiven_A *impellerDiam_A) )^2* polyval(headCurveCoefficeints_A, (  (wGiven_A *impellerDiam_A^3)/(wInputs(ii)*newImpellerDiams(i)^3) ) * Q);
        pumpCombos(count,1) = wInputs(ii);
        pumpCombos(count,2) = newImpellerDiams(i);
        count = count + 1;
    end
end
newPumpHeadMatrix;

BnewPumpHeadMatrix = newPumpHeadMatrix;%getting rid of all negatives to make ymax the highest positive pump head value * 1.1
BnewPumpHeadMatrix(BnewPumpHeadMatrix<0) = 0;

%plotting headloss curves, one curve for each pipe diameter (which is each its own system) 
figure;
for i = 1:length(Dguess)
    plot(Q,totalHeadLossMatrix(:,i),"r")
    hold on;
end
xlabel('Volumetric Flow Rate (ft^{3}/s)')
ylabel('Head (ft)')
title("System Head Curves")
hold off;

%plotting both head loss and pump head curves
figure;
for i = 1:height(newPumpHeadMatrix)
    plot( Q , newPumpHeadMatrix(i,:) ,"b")
    hold on;
end
%yieldYL = get(gca, 'YLim');
%ylim([0,yieldYL(2)]);
hold on;
[ ~ , columns] = size(totalHeadLossMatrix);

for i = 1:columns
    plot(Q,totalHeadLossMatrix(:,i), "r")
    hold on;
end
xlim([0,flowRateMax]);
ylim([0, max(BnewPumpHeadMatrix, [], 'all')*1.1]);
xlabel('Volumetric Flow Rate (ft^{3}/s)')
ylabel('Head (ft)')
title("System Head Curves vs. Single Pump Head Curves")
hold off;


%%Plotting only pump heads at max efficiency
%I cannot just guess diameter anymore. Previously, to scale the curves, I
%just scaled the Q values. However, because I am going point by point, and
%my inputs are Q (the x axis),and w (given range of RPM),I must solve for D
%Since we were already given a range for w, its probably easier to just
%solve for D instead of solving for w

%testCurvesPerVar = 20;
Q =linspace(0.001,flowRateMax,meshPoints);
wInputs = linspace(900*pi/30 , 1800*pi/30, testCurvesPerVar);
newImpellerDiams = [];
mostEffPumpHeads = [];
for i = 1:length(wInputs)
    for ii = 1:length(Q)
        newImpellerDiams(i,ii) = ( Q(ii) / (flowCoeffAtMaxEff * wInputs(i) ))^(1/3);
        mostEffPumpHeads(i,ii) = newImpellerDiams(i,ii)^2 * wInputs(i)^2 * headRiseCoeffAtMaxEff/g;
    end
end
[ ~ , initLength ] = size(totalHeadLossMatrix);
figure;
lowestSystemCurveIdx = find(totalHeadLossMatrix(end,:) == min(totalHeadLossMatrix(end,:))); 
for i = 1:length(wInputs)
    plot( Q , mostEffPumpHeads(i,:),"g")
    hold on;
end
Dguess=linspace(0.5,1.5,initLength);
for i = 1:length(Dguess)
    plot(Q,totalHeadLossMatrix(:,i),"r")
    hold on;
end

%finding intersections
intersections = [];%[Q,h,w,D,pipeDiam,power,cost] ---> [ft^3/s , ft , rad/s , ft , ft , kW , $]
for i = 1:height(mostEffPumpHeads)%looping downwards
    for ii = 1:initLength
        [xx,yy] = polyxpoly(Q,mostEffPumpHeads(i,:),Q,totalHeadLossMatrix(:,ii));
        %possibly 2 intersections,so will just append the row to the matrix
        if isempty(xx) == false
            for iii = 1:height(xx) %the Q where it intersects may not have an entry for a specific D. w is an imput here
                %again, D = ( Q / (flowCoeffAtMaxEff * wInputs) ))^(1/3);
                newD = ( xx(iii) / (flowCoeffAtMaxEff * wInputs(i) ))^(1/3);
                [powUsed,costToRun] = operationCost(Dguess(ii),L,newD,wInputs(i),maxTime,density, powerCoeffAtMaxEff , powerCostRate, pipeSizeCostRate ,[pumpInitialCost, impellerSizeCostRate], [valveInitialCost,valveSizeCostRate,minorLossComponentsQuantities(2)], [elbowInitialCost,elbowSizeCostRate,minorLossComponentsQuantities(1)]);
                intersections = [ intersections ;xx(iii) , yy(iii), wInputs(i), newD , Dguess(ii) ,powUsed, costToRun];
            end
        end
    end
end
%would have been easier to limit the plot minimum to flowrate minimum
%but would have made an ugly graph imo
%removing all intersection points before the minimum flowrate
indecies = find(intersections(:,1)<flowRateMin);
intersections(indecies,:) = [];
xline(flowRateMin , "--");
plot(intersections(:,1) , intersections(:,2),  "." , "MarkerSize" , 10,'color' ,[63, 5, 156]/256 )
hold on;
%plotting cheapest solution
bestCost = min(intersections(:,7))
bestIndex = find(intersections(:,7) == min(intersections(:,7)))
plot(intersections(bestIndex,1) , intersections(bestIndex,2),  "." , "MarkerSize" , 25,'color' ,[3, 152, 252]/256 )
plot(intersections(bestIndex,1) , intersections(bestIndex,2), "hexagram", 'MarkerSize', 20,'color' ,[3, 152, 252]/256,'linewidth',1.5)

xlim([0,flowRateMax]);
ylim([0, max(intersections(:,2))*1.1]);
xlabel('Volumetric Flow Rate (ft^{3}/s)')
ylabel('Head (ft)')
title("System Head Curves vs. Most Efficient Pump Heads (Single Pump System)")
hold off;

%Pipe Cost Only
[~,pipeLineCost] = operationCost(intersections(bestIndex,5),L,newD,0,0,0,0,0, pipeSizeCostRate,[0,0],[0,0,0],[0,0,0])
%Impeller Cost Only
[~,impellerCost] = operationCost(0,0,intersections(bestIndex,4),0,0,0,0,0,0,[pumpInitialCost, impellerSizeCostRate], [0,0,0], [0,0,0])
%Valve Cost Only
[~,valveTotalCost] = operationCost(intersections(bestIndex,5),0,0,0,0,0,0,0,0,[0,0],[valveInitialCost,valveSizeCostRate,minorLossComponentsQuantities(2)],[0,0,0])
%Elbow Cost Only
[~,elbowTotalCost] = operationCost(intersections(bestIndex,5),0,0,0,0,0,0,0,0,[0,0],[0,0,0],[elbowInitialCost,elbowSizeCostRate,minorLossComponentsQuantities(1)])
%Power Cost Only
[~,powerTotalCost] = operationCost(0,0,intersections(bestIndex,4),intersections(bestIndex,3),maxTime,density, powerCoeffAtMaxEff , powerCostRate,0,[0,0],[0,0,0],[0,0,0])

%Best Specs:
disp("Q = "+intersections(bestIndex,1) + " , h = " + intersections(bestIndex,2) + " , w = " + intersections(bestIndex,3) + " , impD = " +intersections(bestIndex,4) + " , pipeD = " + intersections(bestIndex,5) + " , power = " + intersections(bestIndex,6) + " , $ = " + intersections(bestIndex,7))
holdUpWait = pipeLineCost+impellerCost+valveTotalCost+elbowTotalCost+powerTotalCost

%%Functions
function [frictionFactorMatrix] = coleBrookMatrix(rough,diamArray,reynoldsMatrix,resol)
    fFricMatrix = [];
    for i = 1:height(reynoldsMatrix)
        for ii = 1:length(diamArray)
            closeEnough = false;
            fFric = 0;
            tol = 1e-12;%for floating point numbers
            %fprintf('(%f,%f) \n',i,ii)
            %fprintf('Re = %f \n',reynoldsMatrix(i,ii))
            %fprintf('D = %f \n',diamArray(ii))
            while closeEnough == false
                fFricPrev = fFric;
                if abs(1/sqrt(fFric) + 2*log10((rough/diamArray(ii))/3.7+2.51/(reynoldsMatrix(i,ii)*sqrt(fFric)))) > tol %checking inequality with floating point numbers
                    fFric = fFric + resol;
                end
                fFric;%delete the semicolon to see in real time
                %check if leftside - rightside < tolerance, OR 
                %if it overshot, to use the previous value
                if abs(1/sqrt(fFric) + 2*log10((rough/diamArray(ii))/3.7+2.51/(reynoldsMatrix(i,ii)*sqrt(fFric)))) < tol || abs(1/sqrt(fFric) + 2*log10((rough/diamArray(ii))/3.7+2.51/(reynoldsMatrix(i,ii)*sqrt(fFric)))) > abs(1/sqrt(fFricPrev) + 2*log10((rough/diamArray(ii))/3.7+2.51/(reynoldsMatrix(i,ii)*sqrt(fFricPrev)))) 
                    closeEnough = true;
                    if abs(1/sqrt(fFric) + 2*log10((rough/diamArray(ii))/3.7+2.51/(reynoldsMatrix(i,ii)*sqrt(fFric)))) > abs(1/sqrt(fFricPrev) + 2*log10((rough/diamArray(ii))/3.7+2.51/(reynoldsMatrix(i,ii)*sqrt(fFricPrev))))
                        fFric = fFricPrev;
                    end
                end   
            end
            fFric;
            fFricMatrix(i,ii) = fFric;
        end
    end
    frictionFactorMatrix  = fFricMatrix;
end

function [ReNum] = reynoldsNumberMatrix(dens, dynVis, diamArray, volFlowArray)
    for i = 1:length(volFlowArray)
        ReNum(i,:) = 4*dens/(pi * dynVis) * volFlowArray(i)./diamArray ;
        %Vel = 4*volFlowArray(i)./(pi*diamArray.^2);
        %fprintf('V = %f) \n', Vel)
    end   
end

function [hMajorMatrix] = majorLossMatrixFunc(Length, Gravity,fFricMatrix, diameterArray,QArray)
    for i = 1:length(QArray)
        for ii = 1:length(diameterArray)
            hMajorMatrix(i,ii) = (8*Length/(pi^2 * Gravity)) * (fFricMatrix(i,ii)/diameterArray(ii)) * (QArray(i)/(diameterArray(ii)^2))^2;
        end
    end
end

function [hMinorMatrix] = minorLossMatrixFunc(K_m_list,Gravity, diameterArray,QArray)
    for i = 1:length(QArray)
        for ii = 1:length(diameterArray)
            for iii = 1:length(K_m_list)
                hMinorMatrix(i,ii) = 8*K_m_list(iii)/(pi*Gravity) * (QArray(i)/(diameterArray(ii)^2))^2 ;
            end
        end
    end
end

function [pow , cost] = operationCost(pipeDiam,pipeLength,impellerDiam,shaftSpeed,time,dens, powerCoeff , powCostRate, pipeCost ,impellerCostStuff, valveCostStuff, elbowCostStuff)
    %Start Up Costs
    pipeTotal = pipeDiam * pipeLength * pipeCost;
    impellerTotal = impellerCostStuff(1) + impellerCostStuff(2) * impellerDiam;
    valveTotal = (valveCostStuff(1) + valveCostStuff(2) * pipeDiam)*valveCostStuff(3);
    elbowTotal = (elbowCostStuff(1) + elbowCostStuff(2) * pipeDiam)*elbowCostStuff(3);
    startUpCost = pipeTotal + impellerTotal + valveTotal + elbowTotal;
    
    %Calculating Power Used in kW
    power = (dens * powerCoeff * shaftSpeed^3 * impellerDiam^5); % is in lb ft/s. need to convert to kW
    power = power * 0.0013558179483314;
    pow = power;
    %power = density * 32.2 
    
    %Calculating cost of Power
    toHours = time/3600;
    energyCost = powCostRate*toHours*power;
    
    cost = startUpCost + energyCost;
end
