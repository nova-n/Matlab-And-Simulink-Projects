clear;
clc;
close all;

sensorData = table2array(readtable("Ducky Lab 2 Data.xlsx"','Range','A2:B18'));
setPoint = sensorData(1,:);%( distance in inches x_0 , voltage reading v_so )
x = linspace(sensorData(1,1),sensorData(end,1),100);

%%Linear Approximation
linearishStart = 0.75;
startIndex = floor( height(sensorData) * linearishStart );
if startIndex == 0
    startIndex = 1;
end
K_s = polyfit(sensorData(startIndex:end,1),sensorData(startIndex:end,2),1);
linearApprox = setPoint(2) + K_s(1) * ( (x - sensorData(startIndex,1) ) - setPoint(1) ) - sensorData(startIndex,2);

%%Non-Linear Curve Fit
%is actually done by linear curve fitting x vs 1/V_s
%(sensorData(:,1),1./sensorData(:,2))
%if you plot this ^^^, you can see it actually becomes quite linear...
inv = polyfit(sensorData(:,1),1./sensorData(:,2),1);
curveApprox = 1./( inv(1)*x + inv(2) );

metricInv = polyfit(sensorData(:,1),1./(0.0254*sensorData(:,2)),1);
metricCurveApprox = 1./( metricInv(1)*x + metricInv(2) );

plot(sensorData(:,1),sensorData(:,2),".")
hold on;
plot(x,linearApprox)
hold on;
plot(x,curveApprox)

title("$Sensor \ Voltage \ V_{s}  \ vs \ Distance \ (in)$",'Interpreter','latex')
xlabel("$ Distance \ (in) $",'Interpreter','latex') 
ylabel("$Sensor \ Voltage \ V_{s}$",'Interpreter','latex')
legendStuff = ["$Raw \ Data$", strcat("$Linear Approximation, K_{s} = \ " , num2str(K_s(1)) , "( \frac{V}{in} ) $")...
    strcat("$Rational Approximation, V_{s}(x) = \ \frac{1}{ " , num2str(inv(1)) , " x + " , num2str(inv(1)), " } $")];
lgnd= legend(legendStuff,'Location','best');
set(lgnd, 'Interpreter','latex')

%hold on
%plot(sensorData(:,1),1./sensorData(:,2),".") %if you plot this, you can
%see it actually becomes quite linear...

plot1 = strcat("Sensor Voltage vs Distance");
print('-r600','-dpng',plot1);
hold off;