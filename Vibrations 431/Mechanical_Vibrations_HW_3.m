clc;
clear;
close all;

res = 1000;
t = linspace(0,100,res);

%%Problem 5
y = -1250 + 25*t + 1255*exp(1).^(-0.02*t);
y_h = 1255*exp(1).^(-0.02*t);
y_p = -1250 + 25*t;
y_plots = [y;y_h;y_p];
names = {'$ y = -1250 + 25t + 1255e^{-0.02t} $','$ y_{h} = 1255e^{-0.02t} $','$ y_{p} = -1250 + 25t $'};
colors = ["b-" , "r--" , "g-."];
titlesAndLabels = ["$ Homework \ \#3 \ Problem \ \# 5: \ y(t) \ vs \ t $","$ t $","$ y(t) $"]; %"title", "x label", "y label"
figure;
graphsPlotter2D(t,y_plots,names,colors,titlesAndLabels,NaN)
hold off;

%%Problem 7
y = 0.2693*exp(1).^((-4+2*sqrt(3))*t) + 0.0193*exp(1).^((-4-2*sqrt(3))*t) + 1/4;
y_h = 0.2693*exp(1).^((-4+2*sqrt(3))*t) + 0.0193*exp(1).^((-4-2*sqrt(3))*t);
y_p = 1/4 * ones(1,length(t)); %need to have an entire vector of the same constant
y_plots = [y;y_h;y_p];
names = {'$ y \approx \frac{1}{4} + 0.2693e^{(-4+2\sqrt(3))t} - 0.0193e^{(-4+2\sqrt(3))t} $','$ y_{p} \approx 0.2693e^{(-4+2\sqrt(3))t} - 0.0193e^{(-4+2\sqrt(3))t} $','$ y_{p} = \frac{1}{4} $'};
colors = ["b-" , "r--" , "g-."];
titlesAndLabels = ["$ Homework \ \#3 \ Problem \ \# 7 \ y(t) \ vs \ t $","$ t $","$ y(t) $"]; %"title", "x label", "y label"
figure;
graphsPlotter2D(t,y_plots,names,colors,titlesAndLabels,NaN)
hold off;


%%Plotting Function
function plotReturn =  graphsPlotter2D(x,y_graphs,legendStuff,lineColors,graphTitles,limits)
    %figure;
    for i = 1:height(y_graphs)
        plot(x,y_graphs(i,:),lineColors(i),LineWidth=2)
        hold on;
    end
    lgnd= legend(legendStuff);
    set(lgnd, 'Interpreter','latex')
    lgnd.FontSize = 12;
    lgnd.Location = 'best';
    title(graphTitles(1),'Interpreter','latex')
    xlabel(graphTitles(2),'Interpreter','latex') 
    ylabel(graphTitles(3),'Interpreter','latex') 
    yline(0,'LineWidth',3,'HandleVisibility','off')
    if isnan(limits) == false
        xlim(limits(1,:));
        ylim(limits(2,:));
    end
  
    %hold off;
end