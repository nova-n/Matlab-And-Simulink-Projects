clc;
clear;
close all;

res = 100;
t = linspace(0,10,res);

%%Problem 1
y = -5/9 + (5/3)*t + (95/9)*exp(1).^(-3*t);
y_h = (95/9)*exp(1).^(-3*t);
y_p = -5/9 + (5/3)*t;
y_plots = [y;y_h;y_p];
%y_plots(1,:)
names = {'$ y = -\frac{5}{9} + \frac{5}{3}t + \frac{95}{9}e^{-3t} $','$ y_{h} = \frac{95}{9}e^{-3t} $','$ y_{p} = -\frac{5}{9} + \frac{5}{3}t $'};
colors = ["b-" , "r--" , "g-."];
titlesAndLabels = ["$ Homework \ \#2 \ Problem \ \# 1: \ y(t) \ vs \ t $","$ t $","$ y(t) $"]; %"title", "x label", "y label"
figure;
graphsPlotter2D(t,y_plots,names,colors,titlesAndLabels,NaN)
hold off;

%%Problem 2
y = 1/4 + (1/4)*exp(1).^(-2*t) - (1/2)*exp(1).^(-2*t);
y_h = (1/4)*exp(1).^(-2*t) - (1/2)*exp(1).^(-2*t);
y_p = 1/4 * ones(1,length(t)); %need to have an entire vector of the same constant
y_plots = [y;y_h;y_p];
names = {'$ y = \frac{1}{4} + \frac{1}{4}e^{-2t} - \frac{1}{2}e^{-2t} $','$ y_{h} = \frac{1}{4}e^{-2t} - \frac{1}{2}e^{-2t} $','$ y_{p} = \frac{1}{4} $'};
colors = ["b-" , "r--" , "g-."];
titlesAndLabels = ["$ Homework \ \#2 \ Problem \ \# 2: \ y(t) \ vs \ t $","$ t $","$ y(t) $"]; %"title", "x label", "y label"
figure;
graphsPlotter2D(t,y_plots,names,colors,titlesAndLabels,NaN)
hold off;

%%Problem 3
y = -5 + (5/4)*exp(1).^(-1*t) + (5/4)*exp(1).^(t) + (5/2).*cos(t);
y_h = (5/4)*exp(1).^(-1*t) + (5/4)*exp(1).^(t) + (5/2).*cosd(t);
y_p = -5 * ones(1,length(t)); %need to have an entire vector of the same constant
y_plots = [y;y_h;y_p];
names = {'$ y = -5 + \frac{5}{4}e^{-t} - \frac{5}{4}e^{t} +\frac{5}{2}cos(t) $','$ y_{h} = \frac{5}{4}e^{-t} - \frac{5}{4}e^{t} +\frac{5}{2}cos(t) $','$ y_{p} = -5 $'};
colors = ["b-" , "r--" , "g-."];
titlesAndLabels = ["$ Homework \ \#2 \ Problem \ \# 3: \ y(t) \ vs \ t $","$ t $","$ y(t) $"]; %"title", "x label", "y label"
figure;
graphsPlotter2D(t,y_plots,names,colors,titlesAndLabels,[0,3;-10,20])
hold off;

%%Problem 4
y = 1 + 3*exp(1).^((-1/2)*t).*sin((1/2)*t) + - exp(1).^((1/2)*t).*cos((1/2)*t);
y_h = 3*exp(1).^((-1/2)*t).*sin((1/2)*t) + - exp(1).^((1/2)*t).*cos((1/2)*t);
y_p = ones(1,length(t)); %need to have an entire vector of the same constant
y_plots = [y;y_h;y_p];
names = {'$ y = 1 + 3e^{-\frac{1}{2}t}sin(\frac{1}{2}t) - e^{\frac{1}{2}t}cos(\frac{1}{2}t) $','$ y_{h} = 3e^{-\frac{1}{2}t}sin(\frac{1}{2}t) - e^{\frac{1}{2}t}cos(\frac{1}{2}t) $','$ y_{p} = 1 $'};
colors = ["b-" , "r--" , "g-."];
titlesAndLabels = ["$ Homework \ \#2 \ Problem \ \# 4: \ y(t) \ vs \ t $","$ t $","$ y(t) $"]; %"title", "x label", "y label"
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