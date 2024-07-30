clc;
clear;
close all;

fileName = 'C:/Users/n8dsa/Documents/OpenFOAM/bl/postProcessing/residuals/0/residuals.txt';
A_resid = readtable(fileName,'filetype', 'text');


%plotting residuals for y vs u_x along centerline x = 0.5
y = table2array( A_resid(:,1) );
for ii = 2:width(A_resid) -1        
    plot(y, table2array( A_resid(:,ii) ) )
    hold on;
end

legendStuff =["$p$", "$Ux$", "$Uz$","k","$ \omega $"];
title(strcat("$Residuals \ of \ Turbulent \ Boundary \ Layer \ \Re = ",num2str( 8*10^6 ), " \ at \ Center \ Line \ x = 3 $"),'Interpreter','latex')
xlabel("$Itterations$",'Interpreter','latex')
lgnd= legend(legendStuff,'Location','best');
set(lgnd, 'Interpreter','latex')
hold off;


fileName = 'C:/Users/n8dsa/Documents/OpenFOAM/bl/postProcessing/graphUniform/10000/line.xy';
A_velocity = readtable(fileName,'filetype', 'text');
y = table2array( A_velocity(:,1) );
figure;
plot(table2array( A_velocity(:,3) ) , y )
xlabel("$U_x$",'Interpreter','latex')
ylabel("$Y \ Height$",'Interpreter','latex')

