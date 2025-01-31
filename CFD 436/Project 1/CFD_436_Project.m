clc;
clear;
close all;

Re = [400,1000,3200];
meshSizes = 20:20:200;

convergenceCases = [[177,409,697,1078,1078,2091,2707,3391,4135,4936];
                    [314,785,984,1075,1075,1884,2354,2873,3442,4063];
                    [269,960,1536,1940,1940,2046,2322,2882,3492,4125]] ;

legendStuff = [];
for i = 1:length(meshSizes)
    meshSizes(i);
    legendStuff = [legendStuff,strcat("$ Mesh \ Size \ = \ ", num2str( meshSizes(i) ) , " \ \times \ ", num2str( meshSizes(i) ), " $" ) ];
end


%%Problem 1
A_xy = {}; % u(y) along center X
B_xy = {}; % u(x) along center Y
fileNamesXY = {};
%fileNames{1} ='C:/Users/n8dsa/Documents/OpenFOAM/Project_1_Cavity/Re_400/X_CenterLine/Mesh_20/postProcessing/graphUniform/177/line.xy';
%A = importdata(fileNames{1}); %new way to read files in matlab, imports it as an object

for i = 1:length(Re)
    for ii = 1:length(meshSizes)
        fileNamesXY{1,i*ii} = strcat(['C:/Users/n8dsa/Documents/OpenFOAM/Project_1_Cavity/' ...
            'Re_', num2str( Re(i) ) ,'/X_CenterLine/Mesh_' , num2str(meshSizes(ii)) , '/postProcessing/graphUniform/' ...
            num2str( convergenceCases(i,ii) ),'/line.xy']);
        %A{i,ii} = importdata(fileNames{1,i*ii}); %%new way to read files in matlab, imports it as an object
        A_xy{i,ii} = readtable(fileNamesXY{1,i*ii},'filetype', 'text');

        fileNamesXY{2,i*ii} = strcat(['C:/Users/n8dsa/Documents/OpenFOAM/Project_1_Cavity/' ...
            'Re_', num2str( Re(i) ) ,'/Y_CenterLine/Mesh_' , num2str(meshSizes(ii)) , '/postProcessing/graphUniform/' ...
            num2str( convergenceCases(i,ii) ),'/line.xy']);
        %B{i,ii} = importdata(fileNames{2,i*ii}); %%new way to read files in matlab, imports it as an object
        B_xy{i,ii} = readtable(fileNamesXY{2,i*ii},'filetype', 'text');
    end
end


for i = 1:length(Re)
    FigH = figure('Position', get(0, 'Screensize'),'Color',[1 1 1]);
    tiledlayout(1,2);

    nexttile;
    %plotting y vs u_x along centerline x = 0.5
    for ii = 1:length(meshSizes)
        % y = A{i,ii}.data (:,1);
        % w = A{i,ii}.data (:,3);
        y = table2array( A_xy{i,ii}(:,1) );
        w = table2array( A_xy{i,ii}(:,3) );

        plot(w,y)
        hold on;
    end

    title(strcat("$Velocity \ Profiles \ of \ Open \ Lid \ Cavity \ Flow \ \Re = ",num2str( Re(i) ), " \ at \ Center \ Line \ x = 0.5 $"),'Interpreter','latex')
    xlabel("$X \  Velocity \ (m/s) $",'Interpreter','latex')
    ylabel("$Y \  Distance \ (m)$",'Interpreter','latex') 
    lgnd= legend(legendStuff,'Location','best');
    set(lgnd, 'Interpreter','latex')
    hold off;

    nexttile;
    %plotting u_y vs x along centerline y = 0.5
    for ii = 1:length(meshSizes)
        %y = B{i,ii}.data (:,1);
        %w = B{i,ii}.data (:,4);
        y = table2array( B_xy{i,ii}(:,1) );
        w = table2array( B_xy{i,ii}(:,4) );

        plot(y,w)
        hold on;
    end
   
    ylim([-1,1]);
    title(strcat("$Velocity \ Profiles \ of \ Open \ Lid \ Cavity \ Flow \ \Re = ",num2str( Re(i) ), " \ at \ Center \ Line \ y = 0.5 $"),'Interpreter','latex')
    ylabel("$Y \  Velocity \ (m/s) $",'Interpreter','latex')
    xlabel("$X \  Distance \ (m)$",'Interpreter','latex') 
    lgnd= legend(legendStuff,'Location','best');
    set(lgnd, 'Interpreter','latex')
    hold off;

    set(gcf, 'InvertHardCopy', 'off');
    F = getframe(FigH);
    imwrite(F.cdata, strcat('centerline_Re_' , num2str( Re(i) ) , '.png'), 'png')
end


%%Problem 2
legendStuff = ["$Pressure, \  P $", "$ U_x $", "$ U_y $"];

A_resid = {}; % u(y) along center X
B_resid = {}; % u(x) along center Y
for i = 1:length(Re)  
    fileNamesResid{1,i} = strcat(['C:/Users/n8dsa/Documents/OpenFOAM/Project_1_Cavity/' ...
        'Re_', num2str( Re(i) ) ,'/X_CenterLine/Mesh_' , num2str(meshSizes(end)) , '/postProcessing/residuals/' ...
        num2str( 0 ),'/residuals.dat']);
    %A{i,ii} = importdata(fileNames{1,i*ii}); %%new way to read files in matlab, imports it as an object
    A_resid{i} = readtable(fileNamesResid{1,i},'filetype', 'text');

    fileNamesResid{2,i} = strcat(['C:/Users/n8dsa/Documents/OpenFOAM/Project_1_Cavity/' ...
        'Re_', num2str( Re(i) ) ,'/Y_CenterLine/Mesh_' , num2str(meshSizes(end)) ,'/postProcessing/residuals/' ...
        num2str( 0 ),'/residuals.dat']);
    %B{i,ii} = importdata(fileNames{2,i*ii}); %%new way to read files in matlab, imports it as an object
    B_resid{i} = readtable(fileNamesResid{2,i},'filetype', 'text');
end

for i = 1:length(Re)
    FigH = figure('Position', get(0, 'Screensize'),'Color',[1 1 1]);
    tiledlayout(2,1);

    nexttile;
    %plotting residuals for y vs u_x along centerline x = 0.5
    y = table2array( A_resid{i}(:,1) );
    for ii = 2:width(A_resid{1,1})-1 %trying to not plot the N/A column from the vorticity residual.        
        plot(y, table2array( A_resid{i}(:,ii) ) )
        hold on;
    end

    title(strcat("$Residuals \ of \ Open \ Lid \ Cavity \ Flow \ \Re = ",num2str( Re(i) ), " \ at \ Center \ Line \ x = 0.5 $"),'Interpreter','latex')
    xlabel("$Itterations$",'Interpreter','latex')
    lgnd= legend(legendStuff,'Location','best');
    set(lgnd, 'Interpreter','latex')
    hold off;

    nexttile;
    %plotting u_y vs x along centerline y = 0.5
    y = table2array( B_resid{i}(:,1) );
    for ii = 2:width(A_resid{1,1})-1 %trying to not plot the N/A column from the vorticity residual.
        plot(y, table2array( B_resid{i}(:,ii) ) )
        hold on;
    end
   
    title(strcat("$Residuals \ of \ Open \ Lid \ Cavity \ Flow \ \Re = ",num2str( Re(i) ), " \ at \ Center \ Line \ y = 0.5 $"),'Interpreter','latex')
    xlabel("$Itterations$",'Interpreter','latex')
    lgnd= legend(legendStuff,'Location','best');
    set(lgnd, 'Interpreter','latex')
    hold off;

    set(gcf, 'InvertHardCopy', 'off');
    F = getframe(FigH);
    imwrite(F.cdata, strcat('residuals_Re_' , num2str( Re(i) ) , '.png'), 'png')
end

%%Problem 3
legendStuff = ["$Paper$", "$OpenFOAM$"];
paperData_uy = readtable('C:/Users/n8dsa/OneDrive/Documents/MATLAB/TestCaseOpenLidCavity.xlsx', 'Range', 'A3:D19');
paperData_ux = readtable('C:/Users/n8dsa/OneDrive/Documents/MATLAB/TestCaseOpenLidCavity.xlsx', 'Range', 'F3:I19');
for i = 1:length(Re)  
    FigH = figure('Position', get(0, 'Screensize'),'Color',[1 1 1]);
    tiledlayout(1,2);

    nexttile;

    plot( table2array( A_xy{i,end}(:,3) ) ,table2array( A_xy{i,end}(:,1) ) )
    hold on;
    plot(table2array( paperData_uy(:,1+i) ),table2array( paperData_uy(:,1) ) )

    title(strcat("$U(y) \ Profiles \ OpenFOAM \ vs \ Paper \Re = ",num2str( Re(i) ), " \ at \ Center \ Line \ y = 0.5 $"),'Interpreter','latex')
    xlabel("$X \  Velocity \ (m/s) $",'Interpreter','latex')
    ylabel("$Y \  Distance \ (m)$",'Interpreter','latex') 
    lgnd= legend(legendStuff,'Location','best');
    set(lgnd, 'Interpreter','latex')
    hold off;

    nexttile;
    
    plot( table2array( B_xy{i,end}(:,1) ) ,table2array( B_xy{i,end}(:,4) ) )
    hold on;
    plot(table2array( paperData_ux(:,1) ),table2array( paperData_ux(:,1+1) ) )

    title(strcat("$V(x) \ Profiles \ OpenFOAM \ vs \ Paper \Re = ",num2str( Re(i) ), " \ at \ Center \ Line \ x = 0.5 $"),'Interpreter','latex')
    xlabel("$X \  Distance \ (m) $",'Interpreter','latex')
    ylabel("$Y \  Velocity \ (m/s)$",'Interpreter','latex') 
    lgnd= legend(legendStuff,'Location','best');
    set(lgnd, 'Interpreter','latex')
    hold off;

    set(gcf, 'InvertHardCopy', 'off');
    F = getframe(FigH);
    imwrite(F.cdata, strcat('schemesComparison_Re_' , num2str( Re(i) ) , '.png'), 'png')
end

%%Probelm 4
crankNicConvergenceCases = [4936,4063,4812];
eulerConvergenceCases = [4936,4063,4812];
A_xy_otherCases = {}; % u(y) along center X ,top row is crank, bottom is euler
B_xy_otherCases = {}; % u(x) along center Y
fileNamesXY_crank = {};
fileNamesXY_euler = {};
for i = 1:length(Re)
    fileNamesXY_crank{1,i} = strcat(['C:/Users/n8dsa/Documents/OpenFOAM/Project_1_Cavity/' ...
        'Re_', num2str( Re(i) ) ,'/X_CenterLine/Mesh_' , num2str(meshSizes(end)) , '-Crank_Nicolson/postProcessing/graphUniform/' ...
        num2str( crankNicConvergenceCases(i) ),'/line.xy']);
    fileNamesXY_euler{1,i} = strcat(['C:/Users/n8dsa/Documents/OpenFOAM/Project_1_Cavity/' ...
        'Re_', num2str( Re(i) ) ,'/X_CenterLine/Mesh_' , num2str(meshSizes(end)) , '-Euler/postProcessing/graphUniform/' ...
        num2str( eulerConvergenceCases(i) ),'/line.xy']);
    A_xy_otherCases{1,i} = readtable(fileNamesXY_crank{1,i},'filetype', 'text');
    A_xy_otherCases{2,i} = readtable(fileNamesXY_euler{1,i},'filetype', 'text');

    fileNamesXY_crank{2,i} = strcat(['C:/Users/n8dsa/Documents/OpenFOAM/Project_1_Cavity/' ...
        'Re_', num2str( Re(i) ) ,'/Y_CenterLine/Mesh_' , num2str(meshSizes(end)) , '-Crank_Nicolson/postProcessing/graphUniform/' ...
        num2str( crankNicConvergenceCases(i) ),'/line.xy']);
    fileNamesXY_euler{2,i} = strcat(['C:/Users/n8dsa/Documents/OpenFOAM/Project_1_Cavity/' ...
        'Re_', num2str( Re(i) ) ,'/Y_CenterLine/Mesh_' , num2str(meshSizes(end)) , '-Euler/postProcessing/graphUniform/' ...
        num2str( eulerConvergenceCases(i) ),'/line.xy']);
    B_xy_otherCases{1,i} = readtable(fileNamesXY_crank{2,i},'filetype', 'text');
    B_xy_otherCases{2,i} = readtable(fileNamesXY_euler{2,i},'filetype', 'text');
end

legendStuff = ["$Steady$","$Crank-Nicolson$","$Euler$","$Paper$"];

for i = 1:length(Re)  
    FigH = figure('Position', get(0, 'Screensize'),'Color',[1 1 1]);
    tiledlayout(1,2);

    nexttile;

    plot( table2array( A_xy{i,end}(:,3) ) ,table2array( A_xy{i,end}(:,1) ) )
    hold on;
    plot( table2array( A_xy_otherCases{1,i}(:,3) ) ,table2array( A_xy_otherCases{1,i}(:,1) ) )
    hold on;
    plot( table2array( A_xy_otherCases{2,i}(:,3) ) ,table2array( A_xy_otherCases{2,i}(:,1) ) )
    hold on;
    plot(table2array( paperData_uy(:,1+i) ),table2array( paperData_uy(:,1) ) )

    title(strcat("$U(y) \ Profiles \ of \ Different \ \frac{\partial}{\partial t} \ Discretization \ Schemes  \Re = ",num2str( Re(i) ), " \ at \ Center \ Line \ y = 0.5 $"),'Interpreter','latex')
    xlabel("$X \  Velocity \ (m/s) $",'Interpreter','latex')
    ylabel("$Y \  Distance \ (m)$",'Interpreter','latex') 
    lgnd= legend(legendStuff,'Location','best');
    set(lgnd, 'Interpreter','latex')
    hold off;

    nexttile;
    
    plot( table2array( B_xy{i,end}(:,1) ) ,table2array( B_xy{i,end}(:,4) ) )
     hold on;
    plot( table2array( B_xy_otherCases{1,i}(:,1) ) ,table2array( B_xy_otherCases{1,i}(:,4) ) )
    hold on;
    plot( table2array( B_xy_otherCases{2,i}(:,1) ) ,table2array( B_xy_otherCases{2,i}(:,4) ) )
    hold on;
    plot(table2array( paperData_ux(:,1) ),table2array( paperData_ux(:,1+1) ) )

    title(strcat("$V(x) \ Profiles \ of \ Different \ \frac{\partial}{\partial t} \ Discretization \ Schemes \Re = ",num2str( Re(i) ), " \ at \ Center \ Line \ x = 0.5 $"),'Interpreter','latex')
    xlabel("$X \  Distance \ (m) $",'Interpreter','latex')
    ylabel("$Y \  Velocity \ (m/s)$",'Interpreter','latex') 
    lgnd= legend(legendStuff,'Location','best');
    set(lgnd, 'Interpreter','latex')
    hold off;

    set(gcf, 'InvertHardCopy', 'off');
    F = getframe(FigH);
    imwrite(F.cdata, strcat('paperComparisons_Re_' , num2str( Re(i) ) , '.png'), 'png')
end

%%Problem 5
A_cutPlane= {}; % top row is steady, middle is crank, bottom is euler, each column is per Re
fileNames_cutPlane = {}; % top row is steady, middle is crank, bottom is euler

for i = 1:length(Re)
    fileNames_cutPlane{1,i} = strcat(['C:/Users/n8dsa/Documents/OpenFOAM/Project_1_Cavity/' ...
        'Re_', num2str( Re(i) ) ,'/X_CenterLine/Mesh_' , num2str(meshSizes(end)) , '/postProcessing/cutPlaneSurface/' ...
        num2str( convergenceCases(i,end) ),'/cutPlane.xy']);
    fileNames_cutPlane{2,i} = strcat(['C:/Users/n8dsa/Documents/OpenFOAM/Project_1_Cavity/' ...
        'Re_', num2str( Re(i) ) ,'/X_CenterLine/Mesh_' , num2str(meshSizes(end)) , '-Crank_Nicolson/postProcessing/cutPlaneSurface/' ...
        num2str( crankNicConvergenceCases(i) ),'/cutPlane.xy']);
    fileNames_cutPlane{3,i} = strcat(['C:/Users/n8dsa/Documents/OpenFOAM/Project_1_Cavity/' ...
        'Re_', num2str( Re(i) ) ,'/X_CenterLine/Mesh_' , num2str(meshSizes(end)) , '-Euler/postProcessing/cutPlaneSurface/' ...
        num2str( crankNicConvergenceCases(i) ),'/cutPlane.xy']);
    A_cutPlane{1,i} = readtable(fileNames_cutPlane{1,i},'filetype', 'text');
    A_cutPlane{2,i} = readtable(fileNames_cutPlane{2,i},'filetype', 'text');
    A_cutPlane{3,i} = readtable(fileNames_cutPlane{3,i},'filetype', 'text');
end

% Map the 2D data onto a meshgrid x-y plane
Nx = 2000;
Ny = 2000;
titles = ["Steady \ State" , "Crank-Nicolson", "Euler"];
for i = 1:3
    FigH = figure('Position', get(0, 'Screensize'),'Color',[1 1 1]);
    T = tiledlayout(1,3);
    for ii = 1:length(Re)
        x = table2array( A_cutPlane{i,ii}(:,1) );
        y = table2array( A_cutPlane{i,ii}(:,2) );
        vorticity = table2array( A_cutPlane{i,ii}(:,9) );

        xMin = min(x);
        xMax = max(x);
        yMin = min(y);
        yMax = max(y);
        vorticityLowest = -1 * max( abs( rmoutliers( vorticity( find(vorticity < 0) ) ) ) ); 
        vorticityHighest = max(  rmoutliers( vorticity( find(vorticity > 0) ) )  );
        vorticityContourLevels = linspace(vorticityLowest,vorticityHighest,50);

        X = linspace(xMin,xMax,Nx);
        Y = linspace(yMin,yMax,Ny);

        [xq,yq] = meshgrid(X,Y);
        vorticity_q = griddata(x,y,vorticity,xq,yq);
        nexttile;
        contourf(X,Y,vorticity_q, 'LevelList', vorticityContourLevels) %isocontour of vorticity
        hold off;
        title(strcat("$ \Re = ",num2str( Re(ii) ) ," $"),'Interpreter','latex')
        xlabel("$X \  Distance \ (m) $",'Interpreter','latex')
        ylabel("$Y \  Distance \ (m) $",'Interpreter','latex')
    end
    title(T,strcat("$Vorticity \ Isocontours \ Using \ " , titles(i) , " \ \frac{\partial}{\partial t} \ Discretization $"),'Interpreter','latex')
    col_bar = colorbar('Orientation', 'Horizontal');
    col_bar.Layout.Tile = 'north';

    F = getframe(FigH);
    if i == 1
        titles(i) = "Steady  State";
    end
    set(gcf, 'InvertHardCopy', 'off');
    imwrite(F.cdata, strcat('vorticities_Re_' , titles(i) , '.png'), 'png');
end
