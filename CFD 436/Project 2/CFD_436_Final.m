clc;
clear;
close all;

fileName = 'C:/Users/n8dsa/Documents/OpenFOAM/drivaer/drivaer_Re05_xe05_CaseA1/postProcessing/residuals/0/residuals.dat';
A_resid = readtable(fileName,'filetype', 'text');


%plotting residuals for y vs u_x along centerline x = 0.5
y = table2array( A_resid(:,1) );
for ii = 2:width(A_resid) -1        
    %plot(y, table2array( A_resid(:,ii) ) )
    %hold on;
end

%%Problem 1
Re = 5:5:35;
indecies_forOneDigit_Re = find(Re<10);
Re_string = string(Re);
Re_string(indecies_forOneDigit_Re) = "0" + Re_string(indecies_forOneDigit_Re);
cases = ["A1","B1","B2","B3","C2","C3"];%,"B2","B3","C2","C3"];
Re = Re*10^5;
figure;
for i = 1:length(cases)
    Cd_Cases = zeros(1,length(Re));
    for ii = 1:length(Re)
        fileName = strcat('C:/Users/n8dsa/Documents/OpenFOAM/drivaer/drivaer_Re', Re_string(ii) ,...
            '_xe05_Case',cases(i),'/postProcessing/forceCoeffsIncompressible/0/forceCoeffs.dat');
        forceCoeffsTable = readtable(fileName,'filetype', 'text');
        forceCoeffsTable = forceCoeffsTable(3:end,1:end-1); %Gets rid of the titles, which are NaN
        if cases(i) ~= "B3"
            forceCoeffsTable = 0.5 .* forceCoeffsTable; %Need to double, since did half area only
        end
        %forceCoeffsTable = importdata(fileName).data;
        Cd_Cases(1,ii) = forceCoeffsTable{end,3};
        
    end
    
    plot(Re,Cd_Cases, '.-','MarkerSize',12)
    hold on;
end
title("$ C_d \ For \ All \ Cases \ vs \ \Re $",'Interpreter','latex')
xlabel("$ \Re $",'Interpreter','latex')
ylabel("$ C_d $",'Interpreter','latex') 
xlim([0,Re(end)+5*10^5]);
lgnd= legend("$Case \ " + cases + " $",'Location','best');
set(lgnd, 'Interpreter','latex')
hold off;




%%Problem 2
%{
fileName = 'C:/Users/n8dsa/Documents/OpenFOAM/drivaer/drivaer_Re35_xe05_CaseA1/postProcessing/residuals/0/residuals.dat';
lastItteration = readtable(fileName,'filetype', 'text');
lastItteration = table2array( lastItteration(end,1) ); %returns a number. Just gets the last itteration the sim stopped on
fileName = strcat('C:/Users/n8dsa/Documents/OpenFOAM/drivaer/drivaer_Re35_xe05_CaseA1/postProcessing/cuttingPlane/' ...
    , num2str(lastItteration) ,'/yNormal.xy');
cutPlane = table2array( readtable(fileName,'filetype', 'text') );
n = 6; %every nth row of matrix, remove nn rows just to save space on mesh when plotting
nn = 4;
indeciesToRemove = [];
for i = n:n:height(cutPlane)
    indeciesToRemove = [indeciesToRemove,i:i+nn]; 
end
indeciesToRemove = indeciesToRemove(indeciesToRemove < height(cutPlane) ); %incase it is above the length of cutplane
cutPlane(indeciesToRemove,:) = [];
x = cutPlane(:,1);
z = cutPlane(:,3);
p = cutPlane(:,4);
[X,Z] = meshgrid(x,z);
pField = griddata(x,z,p,X,Z);
clearvars indeciesToRemove x z;
pFieldMin = -1 * max( abs( rmoutliers( pField( find(pField < 0) ) ) ) ); 
pFieldMax = max(pField,[],"all");
%pFieldMin = min(pField,[],"all");
pFieldContourLevels = linspace(pFieldMin,pFieldMax,50);
figure;
%contourf(X,Z,pField , 'LevelList', pFieldContourLevels,'LineColor','none');
%plot3(x,z,p,".")
%hold on
%surf(X,Z,pField,'EdgeColor','none','LineStyle','none')
%shading interp
title("$ Case \ A1 \ Pressure \ Field \ at \ Center \ Y \ Plane $",'Interpreter','latex')
xlabel("$X \  Distance \ (m) $",'Interpreter','latex')
ylabel("$Z \  Distance \ (m) $",'Interpreter','latex')
col_bar = colorbar('southoutside','Orientation', 'Horizontal');
hold off;
%}

figure;
tiledlayout(2,1);
nexttile;
meshLevels = [4,5,6,7,9];
for m = 1:length(meshLevels)
    if meshLevels(m) ~= 9
        fileName = strcat('C:/Users/n8dsa/Documents/OpenFOAM/drivaer/Case_A1_Re05_LOWER_MESHES/drivaer_Re35_xe05_CaseA1_Mesh'...
            ,num2str(meshLevels(m)),'/postProcessing/residuals/0/residuals.dat');
        lastItteration = readtable(fileName,'filetype', 'text');
        lastItteration = table2array( lastItteration(end,1) ); %returns a number. Just gets the last itteration the sim stopped on
        fileName = strcat('C:/Users/n8dsa/Documents/OpenFOAM/drivaer/Case_A1_Re05_LOWER_MESHES/drivaer_Re35_xe05_CaseA1_Mesh'...
            ,num2str(meshLevels(m)),'/postProcessing/wallSampling/', num2str(lastItteration) ,'/xWall.xy');
    else
        fileName = 'C:/Users/n8dsa/Documents/OpenFOAM/drivaer/drivaer_Re35_xe05_CaseA1/postProcessing/residuals/0/residuals.dat';
        lastItteration = readtable(fileName,'filetype', 'text');
        lastItteration = table2array( lastItteration(end,1) ); %returns a number. Just gets the last itteration the sim stopped on
        fileName = strcat('C:/Users/n8dsa/Documents/OpenFOAM/drivaer/drivaer_Re35_xe05_CaseA1/postProcessing/wallSampling/' ...
            , num2str(lastItteration) ,'/xWall.xy');
    end
    cutPlane = table2array( readtable(fileName,'filetype', 'text') );
    
    sortedCut = sort(cutPlane,1);
    %plot(sortedCut(:,1),sortedCut(:,4))
    clearvars sortedCut
    %sort data by type
    U = 16;
        xa = cutPlane(:,1);
        ya = cutPlane(:,2);
        za = cutPlane(:,3);
        pa = cutPlane(:,4)/(0.5*U*U);
    
    %   firter data so that it is within the center (x-z) plane +/- 0.01
        dy = 0.01;
        index = ya<dy;
        xb = xa(index);
        yb = ya(index);
        zb = za(index);
        pb = pa(index);
    
        index = yb>(-1*dy);
    
        xa = xb(index);
        ya = yb(index);
        za = zb(index);
        pa = pb(index);
    
    %   firter data so that it is on the upper surface within the center plance
        index = za>0;
        xb = xa(index);
        yb = ya(index);
        zb = za(index);
        pb = pa(index);
    
    %   sort data by x-axis
        [~,index] = sort(xb);
    
        x = xb(index);
        y = yb(index);
        z = zb(index);
        p = pb(index);
        plot(x,p)
        hold on;
        length(p)
end
clearvars pa pb xa xb ya yb za zb
title("$ Case \ A1 \ \Re \ 35 \times 10^5 \ Pressure \ Along \ Car \ Surface \ per \ X \ Distance \ Per \ Different \ Mesh \ Refinement$",'Interpreter','latex')
xlabel("$X \  Distance \ (m) $",'Interpreter','latex')
ylabel("$P \  Pressure \ (Pa) $",'Interpreter','latex')
lgnd= legend("$Mesh \ Refinement \ Level \  \ " + meshLevels + " $",'Location','best');
set(lgnd, 'Interpreter','latex')
hold off;

nexttile;
height(cutPlane)
n = 40; %every nth row of matrix, remove nn rows just to save space on mesh when plotting
nn = 38;
indeciesToRemove = [];
for i = n:n:height(cutPlane)
    indeciesToRemove = [indeciesToRemove,i:i+nn]; 
end
indeciesToRemove = indeciesToRemove(indeciesToRemove < height(cutPlane) ); %incase it is above the length of cutplane
cutPlane(indeciesToRemove,:) = [];
height(cutPlane)
x = cutPlane(:,1);
z = cutPlane(:,3);
p = cutPlane(:,4);
[X,Z] = meshgrid(x,z);
pField = griddata(x,z,p,X,Z);
clearvars indeciesToRemove x z;
pFieldMin = -1 * max( abs( rmoutliers( pField( find(pField < 0) ) ) ) ); 
pFieldMax = max(pField,[],"all");
%pFieldMin = min(pField,[],"all");
pFieldContourLevels = linspace(pFieldMin,pFieldMax,50);

contourf(X,Z,pField , 'LevelList', pFieldContourLevels,'LineColor','none');
%plot3(x,z,p,".")
%hold on
%surf(X,Z,pField,'EdgeColor','none','LineStyle','none')
%shading interp
title("$ Case \ A1 \  \Re \ 35*10^5 \ Pressure \ Along \ Car \ Surface \ at \ View \ Normal \ to \ Y \ Plane $",'Interpreter','latex')
xlabel("$X \  Distance \ (m) $",'Interpreter','latex')
ylabel("$Z \  Distance \ (m) $",'Interpreter','latex')
col_bar = colorbar('southoutside','Orientation', 'Horizontal');
hold off;

%%Problem 3
Re = 5:5:25;
indecies_forOneDigit_Re = find(Re<10);
Re_string = string(Re);
Re_string(indecies_forOneDigit_Re) = "0" + Re_string(indecies_forOneDigit_Re);
Cd_A_Cases = zeros(1,length(Re));
Re = Re*10^5;

for i = 1:length(Re)
    fileName = strcat('C:/Users/n8dsa/Documents/OpenFOAM/drivaer/drivaer_Re', Re_string(i) ,...
        '_xe05_CaseA1/postProcessing/forceCoeffsIncompressible/0/forceCoeffs.dat');
    forceCoeffsTable = readtable(fileName,'filetype', 'text'); %Need to double, since did half area only
    forceCoeffsTable = forceCoeffsTable(3:end,1:end-1); %Gets rid of the titles, which are NaN
    forceCoeffsTable = 2 .* forceCoeffsTable; %Need to double, since did half area only
    %forceCoeffsTable = importdata(fileName).data;
    Cd_A_Cases(1,i) = forceCoeffsTable{end,3};
    
end

figure;
plot(Re,Cd_A_Cases, '.-','MarkerSize',12)
title(strcat("$ C_d \ For \ Case \ A1 \ vs \ \Re $"),'Interpreter','latex')
xlabel("$ \Re $",'Interpreter','latex')
ylabel("$ C_d $",'Interpreter','latex') 
xlim([0,Re(end)]);
hold off;

%%Problem 4
cases = ["B1","B2","B3"];
Cd_B_Cases = array2table(zeros(1,3),'RowNames',["C_d"]);
Cd_B_Cases.Properties.VariableNames = "Case " + cases + " Re_35xe05";
for i = 1:length(cases)
    fileName = strcat('C:/Users/n8dsa/Documents/OpenFOAM/drivaer/drivaer_Re35_xe05_Case',cases(i),...
        '/postProcessing/forceCoeffsIncompressible/0/forceCoeffs.dat');
    forceCoeffsTable = readtable(fileName,'filetype', 'text');
    forceCoeffsTable = forceCoeffsTable(3:end,:); %Gets rid of the titles, which are NaN
    %forceCoeffsTable = importdata(fileName).data;
    Cd_B_Cases(1,i) = forceCoeffsTable(end,3);
end

%%Problem 5
cases = ["B1","C2","C3"];
Cd_C_Cases = array2table(zeros(1,3),'RowNames',["C_d"]);
Cd_C_Cases.Properties.VariableNames = "Case " + cases + " Re_35xe05";
for i = 1:length(cases)
    fileName = strcat('C:/Users/n8dsa/Documents/OpenFOAM/drivaer/drivaer_Re35_xe05_Case',cases(i),...
        '/postProcessing/forceCoeffsIncompressible/0/forceCoeffs.dat');
    forceCoeffsTable = readtable(fileName,'filetype', 'text');
    forceCoeffsTable = forceCoeffsTable(3:end,:); %Gets rid of the titles, which are NaN
    %forceCoeffsTable = importdata(fileName).data;
    Cd_C_Cases(1,i) = forceCoeffsTable(end,3);
end

