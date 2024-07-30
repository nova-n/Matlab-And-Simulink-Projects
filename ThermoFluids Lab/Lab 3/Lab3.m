clc;
clear;
close all;
densityWater = 998;%kg/m^3
dynVisWater = 0.0010016	;%at room temp!
densityAir = 1.204;%kg/m^3 at room temp!!!
dynVisAir = 1.81*10^-5;%at room temp!
gravity = 9.81;%m/s^2
Lmintom3 = 1/60000;
mmtom = 1/1000;
bartoPa = 100000;
resolution = 100;

smallDiam = 13.7*mmtom;%in m
largeDiam = 26.4*mmtom;%in m
straightLengthA = 914.4*mmtom;%in m
straightLengthL = 901.7*mmtom;%in m

DBRange = 'A4:I12';
LBRange = 'A4:M12';

labData = {readtable("EGME306B_Ex03_DataSheet_NEW",'Sheet',"DB_heights",'Range',DBRange),readtable("EGME306B_Ex03_DataSheet_NEW",'Sheet',"LB_heights",'Range',LBRange)};
labData{1,1}{:,2:7} = labData{1,1}{:,2:7}*mmtom;%converting manometer heights to m
labData{1,2}{:,2:11} = labData{1,2}{:,2:11}*mmtom;
labData{1,1}{:,8} = labData{1,1}{:,8}*bartoPa;%converting pressure to Pa, from pressure transducer in tap 19
labData{1,2}{:,12} = labData{1,2}{:,12}*bartoPa;%converting pressure to Pa, from pressure transducer in tap 19
labData{1,1}{:,9} = labData{1,1}{:,9}*Lmintom3;%converted flowrates to m^3/s
labData{1,2}{:,13} = labData{1,2}{:,13}*Lmintom3;%converted flowrates to m^3/s\

flowRatesDB = labData{1,1}{:,9}; %dark blue flowrates, converted to m^3/s
flowRatesLB = labData{1,2}{:,13}; %light blue flowrates, converted to m^3/s
flowRateCells = [[flowRatesDB]' ; [flowRatesLB]'];

Re_DB = (4*densityWater)/(pi*dynVisWater*smallDiam) * flowRatesDB; %thru straight section of intereset, section A, had small diam
Re_LB = (4*densityWater)/(pi*dynVisWater*largeDiam) * flowRatesLB; %thru straight section of intereset, section L, had large diam
%Seems like all flows are laminar

%%Getting Friction Factor of last stretch of pipe for each system
[frictionFactorsDB] = frictionFactorSolver(labData{1,1},smallDiam,flowRatesDB,straightLengthA,labData{1,1}{:,4},labData{1,1}{:,5},gravity);
[frictionFactorsLB] = frictionFactorSolver(labData{1,2},largeDiam,flowRatesLB,straightLengthL,labData{1,2}{:,3},labData{1,2}{:,4},gravity);
fricFactorCells = [frictionFactorsDB ; frictionFactorsLB];

%combining both data to get a best fit line
Re_combined = cat(1,Re_LB,Re_DB);
frictionFactors_combined = cat(1,[frictionFactorsLB]',[frictionFactorsDB]');
guessesC0 = [0.001,1.75,-0.4];
[c,~] = lsqcurvefit(@Nathan_Delos_Santos_power_func, guessesC0,Re_combined,frictionFactors_combined);
Re = linspace(min(Re_combined),max(Re_combined),resolution);
Re_forBFL =linspace(1500,35000,resolution);
frictionFacs = c(1)+c(2)*Re_forBFL.^c(3);


legendStuff = ["$Dark \ Blue \ Pipe$","$Light \ Blue \ Pipe$","$Best \ Fit \ Line$"];
figure;
plot(Re_DB,frictionFactorsDB, 'ob-');
hold on;
plot(Re_LB,frictionFactorsLB, 'or-');
hold on;
plot(Re_forBFL,frictionFacs,'k--');
hold off;
lgnd= legend(legendStuff);
set(lgnd, 'Interpreter','latex')
xlabel("$ \Re$",'Interpreter','latex') 
ylabel("$Friction \ Factor, f$",'Interpreter','latex') 
title("$ \Re \ vs \ Friction \ Factor \ For \ Both \ Pipes $",'Interpreter','latex')

plot1 = "Friction Factor vs Re for last stretch";
print('-r600','-dpng',plot1);

%%Finding minor loss coefficients (for each flow rate) of...sudden expansion, K_exp, in the light blue pipe
%sudden expansion between taps 7 & 8, flowing from 7 to 8
K_exp = [];
for i = 1:height(labData{1,2})
    K_exp(i) = (gravity*pi^2*smallDiam^4)*(labData{1,2}{i,2} - labData{1,2}{i,3})/(8*flowRatesLB(i)^2) +1 - (smallDiam/largeDiam)^4;
end

%%Finding minor loss coefficients (for each flow rate) of...sudden contraction, K_cont, in the light blue pipe
%sudden contraction between taps 9 & 10, flowing from 9 to 10
K_cont = [];
for i = 1:height(labData{1,2})
    K_cont(i) = (gravity*pi^2*smallDiam^4)*(labData{1,2}{i,4} - labData{1,2}{i,5})/(8*flowRatesLB(i)^2) - 1 + (smallDiam/largeDiam)^4;
end

%%Plotting Loss Coefficients vs Reynolds Numbers
legendStuff = ["$K_{expansion}$","$K_{contraction}$"];
figure;
plot(Re_LB,K_exp, 'ob-');
hold on;
plot(Re_LB,K_cont, 'or-');
hold off;
lgnd= legend(legendStuff);
set(lgnd, 'Interpreter','latex')
xlabel("$ \Re$",'Interpreter','latex') 
ylabel("$Minor \ Loss \ Coefficients$",'Interpreter','latex') 
title("$ \Re \ vs \ Minor \ Loss \ Coefficients \ For \ Both \ Pipes$",'Interpreter','latex')
plot1 = "Minor Loss Coefficients vs Re";
print('-r600','-dpng',plot1);

figure;
plot(Re_LB,K_exp, 'ob-');
hold off;
xlabel("$ \Re$",'Interpreter','latex') 
ylabel("$Minor \ Loss \ Coefficient \ of \ Sudden \ Expansion$",'Interpreter','latex') 
title("$ \Re \ vs \ Minor \ Loss \ Coefficients \ of \ Sudden \ Expansion$",'Interpreter','latex')
plot1 = "Minor Loss Coefficient of Sudden Expansion vs Re";
print('-r600','-dpng',plot1);

figure;
plot(Re_LB,K_cont, 'or-');
hold off;
xlabel("$ \Re$",'Interpreter','latex') 
ylabel("$Minor \ Loss \ Coefficient \ of \ Sudden \ Contraction$",'Interpreter','latex') 
title("$ \Re \ vs \ Minor \ Loss \ Coefficients \ of \ Sudden \ Contraction$",'Interpreter','latex')
plot1 = "Minor Loss Coefficient of Sudden Contraction vs Re";
print('-r600','-dpng',plot1);


%%Finding average minor loss coefficients, avg_K_bend (per flowrate) of bend vs R/D for both pipoes
R_B = 0; R_C = 1.5*smallDiam; R_G = 50.8; R_H = 101.6; R_J = 152.5; %Radii of bends in mm
%R_D_ratio_DB = mmtom*[R_B,R_C]/smallDiam;
%R_D_ratio_LB = mmtom*[R_G,R_H,R_J]/smallDiam;
R_D_ratio = mmtom*[R_B,R_C,R_G,R_H,R_J]/smallDiam;

%lengths of straight pipe between manometers
L_DB = [0.94,0.93];
DB_bend_tapNums = [5,6 ; 1,2];% [B,C]
DB_tap_heightDiffs = [R_B + 0, R_C + 812.8]*mmtom; %height between 5&6 is same, since is 0 bend radius

L_LB = [0.92,0.935,0.88;];
LB_bend_tapNums = [9,10 ; 5,6 ; 7,8];% [G,H,J], manometerNumbers are actually [15,16 ; 11,12 ; 13,14], but is offset in excel sheet
LB_tap_heightDiffs = [-1*(R_G + 21.59), R_H + 723.9, R_J + 21.59]*mmtom;%any height going up is negative

DB_K_b = bendCoeffSolver(L_DB, labData{1,1}, DB_bend_tapNums, ...
    DB_tap_heightDiffs, flowRatesDB, densityWater, dynVisWater, smallDiam, gravity);

LB_K_b = bendCoeffSolver(L_LB, labData{1,2}, LB_bend_tapNums, ...
    LB_tap_heightDiffs, flowRatesLB, densityWater, dynVisWater, smallDiam, gravity);

avgK_b = [mean(DB_K_b(:,1)),mean(DB_K_b(:,2)),mean(LB_K_b(:,1)),mean(LB_K_b(:,2)),mean(LB_K_b(:,3))]; %[dark blue averages: 2 columns, light blue averages: 3 columns]

figure;
plot(R_D_ratio,avgK_b);
hold off;
xlabel("$ Bend Radius/Pipe Diameter, \ r/D$",'Interpreter','latex') 
ylabel("$Average \ Bend \ Loss \ Coefficients, K_{bend_{average}}$",'Interpreter','latex') 
title("$ r/D \ vs \ Average \ Bend \ Loss \ Coefficients \ For \ Both \ Pipes $",'Interpreter','latex')
plot1 = "Bend Loss Coefficients vs r to D ratio";
print('-r600','-dpng',plot1);

%%Valve Losses
Re_DB_Valves = (4*densityWater)/(pi*dynVisWater*smallDiam) * flowRatesDB; %thru straight section of intereset, section A, had small diam
Re_LB_Valves = (4*densityWater)/(pi*dynVisWater*smallDiam) * flowRatesLB; %thru straight section of intereset, section L, had large diam
%taps 19 and 20 measured pressure

%Dark Blue had gate valve, Light Blue had globe valve
K_Gate_Valve = ((pi^2 * smallDiam^4)./(8*densityWater*flowRatesDB.^2)).*labData{1,1}{:,8};
K_Globe_Valve = ((pi^2 * smallDiam^4)./(8*densityWater*flowRatesLB.^2)).*labData{1,2}{:,12}; 

figure;
loglog(Re_DB_Valves,K_Gate_Valve,'ob-')
hold on;
loglog(Re_LB_Valves,K_Globe_Valve,'rv-')
hold off;
legendStuff = ["$Gate \ Valve \ Losses, K_{Gate}$","$Globe \ Valve \ Losses, K_{Globe}$"];
lgnd= legend(legendStuff);
set(lgnd, 'Interpreter','latex')
xlabel("$ \Re$",'Interpreter','latex') 
ylabel("$Valve \ Loss \ Coefficients$",'Interpreter','latex') 
title("$ \Re \ vs \ Valve \ Loss \ Coefficients \ For \ Both \ Pipes$",'Interpreter','latex')
plot1 = "Valve Loss Coefficients vs Re";
print('-r600','-dpng',plot1);

mean(K_exp(3:9))
mean(K_cont(3:9))

function [frictionFactorArray] = frictionFactorSolver(labSheet,diam,flowRates,sectionLength,manometerHeight1,manometerHeight2,grav)
    for i = 1:height(labSheet)
     frictionFactorArray(i) = ( (pi^2*grav * diam^5)/(8*flowRates(i).^2 * sectionLength) ) * (manometerHeight1(i) - manometerHeight2(i));   
    end
end

function [K_bend] = bendCoeffSolver(straightLengths, labSheet, tapNums, ...
    heightDiffs, flowRates, fluidDensity, fluidDynVis, pipeDiam, grav)

    K_bend = [];
    
    %note, the eqn takes into account the straight lengths (via the majorlosses)
    %because we didn't actually measure purely across the bend. We used
    %the manometers before and after the bends, and in between those two
    %readings, were the straight sections, as well as the bends
    
    
    %the (h_1 - H-2 + z_1 - z_2) part has to be looped. Will store in a variable
    %remember, z_1-z_2 is already stored as the tap_heightDiffs variable
    pressureAndHeights = [];
    for ii = 1:length(straightLengths)
            %each column is the individual elblows
            %each row is per flowrate

            pressureAndHeights(:,ii) = (labSheet{:,tapNums(ii,1) + 1} - labSheet{:,tapNums(ii,2) + 1}...
                + heightDiffs(ii));
            %put +1, since first column is for flowrates. All manometer columns were shifted right
    end
    pressureAndHeights;
    
    for i = 1:length(flowRates)
        %manometer readings are h_1 and h_2 from eqn in worksheet
        %z_1 - z_2 come from the fact that monometers measure at different
        %heights, and is from tap_heightDiffs variables
    
        systemRe=4*fluidDensity*flowRates(i)/(pi*fluidDynVis*pipeDiam);
        straightSection_f = 0.0785/(systemRe^0.25); % is from blasius approximation
        
        %Each row is per flowrate, and each column is k_b for all elbows measured
        K_bend(i,:) = (grav*pi^2*pipeDiam^4)/(8*flowRates(i)^2)*pressureAndHeights(i)*straightSection_f.*straightLengths./pipeDiam;
    end
end

