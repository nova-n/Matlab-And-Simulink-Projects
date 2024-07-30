clear;
clc;
close all;

resolution = 100;

%%Lab Data
L_beam = 32/12;%inches to ft
b_Beam = 1/12;%feet
h_Beam = 0.5/12;%feet
E_beam = 30*10^6 * 144;%psi to psf
I_beam = (b_Beam*h_Beam^3)/12; %ft^4
m_beam = (4 + 9.94/16)/32.2;%slugs
m_motorAssem = (10 + 4.34/16)/32.2;%slugs


range = 'A2:C10';
labData = table2array(readtable("Lab Data.xlsx",'Range',range));
%columns:
% Additional Weights | Tachometer Reading (of entire system)  | FFT (of entire system)

labData(:,1) = labData(:,1)/32.2 + m_motorAssem;%converting lbm to slugs, and adding the motor assembly mass
labData(:,2) = labData(:,2)*(pi/30);%converting rpm to rad/s
labData(:,3) = labData(:,3)*(2*pi);%converting Hz to rad/s

%%From Theoretical Calculations
%See derivation in report
m_beamConcentrated_theoretical = 0.4857*m_beam
k_equivalent_theoretical = 48*E_beam*I_beam/L_beam^3 %lbf/ft
%w_n_theoretical = sqrt(k_equivalent_theoretical/m_beamConcentrated) %rad/s
w_n_theoretical = (pi^2)*sqrt((E_beam*I_beam)/(m_beam*L_beam^3))

%%From Tachometer
[k_equivalent_tach,m_beamConcentrated_tach,w_n_tach] = ...
    dunkerleySolver([labData(:,1),labData(:,2)],resolution,"Tachometer")

%%From FFT
[k_equivalent_FFT,m_beamConcentrated_FFT,w_n_FFT] = ...
    dunkerleySolver([labData(:,1),labData(:,3)],resolution,"FFT")

%%Tabulating Results Data
results = [ [k_equivalent_theoretical;k_equivalent_tach;k_equivalent_FFT] , ...
    [m_beamConcentrated_theoretical;m_beamConcentrated_tach;m_beamConcentrated_FFT] , ...
    [w_n_theoretical;w_n_tach;w_n_FFT]];
array2table(results, 'VariableNames', {'k_equivalent', 'm_concentratedBeam', 'w_nBeam'} , ...
    'RowNames',["Theoretical","Tachometer","FFT"])

percentError = [ [abs( (k_equivalent_tach - k_equivalent_theoretical)/k_equivalent_theoretical) ; ...
    abs( (k_equivalent_FFT - k_equivalent_theoretical)/k_equivalent_theoretical)] ...
    , [abs( (m_beamConcentrated_tach - m_beamConcentrated_theoretical)/m_beamConcentrated_theoretical) ; ...
    abs( (m_beamConcentrated_FFT - m_beamConcentrated_theoretical)/m_beamConcentrated_theoretical)] ...
    , [abs( (w_n_tach - w_n_theoretical)/w_n_theoretical) ; abs( (w_n_FFT - w_n_theoretical)/w_n_theoretical)] ];
percentError = percentError * 100;

array2table(percentError, 'VariableNames', {'k_equivalent', 'm_concentratedBeam', 'w_nBeam'} , ...
    'RowNames',["Tachometer","FFT"])

function [k_equivalent,m_beamConcentrated,w_n] = dunkerleySolver(dataIn,res,titleStuff)    
    lineEqn = polyfit(dataIn(:,1), 1./(dataIn(:,2).^2),1);    
    M = linspace(0,dataIn(end,1),res);%from 0lbs to 16lbs 
    %solving for x intercept (will be negative), and is negative m_beamConcentrated
    m_beamConcentrated = lineEqn(2)/lineEqn(1);
    M = linspace(-m_beamConcentrated,0,res) + M;
    %solving for spring constant
    k_equivalent = 1/lineEqn(1);
    %solving for natural frequncy of beam itself
    w_n = sqrt(k_equivalent/m_beamConcentrated);
    
    %plotting
    figure;
    pointPlot1 = plot(dataIn(:,1), 1./(dataIn(:,2).^2),"." );
    hold on;
    plot(-m_beamConcentrated,0,".",'color',get(pointPlot1,'color'),'HandleVisibility','off','MarkerSize',16);
    plot(-m_beamConcentrated,0,"hexagram",'color',get(pointPlot1,'color'),'MarkerSize',12,'linewidth',2);
    hold on;
    pointPlot2 = plot(M,lineEqn(1)*M + lineEqn(2));
    hold on;
    plot(0,lineEqn(2),".",'color',get(pointPlot2,'color'),'HandleVisibility','off','MarkerSize',16);
    plot(0,lineEqn(2),"hexagram",'color',get(pointPlot2,'color'),'MarkerSize',12,'linewidth',2);
    hold on
    xline(0,"--")
    title(strcat("$Dunkerly \ Line \ Plot \ From \ ", titleStuff , " \ Data $"),'Interpreter','latex')
    xlabel("$ M = m_{MotorAssembly} + m_{AddedWeights}  \ ( slugs ) $",'Interpreter','latex') 
    ylabel("$ \frac{1}{(w_{n,system}) ^2} \ ( \frac{s^2}{rad^2} ) $",'Interpreter','latex',FontSize= 18)
    legendStuff = ["$Lab \ Data$" ,strcat("$ m_{beamConcentrated} = " , num2str(m_beamConcentrated) , " \  (slugs) $") ...
        strcat("$Best \  Fit \ Line, \ k_{equivalent} = " , num2str(k_equivalent) , " \ ( \frac{lb_f}{ft} ) $")...
        strcat("$ w_{n,beam} = " , num2str(w_n) , " \ ( \frac{rad}{s} ) $")];
    lgnd= legend(legendStuff,'Location','best');
    set(lgnd, 'Interpreter','latex')
    hold off;   
    plot1 = strcat("Dunkerly Line Plots For", titleStuff ,  " data");
    print('-r600','-dpng',plot1);
end



