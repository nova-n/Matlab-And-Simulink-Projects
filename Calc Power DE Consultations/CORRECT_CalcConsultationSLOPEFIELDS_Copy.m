clc;
clear;
%%Defining Mesh
    increments = 200;
    minutes = 60;
    Start = 8;
    meshy = linspace(0,minutes,increments);
    [t,P] = meshgrid(meshy,meshy);
    Graphs = tiledlayout(3,1);
    Graphs.Title.String = sprintf("%f",Start)+"kW/min Start";
    annotation('textbox', [0.75, 0.875, 0.1, 0.1], 'String', "Mesh Size: " + increments + " X " + increments)
    maxY = 20;
    textOffSet = 0.3;
%%Explanation:
    %We are only given the rate of the generated power.
    %Integrating once would give you the power output at each point
    %Integrating twice would give the total power generated after an hour
    
%%Rate at which power is generated:
    %Calculating
        PowerRate = (P.*(P-5).^2)./(50+P.*t);
        L = 1./(sqrt(1+PowerRate.^2));
        M = PowerRate./(sqrt(1+PowerRate.^2));    
    %Plotting
        nexttile
        quiver(t,P,L,M);
        axis on
        axis([0,minutes 0,maxY]);
        xlabel("Time in Minutes")
        ylabel("Power per Minute in kW/min")
        hold on
        %SlopeRate = @(t,P)(P*(P-5)^2)/(50+P*t);
        %[t,P] = ode45(SlopeRate,meshy,start);
        %plot(t,P,'b','LineWidth',1.5);
        hold off

%%Power at each instance
    %Plotting
        nexttile
        quiver(t,P,L,M);
        axis on
        axis([0,minutes 0,maxY]);
        xlabel("Time in Minutes")
        ylabel("Power at Each Instant in kW")
        hold on
        SlopeRate = @(t,P)(P*(P-5)^2)/(50+P*t);
        [t,P] = ode45(SlopeRate,meshy,Start);
        plot(t,P,'r','LineWidth',1.5);
        hold off
    %Calculating
        sol = ode45(SlopeRate,meshy,Start);
        PowerSlope = deval(sol,t);
        PowerInstance = [PowerSlope(1)];
        for i = 1:length(t)-1
            PowerInstance(i+1) = PowerInstance(i) + PowerSlope(i+1);
        end
        PowerInstance = PowerInstance/increments;      
       
%%Total Power Up To Time
    %Calculating
        PowerSum = [PowerInstance(1)];
        for ii = 1:length(t)-1
            PowerSum(ii+1) = PowerSum(ii) + PowerInstance(ii+1);
        end
        PowerSum = PowerSum/increments;   
    %Plotting 
        nexttile
        plot(t,PowerSum,'g','LineWidth',1.5)
        axis([0,minutes 0,max(PowerSum, [], 'all')]);
        xlabel("Time in Minutes")
        ylabel("Energy in kWxMin")
        
        text(0.75,max(PowerSum, [], 'all')-textOffSet,"$$ \int$$P(t) = E(t)",'Interpreter','latex','FontSize',12);
        text(0.75,(max(PowerSum, [], 'all')-textOffSet)*0.75,sprintf("= E(%.2f) - E(0)",minutes),'Interpreter','latex','FontSize',12);
        text(0.75,(max(PowerSum, [], 'all')-textOffSet)*0.5,sprintf("= %f - %f",PowerSum(end),PowerSum(1)),'Interpreter','latex','FontSize',12);
        text(0.75,(max(PowerSum, [], 'all')-textOffSet)*0.25,sprintf("= %fkW Minutes = %fkWH | Total Energy Released In %.2f minutes",(PowerSum(end) - PowerSum(1)),(PowerSum(end) - PowerSum(1))/60,minutes),'Interpreter','latex','FontSize',12);
        hold off;
 
        
        
        
        
        
        
        
        
        
        
    
