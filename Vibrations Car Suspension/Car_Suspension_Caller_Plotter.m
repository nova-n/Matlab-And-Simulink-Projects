clear ; clc ; close all ;

m = [320,44];
k = [3.2*10^4,1.8*10^5];
c = [3430];

v = [5,10,35] * 5280*12*0.0254/3600;

h = 0.1;
L = 3.7;

resolution = 1500;
tStart =0; % time starts at 0 seconds
tEnd = 10; % time ends at this amount of time
% specify the time vector
time = linspace(tStart, tEnd, resolution);
yTime = linspace(tStart,tEnd+1,resolution);
% initial conditions
x_0 =[0,0]; % in m
x_dot_0 =[0,0]; % in m / s
IC =[ x_0(1) , x_0(2) , x_dot_0(1), x_dot_0(2) ]';


for i = 1:length(v)
    d = 1*v(i);
    y = zeros([1,length(time)]);
    [ ~, idxStart] = min(abs(time - d/v(i) ));
    [ ~, idxEnd] = min(abs(time - (d+L)/v(i) ));
    y(idxStart:idxEnd) = h*sin( (pi*v(i)/L) * (time(idxStart:idxEnd) - d/v(i) ));
    
    % using the ode solver
    [t , z]= ode45 ( @ (t , z ) multipleDOF_Car_Suspension_Function(t,z,y,yTime,m,c,k) , time , IC );
    %only passing time into function, and not just ode45, so that it can
    %interpolate the input y for any time t

    figure; 
    plot(t , z(: ,1)) %plots all of z_1, aka, all of x_1
    hold on
    plot(t , z(: ,2)) %plots all of z_2, aka, all of x_2
    hold on
    plot(t,y,'k--')
    legendStuff = ["$Mass \ 1, \ x_{1}$","$Mass \ 2, \ x_{2}$","$Road \ Surface, \ y(t)$"];
    lgnd= legend(legendStuff);
    set(lgnd, 'Interpreter','latex')
    lgnd.Location = 'best';
    title( strcat("Positions of Masses on Multiple DOF Suspension at ", num2str(v(i)/(5280*12*0.0254/3600)) ,"mph"),...
        'Interpreter','latex')
    xlabel("$ t, \ time \ (s) $",'Interpreter','latex') 
    ylabel("$X, \ position \ (m)$",'Interpreter','latex')
    hold off;

    plot1 = strcat("Positions of Masses on Multiple DOF Suspension at ", num2str(v(i)/(5280*12*0.0254/3600)) ,"mph");
    print('-r600','-dpng',plot1);

    %Finding maximum acceleration of m1, the chassis
    accel = [ (1/m(1)) * ( c(1)*z(:,4) + k(1)*z(:,2) - c(1)*z(:,3) - k(1)*z(:,1) )  , ...
        (1/m(2)) * ( k(2)*y + c(1)*z(:,3) + k(1)*z(:,1) - c(1)*z(:,4) - (k(1)+k(2))*z(:,2) ) ];
    [maxAccel1,accIdx1] = max(abs(accel(:,1)));%note, will need abs to get max negative, but will only return pos values
    [maxAccel2,accIdx2] = max(abs(accel(:,2)));
    maxAccel1 = accel(accIdx1,1);
    maxAccel2 = accel(accIdx2,2);

    figure; 
    pointPlot1 = plot(t , accel(:,1) ); %plots all of dot_dot_z_1, check formula on state space work
    hold on
    %plot(time(accIdx1),maxAccel1,'HandleVisibility','off','MarkerSize',25);
    %(1/m(2))*( k(2)*y'+c(1)*z(:,3)+k(1)*z(:,1)-c(1)*z(:,4)-(k(1)+k(2))*z(:,2) )
    %(1/m(2))*( -k(2)*(z(:,2)-y') + k(1)*(z(:,1)-z(:,2)) + c(1)*(z(:,3)-z(:,4)) )
    pointPlot2 = plot(t , accel(:,2)  ); %plots all of dot_dot_z_2, check formula on state space work
    %plot(time(accIdx2),maxAccel2,'color',get(pointPlot,'color'),'HandleVisibility','off','MarkerSize',10);

    hold on
    plot(time(accIdx1),maxAccel1,".",'color',get(pointPlot1,'color'),'HandleVisibility','off','MarkerSize',8);
    hold on
    plot(time(accIdx1),maxAccel1,"hexagram",'color',get(pointPlot1,'color'),'MarkerSize',6,'linewidth',1);
    hold on
    plot(time(accIdx2),maxAccel2,".",'color',get(pointPlot2,'color'),'HandleVisibility','off','MarkerSize',8);
    hold on
    plot(time(accIdx2),maxAccel2,"hexagram",'color',get(pointPlot2,'color'),'MarkerSize',6,'linewidth',1);

    legendStuff = ["$Mass \ 1, \ \ddot{x}_{1}$","$Mass \ 2, \ \ddot{x}_{2}$"...
        ,strcat("$Max \ Acceleration \ of \ (m_{1}) \ = \ "...
        ,num2str(maxAccel1) ," \ \frac{m}{s^2} \ $"),...
        strcat("$Max \ Acceleration \ of \ (m_{2}) \ = \ ",num2str(maxAccel2) ," \ \frac{m}{s^2} \ $")];
    lgnd= legend(legendStuff);
    set(lgnd, 'Interpreter','latex')
    lgnd.Location = 'best';
    title( strcat("Accelerations of Masses on Multiple DOF Suspension at ", ...
        num2str(v(i)/(5280*12*0.0254/3600)) ,"mph"),'Interpreter','latex')
    xlabel("$ t, \ time \ (s) $",'Interpreter','latex') 
    ylabel("$ \ddot{X}, \ acceleration \ ( \frac{m}{s^2} )$",'Interpreter','latex')

    
    hold off;

    plot1 = strcat("Accelerations of Masses on Multiple DOF Suspension at ", num2str(v(i)/(5280*12*0.0254/3600)) ,"mph");
    print('-r600','-dpng',plot1);
end



