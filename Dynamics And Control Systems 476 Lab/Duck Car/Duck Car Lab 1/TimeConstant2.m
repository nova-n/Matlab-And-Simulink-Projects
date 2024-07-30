% Ian Mark Bautista
% Lab 1 - 2nd method for Time Constant

clear;
clc;

ThreeV_PosData = readtable("3V_Ducky_Position.txt");
FourV_PosData = readtable("4v_Ducky_position.txt");

ThreeVData = table2array([readtable("3V_Ducky_Velocity.txt"),ThreeV_PosData(:,2)]); %both pos and vel have same timestamps
FourVData = table2array([readtable("4v_Ducky_velocity.txt"), FourV_PosData(:,2) ]);
FiveVData = table2array(readtable("Duck Car5V.txt"));

vel_stdy_3v = .9893;
vel_stdy_4v = 1.5156;
vel_stdy_5v = 2.2347;

theta_time_3v = log((vel_stdy_3v - ThreeVData(:,2))/vel_stdy_3v);

theta_poly_3v = polyfit(ThreeVData(:,1),theta_time_3v, 2);

plot(ThreeVData(:,1), theta_poly_3v)