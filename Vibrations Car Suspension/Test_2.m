clear ; clc ; close all ;
tstart =0; % time starts at 0 seconds
tend =10; % time ends at this amount of time
% specify the time vector
time = tstart :.1: tend ;
% initial conditions
x_0 =0; % in m
x_dot_0 =0; % in m / s
IC =[ x_0 x_dot_0 ]';
% Second version - passing the function and variables in as an ←↩ inputs
% define the constants
m =1; % in kg
c =3; % in N * s / m
k =5; % in N / m
% create a time vector for the forcing function
ft = linspace( tstart , tend ,1e4 ) ;
% define the focing function
fval = ft *0+1;
% using the ode solver
[ tex2 , zex2 ]= ode45 ( @ (t , z ) ExSDOF2 (t ,z ,m ,c ,k , ft , fval ) , time , IC ) ;
hold on ;
plot ( tex2 , zex2 (: ,1) , 'r ')