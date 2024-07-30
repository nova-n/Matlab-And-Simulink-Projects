clear ; clc ; close all ;
tstart =0; % time starts at 0 seconds
tend =10; % time ends at this amount of time
% specify the time vector
time = tstart :.1: tend ;
% initial conditions
x_0 =0; % in m
x_dot_0 =0; % in m / s
IC =[ x_0 x_dot_0 ]';
% using the ode solver
[t , z ]= ode45 ( 'ExSDOF' , time , IC );
figure; 
plot(t , z (: ,1) , 'k ')
xlabel ( ' time ( s ) ')
ylabel ( ' position ( m ) ')