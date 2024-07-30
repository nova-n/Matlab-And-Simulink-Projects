% This program is used to run ode45 with the simple pendulum (linearized 
% and nonlinear) system from lecture
% Written by: Dr H L Weiss
% Written on: 4/25/19

clear; clc; close all;

tstart=0; % time starts at 0 seconds
tend=10; % time ends at this amount of time
g=9.81; % m/s^2
L=0.2;  % m

% specify the time vector
time=linspace(tstart,tend,1e4);

% initial conditions
th_0=45*pi/180; % in rad
th_dot_0=0; % in rad/s

IC=[th_0 th_dot_0]';

% using the ode solver
[t,z]=ode45(@(t,z) ExPenLin(t,z,g,L),time,IC);
[tnon,znon]=ode45(@(t,z) ExPenNonlin(t,z,g,L),time,IC);

figure; hold on
plot(t,z(:,1), 'k')
plot(tnon,znon(:,1), 'r')
xlabel('time (s)')
ylabel('position (rad)')
legend('Linearized system', 'Nonlinear system')

figure; hold on
plot(t,z(:,1)*180/pi, 'k')
plot(tnon,znon(:,1)*180/pi, 'r')
xlabel('time (s)')
ylabel('position (deg)')
legend('Linearized system', 'Nonlinear system')
