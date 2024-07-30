clc;
clear;
close all;
%example 4.02
M = 41;
L = 10;
dx = 2*L/(M-1);

%swap his N for M, and his M for N

%note, there DOES exist an exact solution to the logistic eqn
f = @(x)( 1./(1+exp(-x)) ); %is the function of logistic eqn
F = @(u)( u*(1-u) ); %is the derivative of logistic eqn

%define initial condition
u0 = f(-L); %computing actual value at starting point

%define arrays
x = linspace(-L,L,M); %for discretized
z = linspace(-L,L,101); %for exact
u = zeros(M,3); %storing all three methods in the U vector, each column being for forward euler, midpoint rule, and 4th order runge-kutta

%defining initial condition
u(1,:) = u0; %ALL methods start with the same value

%implement Runge-Kutta methods
for n = 1:M-1
    %forward euler in first column
    u(n+1,1) = u(n,1) + F( u(n,1) )*dx;

    %midpoint rule in second column
    us = u(n,2) + 0.5*F( u(n,2) )*dx; %intermediate step
    u(n+1,2) = u(n,2) + F( us )*dx;

    %4th order runge kutta in third column
    k1 = F( u(n,3) );
    k2 = F( u(n,3) + 0.5*k1*dx );
    k3 = F( u(n,3) + 0.5*k2*dx );
    k4 = F( u(n,3) + k3*dx );

    u(n+1,3) = u(n,3) + dx*(k1+2*k2+2*k3+k4)/6;
end

figure;
plot(x,u(:,1),'-b')
hold on;
plot(x,u(:,2),'-r')
hold on;
plot(x,u(:,3),'-g')
hold on;
plot(z,f(z),':k')
hold on;
