clear;
clc;

%Excersise 2.01

%%Constants
N = 21; %number of points
L=1;
M = 100; %number of time steps

alpha = 1;
dt = 0.001;
dx = L /(N-1);
Fo = alpha*dt/(dx)^2;

%%Vector u
u = zeros(N,M);
x = linspace(0,L,N);

%%initial and Boundary Conditions
u(:,1) = [sin(pi*x)]';

%%Marcing time forward
TO = 1;
TL = 2;
for t_step=1:M-1 %time loop
    t = t_step*dt;
    for i=2:N-1 %space loop
        u(i,t_step+1) = u(i,t_step) + Fo*( u(i+1,t_step) - 2*u(i,t_step) + u(i-1,t_step) );
    end
    %update boundary conditions
    %u(1,t_step+1) = TO;
    %u(N,t_step+1) = TL;
end

%u(:,i+1) = u(i,n) + Fo * ( u(i));