clear;
clc;

%Excersise 2.01

%%Constants
N = 201; %number of points
L=1;
M = 100; %number of time steps

c = 1;
dt = 0.001;
dx = L /(N-1);
Co = c*dt/dx;

%%Vector u
u = zeros(N,M);
x = linspace(0,L,N);

%%initial and Boundary Conditions
u(:,1) = [exp(-200*(x-0.25).^2)]';

%%Marcing time forward
for t_step=1:M-1
    t = t_step*dt;
    for i=2:N-1
        u(i,t_step+1) = u(i,t_step) - Co*( u(i,t_step) - u(i-1,t_step) );
    end
end

%u(:,i+1) = u(i,n) + Fo * ( u(i));