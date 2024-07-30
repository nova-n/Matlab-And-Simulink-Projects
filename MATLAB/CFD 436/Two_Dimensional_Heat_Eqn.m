%ex5.01 2D Heat Eqn
clear;
clc;
close all;

M = 81;
L = 1;
N = 1500;
alpha = 0.1;
dt = 0.0005;
dx = 2*L/(M-1);
Fo = alpha*dt/(dx^2); %same Fo in bith x and y, since same spacing in both x and y

%dirichilet boundary conditions (all set to 0, so no action needed)

%define arrays
x = linspace(-L,L,M);
y = x; %same mesh size in y
u = zeros(M,M,N); %is a 3d matrix, x,y,t, is u(i,j) at timestamp n

%defining source term, this time, using meshgrid, is needed anyway since plotting in 2D
[X,Y] = meshgrid(x,y); %plotting 2D array onto a contour plot 
S = 2 * (2 - X.^2 - Y.^2);
%contour(X,Y,S);
%contourf(X,Y,S);

%%Implement FTCS for 2d heat eqn
for n = 1:N-1
    for j = 2:M-1 %remember, position 1 is the boundary condition
        for i = 2:M-1
            u(i,j,n+1) = u(i,j,n) + Fo * ( u(i+1,j,n) - 2*(u(i,j,n)) + u(i-1,j,n) ...
            + u(i,j+1,n) - 2*(u(i,j,n)) + u(i,j-1,n) + dt*S(i,j) );
        end
    end
end

figure;
contourf(X,Y, u(:,:,M)) %at final timestep
hold off;