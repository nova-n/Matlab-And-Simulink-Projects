clear;
clc;

%Excersise 3.01

%%Constants
N = 21; %number of points
L=1;
M = 1001; %number of time steps

alpha = 1;
dt = 0.002;
dx = L /(N-1);
Fo = alpha*dt/(dx)^2;

%%Boundary Conditions
uo = 0;
uL = 0;

%%Vector u
u = zeros(N,M);
x = linspace(0,L,N);
%need array for system of equations
A = zeros(N,N); %square matrix A
%note, matrix A was almost a diagonal matrix, just add row vector [-Fo,1+2*Fo , -Fo] to each "diagonal of the matrix"
a = - [-Fo,1+2*Fo,-Fo];
for i = 2:N-1
    addRow = zeros(1,N);
    addRow(i-1:i+ (-1+length(a)) -1) = addRow(i-1:i+ (-1+length(a)) -1) + a; %starts at i-1, since i starts at 2
    A(i,:) = A(i,:) + addRow;
end
A(1,1) = 1;
A(N,N) = 1;
A
det(A)
%Take inverse of A
A_inv = inv(A)

%define the source term
s = @(x,t)( (pi^2 -1)*exp(-t)*sin(pi*x) ); % @ is a handle function, define the function in the second set of parenthasees
%x and t are redefined, within the scope of the function
uo = @(x)(sin( pi*x ));

%write the initial condition
u(:,1) = uo( x' ); %inputting the transpose of the already defined (global) x
% figure;
% plot(x,u(:,1))
% hold off;

%implement backward euler
for timeStep = 1:M-1
    t = dt*timeStep;
    u(:,timeStep+1) = A_inv*(u(:,timeStep)) + dt*s(x',t);
    %Adding source term in discretization with the " +dt*s(x',t) "
end
figure;
plot(x,u(:,:))
hold off;
