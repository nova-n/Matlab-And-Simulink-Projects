clear;
clc
%example: invicid burgers eqn

M = 201; %space steps
N = 101; %time steps
L = 1;
dt = 0.0005;
dx = L/(N-1);

beta = 0.5*dt/dx;

u0 = 0.5;

gamma = @(q)( 1 + beta*(2*q(2) - q(1) ) ); %here, q takes u, and u-1

%Spatial discretization
F = @(q)( -0.5*q(2)*( q(3) - q(1) )/dx );%central differencing

%Defining arrays
x = linspace(0,L,M);
u = zeros(M,N);
A = zeros(M,M);
a = [-beta,beta+2,-1]; %is the "triplet" in each diagonal of the matrix

%Initial conditions
uo = @(x)( 0.5 + 0.5*exp(-200*(x-0.5).^2) );
u(:,1) = uo(x');

A(1,1) = 1;
b(1) = u0;
%is implicit solution
%begin time stepping
for n= 1:(N-1)
    %crank nicholson
    %rebuild the matrix at each step
    for i =2:M %note, dont go to M-1, since boundary at end is not defined, so not really a boundary
        j = (i-1):i;
        A(i,j) = [-beta*u(i,n), gamma( u(j,n) )];
        b(i) = u(i,n);
    end
    u(:,n+1) = A\b';
end

figure;
plot(x,u)
hold off;