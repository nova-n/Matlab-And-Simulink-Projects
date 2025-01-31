clear;
clc;

%Excersise 3.02

%%Constants
N = 21; %number of points
L=1;
M = 1001; %number of time steps

alpha = 1;
dt = 0.002;
dx = L /(N-1);
beta = 0.5*alpha*dt/(dx^2);

%%Vector u
u = zeros(N,M);
x = linspace(0,L,N);

beta = 0.5*alpha*dt/(dx^2);
% Dirichlet boundary conditions
u0 = 0; % BC at x = 0
uL = 0; % BC at x = L
% Define arrays
x = linspace(0,L,N); % domain in x



%need array for system of equations
A = zeros(N,N); %square matrix A
b = zeros(N,1); % square matrix A
a = [-beta,1+2*beta,-beta];
b = zeros(N,1);% Define the initial condition

uo = @(x)( exp( -200*(x-0.5).^2 ) );
u(:,1) = uo( x' );
u
% Build matrix A
for i=2:(N-1)
    j = (i-1):(i+1);
    A(i,j) = a;
end
% Impose boundary conditions
%have not yet implemented the boundary conditions, where we force the derivatives to equal 0
%so must use forward or backwards differencing, cannot use central
%at left boundary, use forward differencing, at right boundary, use backwards differencing
A(1,1:3) = [-3,4,-1]; %forward diff at i = 1 (left boundary)
A(N,(N-2):N) = [1,-4,3];%backward diff at i = N (right boundary)
% Compute the inverse of A
A_inv = inv(A);
%note, matrix A was almost a diagonal matrix, just add row vector [-Fo,1+2*Fo , -Fo] to each "diagonal of the matrix"

% for i = 2:N-1
%     addRow = zeros(1,N);
%     addRow(i-1:i+ (-1+length(a)) -1) = addRow(i-1:i+ (-1+length(a)) -1) + a; %starts at i-1, since i starts at 2
%     A(i,:) = A(i,:) + addRow;
% end
%implement crank-nicholson
for timeStep = 1:M-1
    %update known vector b
    for i = 2:N-1
        b(i) = beta*u(i-1,timeStep) + (1-2*beta)*u(i,timeStep) + beta*u(i+1,timeStep);
    end
    b(1) = 0; %derichelet boundary conditions, like in the problem
    b(N) = 0;
    b;
    u(:,timeStep+1) = A_inv * b;
end
figure;
plot(x,u(:,:))