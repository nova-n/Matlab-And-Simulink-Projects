
% -------------------------------------------------------------------------
% EGME 436, Exercise 3.02 - Crank-Nicholson of the 1D heat equation
%
% Output files :
% Plots
%
% Code written by : S.Mayoral
% Last update : 2-Feb-2024 by S.Mayoral (
% ----------------------------------------------------------------------------
%
% INITIALIZE CODE
% ----------------------------------------------------------------------------
% Clear memory
clc
format short
% Constants
caseName = 'ex302'; % name of files
N = 21; % No. points in x
L = 1; % size of the domain
M = 1001; % No. time steps Dt
alpha = 1; % thermal diffusivity
dt = 0.002; % time-step
dx = L /(N-1); % cell size
% Compute the Fourier No.
beta = 0.5*alpha*dt/(dx^2);
% Dirichlet boundary conditions
u0 = 0; % BC at x = 0
uL = 0; % BC at x = L
% Define arrays
x = linspace(0,L,N); % domain in x
u = zeros(N,M); % array to store the solution u
A = zeros(N,N); % square matrix A
b = zeros(N,1); % square matrix A
a = [-beta, 1+2*beta, -beta]; % Coefficients of matrix A
% Define the initial condition
uo = @(x)( exp( -200*(x-0.5).^2 ) );
u(:,1) = uo( x' );
u
% Build matrix A
for i=2:(N-1)
    j = (i-1):(i+1);
    A(i,j) = a;
end
% Impose boundary conditions
A(1,1:3) = [-3, 4,-1];
A(N,(N-2):N) = [ 1,-4, 3];
% Compute the inverse of A
Ai = inv(A);
% Implement Backward Euler with CS for heat equation
for n = 1:(M-1)
    % Build the known vector b
    for i = 2:(N-1)
        b(i) = beta*u(i-1,n) + (1-2*beta)*u(i,n) + beta*u(i+1,n);
    end
    b(1) = 0;
    b(N) = 0;
    b;
    u(:,n+1) = Ai*b;
end
% Plot solution
figHandle = figure('Position', [100, 150, 350, 290]);
plot(x,u(:,1),'k-')
hold on
plot(x,u(:,251),'b-.')
plot(x,u(:,501),'r--')
plot(x,u(:,751),'g-')
plot(x,u(:,1001),'k:')
hold off
xlabel('x','FontSize',9)
ylabel('u(x,t)','FontSize',9)
axis([0 1 -0.2 1.2])
set(gca,'FontSize',9);
legend1 = legend('t=0.0','t=0.5','t=1.0','t=1.5','t=2.0');
set(legend1,'EdgeColor',[1 1 1],'FontSize',9,'Location','northwest');
% Export as a PNG file or PDF
%pngFile = strcat(caseName,'.png');
%pdfFile = strcat(caseName,'.pdf');
%exportgraphics(figHandle,pngFile,'Resolution',300)
%exportgraphics(figHandle,pdfFile,'ContentType','vector')
% Generate GIF file
gifHandle = figure('Position', [500, 150, 350, 290]);
plot(x,u(:,1),'b-')
xlabel('x','FontSize',9)
ylabel('u(x,t)','FontSize',9)
axis([0 1 -0.2 1.2])
set(gca,'FontSize',9)
%gifFile = strcat(caseName,'.gif');
%exportgraphics(gifHandle,gifFile);
% for n = 2:10:M % Note: the 10 speeds-pd the GIF
%     plot(x,u(:,n),'b-')
%     xlabel('x','FontSize',9)
%     ylabel('u(x,t)','FontSize',9)
%     axis([0 1 -0.2 1.2])
%     set(gca,'FontSize',9)
%     exportgraphics(gifHandle,gifFile,Append=true);
% end
