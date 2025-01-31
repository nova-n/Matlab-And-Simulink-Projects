clc;
clear;
close all;

%%Problem 5a 
dt = 0.01;
t = 0:dt:0.06;
N = length(t); %time steps
Fo = 0.25; %can choose any dx and alpha, as long as Fo = 0.25 for this problem
L = 1;
M = 51;%grid points on x axis 
dx = L/(N-1);
alpha = (Fo*dx^2)/dt;
x = [linspace(0,L,M)]';
u = zeros(M,N);


%Problem 5A
u = initConditions(N,M,x,u);
u = FTCS_time_march(N,M,u,Fo);
graphPlotter(x,u,t,"$Problem 5A$")
%Problem 5B
u = initConditions(N,M,x,u);
Fo = 1;
alpha = (Fo*dx^2)/dt;
u = FTCS_time_march(N,M,u,Fo);
graphPlotter(x,u,t,"$Problem 5B$")

%Problem 5C
Fo = 0.25;
alpha = (Fo*dx^2)/dt;
u = initConditions(N,M,x,u);

%need array for system of equations
A = zeros(M,M); %square matrix A
a = [-Fo,1+2*Fo,-Fo];
b = zeros(M,1);% Define the initial condition
%uo = @(x)( exp( -200*(x-0.5).^2 ) );
%u(:,1) = uo( x' );
% Build matrix A
for i=2:(M-1)
    j = (i-1):(i+1);
    A(i,j) = a;
end
% Impose boundary conditions
%have not yet implemented the boundary conditions, where we force the derivatives to equal 0
%so must use forward or backwards differencing, cannot use central
%at left boundary, use forward differencing, at right boundary, use backwards differencing
A(1,1:3) = [-3,4,-1]; %forward diff at i = 1 (left boundary)
A(M,(M-2):M) = [1,-4,3];%backward diff at i = N (right boundary)
% Compute the inverse of A
A_inv = inv(A);
%note, matrix A was almost a diagonal matrix, just add row vector [-Fo,1+2*Fo , -Fo] to each "diagonal of the matrix"

%implement crank-nicholson
for timeStep = 1:N-1
    %update known vector b
    for i = 2:M-1
        b(i) = Fo*u(i-1,timeStep) + (1-2*Fo)*u(i,timeStep) + Fo*u(i+1,timeStep);
    end
    b(1) = 0; %derichelet boundary conditions, like in the problem
    b(N) = 0;
    b;
    u(:,timeStep+1) = A_inv * b;
end
graphPlotter(x,u,t,"$Problem 5C$")
plot1 = "CFDHW1_PROBLEM3";
print('-r600','-dpng',plot1);
hold off;


function [u_return] = initConditions(timeSteps,spaceSteps,x_vals,u_vals)
    u_vals = zeros(spaceSteps,timeSteps); %each column is a spatial point, and going right gives you the time step;
    %defining initial conditions, where u is 0 from x<=0.5, and is 1-4*abs(x-0.75) for x>0.5
    firstSlopeIndex = find( abs(x_vals - 0.5) <= 0.021); %actually don't want abs, since want <= 0.5
    firstSlopeIndex = firstSlopeIndex(end);
    u_vals(firstSlopeIndex:end,1) = 1- 4.* abs(x_vals(firstSlopeIndex:end) - 0.75);
    u_return = u_vals;
end

function [u_return] = FTCS_time_march(timeSteps,spaceSteps,u_vals,Fo_val)
%since have 1st spatial distribution at time = 0, then can start at second time index
    for n = 1:timeSteps-1 %because doing n+1 each time, then cannot go to end of time boundary
        u_vals(1,:) = 0;%boundary at start is always set to 0 for this problem
        u_vals(spaceSteps,:) = 0;%boundary at end is always set to 0 for this problem
        for i = 2:spaceSteps-1 %same idea here, but in space
            u_vals(i,n+1) = Fo_val* ( u_vals(i+1,n) - 2*u_vals(i,n) + u_vals(i-1,n) ) + u_vals(i,n);
        end
    end
    u_return = u_vals;
end

function graphPlotter(x_vals,u_vals,t_vals,titleStuff)
    figure;
    plot(x_vals,u_vals(:,:))
    title(titleStuff,'Interpreter','latex')
    xlabel("$ Space, \ x $",'Interpreter','latex')
    ylabel("$u(x,t_?)$",'Interpreter','latex')
    legendStuff = "$t = " + string(t_vals) + " s $";
    lgnd= legend(legendStuff,'Location','best');
    set(lgnd, 'Interpreter','latex')
    plot1 = strcat("CFD_HW1_",titleStuff);
    print('-r600','-dpng',plot1);
    hold off;
end