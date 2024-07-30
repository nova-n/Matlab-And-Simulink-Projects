function z_dot=ExSDOF(t,z)
% This is the function for the ode solvers for a MSD system with forcing
% f(t)

% define constants
m=2; % in kg
c=1; % in N*s/m
k=5; % in N/m

% define forcing function
f=1; % unit step
%f=2*sin(t)+5;
%f=0; % no force

% state space model
z_dot(1)=z(2);
z_dot(2)=1/m*(f-c*z(2)-k*z(1));

z_dot=z_dot'; % output needs to be a column vector

