function z_dot=ExSDOF2(t,z,m,c,k,ft,fval)
% This is the function for the ode solvers for a MSD system with forcing
% f(t) where the variables m,c,k, and function f are passed into this file

% forcing function
f=interp1(ft,fval,t); % interpolate the data (ft, fval) at time t

% state space model
z_dot(1)=z(2);
z_dot(2)=1/m*(f-c*z(2)-k*z(1));

z_dot=z_dot'; % output needs to be a column vector
