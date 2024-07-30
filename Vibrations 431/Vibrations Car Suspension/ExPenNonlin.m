function z_dot=ExPenNonlin(t,z,g,L)
% This is the function for the ode solvers for a simple pendulum system
% that has NOT been linearized


z_dot(1)=z(2); 
z_dot(2)=-g/L*sin(z(1));

z_dot=z_dot';