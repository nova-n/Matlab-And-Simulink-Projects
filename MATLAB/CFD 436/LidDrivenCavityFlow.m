%ex5.03 Lid Driven Cavity Flow
clear;
clc;
close all;

Mx = 41;
My = 41;
N = 700;
maxItterations = 50;
L = 1;
U_lid = 1;
omega = 1.8; %Using SOR
dx = L/(Mx-1);
dy = L/(My-1);
x = linspace(0,L,maxItterations);
y = x; %same mesh size in y
[X,Y] = meshgrid(x,y);

%fluid properties
rho = 1;
nu = 0.1;
dt = 0.001;

%building arrays
u = zeros(Mx,My,N); %velocity in x
v = u; %velocity in y
p = u; %pressure
U = u; %velocity magniture (for post processing)
S = zeros(Mx,My); %source term for pressure Poisson eqn
    %note, is not needed for every time step as of now, since will be computed and
    %overwritten at each timestep

%Boundary Condition
u(:,My,:) = U_lid;

%Implement CFD (navier stokes, poisson eqn)
for n = 1:N-1
    %Building source term at each timestep
    for j=2:My-1
        for i=2:Mx-1
            du_dx = 0.5*( u(i+1,j,n) - u(i-1,j,n) )/dx;%partial u wrt x
            du_dy = 0.5*( u(i,j+1,n) - u(i,j-1,n) )/dy;%partial u wrt y
            dv_dx = 0.5*( v(i+1,j,n) - v(i-1,j,n) )/dx;%partial v wrt x
            dv_dy = 0.5*( v(i,j+1,n) - v(i,j-1,n) )/dy;%partial v wrt y

            S(i,j) = rho*( du_dx + dv_dy )/dt - rho*( du_dx^2 + 2*du_dy*dv_dx + dv_dy^2 );
        end
    end

    %%Solving pressure poisson eqn
    %%Implementing Gauss-Sidel SORol) && (itt < maxItterations)
    err = 1;
    while (err > tol) && (itt < maxItterations)
        itt = itt+1;
        for j = 2:My-1
            for i = 2:Mx-1 %remember, is steady, so instead of timestep march, is ittereation step march
                p0 = 2*(dx*dx + dy*dy);
                p1 = dy^2 * ( p(i+1,j,n) + p(i-1,j,n) + dx^2.*(p(i,j+1,n)+p(i,j-1,n)) ); %can overwrite current p space at time n by using (i+1,j,n)
                p2 = rho*dx^2*dy^2 * S(i,j);
                %intermediate step
                ps = p1/p0 - p2/p0;

                pnew = p(i,j,n) + omega*(ps - p(i,j,n) );
            end
        end
        %update error
        err = max( max( abs( p(i,j,n) - pnew )))/ max(max(abs(f))); %is 2d matrix, so take max of the max
    end
end

