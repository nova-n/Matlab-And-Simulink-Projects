%ex5.02 2D Heat Eqn Point Jacopi
clear;
clc;
close all;

maxItterations = 1500;
tol = 1e-4;
omega = 1.8;

M = 21;
L = 1;
alpha = 1;
dx = 2*L/(M-1);

%dirichilet boundary conditions (all set to 0, so no action needed)

%define arrays
x = linspace(-L,L,M);
y = x; %same mesh size in y
u_jacopi = zeros(M,M,maxItterations); %is a 3d matrix, x,y,t, is u(i,j) at timestamp n
u_gausssidel = u_jacopi;
u_SOR = u_jacopi;

%defining exact solution, to use as reference to terminate loop
[X,Y] = meshgrid(x,y); %plotting 2D array onto a contour plot 
f = (X.^2 -1).*(Y.^2 -1);
figure;
contour(X,Y,f);
%contourf(X,Y,f);
hold off;

%building source term
S = 2.* (X.^2 + Y.^2 - 2);

%%Implementing Point-Jacopi
err = 1;% 0.01%percent error, shown by tol. Starts as something that doesn't meet condition
itt = 0;

while (err > tol) && (itt < maxItterations)
    itt = itt+1;
    for j = 2:M-1
        for i = 2:M-1 %remember, is steady, so instead of timestep march, is ittereation step march
            u_jacopi(i,j,itt+1) = 0.25* ( u_jacopi(i+1,j,itt) + u_jacopi(i-1,j,itt) + u_jacopi(i,j+1,itt) + u_jacopi(i,j-1,itt) )...
            -0.25*dx^2 * S(i,j);
        end
    end
    %update error
    err = max( max( abs( f - u_jacopi(:,:,itt+1) ))) / max(max(abs(f))); %is 2d matrix, so take max of the max
end
err
itt_jacopi = itt
figure;
contour(X,Y,u_jacopi(:,:,itt_jacopi))
hold off;

%%FIX BELOW LATER!!!!!!!!

%%Implementing Gauss-Sidel
while (err > tol) && (itt < maxItterations)
    itt = itt+1;
    for j = 2:M-1
        for i = 2:M-1 %remember, is steady, so instead of timestep march, is ittereation step march
            u_gausssidel(i,j,itt+1) = 0.25* ( u_gausssidel(i+1,j,itt) + u_gausssidel(i-1,j,itt+1) + u_gausssidel(i,j+1,itt) + u_gausssidel(i,j-1,itt+1))...
            -0.25*dx^2 * S(i,j);
        end
    end
    %update error
    err = max( max( abs( f - u_gausssidel(:,:,itt+1) )))/ max(max(abs(f))); %is 2d matrix, so take max of the max
end
err
itt_gausssidel = itt
figure;
contour(X,Y,u_jacopi(:,:,itt_gausssidel))
hold off;


%%Implementing Gauss-Sidel SORol) && (itt < maxItterations)
while (err > tol) && (itt < maxItterations)
    itt = itt+1;
    for j = 2:M-1
        for i = 2:M-1 %remember, is steady, so instead of timestep march, is ittereation step march
            %intermediate step
            us = 0.25* ( u_SOR(i+1,j,itt) + u_SOR(i-1,j,itt+1) + u_SOR(i,j+1,itt) + u_SOR(i,j-1,itt+1))...
            -0.25*dx^2 * S(i,j);

            u_SOR(i,j,itt+1) = u_SOR(i,j,itt) + omega*(us-u_SOR(i,j,itt))
        end
    end
    %update error
    err = max( max( abs( f - u_SOR(:,:,itt+1) )))/ max(max(abs(f))); %is 2d matrix, so take max of the max
end
err
itt_SOR = itt
figure;
contour(X,Y,u_jacopi(:,:,itt_SOR))
hold off;

