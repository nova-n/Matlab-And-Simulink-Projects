clc;
clear;
close all;

%% Import Data
fileName = 'C:/Users/n8dsa/Documents/OpenFOAM/cavity_project/cavity/postProcessing/cutPlaneSurface/459/cutPlane.xy';
A = importdata(fileName) %%new way to read files in matlab, imports it as an object

x = A.data (:,1);
y = A.data (:,2);
u_x = A.data (:,4);
u_y = A.data (:,5);
vor= A.data (:,9);

%% Now, creating streamlines the hard way, integrating numerically

%% Map the 2D data onto a meshgrid x-y plane
Nx = 100;
Ny = 100;
xMin = min(x);
xMax = max(x);
yMin = min(y);
yMax = max(y);

X = linspace(xMin,xMax,Nx);
Y = linspace(yMin,yMax,Ny);

[xq,yq] = meshgrid(X,Y);

u_x_q = griddata(x,y,u_x,xq,yq);
u_y_q = griddata(x,y,u_y,xq,yq);
vor_q = griddata(x,y,vor,xq,yq);

figure;
contourf(X,Y,vor_q) %velocity in x dir
hold off;
figure;
contourf(X,Y,u_y_q) %velocity in y dir
hold off;

%%Compute the stream function
psi = zeros(Nx,Ny);
for i = 2:Nx-1
    for j = 2:Ny-1
        psi(i,j) = u_x_q(i,j)* (Y(j+1) -Y(j) ) + psi(i,j-1);
    end
end

figure;
contour(X,Y,psi) %plotting streamlines
hold off;