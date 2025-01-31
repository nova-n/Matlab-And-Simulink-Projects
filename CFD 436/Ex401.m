clear;
clc;
%Ex 4.01
c = 1;
alpha = 0.1;
M = 21;
L = 1;
dx = L/(M-1);

%swap his N for M, and his M for N

%Boundary Conditions (Dirichilet, since fixed values)
u0 = 0;
uL = 1;

%beta
beta = c*dx/alpha;

%define arrays;
x = linspace(0,L,M);
A = zeros(M,M); %for upwind scheme, UDS
D = zeros(M,M); %central differencing
b = zeros(M,1);

%for each "triplet" in each row of the matrix, except for the corners
a = [-(beta + 1), beta+2, -1]; %is upwind
d = [-(beta + 1), 2 , beta-1]; %is central differencing

%build matricies
for i =2:(M-1)
    j = (i-1):(i+1); %so can shhift the "center" of the triplet in the matrix
    A(i,j) = a; %upwind
    D(i,j) = d;%central differencing
end

%Adding Boundary Conditions
A(1,1) = 1;
A(M,M) = 1;
D(1,1) = 1;
D(M,M) = 1;

%update b vector;
b(1) = u0;
b(M) = uL;

%solve the system of equations in the matrix, where each u value in the u matrix is the value at each point of x
u = A\b; %upwind
v = D\b; %central differencing, using v for vector instead of u, but is same thing

%exact solution, to be compared to u and v
Pe = c*L/alpha;
z = linspace(0,L,100);
w = ( exp(z*Pe/L) - 1 )/( exp(Pe) - 1 );

figure;
plot(x,u,'--b')
hold on;
plot(x,v,'--r')
hold on;
plot(z,w,'k:')
hold off;
