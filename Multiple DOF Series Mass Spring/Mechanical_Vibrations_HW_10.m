clc;
clear;
close all;

res = 1000;
t_max = 60;
t = linspace(0,t_max,res);

w_n_Array = [sqrt(2)/2,2]; %[w_n1, w_n2, ...]

u_vectors_Matrix = [[1,3];[-1,3]]; %[u_1 ; u_2, ...]
%u_1 = [1,3];
%u_2 = [-1,3];

%PHI IS IN DEGREES

%%Part A
A_Array = [sqrt(2)/6,0]; %[A_1,A_2, ...]
phi_Array = [90,0]; %[phi_1, phi_2, ...] in degrees
%NOTE: phi_2 doesn't matter for part A

X_vals = coupledDoFPlotter(A_Array, phi_Array,w_n_Array,u_vectors_Matrix,t,"$Problem \ 2A$");

%%Part C
A_Array = [0/2,-sqrt(2)/12]; %[A_1,A_2, ...]
phi_Array = [0,0]; %[phi_1, phi_2, ...] in degrees
%NOTE: phi_1 doesn't matter for part C

X_vals = coupledDoFPlotter(A_Array, phi_Array,w_n_Array,u_vectors_Matrix,t,"$Problem \ 2C$");

%%Coupled 3 DoF Test
A_Array = zeros(3);
w_n_Array = zeros(3);
u_vectors_Matrix = [];
A_Array = [1,-sqrt(2)/12,1]; %[A_1,A_2, ...]
phi_Array = [40,70,20]; %[phi_1, phi_2, ...] in degrees
%NOTE: phi_1 doesn't matter for part C
w_n_Array = [sqrt(5)/5,2,4]; %[w_n1, w_n2, ...]

u_vectors_Matrix = [[1,3,4];[-1,3,5];[5,4,6]]; %[u_1 ; u_2, ...]


X_vals = coupledDoFPlotter(A_Array, phi_Array,w_n_Array,u_vectors_Matrix,t,"$Coupled \ 3 \ DoF \ TEST$");

%%Coupled 4 DoF Test
A_Array = zeros(4);
w_n_Array = zeros(4);
u_vectors_Matrix = [];
A_Array = [6,-sqrt(2)/12,1,4]; %[A_1,A_2, ...]
phi_Array = [40,70,20,10]; %[phi_1, phi_2, ...] in degrees
%NOTE: phi_1 doesn't matter for part C
w_n_Array = [sqrt(5)/5,2,4,5]; %[w_n1, w_n2, ...]

u_vectors_Matrix = [[1,3,4,1];[-1,3,5,8];[5,4,6,7];[1,2,3,4]]; %[u_1 ; u_2, ...]


X_vals = coupledDoFPlotter(A_Array, phi_Array,w_n_Array,u_vectors_Matrix,t,"$Coupled \ 4 \ DoF \ TEST$");

%%Coupled 5 DoF Test
A_Array = zeros(5);
w_n_Array = zeros(5);
u_vectors_Matrix = [];
A_Array = [1,6,-sqrt(2)/12,1,4]; %[A_1,A_2, ...]
phi_Array = [0,40,70,20,10]; %[phi_1, phi_2, ...] in degrees
%NOTE: phi_1 doesn't matter for part C
w_n_Array = [sqrt(5)/5,2,4,5,8]; %[w_n1, w_n2, ...]

u_vectors_Matrix = [[1,3,4,1,0];[-1,3,5,8,1/2];[5,4,6,7,sqrt(2)];[1,2,3,4,5];[5,4,3,2,1]]; %[u_1 ; u_2, ...]


X_vals = coupledDoFPlotter(A_Array, phi_Array,w_n_Array,u_vectors_Matrix,t,"$Coupled \ 5 \ DoF \ TEST$");


%%Coupled n-DoF Test
DoF = 21;
% randMin = 0;
% randMax = 10;
A_Array = zeros(DoF);
w_n_Array = zeros(DoF);
u_vectors_Matrix = [];

for i = 1:DoF
    randMin = -10;
    randMax = 10;
    randomNum = ((randMax-randMin)*rand + randMin)/(randMax-randMin);
    A_Array(i) = randomNum;
    randMin = 0;
    randMax = 7;
    randomNum = (randMax-randMin)*rand + randMin;
    w_n_Array(i) = randomNum;
    for ii = 1:DoF
        randomNum = ((randMax-randMin)*rand + randMin)/(randMax-randMin);
        u_vectors_Matrix(i,ii) = randomNum;
    end
end

w_n_Array = sort(w_n_Array);

randMin = -360;
randMax = 360;

for i = 1:DoF
    randomNum = (randMax-randMin)*rand + randMin;
    phi_Array(i) = randomNum;
end

X_vals = coupledDoFPlotter(A_Array, phi_Array,w_n_Array,u_vectors_Matrix,t,strcat("$Coupled \ ", num2str(DoF), " \ DoF \ TEST$"));



function [coupledDoFValues] = coupledDoFPlotter(A_Values, phi_Values,natFreqs,u_Vectors,t_vals,titled)
    X = zeros([length(A_Values),length(t_vals)]);
    for i = 1:length(A_Values) %gets X mass1, X mass2, ... and selects which u vector to use
        %X(i,:) = X(i,:) + A_Values(1)*sin(natFreqs(1)*t_vals + phi_Values(1)*180/pi)*u_Vectors(1,i);
        for ii = 1:length(A_Values) %within X mass1, goes thru A_1, A2..., phi_1, phi_2..., selects
            %which value in the specified u vector to use
            X(i,:) = X(i,:) + A_Values(ii)*sin(natFreqs(ii)*t_vals + phi_Values(ii)*pi/180)*u_Vectors(ii,i);
        end
    end
    coupledDoFValues = X;
    figure;
    for i = 1:length(A_Values)
        plot(t_vals,X(i,:))
        hold on;
    end
    hold off;
    a = strcat("$Mass \ ", num2str(1)," \ X \ Position$");
    legendStuff = strings(1,length(A_Values)); %creates aray of empty string
    %legendStuff(1) = strcat("$Mass \ ", num2str(1)," \ X \ Position$");
    %legendStuff(2) = strcat("$Mass \ ", num2str(2)," \ X \ Position$");    
    %legendStuff = [strcat("$Mass \ ",num2str(1)," \ X \ Position$"),strcat("$Mass \ ", num2str(2)," \ X \ Position$")];
    for i = 1:length(A_Values)
        a = strcat("$Mass \ ", num2str(i)," \ X \ Position : \ ");
        b = "";
        for ii = 1:length(A_Values)
            if ii > 1
              b = strcat(b, " \ + \ ") ;
            end
            b = strcat(b, num2str(u_Vectors(ii,i)*A_Values(ii)),"sin(",num2str(natFreqs(ii)),"t \ + \ ", num2str(phi_Values(ii)), "^{\circ}  )");
        end
        c = strcat(a,b," $");
        legendStuff(i) = c;
    end
    lgnd= legend(legendStuff,'Location','northoutside');
    set(lgnd, 'Interpreter','latex')
    %lgnd.Location = 'northoutside';	
    title(titled,'Interpreter','latex')
    xlabel("$ t, \ time \ (s) $",'Interpreter','latex') 
    ylabel("$X, \ position \ (m)$",'Interpreter','latex') 
end

%%This is the model
% A_1 = sqrt(2)/2;
% A_2 = 0;
% phi_1 = 90;
% phi_2 = 0; %doesn't matter for part A
% 
% u_1 = [1;3];
% u_2 = [-1;3];
% 
% X= [;];
% X(1,:) = A_1*sin(w_n1*t + phi_1*180/pi)*u_1(1) + A_2*sin(w_n2*t + phi_2*180/pi)*u_2(1);
% X(2,:) = A_1*sin(w_n1*t + phi_1*180/pi)*u_1(2) + A_2*sin(w_n2*t + phi_2*180/pi)*u_2(2);


