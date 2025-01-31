clc;
clear;
close all;

res = 1000;
t_max = 60;
t = linspace(0,t_max,res);

m = 10; k = 100;

DOF = 4;

m_matrix = m*eye(DOF);
k_matrix = [[5*k,0,-k,0];[0,4*k,0,0];[-k,0,5*k,0];[0,0,0,6*k]];

[u_vectors,w_n_squared] = eig(k_matrix,m_matrix);

w_n = sqrt(w_n_squared);

%normalizing u vectors
for i = 1:DOF
    u_vec = u_vectors(:,i)
    lowestNum = abs(min(u_vec(u_vec ~=0)))
    u_vectors(:,i) =  u_vectors(:,i)/lowestNum;
end

for i = 1:DOF
    disp("u_"+num2str(i) + "=")
    disp( "u_" + num2str(i) + "," + num2str([1:DOF]')+" = [" + num2str(u_vectors(:,i)) + "]" )
end

for i = 1:DOF
    disp("w_n"+num2str(i) + "=" + num2str(w_n(i,i)))
end
