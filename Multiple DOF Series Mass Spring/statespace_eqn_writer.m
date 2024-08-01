clc;
clear;

syms x
disp("\bullet  \ z_{1} = x_{1} \\")
disp("\bullet  \ z_{2} = x_{2} \\")
disp("\bullet  \  z_{3} = \dot{z}_{1} \ = \dot{x}_{1} \\")
disp("\bullet  \ z_{4} = \dot{z}_{2} \ = \dot{x}_{2} \\")
disp("\bullet  \  z_{5} = \dot{z}_{3} \ = \ddot{z}_{1} \ = \ddot{x}_{1} \\")
disp("\bullet  \ z_{6} = \dot{z}_{4} \ = \ddot{z}_{2} \ = \ddot{x}_{2} \\")
disp("\\")
disp("\vec{z} = ")

sI_A = sym(  [[x,-1,0];[0,x-1,-2];[5,4,x+3]]  );
Ident = sym(  [[1,0,0];[0,1,0];[0,0,1]]  );