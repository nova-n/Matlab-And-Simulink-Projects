%ex5.02 2D Heat Eqn Point Jacopi
clear;
clc;
close all;

maxItterations = 21;
tol = 1e-4;
omega = 1.8

M = 21;
L = 1;
N = 1500;
alpha = 1;
dt = 0.0005;
dx = 2*L/(M-1);
Fo = alpha*dt/(dx^2); %same Fo in bith x and y, since same spacing in both x and y