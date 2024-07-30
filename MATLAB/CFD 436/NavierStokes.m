%ex5.03 2D NAvier Stokes
clear;
clc;
close all;

M = 41;
L = 1;
N = 1500;
alpha = 1;
dt = 0.0005;
dx = 2*L/(M-1);
Fo = alpha*dt/(dx^2); %same Fo in bith x and y, since same spacing in both x and y