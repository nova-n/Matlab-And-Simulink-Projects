clc;
clear;
close all;

fileName = 'C:/Users/n8dsa/Documents/OpenFOAM/cavity_project/cavity/postProcessing/graphUniform/459/line.xy'
A = importdata(fileName) %%new way to read files in matlab, imports it as an object

y = A.data (:,1);
w = A.data (:,8);

plot(y,w)