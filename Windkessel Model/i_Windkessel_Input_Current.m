clc;
clear;
%close all;

I_0 = 500;
t_s = 0.3;
t_end_cyc = 0.8;
resolution = 1000;
sampleTime = [t_end_cyc/resolution , 0]; %period,offset

t = linspace(0,t_end_cyc,resolution);
%t_flat = linspace(t_s,t_end_cyc,resolution);
t_humplast_index = max(find(t <= t_s));

i = zeros(1,resolution);
i(1:t_humplast_index) = I_0 * sin( pi* t(1:t_humplast_index) / t_s ) .^2;
i(t_humplast_index + 1 : resolution) = 0;

i_of_t = [t',i'];

R = [1,NaN,NaN;1,0.79,0.63;1,0.79,0.63];
C = [1,NaN,NaN;1,1.75,5.16;1,1.22,2.53];
r = [NaN,NaN,NaN;0.05,0.033,0.03;0.05,0.056,0.045];
L = [NaN,NaN,NaN;NaN,NaN,NaN;0.005,0.0051,0.0054];

%plot(t,i)