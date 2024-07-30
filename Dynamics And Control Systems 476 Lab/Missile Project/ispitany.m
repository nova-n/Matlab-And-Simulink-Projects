clc;
close all;

resol = 1000;
t = linspace(0,10,resol);
t = t';

%Original Comparison Inputs
%alphaDesired = (5*pi/180)*ones(1,resol);
%betaDesired = (0*pi/180)*ones(1,resol);

%Step Inputs
alphaDesired = (5*pi/180)*ones(1,resol);
betaDesired = (0*pi/180)*ones(1,resol);

%Stronger Step Inputs
%alphaDesired = (25*pi/180)*ones(1,resol);
%betaDesired = (25*pi/180)*ones(1,resol);

%Ramp Inputs
%alphaDesired = (5*pi/180)*t;
%betaDesired = (5*pi/180)*t;

%Stronger Ramp Inputs
%alphaDesired = (45*pi/180)*t;
%betaDesired = (45*pi/180)*t;

%Sinusoidal Inputs
%alphaDesired = (45*pi/180)*sin(2*pi*t);
%betaDesired = (45*pi/180)*cos(2*pi*t);

%Faster Sinusoidal Inputs
%alphaDesired = (45*pi/180)*sin(8*pi*t);
%betaDesired = (45*pi/180)*cos(8*pi*t);


alphaDesiredSimIn = timeseries(alphaDesired,t);
betaDesiredSimIn = timeseries(betaDesired,t);

%secs = seconds(t);
