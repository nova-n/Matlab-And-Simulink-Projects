% -------------------------------------------------------------------------
%   EGME 436, Exercise 9.01 - Steady convection-diffusion equation with the
%   finite volume method
%
%   Output files :
%    Plot
%
%   Code written by :  S.Mayoral
%   Last update     :  18-Mar-2024 by S.Mayoral
%   ----------------------------------------------------------------------------
%
%   INITIALIZE CODE
%   ----------------------------------------------------------------------------
%   Clear memory
    clc
    format short

%   Constants
    caseName = 'ex901';         %  name of files
    N = 101;                    %  No. points in x
    L = 1;                      %  size of the domain
    Gamma = 0.2;                %  thermal diffusivity
    rho = 1;                    %  fluid density
    c =  1.0;                   %  advection speed
    dx = L /(N-1);              %  cell size
    
%%  -----------------
%   START CODING FROM HERE

    %   Compute D and F
    F = rho*c;
    D = Gamma/dx;
    
    a_E = D-0.5*F;
    a_W = D+0.5*F;
    a_P = a_E + a_W; %seems like no source term
    
    %   Dirichelet Boundary Conditions
    u_0 = 0;
    u_L = 1;
    
    %Defining Arrays
    x = linspace(0,L,N);
    A = zeros(N,N);
    b = zeros(N,1);
    
    a = [-a_W,a_P,-a_E]; %that is the mini row you see patterned accross the diagonal
    
    %Building Matrix A Without Boundary Conditions
    for i = 2:N-1
        j = (i-1):(i+1); %this starts at the indecie left of the current a_P, and ends at the current a_E
        A(i,j) = a;
    end
    
    %Impose Boundary Conditions on Matrix A
    A(1,1:2) = [a_W+2*D, -a_E];
    A(N,N-1:N) = [-a_W a_E + 2*D];
    
    %Building Known Vector, b
    b(1) = u_0*(2*D+F); %No source term, so almost all of b is 0, but boundary conditions may or may not be 0;
    b(N) = u_L*(2*D-F);
    
    %Now, Solve System Of Equations
    u = A\b;

%   
%%  -----------------
%   Compute the exat solution
    Pe = c*L/Gamma;
    z = linspace(0,L,100); 
    w = ( exp(z*Pe/L)-1 ) / ( exp(Pe)-1 );


%   Plot solution 
    figHandle = figure('Position', [100, 150, 350, 290]);
    plot(z,w,'k-.')
    hold on
    plot(x,u,'r-')
    hold off
    xlabel('x','FontSize',9)
    ylabel('u(x,t)','FontSize',9)
    axis([0 1 -0.2 1.2])
    set(gca,'FontSize',9);
    legend1 = legend('Exact','FVM');
    set(legend1,'EdgeColor',[1 1 1],'FontSize',9,'Location','northwest');

%   Export as a PNG file or PDF
    pngFile = strcat(caseName,'.png');
    pdfFile = strcat(caseName,'.pdf');
    exportgraphics(figHandle,pngFile,'Resolution',300)
    exportgraphics(figHandle,pdfFile,'ContentType','vector')





