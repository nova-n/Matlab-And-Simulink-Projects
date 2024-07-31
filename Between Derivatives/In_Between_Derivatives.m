clc;
clear;
load census.mat
lineThick = 5;
maxGridX=5;
maxGridY=5;
derivativeCount=5;
integralCount = 5;
gridMesh = 100;
syms X;
%Polynomials/Linear Functions
inBetweenGraph(((X-1)^4)/8+X,maxGridX,maxGridY,gridMesh,derivativeCount,integralCount,lineThick,"((x-1)^4)/8+x")
inBetweenGraph(X,maxGridX,maxGridY,gridMesh,derivativeCount,integralCount,lineThick,"x")
inBetweenGraph(abs(X)^(2/3),maxGridX,maxGridY,gridMesh,derivativeCount,integralCount,lineThick,"x^{2/3}")
%Trigonometric Functions
inBetweenGraph(sin(X),maxGridX,maxGridY,gridMesh,derivativeCount,integralCount,lineThick,"sin(x)")
inBetweenGraph(asin(X),maxGridX,maxGridY,gridMesh,derivativeCount,2,lineThick,"arcsin(x)")
inBetweenGraph(acos(X),maxGridX,maxGridY,gridMesh,derivativeCount,2,lineThick,"arccos(x)")
inBetweenGraph(atan(X),maxGridX,maxGridY,gridMesh,derivativeCount,integralCount,lineThick,"arctan(x)")
inBetweenGraph(sinh(X),maxGridX,maxGridY,gridMesh,derivativeCount,integralCount,lineThick,"sinh(x)")
inBetweenGraph(tanh(X),maxGridX,maxGridY,gridMesh,derivativeCount,0,lineThick,"tanh(x)")
%Exponential And Logarithmic Functions
inBetweenGraph(log(X),maxGridX,maxGridY,gridMesh,derivativeCount,integralCount,lineThick,"ln(x)")
inBetweenGraph(exp(-1*(X)^2),maxGridX,maxGridY,gridMesh,derivativeCount,integralCount,lineThick,"e^{-x^2}")
inBetweenGraph(exp(X),maxGridX,maxGridY,gridMesh,derivativeCount,integralCount,lineThick,"e^x")
%{
inBetweenGraph(abs(X)^X,maxGridX,maxGridY,gridMesh,derivativeCount,-1,lineThick,"|x|^x")
inBetweenGraph((abs(X)^X)/X,maxGridX,maxGridY,gridMesh,derivativeCount,-1,lineThick,"|x|^x/x")
inBetweenGraph(X/(abs(X)^X),maxGridX,maxGridY,gridMesh,derivativeCount,-1,lineThick,"x/{|x|^x}")
%}
%Rational Functions
inBetweenGraph(1/X+X,maxGridX,maxGridY,gridMesh,derivativeCount,integralCount,lineThick,"1/x + x")
%inBetweenGraph(factorial(abs(X)),maxGridX,maxGridY,gridMesh,-1,-1,lineThick,"x!")

function inBetweenGraph(inputFunction, xMax, yMax, meshResolution, derivatives,integrals,lineThickness,titleFunction)
    x = linspace(-xMax,xMax,meshResolution);
    x = double(x);
    figure;
    derivativePanel = ones(length(x),1);
    FUNC = inputFunction;
    y=NaN(length(x),length(x));
    xHist = [];
    yHist = [];
    derivativeHist = [];

    for i = 0:derivatives+1
        if i > 0
            FUNC = diff(FUNC);
        end
        derivativePanel = ones(length(x),1);
        derivativePanel = i.*derivativePanel;
        xHist = [xHist, x];
        yHist = [yHist, subs(FUNC,x)];
        derivativeHist = [derivativeHist, [derivativePanel]'];
        plotThick = lineThickness;
        if i == 0 
            plotThick = 2*lineThickness;
        end
        plot3(x,derivativePanel,subs(FUNC,x),'LineWidth',plotThick);
        hold on;
    end
    FUNC = inputFunction;
    for i = 1:integrals+1
        FUNC = int(FUNC);
        derivativePanel = ones(length(x),1);
        derivativePanel = -1.*i.*derivativePanel;
        xHist = [x,xHist];
        yHist = [subs(FUNC,double(x)),yHist];
        derivativeHist = [[derivativePanel]',derivativeHist];
        plot3(x,derivativePanel,subs(FUNC,x),'LineWidth',lineThickness);
        hold on;
    end
    planeX = [x(1),x(length(x))];
    planeY = [[-yMax/2,yMax]',[-yMax/2,yMax]'];
    for i = -1.*integrals:derivatives
        planeDerivative = ones(2,1);
        planeDerivative = i.*planeDerivative;
        if i == 0
            surf(planeX,planeDerivative,planeY,'MarkerSize',0.01,'FaceAlpha',0.375)
        else
            surf(planeX,planeDerivative,planeY,'MarkerSize',0.01,'FaceAlpha',0.0625)
        end
    end
    view(60,30)
    
    [xData, yData, zData] = prepareSurfaceData( xHist,derivativeHist,yHist);

    % Set up fittype and options.
    ft = 'cubicinterp';

    % Fit model to data.
    [fitresult, gof] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );

    % Plot fit with data.
    h = plot(fitresult);
    h.FaceAlpha = 0.5;
    %legend( h, 'Between Derivatives', 'yHist vs. xHist, derivativeHist', 'Location', 'NorthEast', 'Interpreter', 'none' );
    grid off;
    rotate3d on;
    xlim([-xMax,xMax]);
    zlim([-yMax/2,yMax]);
    ylim([-1*integrals,derivatives]);
    xlabel("X")
    ylabel("Derivative Count")
    zlabel("Y")
    title("In Between The Integrals And Derivatives Of: f(x) = " + titleFunction)
    hold off;
end
function [fitresult, gof] = createFit(x, derivative, y)
end