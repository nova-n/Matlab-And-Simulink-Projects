clear;
clc;
E=2; 
D=3;
Re = 5;
resolution = 0.011;

fFriction = coleBrook( E, D, Re, resolution)

function [frictionFactor] = coleBrook(rough,diam,reynolds,resol)
    closeEnough = false;
    fFric = 0;
    tol = 1e-12;%for floating point numbers
    while closeEnough == false
        fFricPrev = fFric;
        if abs(1/sqrt(fFric) + 2*log10((rough/diam)/3.7+2.51/(reynolds*sqrt(fFric)))) > tol %checking inequality with floating point numbers
            fFric = fFric + resol;
        end
        fFric;%delete the semicolon to see in real time
        %check if leftside - rightside < tolerance, OR 
        %if it overshot, to use the previous value
        if abs(1/sqrt(fFric) + 2*log10((rough/diam)/3.7+2.51/(reynolds*sqrt(fFric)))) < tol || abs(1/sqrt(fFric)...
                + 2*log10((rough/diam)/3.7+2.51/(reynolds*sqrt(fFric)))) > abs(1/sqrt(fFricPrev)...
                + 2*log10((rough/diam)/3.7+2.51/(reynolds*sqrt(fFricPrev)))) 
            closeEnough = true;
            if abs(1/sqrt(fFric) + 2*log10((rough/diam)/3.7+2.51/(reynolds*sqrt(fFric)))) > abs(1/sqrt(fFricPrev)...
                    + 2*log10((rough/diam)/3.7+2.51/(reynolds*sqrt(fFricPrev))))
                fFric = fFricPrev;
            end
        end   
    end
     frictionFactor = fFric;
end