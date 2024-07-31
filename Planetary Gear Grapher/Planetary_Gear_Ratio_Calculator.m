clc;
clear;
%%Planetary Gear Ratio Calculator
meshSize = 50;
minimumTeeth = 5;
Sun = 15;
Planet = 21;
Ring = 2*Planet + Sun;
S = minimumTeeth:1:meshSize;
P = minimumTeeth:1:meshSize;
Graphs = tiledlayout(1,3);
Graphs.Title.String = "All Possible Planetary Gear Ratios Graphed";
dotSize = 8.75;
desiredRatio = 4;

nexttile
heldPlanet(S,P,minimumTeeth,meshSize,dotSize,Sun,Planet,desiredRatio)
nexttile
heldSun(S,P,minimumTeeth,meshSize,dotSize,Sun,Planet)
nexttile
heldCarrier(S,P,minimumTeeth,meshSize,dotSize,Sun,Planet)
figure;
heldPlanet(S,P,minimumTeeth,meshSize,dotSize,Sun,Planet,desiredRatio)
figure;
heldSun(S,P,minimumTeeth,meshSize,dotSize,Sun,Planet)
figure;
heldCarrier(S,P,minimumTeeth,meshSize,dotSize,Sun,Planet)



%Sliders And Inputs
    
%Ring held in place
function [heldPlanetGraph] = heldPlanet(sun,planet,minimum,size,dot,input_sun,input_planet,desired)
    ratio = 0;
    ratioList = [];
    title("Ring Gear Held Stationary");
    xlabel("Sun Teeth")
    ylabel("Planet Teeth")
    zlabel("Output Ratio")
    rotate3d on;
    view(45,22.5)
    for i = 1:size-minimum
        hold on;
        for ii = 1:size-minimum
            ratio = (2*planet(i)+2*sun(ii))/sun(ii);
            ratioList = [ratioList,ratio];
            plot3(sun(ii),planet(i),ratio,"b--.",'MarkerSize',dot)
            if ratio == desired
                fprintf('ratio of: %f Planet: %f Sun: %f \n',ratio,planet(i),planet(ii))
            end
        end
    end
    [s,p] = meshgrid(sun,planet);
    ratioMesh = (2.*s+2.*p)./s;
    xlim([minimum,size])
    ylim([minimum,size])
    zlim([min(ratioList),max(ratioList)])
    surf(s,p,ratioMesh)
    alpha 0.125
    hold off;
    text = "Ring Held"
    ratio = (2*input_planet+2*input_sun)/input_sun
    
end

%Sun held in place
function [heldSunGraph] = heldSun(sun,planet,minimum,size,dot,input_sun,input_planet)
    ratio = 0;
    ratioList = [];
    title("Sun Gear Held Stationary");
    xlabel("Sun Teeth")
    ylabel("Planet Teeth")
    zlabel("Output Ratio")
    rotate3d on;
    view(60,2.5)
    for i = 1:size-minimum
        hold on;
        for ii = 1:size-minimum
            ratio = (2*planet(i)+2*sun(ii))/(2*planet(i)+sun(ii));
            ratioList = [ratioList,ratio];
            plot3(sun(ii),planet(i),ratio,"r--.",'MarkerSize',dot)
        end
    end
    [s,p] = meshgrid(sun,planet);
    r = 2.*p+s;
    ratioMesh = (r+s)./r;
    xlim([minimum,size])
    ylim([minimum,size])
    zlim([min(ratioList),max(ratioList)])
    surf(s,p,ratioMesh)
    alpha 0.125
    hold off;
    text = "Sun Held"
    ratio = (2*input_planet+2*input_sun)/(2*input_planet+input_sun)
end
%Carrier held in place
function [heldCarrierGraph] = heldCarrier(sun,planet,minimum,size,dot,input_sun,input_planet)
    ratio = 0;
    ratioList = [];
    title("Planet Carrier Held Stationary");
    xlabel("Sun Teeth")
    ylabel("Planet Teeth")
    zlabel("Output Ratio")
    rotate3d on;
    view(45,-22.5)
    for i = 1:size-minimum
        hold on;
        for ii = 1:size-minimum
            ratio = -1*(2*planet(i)+sun(ii))/sun(ii);
            ratioList = [ratioList,ratio];
            plot3(sun(ii),planet(i),ratio,"g--.",'MarkerSize',dot)
        end
    end
    [s,p] = meshgrid(sun,planet);
    r = 2.*p+s;
    ratioMesh = -r./s;
    xlim([minimum,size])
    ylim([minimum,size])
    zlim([min(ratioList),max(ratioList)])
    surf(s,p,ratioMesh)
    alpha 0.125
    hold off;
    text = "Carrier Held"
    ratio = -1*(2*input_planet+input_sun)/input_sun
end

