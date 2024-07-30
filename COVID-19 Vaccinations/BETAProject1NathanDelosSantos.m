clc;
clear;

%%Project 1 Nathan Delos Santos
    population = readtable("cali_county_pop.csv");
    vaccines = readtable("covid19vaccinesbycounty.csv");
    partsPer = 1000000;%Parts per Million in this case
    mov = 7;%days
    vaccines.Properties.VariableNames(15) = "New_People_With At_Least_One_Dose";
    vaccines.Properties.VariableNames(13) = "New_People_Fully_Vaccinated";
    vaccines.Properties.VariableNames(3) = "Total_Doses";
    vaccines.Properties.VariableNames(4) = "Cumulative_Total_Doses";
    vaccines.Properties.VariableNames(6) = "Cumulative_Pfizer_Doses";
    vaccines.Properties.VariableNames(8) = "Cumulative_Moderna_Doses";
    vaccines.Properties.VariableNames(10) = "Cumulative_J&J_Doses";

%%Top Populous Counties
    topNumber = 5;                                                          %Sets how many counties to look at based on population
    topPop = sortrows(population, "Population", "descend");                 %Puts the most populated counties first
    topPop([topNumber+1:height(topPop)],:) = [];                            %Everything after the top few counties is erased
    topPopNames = topPop.County;                                            %Gathers a list of names of those counties. The previous line gathered ALL information about the county.
    topPopNames = categorical(topPopNames);                                 %makes the array usable elsewhere
    index = [];                                                             %Creating an array of the indexes of the counties
    for i = 1:topNumber
        index(i) = max(find(vaccines.county == topPopNames(i)));            %Finds the most recent index of the county in the vaccine table
    end
    
%%Calling Chart Functions
    barGraphOutput(index,topPopNames,["Cumulative_Total_Doses"],topNumber,vaccines,"Total Vaccinations")
    barGraphOutput(index,topPopNames,["Cumulative_Pfizer_Doses","Cumulative_Moderna_Doses","Cumulative_J&J_Doses"],topNumber,vaccines,"Total Vaccinations By Manufacturer")
    percentVax(index,topPopNames,topPop,["Cumulative_Total_Doses"],topNumber,vaccines,"")
    rollingAvgChart(["New_People_With At_Least_One_Dose","New_People_Fully_Vaccinated","Total_Doses"],vaccines,topPopNames,topPop,topNumber,partsPer,mov)
    
%%Doses Bar Charts
    function [barOutput] = barGraphOutput(countiesIndex,countyNames,manufacturers,number,vaxTable,titles) 
    %Function that returns a bar graph for a "manufacturer" for each county. Takes the inputs: index of counties, county names, manufacturers, number of counties, vaccine table, and a title
        cumulative = [];
        vaxName = [];
        colors = [];
        for i = 1:length(manufacturers)                                     
            vaxName=[vaxName,strrep(manufacturers(i),"_"," ")];             %Removes the underscores "_" from the manufacturer names for the legend
            colors=[colors,[rand,rand,rand]'];                              %Randomly generates colors for the many possible manufacturers
        end
        figure;
        for i = 1:number                                                    %Plots Bars in groups of (how many manufacturers), for each county
            for ii = 1:length(manufacturers)
                cumulative(ii) = vaxTable{countiesIndex(i),manufacturers(ii)}; %For this one specific county, documents the cumulative for each manufacturer
            end
            b=bar(countyNames(i),cumulative);                               %Adds (how many manufacturers) bars at a time
            for c = 1:length(manufacturers)                                 %Assigns each bar one of (how many manufacturers) colors.
                b(c).FaceColor = colors(1:3,c);
            end
            legend(vaxName)                                                 %Adds the legend without the underscores
            hold on                                                         %Allows for the other counties' bars to be displayed on the same plot
            legend show
        end
        title(titles + " By County" + " Up To " + datestr(max(vaxTable.administered_date))+ " (Top" + sprintf(" %.0f ",number)+ "Most Populous)") %Adds the title. Lets you know how many counties, and most recent date.
        ax=gca;
        ax.YGrid = 'on';
        ax.YAxis.Exponent = 0;
        ytickformat("%.0f");
        hold off;
    end                                                                     %Repeats for all counties
    
%%Percentage Vaccination
    function [barPercent] = percentVax(countiesIndex,countyNames,popul,manufacturers,number,vaxTable,titles)
    %Function that returns a bar graph for the percentage vaccinated by a "manufacturer" for each county. Takes the inputs: index of counties, county names, county population, manufacturers, number of counties, vaccine table, and a title
        colors = [];
        vaxName = [];
        figure;
        for i = 1:length(manufacturers)
            vaxName=[vaxName,strrep(manufacturers(i),"_"," ")];             %Removes the underscores "_" from the manufacturer names for the legend
            colors=[colors,[rand,rand,rand]'];                              %Randomly generates colors for the many possible manufacturers
        end
        for i = 1:number                                                    %Plots Bars in groups of (how many manufacturers), for each county
            for ii = 1:length(manufacturers)
                cumulative(ii) = 100*vaxTable{countiesIndex(i),manufacturers(ii)}/popul.Population(i); %For this one specific county, documents the cumulative percentages for each manufacturer
            end
            b=bar(countyNames(i),cumulative);                               %Adds (how many manufacturers) bars at a time
            for c = 1:length(manufacturers)                                 %Assigns each bar one of (how many manufacturers) colors.
                b(c).FaceColor = colors(1:3,c);                             
            end
            legend(vaxName)                                                 %Adds the legend without the underscores
            hold on                                                         %Allows for the other counties' bars to be displayed on the same plot
            legend show
        end
        title("Vaccination Percentage By County" + titles + " Up To " + datestr(max(vaxTable.administered_date))+ " (Top" + sprintf(" %.0f ",number)+ "Most Populous)") %Adds the title. Lets you know how many counties, and most recent date.
        ax=gca;
        ax.YGrid = 'on';
        ylim([0,100])
        ytickformat('%g\%')
        hold off;
    end                                                                     %Repeats for all counties
    
%%Rolling Average Line Charts
    function [roll] = rollingAvgChart(factors,vaxTable,countyNames,countyPop,number,parts,moving)
    %Function that returns (how many factors) line plots, displaying the factors on a 7 day moving average for each county
        for i = 1:factors.length                                            %Loops through one factor at a time
            vaxNumbers = [];                                                %Resets the factor numbers after each factor
            for ii = 1:number
                vaxNumbersIndex = [];                                       %Resets the array of indecies after each county
                vaxNumbersIndex = find(vaxTable.county == countyNames(ii)); %Finds all of the indecies for that county for that specific factor. Each index is a day.
                temp = [];                                                  %Resets the temporary variable for each county
                for iii = 1:length(vaxNumbersIndex)                         %Loops through each day for that county for that factor
                    temp(iii) =  [vaxTable{vaxNumbersIndex(iii),factors(i)}]'; %Adds the factor number to the temporary array AS A COLUMN. This way, it is easier for humans to see the divisions for each county.
                    temp(iii) = parts * temp(iii)/countyPop{ii,"Population"}; %Divides that newly added factor by the population, and multiplies it by the "Parts Per" variable
                end
                vaxNumbers(:,ii) = temp;                                    %Adds the column of factor numbers
            end
            strFactor = strrep(factors(i),"_"," ");                         %Removes the underscores "_" from the factors for the legend
            figure;
            title(strFactor + " Per Day Per "+ sprintf("%.0f",parts) + " Up To " + datestr(max(vaxTable.administered_date))+ "(Top" + sprintf(" %.0f ",number)+ "Most Populous)") %Adds the title. Lets you know how many counties, and most recent date.
            ax=gca;
            ax.YGrid = 'on';
            for i = 1:length(countyNames)                                   %For this one factor, plots ALL top populous counties' stats.
                plotVal = movmean(vaxNumbers(:,i),moving);                  %Creates a rolling avereage of the data. You can choose how many days to make a rolling average. I chose 7.
                hold on                                                     %Allows for the other counties' lines to show up on the same plot.
                plot(vaxTable{find(vaxTable.county == countyNames(i)),"administered_date"},plotVal,'DisplayName',string(countyNames(i))) %Plots a county's factor stats per day. 
                legend show
                legend('Location','northwest')
                title(legend,sprintf("%.0f",moving)+' Day Rolling Average','FontSize',12) %Titles the legend, letting you know that it is a rolling average
            end
            xlim([min(vaxTable.administered_date) max(vaxTable.administered_date)]) %Sets the x limits from the minimum date t the latest date
            plotVal = [];                                                   %Resets the rolling average array for the next plot
            hold off                                                        %Ends the plot so that no new lines will be plotted onto the same graph
        end                                                                 %Repeats for all factors
    end