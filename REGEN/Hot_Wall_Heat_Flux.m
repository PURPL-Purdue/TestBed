clear;
clc;
figure(1);
clf;
figure(2);
clf;
figure(3);
clf;
pause(0.1)

%% Inputs and Variable Definitions
config = (readmatrix("Config_Hot_Wall_Heat_Flux.xlsx")); %Config file for ease of change
config = config(:,2); 
fast = 1; %Specifies fast or slow algorithum

%Config and value imports
safetyFactor = config(1); %Safety Factor
targetMaxTemp = config(2)/safetyFactor; %Target max temp (K)
materialDensity = config(3); %Material density (kg/m^3)
materialSpecificHeat = config(5); %Material specific heat (j/kg*k)

time = config(7);%Time length of heat transfer iterations until convirgence
radialSliceThickness = config(8); %Iteration refinement value (meters)
deltaT = config(6); %Iteration refinement value (seconds)
staticTemp = 250; %Referance value do not change
K = config(4); %(Material thermal conductivity)
QDotStep = config(9); %Qdot step (J/s)
convirganceThreshold = config(11)/100; %Percent change convergance criteria
minimumIterations = config(10); %Minimum iterations until convergence threshold is active

engineContourTemp = readmatrix("Engine Contour Cleaned and Sorted (Metric).csv"); %Imports the engine contour
engineContour = engineContourTemp; 

targetTemps = readmatrix("HotWall Temperature Distribution.xlsx"); %Imports target temperature distribution
targetTemps = targetTemps(:,2);

%targetQDots = readmatrix("HotWall QDot Distribution.xlsx"); %Optional Qdot read in to refine output with smaller Qdot steps
targetQDots = zeros(998); %Initializes an empty initial Qdot distribution
targetQDots = targetQDots(:,2);

wallThickness = 0.002; %Meters !!!To be replaced by stress calcs!!!

%Output instantiation
steadyStateHsubGDistribution = zeros(length(engineContour(:,1))); %Heat transfer coeffecient matrix

%CEA Wrapper output read in
fluidProperties = readmatrix("CEAOutFrozen.xlsx");
%fluidProperties = readmatrix("CEAOutEquilibriumFAC.xlsx");
fluidProperties = fluidProperties(3:end,:);
fluidProperties = double(fluidProperties);

%% Iteration Control
SSSF = 0;
loops = 1;
while SSSF ~= 1
    hotwallTemps = ones(length(engineContour(:,1)));
    e = 1;
    while e <= length(hotwallTemps)
        hotwallTemps(e) = 272;
        e = e + 1;
    end
%% Convirgance Structure
    maxTemps = [];
    axialNum = 1;
    bigArray = zeros(time / deltaT, wallThickness / radialSliceThickness);
    convirgenceNums = [];
    while axialNum <= length(engineContour(:,1))
        if loops == 1 || incompliance(axialNum) == 1 
            currentRadius = engineContour(axialNum,2);
            currentTemperature = hotwallTemps(axialNum);
            axialPosition = engineContour(axialNum);
            currentQDot = targetQDots(axialNum); %To do: Make all input excel files have same length as contour   
            if axialNum == length(engineContour(:,1))
                length_thickness = 0.0002;
            else
                length_thickness = engineContour(axialNum + 1,1) - engineContour(axialNum,1);
            end
            if axialNum == 700
                pause(0.1)
            end
            %Fluid property recall
            hit = 1;
            i=1;
            while i <= length(fluidProperties(:,1))
               if fluidProperties(i,1) >= axialPosition
                   hit = i;
                   break;
               end 
               i = i+1;
            end
            currentProperties = fluidProperties(hit,:)';
            gasTemperature = currentProperties(6);
   
            %Layer volumes initialization and calculation (m^3)
            radialAngle = 2 * 3.14159;
            numLayers = wallThickness / radialSliceThickness; %Calcs Number of layers
            layers = [currentRadius:radialSliceThickness:(wallThickness + currentRadius)]; %Makes vector of layers/diams of layers
            layer_vols = [(radialAngle/2) * (((currentRadius + radialSliceThickness)^2 - (currentRadius)^2)) * length_thickness]; %Calcs first layer volume
            a = 2;
            while a < numLayers + 1
                layer_vols(end+1) = (radialAngle / 2) * length_thickness * ( (layers(a))^2 - (layers(a-1))^2 ); %Layer Volume calculation
                a = a + 1;
            end

            %Calculation of layer masses and surface areas
            layerMasses = layer_vols * materialDensity; %Layer masses initialization
            layer_inner_SAs = layers .* radialAngle .* length_thickness; %SA initialization

            %Qs initialization
            %Note: Qs is an array of thermal energies for each layer in joules.
            Qs = [[layerMasses .* materialSpecificHeat .* staticTemp]]; %Fills first row of Qs with initial Q
            Qs(1,1) = currentTemperature * layerMasses(1) * materialSpecificHeat; %Prepares next row for calculations
            
            %Heat Propigation Control Loop
            time_num = 1;
            convergence = 0;
            %while time > time_num * deltaT
            while convergence == 0  %Convirgance based analysis
                layer_num = 2;
                
                %Convective heat transfer calculation
                Qs(time_num + 1,:) = 0; %Initializes the current row
                h_g = H_g_From_Temperature(currentTemperature, currentProperties);
                init_temp = Q_to_heat((Qs(time_num,1)),layerMasses(1),materialSpecificHeat);
                convection_Qdot = h_g * layer_inner_SAs(1) * (init_temp - gasTemperature); %Calculates the heat transfer rate
                convection_delta_Q = -convection_Qdot * deltaT; %Calculates the change in thermal energy of the innermost layer over one time step
                Qs(time_num+1,1) = Qs(time_num,1) + convection_delta_Q; %Adds energy change to the layer
            
                %Conductive heat transfer calculations
                while numLayers + 1 > layer_num
                    t_hot =  Q_to_heat(Qs(time_num+1, layer_num-1), layerMasses(layer_num-1), materialSpecificHeat); %Temperature of hot layer
                    t_cold = Q_to_heat(Qs(time_num, layer_num), layerMasses(layer_num), materialSpecificHeat); %Temperature of cold layer
                    delta_Q = Qdot(K, layer_inner_SAs(layer_num),t_hot, t_cold, radialSliceThickness) * deltaT; %Calculation of change in thermal energy between two layers
                    Qs(time_num+1,layer_num) = Qs(time_num,layer_num) + delta_Q; %Heat energy gained by cold layer
                    Qs(time_num+1,layer_num-1) = Qs(time_num + 1,layer_num-1) - delta_Q; %Heat energy lost from hot layer
                    layer_num = layer_num + 1; %Increments layer 
                end
            
                %Convection Q subtract
                Qs(time_num+1, end) = Qs(time_num+1,end) - (targetQDots(axialNum) * deltaT);
                
                %Test for convirgance
                la = Q_to_heat_matrix(Qs(time_num,:), layerMasses(:), materialSpecificHeat); 
                previousTemperatures = la(:,1);
                el = Q_to_heat_matrix(Qs(time_num+1,:), layerMasses(:), materialSpecificHeat); 
                newTemperatures = el(:,1);
                if all(abs((previousTemperatures-newTemperatures)./previousTemperatures) < convirganceThreshold) && minimumIterations < time_num
                    convergence = 1;
                    convirgenceNums(end + 1) = time_num;
                end
                time_num = time_num + 1; %Increment time number
            end

            % Conversion of Thermal Energies to Temperatures
            
            %Conversion loop
            col = 1;
            row = 1;
            final_temperatures = [];
            while row < length(Qs(:,1)) + 1
                while col < width(Qs) + 1
                    final_temperatures(row, col) = Q_to_heat(Qs(row,col), layerMasses(col), materialSpecificHeat); %Temperature of single layer
                    col = col + 1;
                end
                col = 1;
                row = row + 1;
            end

            % Final Output
            z = 1;
            while z <= length(final_temperatures(:,1))
                bigArray(z,:,axialNum) = final_temperatures(z,:);
                z = z + 1;
            end
            steadyStateHsubGDistribution(axialNum,loops) = h_g;
        else
            bigArray(:,:,axialNum) = ones(length(bigArray(:,1,1)),length(layers)-1,1);
            steadyStateHsubGDistribution(axialNum,loops) = steadyStateHsubGDistribution(axialNum,loops-1);
            convirgenceNums(end + 1) = 1;
        end
        axialNum = axialNum + 1;
    end
    %Search for convergence temp dist
    s = 1;
    
    while s <= length(bigArray(1,1,:))
        if loops ~= 1 && not(incompliance(s))
            maxTemps = [maxTemps tempMaxtemps(s)];
        else
            maxTemps = [maxTemps bigArray(convirgenceNums(s),1,s)];
        end
        s = s + 1;
    end

    incompliance = (maxTemps > targetTemps');
    incompliance = incompliance(:)';
    if fast == 0
        alterations = (incompliance).* QDotStep;
        targetQDots = targetQDots + alterations;
    else
        fastChanges = (3/1000).*(maxTemps - targetTemps') + 0.1;%8./(1+2.71828.^-(((1/100).*(maxTemps - targetTemps'))-22))+1+(1/1000).*(maxTemps - targetTemps');
        alterations = (incompliance .* fastChanges).* QDotStep;
        targetQDots = targetQDots + alterations';
    end
    %Test for 
    if all(alterations == 0)
        SSSF = 1;
    end
    writematrix([engineContour(:,1), targetQDots], "testtwo.csv")
    tempMaxtemps = maxTemps .* not(incompliance);

    figure(1)
    hold on;
    grid on;
    plot(engineContour(:,1),targetQDots)
    title("Heat Flux")
    xlabel("Axial Position (m)")
    ylabel("QDot (J/s)")

    figure(2)
    plot(engineContour(:,1), maxTemps,"Color","red")
    hold off;
    grid on;
    title("Current Temperature Distribution")
    xlabel("Axial Position (m)")
    ylabel("Temperature (K)")
    
    figure(3)
    plot(engineContour(:,1), steadyStateHsubGDistribution(:,loops))
    hold off;
    grid on;
    title("Heat Transfer Coeffecient")
    xlabel("Axial Position (m)")
    ylabel("Heat Transfer Coeffecient (W/m^2K)")
    
    

    fprintf("\n Num Loops: %f", loops)
    pause(0.1)
    loops = loops + 1;
    fprintf("\n")
end