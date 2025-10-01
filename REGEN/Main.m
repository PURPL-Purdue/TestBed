%% Main Compiled Code

widthArray = linspace(0.001, 0.005, 20); %m %channel width sweep %CHECK WITH LITERATURE
heightArray = linspace(0.001, 0.005, 20); %m %channel height sweep %CHECK WITH LITERATURE

heightStepNumber = 67;

numWidth = length(widthArray);
numHeight = length(heightArray);
numSteps = 10;
%% Initialize all arrays and matrices
flowTempMatrix = zeros(numWidth, numHeight, numSteps); 
geometryMap = zeros(numWidth, numHeight, numSteps);  

flowVelocityMatrix = zeros(numWidth, numHeight, numSteps);
flowPressureMatrix = zeros(numWidth, numHeight, numSteps);
wallTempArray = zeros(1, numSteps);

%% Height Step initialization
syms x;
steps = piecewise(x >= 0 & x <= 0.50777934936 * pi,(-2 * sin(x+(0.192 * pi)))+3.14856, x > 0.50777934936 * pi & x <= pi, 3.14856);
    n = pi/heightStepNumber;
    step = 1;
    for i = 0:n:(pi-n)
        heightStepArray(step) = int(steps,i,i+n);
        step = step +1;
    end

%% Main Loop
for widthValue = 1:length(widthArray) %width value sent to calculateWallTemp from width array
    for heightValue = 1:length(heightArray) %heigth value sent to calculateWallTemp from height array
        width = widthArray(widthValue);
        height = heightArray(heightValue);

        [flowTemp,flowVelocity, flowPressure] = calculateWallTemp(numChannels, heightStepArray, flowTempMatrix, flowVelocityMatrix, flowPressureMatrix, height, width, heightValue, widthValue, newFluidProperties);
        %Flow Temp, Pressure, Velocity are outputted arrays which contain values for *1* channel dimension combination
        
        %create geometry map
        for(x = 1:1:length(flowTemp))
            if flowTemp(x) == 999999999999999999999999999999999999999999
                geometryMap(widthValue, heightValue, x) = 0; %fail
            else
                geometryMap(widthValue, heightValue, x) = 1; %pass
            end
        end
    end
end

