%% Main Compiled Code
updatedHeightValues = readmatrix("wallThicknessesGoated.xlsx");
%widthArray = linspace(0.02/39.37, 0.040/39.37, 10); %m %channel width sweep %CHECK WITH LITERATURE
%heightArray = linspace(0.04/39.37, 0.125/39.37, 10); %m %channel height sweep %CHECK WITH LITERATURE
widthArray = [0.059,0.0394,0.122]/39.37;%updatedHeightValues(7,2:46);
heightArray =  [0.091,0.0394,0.059]/39.37;%updatedHeightValues(6,2:46);
wall_thicknessMatrix = [0.0394,0.0197,0.0394]/39.37;
heightStepNumber = 45;
numChannels = 62;
startConvInd = 28;
throatInd = 36;
%% Initialize all arrays and matrices
flowTempMatrix = zeros(heightStepNumber); %Matrices to store all pressure,velocity and temp data from calculateWallTemp
flowVelocityMatrix = zeros(heightStepNumber);
flowPressureMatrix = zeros(heightStepNumber);

%% Height Step initialization % Not sure if this works, may scrap for even height steps (worked with PSP data)

heightStepArray = linspace(0,0.227382320000000,heightStepNumber);

%% Run NASA CEA and retrieve values
fluidProperties = readmatrix("CEAOutFzPsp_IAC_11-23-25.xlsx"); %pull all nasaCEA values into fluidProperties
% if newFluidProperties errors and cuts off a row, change the middle value
% in the heightStepArray initialization call to be whatever the ACTUAL end
% length is set to.
fluidProperties(1,:) = [];
y = 1; 
r = 1;
axialDist = (fluidProperties(:,1));
newFluidProperties = zeros(length(heightStepArray),10);
newFluidProperties(:,1) = heightStepArray;
chamberDiameter = [];
%chamberPlot = readmatrix("Engine Contour Cleaned and Sorted (Metric).csv");
T_l_reqMatrix = [];
updatedTemps = [];
updatedPressure = [];
updatedVelocity = [];
%heightMatrix = zeros(heightStepNumber);
while y <= length(heightStepArray) % translating CEA outputs to height step number length output by averaging values over height step number
    
    
    a = r; % MAY NEED TO CHANGE BASED ON WHAT GETS READ FROM EXCEL FILE (add 2 or something to accoutn for text)
    
    %sumDiameter = 0;
    sumAEAT = 0;
    sumPrandtl = 0;
    sumMach = 0;
    sumGamma = 0;
    sumT = 0;
    sumVisc = 0;
    sumCp = 0;
    sumP = 0;
    sumCstar = 0;                

    while a <= length(axialDist)
        
        if a-r==0
            %sumDiameter = sumDiameter + chamberPlot(a,2);
            sumAEAT = sumAEAT+fluidProperties(a,2);
            sumPrandtl = sumPrandtl+fluidProperties(a,3);
            sumMach = sumMach+fluidProperties(a,4);
            sumGamma = sumGamma+fluidProperties(a,5);
            sumT = sumT+fluidProperties(a,6);
            sumVisc = sumVisc+fluidProperties(a,7);
            sumCp = sumCp+fluidProperties(a,8);
            sumP = sumP+fluidProperties(a,9);
            sumCstar = sumCstar+fluidProperties(a,10);
            a=a+1;

        elseif axialDist(a) < heightStepArray(y)
            %sumDiameter = sumDiameter + chamberPlot(a,2);
            sumAEAT = sumAEAT+fluidProperties(a,2);
            sumPrandtl = sumPrandtl+fluidProperties(a,3);
            sumMach = sumMach+fluidProperties(a,4);
            sumGamma = sumGamma+fluidProperties(a,5);
            sumT = sumT+fluidProperties(a,6);
            sumVisc = sumVisc+fluidProperties(a,7);
            sumCp = sumCp+fluidProperties(a,8);
            sumP = sumP+fluidProperties(a,9);
            sumCstar = sumCstar+fluidProperties(a,10);


            
            a=a+1;
        else
            divFactor = a-r;
            %chamberDiameter(y,1) = 2*(sumDiameter/divFactor);
            newFluidProperties(y,2) = sumAEAT/divFactor;
            newFluidProperties(y,3) = sumPrandtl/divFactor;
            newFluidProperties(y,4) = sumMach/divFactor;
            newFluidProperties(y,5) = sumGamma/divFactor;
            newFluidProperties(y,6) = sumT/divFactor;
            newFluidProperties(y,7) = sumVisc/divFactor;
            newFluidProperties(y,8) = sumCp/divFactor;
            newFluidProperties(y,9) = sumP/divFactor;
            newFluidProperties(y,10) = sumCstar/divFactor;
            r = a;
            break; 
        
        end
    end
    y=y+1;
end


i = 1;

for a = heightStepArray
    
    if((i)<=startConvInd)
        chamberDiameter(i) = 0.09525; % set diameter in m
        
    elseif((i)<=throatInd)
        chamberDiameter(i) = (((heightStepArray(i)-0.132334)*-1.245) +0.09525);
    else
        chamberDiameter(i) = (((heightStepArray(i)-0.178816)*0.5129) +0.037338);
    end

    i=i+1;
end


chamberDiameter = flip(chamberDiameter)';
newFluidProperties = flip(newFluidProperties,1);
%% Main Loop
% for widthValue = 1:length(widthArray) %width value sent to calculateWallTemp from width array
% 
%         width = widthArray(widthValue);
%   
% 
        [flowTempMatrix,flowVelocityMatrix, flowPressureMatrix,T_l_reqMatrix, wall_thicknessMatrix, updatedTemps, updatedPressure, updatedVelocity] = calculateWallTemp3(updatedTemps, updatedPressure, updatedVelocity, heightArray,T_l_reqMatrix, chamberDiameter,wall_thicknessMatrix,numChannels, heightStepArray, flowTempMatrix, flowVelocityMatrix, flowPressureMatrix, widthArray, newFluidProperties);
        %[flowTemp,flowVelocity,flowPressure, T_l_reqMatrix, wall_thicknesses, updatedTemps,updatedPressure,updatedVelocity, heightMatrix] = calculateWallTemp2(updatedTemps, updatedPressure, updatedVelocity, heightMatrix, heightArray,T_l_reqMatrix, chamberDiameterArray, wall_thicknesses,channelNum, heightStepArray, flowTempMatrix, flowVelocityMatrix, flowPressureMatrix, widthArray, newFluidProperties)
        %Flow Temp, Pressure, Velocity are outputted arrays which contain values for *1* channel dimension combination
        
        % if flowTempMatrix(widthValue,heightValue,length(heightStepArray)) == -1 || flowTempMatrix(heightArray,widthValue,heightValue,length(heightStepArray)) == 0
        %         geometryMap(widthValue, heightValue) = 0;%fail
        %     else
        %         geometryMap(widthValue, heightValue) = 1; %pass
        % end
        %wallthicknessesGoated = wall_thicknessMatrix(1,10,:);
        
% 
% 
% 
% end

%create geometry map
        
           
