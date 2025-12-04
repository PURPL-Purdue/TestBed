%% Main Compiled Code
%widthArray = linspace(0.02/39.37, 0.040/39.37, 10); %m %channel width sweep %CHECK WITH LITERATURE
%heightArray = linspace(0.04/39.37, 0.125/39.37, 10); %m %channel height sweep %CHECK WITH LITERATURE

clc; clear;

m_to_in = 39.3701;

widthArray = [0.02,0.02,0.02] / m_to_in ;% updatedHeightValues(7,2:46);
heightArray =  [0.125,0.03,0.125] / m_to_in;% updatedHeightValues(6,2:46);
wall_thicknessMatrix = [0.05,0.05,0.05] / m_to_in;
heightStepNumber = 45;
numChannels = 60;
converge_index = 23;
throat_index = 14;
generate_new_CEA = false;

p_c = 250; % Main chamber pressure (psi)
OF = 1; % OF Ratio (N/A)
m_dot = 1.2637083; % Total mass flow (kg/s)
contour_name = "Engine Contour Cleaned and Sorted (Metric).csv";
engineContour = readmatrix(contour_name);

T_start= 298; % Flow Initial Temp in degrees K 
P_start = 5102000; % Flow initial Pressure in Pa
rho_start = 810; % Coolant initial density in kg/m^3
m_flow_total = m_dot / (1 + OF); % total coolant mass flow in kg/s
mass_flow = m_flow_total/numChannels; % Precalcuated mass flow based on # of channels in Malestrom
T_target = 400; % target gas-side hotwall temp in degrees K (530 for 7075, 773 for copper)

k_w = 130; % thermal conductivity of the wall (W/m*K) %copper 401, 7075 130
surfaceRoughness = 0.000002; % average height of roughness (chosen from engineering toolbox/elementum) in m
CTE = 0.0000232; % Material's coefficient of thermal expansion in (%change/K)
youngsModulus = 71700000000; %Pa
poissonsRatio = 0.33; %

idx = engineContour(:,3) == 1;        % Look for the point in the chamber where the area ratio is 1
throatDiameter = engineContour(idx,5) * 2 * m_to_in;  % Obtain the throat diameter value (in)

chamberRad = 0.23; % chamber converging radius in (in)
chamberLength = floor((engineContour(end,4) - (engineContour(1,4))) * m_to_in*100) / 100; % Chamber length (in), contencated 

inputValues = [T_start, P_start, rho_start, mass_flow, T_target, k_w, numChannels, surfaceRoughness, CTE, youngsModulus, throatDiameter, chamberRad,poissonsRatio,chamberLength];

%% Initialize all arrays and matrices,
flowTempArray = zeros(1,heightStepNumber); %Matrices to store all pressure,velocity and temp data from calculateWallTemp
flowVelocityArray = zeros(1,heightStepNumber);
flowPressureArray = zeros(1,heightStepNumber);

%% Height Step initialization % Not sure if this works, may scrap for even height steps (worked with PSP data)

heightStepArray = linspace(0,chamberLength/39.37,heightStepNumber);

%% Run NASA CEA and retrieve values

if generate_new_CEA == true
    CEAOut(p_c,OF,contour_name)
end

fluidProperties = readmatrix("CEA_Maelstrom.xlsx"); %pull all nasaCEA values into fluidProperties
% if newFluidProperties errors and cuts off a row, change the middle value
% in the heightStepArray initialization call to be whatever the ACTUAL end
% length is set to.
fluidProperties(1,:) = [];
y = 1; 
r = 1;
axialDist = (fluidProperties(:,1));
newFluidProperties = zeros(length(heightStepArray),10);
newFluidProperties(:,1) = heightStepArray;
%chamberPlot = readmatrix("Engine Contour Cleaned and Sorted (Metric).csv");
T_l_reqMatrix = [];
updatedTemps = [];
updatedPressure = [];
updatedVelocity = [];
%heightMatrix = zeros(heightStepNumber);

%% New fluid Properties
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


throatArea = pi*(throatDiameter/(2*39.37))^2;
chamberDiameter1 = 2*sqrt((throatArea*newFluidProperties(:,2))/pi);
newFluidProperties = flip(newFluidProperties,1);
chamberDiameter = flip(chamberDiameter1);
% diametercontour = engineContour(:,5);
%% Calculate Wall Temp
[flowTempArray,flowVelocityArray, flowPressureArray,T_wgFinal, finEfficiency, Qdot, finQdot, T_wl_Array, h_l_Array, h_g_Array, vonMises,sigma_long, sigma_circ, sigma_rad] = FINAL_CALCWALLTEMP(converge_index, throat_index, heightArray, widthArray, wall_thicknessMatrix, chamberDiameter, heightStepNumber, newFluidProperties, inputValues);