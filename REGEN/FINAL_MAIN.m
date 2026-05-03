%% Main Compiled Code

clc; clear;

repoDir = fileparts(mfilename('fullpath'));
parentDir = fullfile(repoDir, '..');
addpath(parentDir);
outputDir = parentDir;

m_to_in = 39.3701;
lb_to_kg = 0.453592;
psi_to_pa = 6894.76;

file_name = "Maelstrom";

yaml_struct = py.yaml.safe_load(fileread(file_name + ".yaml"));
data = struct(yaml_struct);

p_c             = double(data.chamber_pressure);
OF              = double(data.of_ratio);

chamberDiameter = double(data.chamber_diameter);

totalLength     = double(data.total_length);
throatDiameter  = double(data.throat_diameter);
exitDiameter    = double(data.exit_diameter);
convergingAngle = double(data.converging_angle);
divergingAngle  = double(data.diverging_angle);
convergingFillet= double(data.converging_fillet);
throatFillet    = double(data.throat_fillet);

mdot_total = (p_c * double(data.throat_area)) / (double(data.cstar) * 3.28084) * 32.174;
mdot_coolant = mdot_total / (1 + OF);

numChannels = 60;
widthArray = [0.02,0.02,0.02] / m_to_in;
heightArray = [0.125,0.04,0.125] / m_to_in;
wall_thicknessMatrix = [0.04,0.04,0.04] / m_to_in;

T_start= 298;
P_start = 450;
rho_start = 791.26;
T_target = 530;

heightStepNumber = 50;
contourResolution = 620;

generate_new_CEA = true;
generate_new_Contour = true;

%% Material properties
k_w = 130;
surfaceRoughness = 0.0000032;
CTE = 0.0000212;
youngsModulus = 71700000000;
poissonsRatio = 0.33;

mdot_coolant = mdot_coolant * lb_to_kg;
mdot_channel = mdot_coolant/numChannels;

hotwallGeometry = [chamberDiameter, throatDiameter, exitDiameter, ...
                   convergingAngle, divergingAngle, totalLength, ...
                   convergingFillet, throatFillet];

%% Generate contour (SAVED IN PARENT DIR)
if generate_new_Contour
    generateContour(hotwallGeometry, file_name, contourResolution, outputDir);
end

engineContour = readmatrix(fullfile(outputDir, "Contour_" + file_name + ".xlsx"));
idx = find(engineContour(:,3) == 1, 1);

throatDiameter = engineContour(idx,2) * 2 * m_to_in;
chamberLength = floor((engineContour(end,1)) * m_to_in * 100) / 100;
filletRad = convergingFillet;

%% Arrays
heightStepArray = linspace(0,chamberLength / m_to_in ,heightStepNumber);

%% CEA RUN (SAVED IN PARENT DIR)
if generate_new_CEA
    CEAOut(p_c, OF, file_name, outputDir);
end

fluidProperties = readmatrix(fullfile(outputDir, "CEA_" + file_name + ".xlsx"));
fluidProperties(1,:) = [];
axialDist = fluidProperties(:,1);

%% Bin averaging
newFluidProperties = zeros(length(heightStepArray), 10);
newFluidProperties(:,1) = heightStepArray;

for y = 1:length(heightStepArray)

    if y < length(heightStepArray)
        idx = axialDist >= heightStepArray(y) & axialDist < heightStepArray(y+1);
    else
        idx = axialDist >= heightStepArray(y);
    end

    if any(idx)
        newFluidProperties(y, 2:10) = mean(fluidProperties(idx, 2:10), 1);
    end
end

throatArea = pi*(throatDiameter/(2*39.37))^2;
diameterArray = flip(2*sqrt((throatArea*newFluidProperties(:,2))/pi));
newFluidProperties = flip(newFluidProperties,1);

converge_index = find(newFluidProperties(:,2) == newFluidProperties(end,2), 1, 'first');
[a, throat_index] = min(newFluidProperties(:,2));

inputValues = [T_start, P_start * psi_to_pa, rho_start, mdot_channel, ...
               T_target, k_w, numChannels, surfaceRoughness, CTE, ...
               youngsModulus, throatDiameter, filletRad, poissonsRatio, ...
               chamberLength, throat_index];

%% Solve wall temperature
[flowTempArray,flowVelocityArray, flowPressureArray,T_wgFinal, ...
 finEfficiency, Qdot, finQdot, T_wl_Array, h_l_Array, h_g_Array, ...
 vonMises,sigma_long, sigma_circ, sigma_rad] = ...
 FINAL_CALCWALLTEMP(converge_index, throat_index, heightArray, ...
 widthArray, wall_thicknessMatrix, diameterArray, ...
 heightStepNumber, newFluidProperties, inputValues);