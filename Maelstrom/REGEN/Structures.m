function [vonMises,sigma_long, sigma_circ, sigma_rad, output] = Structures(T_wl, T_wg, wallThicknesses, channelWidths, chamberPressure, chamberDiameter, deltaP, Qdot, inputValues)
%Note: Function returns [Qdot, -1, -1] if a wall thickness could not be
%found that meets the strength requirements

%% 
CTE = inputValues(9);
youngsModulus = inputValues(10);
poissonsRatio = inputValues(13);
thermalConductivity = inputValues(6);

OEffectives = zeros(length(T_wl),1);
thermalStress = zeros(length(T_wl),1);
i = 1;
while i <= length(T_wl)
    
    temperatureDist(:) = linspace(T_wg(i),T_wl(i),length(wallThicknesses)); %Temperature distribution for integration
    
    %% **CHANGE BELOW** Correlations for Yield strength vs temp for material. Only used in comparison, not used for calculations
    %Correlation between copper's yield/ultimate tensile strength in MPa and temperature in K
    %Copper:
    %calcPainArr = 191.31 + 0.65634 .* temperatureDist - 1.85 .* 10.^(-3).* (temperatureDist).^2 +1.0185 .* 10.^(-6) .* (temperatureDist).^3; %
    % 7075
    calcPainArr = 525.35614 ./ (1 + exp(-(-0.0181864 .* (temperatureDist) + 8.7044)));
    %6061 T6 
    %calcPainArr = ;
    
    %Thermal stress component
    thermalStress(i) = sum(CTE.* (temperatureDist-293) .* youngsModulus)/length(T_wl); % Longitudinal thermal stress given as FL/EA = alpha * L *dT, therefore, F/A = E * alpha *deltaT
    %Yield strength integrated across wall thickness
    OEffectives(i) = sum(calcPainArr .* wallThicknesses(i)) .* 1000000;  %Integrates yield strength array to find effective strength for each wall thickness in Pa
    
    i = i + 1;
end

effectiveRatio = (channelWidths.* inputValues(7))./(pi*(chamberDiameter+wallThicknesses));

hoopStress = chamberPressure.*chamberDiameter./(2.*wallThicknesses);
comp_stress = (chamberPressure*2) + deltaP;
compress_circ = effectiveRatio .*deltaP.*chamberDiameter./wallThicknesses;
thermal_circ = (Qdot.* youngsModulus.* CTE.*wallThicknesses)/(2*(1-poissonsRatio)*thermalConductivity);

%Searches and finds coolant wall temperature and Wmin that *just* meets the
%strength requirements.
a = 0;
while a < length(chamberDiameter)
    a = a + 1;

    sigma_rad(a) = comp_stress(a);
    sigma_circ(a) = thermal_circ(a) + compress_circ(a); %hoopStress(a) + 
    sigma_long(a) = thermalStress(a);
    vonMises(a) = sqrt(((sigma_rad(a)+sigma_circ(a))^2)+((sigma_circ(a)-sigma_long(a))^2) + ((sigma_rad(a) + sigma_long(a))^2));
    
    if OEffectives(a) >= vonMises
        
        display(sigma_rad)
        display(sigma_circ)
        display(sigma_long)
        display(vonMises)
        break
    end

    output = [thermal_circ;compress_circ;calcPainArr;temperatureDist;];
   
end



