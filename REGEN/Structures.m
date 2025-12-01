function [vonMises,sigma_long, sigma_circ, sigma_rad] = Structures(T_wl, T_wg, wallThicknesses, channelWidths, chamberPressure, chamberDiameter, deltaP, inputValues)
%Note: Function returns [Qdot, -1, -1] if a wall thickness could not be
%found that meets the strength requirements

%% 
CTE = inputValues(9);
youngsModulus = inputValues(10);

bendMaxs = 3 .* deltaP.* channelWidths .^ 2 ./ (4 .* wallThicknesses.^ 2); %Calcs expected max bending stress in Pa
%TMaxs = 3 .* deltaP .* channelWidth ./ (0.577 .* 4 .* wallThicknesses); %Calcs expected tensile stress in Pa
OEffectives = zeros(length(T_wl),1);
thermalStress = zeros(length(T_wl),1);
hoopStress = zeros(length(T_wl),1);
i = 1;
while i <= length(T_wl)
    temperatureDist(:) = linspace(T_wg(i),T_wl(i),length(wallThicknesses)); %Temperature distribution for integration
    %Correlation between copper's yield/ultimate tensile strength in MPa and temperature in K
    %Copper:
    %calcPainArr = 191.31 + 0.65634 .* temperatureDist - 1.85 .* 10.^(-3).* (temperatureDist).^2 +1.0185 .* 10.^(-6) .* (temperatureDist).^3; %
    % 7075
    calcPainArr = 525.35614 ./ (1 + exp(-(-0.0181864 .* (temperatureDist) + 8.7044)));

    %6061 T6 
    %calcPainArr = ;

    thermalStress(i) = sum(CTE.* (temperatureDist-293) .* youngsModulus)/length(T_wl);
    hoopStress(i) = chamberPressure(i)*chamberDiameter(i)/(2*wallThicknesses(i));
    
    OEffectives(i) = sum(calcPainArr .* wallThicknesses(i)) .* 1000000;  %Integrates yield strength array to find effective strength for each wall thickness in Pa
    OMaxs(i) = sqrt((bendMaxs(i)));

    
    i = i + 1;
end

comp_stress = (chamberPressure*2) + deltaP;



%Searches and finds coolant wall temperature and Wmin that *just* meets the
%strength requirements.
a = 0;
while a < length(OMaxs)
    a = a + 1;

    sigma_rad = comp_stress + OMaxs(a);
    sigma_circ = thermalStress(a) + hoopStress(a);
    sigma_long = thermalStress(a);
    vonMises = sqrt(((sigma_rad+sigma_circ)^2)+((sigma_circ-sigma_long)^2) + ((sigma_rad + sigma_long)^2));
    if OEffectives(a) >= vonMises
        
        display(sigma_rad)
        display(sigma_circ)
        display(sigma_long)
        display(vonMises)
        break
    end
end



