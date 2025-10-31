function [QDot, Tsubc, Wmin] = HeatFluxFunction(heightStepNumber,heightStepArray,chanW,deltaP,thermalConductivity,targetTemp, fluidInfo)
%Note: Function returns [Qdot, -1, -1] if a wall thickness could not be
%found that meets the strength requirements

wallThicknesses = linspace(0.000762,0.00254,length(heightStepArray)); %Prospective wall thicknesses range in meters
gasTemp = fluidInfo(:,6); %Hot gas temperature in K
flux = (H_g_From_Temperature(targetTemp,fluidInfo) .* (gasTemp - targetTemp))'; %Calcs heat flux in W/m^2
Tc = -(wallThicknesses .* flux(heightStepNumber) ./ thermalConductivity - targetTemp); %Calcs prospective coolant wall temperatures ,,flux:
bendMaxs = 3 .* deltaP.* chanW .^ 2 ./ (4 .* wallThicknesses.^ 2); %Calcs expected max bending stress in Pa
TMaxs = 3 .* deltaP .* chanW ./ (0.577 .* 4 .* wallThicknesses); %Calcs expected tensile stress in Pa
OEffectives = zeros(length(Tc),1);
i = 1;
while i <= length(Tc)
    calcPainArr = [];
    temperatureDist(:) = linspace(targetTemp,Tc(i),length(wallThicknesses)); %Temperature distribution for integration
    %Corralation between copper's yield/ultimate tensile strength in MPa and
    %temperature in K
    %Copper
    %calcPainArr = 191.31 + 0.65634 .* temperatureDist - 1.85 .* 10.^(-3).* (temperatureDist).^2 +1.0185 .* 10.^(-6) .* (temperatureDist).^3; %
    
    % 7075
    calcPainArr = 525.35614 ./ (1 + exp(-(-0.0181864 .* (temperatureDist) + 8.7044)));

    OEffectives(i) = sum(calcPainArr .* wallThicknesses(i)) .* 1000000;  %Calcs effective strength for each wall thickness in Pa
    if bendMaxs(i) >= TMaxs(i)
        OMaxs(i) = bendMaxs(i);
    
    else
       OMaxs(i) = TMaxs(i);
      
    end
    i = i + 1;
end

Wmin = -2;
Tsubc = -2;
%Searches and finds coolant wall temperature and Wmin that *just* meets the
%strength requirements.
a = 0;
while a < length(OMaxs)
    a = a + 1;
    if OEffectives(a) >= OMaxs(a)
        Wmin = wallThicknesses(a);
        Tsubc = Tc(a);
        
        break
    end
end

%Outputs
QDot = flux;

