function [QDot, Tsubc, Wmin] = HeatFlux(chanW,deltaP,heightStepLength,thermalConductivity,targetTemp, fluidInfo)
%Note: Function returns [Qdot, -1, -1] if a wall thickness could not be
%found that meets the strength requirements
wallThicknesses = 0.001:0.001:0.02; %Prospective wall thicknesses range in meters
gasTemp = fluidInfo(6); %Hot gas temperature in K
flux = H_g_From_Temperature(targetTemp,fluidInfo) .* (gasTemp - targetTemp); %Calcs heat flux in W/m^2
Tc = -(wallThicknesses .* flux ./ thermalConductivity - targetTemp); %Calcs prospective coolant wall temperatures
bendMaxs = 3 .* deltaP.* chanW .^ 2 ./ (4 .* wallThicknesses.^ 2); %Calcs expected max bending stress in Pa
TMaxs = 3 .* deltaP .* chanW ./ (0.577 .* 4 .* wallThicknesses); %Calcs expected tensile stress in Pa
OEffectives = zeros(length(Tc),1);
i = 1;
while i <= length(Tc)
    calcPainArr = [];
    temperatureDist(:) = linspace(targetTemp,Tc(i),length(wallThicknesses)); %Temperature distribution for integration
    %Corralation between copper's yield/ultimate tensile strength in MPa and
    %temperature in K
    %calcPainArr = 191.31 + 0.65634 .* temperatureDist - 1.85 .* 10.^(-3) .* (temperatureDist).^2 +1.0185 .* 10.^(-6) .* (temperatureDist).^3;
    calcPainArr = 1570.3 - 14.184 .* temperatureDist + 5.6410 .* 10 .^ (-2) .* (temperatureDist) .^ 2 - 1.0592 .* 10 .^ -4 .* (temperatureDist) .^ 3 + 9.2881 .* 10 .^ -8 .* (temperatureDist) .^4 - 3.086 .* 10 .^ -11 .* (temperatureDist) .^ 5;
    OEffectives(i) = sum(calcPainArr .* wallThicknesses(i)) .* 1000000;  %Calcs effective strength for each wall thickness in Pa
    if bendMaxs(i) >= TMaxs(i)
        OMaxs(i) = bendMaxs(i);
    else
       OMaxs(i) = TMaxs(i);
    end
    i = i + 1;
end

Wmin = -1;
Tsubc = -1;
%Searches and finds coolant wall temperature and Wmin that *just* meets the
%strength requirements.
a = 1;
while a <=length(OMaxs) 1
    if OEffectives(a) >= OMaxs(a)
        Wmin = wallThicknesses(a);
        Tsubc = Tc(a);
        break
    end
    a = a + 1;
end

%Outputs
QDot = flux;
Tsubc = Tc(a);

