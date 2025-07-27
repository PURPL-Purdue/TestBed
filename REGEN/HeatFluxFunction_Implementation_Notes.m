%Read in fluid properties upon code startup
%fluidProperties = readmatrix("CEAOut.xlsx");
%fluidProperties = readmatrix("CEAOutFrozen.xlsx");
fluidProperties = readmatrix("CEAOutEquilibriumFAC.xlsx");
%fluidProperties = readmatrix("CEAOutCombo.xlsx")
fluidProperties = fluidProperties(3:end,:);
fluidProperties = double(fluidProperties);

%Read in desired temperature distribution
targetTemps = readmatrix("HotWall Temperature Distribution.xlsx"); %Imports target temperature distribution
targetTemps = targetTemps(:,2);

%Define values in code iteration. Make sure targetTemp, channelW, fluidInfo,
%and deltaP are updated every height iteration
targetTemp = targetTemps(1);
channelW = 0.008; %Width of the coolant channel in meters
thermalCondutivity = 401; %Material thermal conductivity in W/mk
fluidInfo = fluidProperties(500,:); %Select the row corresponding to the axial position being analysed 
                                  %and read in the corresponding fluid properties
heightStepLength = 0.001; %length of the height step being analysed in meters
deltaP = 10000000; %The absolute value of the delta P between coolant channel and the chamber in Pa. 
%               Note: the active fluidInfo row should have chamber pressure
%               in bar


[Qdot,Tc,Wmin] = HeatFlux(channelW, deltaP, heightStepLength, thermalConductivity, targetTemp, fluidInfo)