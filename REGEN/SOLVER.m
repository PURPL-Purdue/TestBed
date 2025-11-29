function [T_wgFinal, flowTemp, flowPressure, flowVelocity, finEfficiency, Qdot, finQdot, T_wl_Array, h_l_Array, h_g_Array] = SOLVER(throat_index, converge_index, chamberDiameterArray, heightStepNumber, newFluidInfo, channelWidth, channelHeight, wallThickness)
  % clc;
  % clear;
%% Inputs (meters, Kelvin, seconds, Pascals, Watts) 
%dimensions
% channelWidths = channelWidth; %Should be fully filled matrices (1 value per heightstep of ACTUAL VALUES)
% channelHeights = channelHeight;
% wallThicknesses = wallThickness;
% converge_index is start of engine contraction (index from nozzle),
% throat_index is point of peak contraction, channelWidth, channelHeight,
% wallThickness are 3 value arrays with [start value, throat value, end
% value]
channelWidths = [linspace(channelWidth(1),channelWidth(2),throat_index),linspace(channelWidth(2),channelWidth(3),converge_index-throat_index),(channelWidth(3)*ones(1,(heightStepNumber - converge_index)))];
channelHeights = [linspace(channelHeight(1),channelHeight(2),throat_index),linspace(channelHeight(2),channelHeight(3),converge_index-throat_index),(channelHeight(3)*ones(1,(heightStepNumber - converge_index)))];
wallThicknesses = [linspace(wallThickness(1),wallThickness(2),throat_index),linspace(wallThickness(2),wallThickness(3),converge_index-throat_index),(wallThickness(3)*ones(1,(heightStepNumber - converge_index)))];


%combustion properties
exhaustTemps = newFluidInfo (:,6);
chamberPressure = newFluidInfo(:,9)*100000; % Chamber Pressure (Pa)

%material properties
thermalConductivity = 161; % 7075 W/m*K

hotwallTemps = [];

% calculateWallTemp Inputs
T_start= 298; % K
P_start = 2551000; % Pa
rho_start = 1000; %kg/m^3, changed coolant density to RP-1 at standard temp

m_flow_total = 2.26796; %kg/s --> Calculated this by multiplying the total water mass flow by the ratio of density of RP-1 to water at standard temp
channel_number = 62;
mass_flow = m_flow_total/channel_number; % Precalcuated mass flow based on # of channels in Malestrom
v_start = mass_flow/(channelWidths(1) * channelHeights(1) * rho_start); %m/s



% initialize/re-define matrices & arrays
flowTemp = [T_start];
flowVelocity = [v_start];
flowPressure = [P_start];
updatedPressure = []; % if on a different height step, initialize as the previous height step's value
updatedTemps = [];
T_wgFinal = zeros(1,heightStepNumber);
xz = 1;
%% SOLVER


 
heightStep = 1; % current height step

while (heightStep < heightStepNumber)
    xz = 0;
    height = channelHeights(heightStep);
    width = channelWidths(heightStep);


    heightStepArray = linspace(0,8.97/39.37,heightStepNumber);
    currentHeightStep = heightStepArray(2)-heightStepArray(1);
    
    
    chamberDiameter = chamberDiameterArray(heightStep);
    if (heightStep==1)
        
        pressure = P_start; %initialize all variables as input variables if on first height step
        temp = T_start;
        velocity = v_start;

    else

        pressure = flowPressure(heightStep-1); % if on a different height step, initialize as the previous height step's value
        temp = flowTemp(heightStep-1);
        velocity = mass_flow/(width*height*density);

    end

    T_wg_InitialGuess = 700; % Set initial guess to be upper threshold for wall temp
    T_wg = T_wg_InitialGuess;
    
    


    notCorrect = true;

    while(notCorrect)
        %hot side
        h_g = H_g_From_Temperature(T_wg, newFluidInfo);
        T_r = exhaustTemps(heightStep);
        Q_dotIN = h_g(heightStep)*(T_r-T_wg);
        
        %cool side
        T_wl = -(Q_dotIN*wallThicknesses(heightStep)/thermalConductivity)+T_wg;
        

            %% H_L finder

            hyd_diam = (2*width*height)/(height+width);
            density = double(py.CoolProp.CoolProp.PropsSI('D','P',pressure,'T',temp,'Water'));
        
            % Dynamic viscosity [Pa·s]
            dyn_visc = double(py.CoolProp.CoolProp.PropsSI('VISCOSITY','P',pressure,'T',temp,'Water'));%dynamic viscosity is normal viscosity, for kinematic viscosity, divide by density
        
        
            % Calculate Reynolds number
            Re = (density * velocity * hyd_diam) / dyn_visc;
           
        
            % Calculate Prandtl Number
            % Thermal conductivity [W/(m·K)]
            kf = double(py.CoolProp.CoolProp.PropsSI('CONDUCTIVITY','P',pressure,'T',temp,'Water'));
            
            %^Use empirical data/curves found from papers in regen channel
        
            % Specific heat capacity [J/(kg·K)]
            cp = double(py.CoolProp.CoolProp.PropsSI('CPMASS','P',pressure,'T',temp,'Water'));
            
            % Prandtl number [-]
            Pr = cp * dyn_visc / kf;
           
          
        
            %% Sieder Tate Nusselt's Number
            % Get viscosity at wall temperature for Sieder-Tate correction
            mu_wall = double(py.CoolProp.CoolProp.PropsSI('VISCOSITY','P',pressure,'T',T_wl,'Water'));
            
            % Calculate Nusselt number using Sieder-Tate correlation
            Nu = 0.027 * Re^(4/5) * Pr^(1/3) * (dyn_visc/mu_wall)^0.14;
            
            surfaceRoughness = 0.000006; % measure for how rough channels are after slit saw (chosen from engineering toolbox)
            frictionFactor = (1/(-1.8*log10(((surfaceRoughness/hyd_diam)/3.7)^(1.11)+(6.9/Re))))^2; %new friction factor which is so cool (from vincent and haaland)
            
            % Calculate convective heat transfer coefficient [W/(m^2·K)]
            h_l = Nu * kf / hyd_diam;
           
            
            % Calculate angle of channel slice
            angle_channel = (width/(pi*(chamberDiameter+2*(wallThicknesses(heightStep)))))*360;
            
            % Area of Fin (m^2)
            A_fin = currentHeightStep*2*height;
        
            % Calculate fin width at start of channels (m)
            fin_width = ((pi*(chamberDiameter+2*(wallThicknesses(heightStep)))) - (channel_number*width))/channel_number;
            
            % Area of Wall on Hotwall side (m^2)
            A_wallG = pi*currentHeightStep *chamberDiameter*angle_channel/360;
        
            % Finning Parameter- measure of convection from fins
            fin_param = sqrt((2*h_l*(fin_width+currentHeightStep))/(thermalConductivity*currentHeightStep*fin_width));
        
            % Fin Efficiency- fin convection efficiency
            fin_efficiency = tanh(fin_param*height)/(fin_param*height);
         

            % Area of Wall on Coolant side (m^2)
            A_wallLiquid = currentHeightStep*width;
            A_wallL = A_wallLiquid+(A_fin*fin_efficiency);
        xz = xz+1;
        %Calculate required heat flux, qdotL_total
        Q_dotOUT = h_l*A_wallL*(T_wl-temp)/(A_wallG);

            if(abs(Q_dotIN-Q_dotOUT) < 0.01*(Q_dotIN))
            T_wgFinal(heightStep) = T_wg;
            notCorrect = false;
            elseif(Q_dotIN<Q_dotOUT)
                T_wg = T_wg * 0.95;
            else
                T_wg = T_wg * 1.05;
            end
        
    end
    display(xz)
        
    % Calculate Coolant Temp increase, delta_T (K)
    delta_T = Q_dotOUT*A_wallL/(mass_flow*cp);
    flowTemp(heightStep) = temp+delta_T; % Previous step flow temp plus deltaT
   
    delta_P = frictionFactor * currentHeightStep * density * (velocity^2)/ (2*hyd_diam); %Frictional static pressure drop across the channel
    flowPressure(heightStep) = pressure - delta_P; %Flow pressure
    
    % Calculate Coolant Velocity Increase via Bernoulli's
    flowVelocity(heightStep) = mass_flow/(width*height*density); % Flow velocity    %% 
       
    % Add any extraneous matrix data here
    finEfficiency(heightStep) = fin_efficiency;
    Qdot(heightStep) = Q_dotOUT;
    finQdot(heightStep) = ((A_fin*fin_efficiency)/A_wallL)*Q_dotOUT;
    T_wl_Array(heightStep) = T_wl;
    h_l_Array(heightStep) = h_l;
    h_g_Array(heightStep) = h_g(heightStep);


    


    heightStep = heightStep+1;
end
    %% Rest of Engine Solver
end