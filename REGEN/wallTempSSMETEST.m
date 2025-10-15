function [flowTemp,flowVelocity,flowPressure, wall_thicknesses] = calculateWallTemp(wall_thicknesses,channelNum, heightStepArray, flowTempMatrix, flowVelocityMatrix, flowPressureMatrix, height, width, heightValue, widthValue, newFluidProperties)
    %% Inlet Condition Values
    T_start= 298; % K
    P_start = 5516000; % Pa
    rho_start = 810; %kg/m^3, changed coolant density to RP-1 at standard temp
    m_flow_total = 0.6109090909; %kg/s --> Calculated this by multiplying the total water mass flow by the ratio of density of RP-1 to water at standard temp
    channel_number = channelNum;
    mass_flow = m_flow_total/channel_number; % Precalcuated mass flow based on # of channels in Malestrom
    v_start = mass_flow/(height * width * rho_start); %m/s
    chamberDiameter = 0.0762; % diameter of chamber (m)
    
    
    chamberPressure = 3447378.6466; % Chamber Pressure (Pa)
    k_w = 130; % thermal conductivity of the wall (W/m*K)
    gravity = 9.83; %m/s^2
    

    % initialize/re-define matrices & arrays
    flowTemp = flowTempMatrix;
    flowVelocity = flowVelocityMatrix;
    flowPressure = flowPressureMatrix;
    height_steps = heightStepArray;
    T_target = 400; % target gas-side hotwall temp

    wInd = widthValue;
    hInd = heightValue;

for heightStepNumber = 1:1:length(height_steps)   
% heightStepNumber = 1;
    currentHeightStep = 0.1796/(length(height_steps)-1);
    if (heightStepNumber==1)

        hotWall_dP = P_start - chamberPressure; %calculate dP for structures (Pa)
        pressure = P_start; %initialize all variables as input variables if on first height step
        temp = T_start;
        velocity = v_start;

    else

        hotWall_dP = flowPressure(wInd, hInd, heightStepNumber-1) - chamberPressure; %calculate dP for structures (Pa)
        pressure = flowPressure(wInd, hInd, heightStepNumber-1); % if on a different height step, initialize as the previous height step's value
        temp = flowTemp(wInd, hInd, heightStepNumber-1);
        velocity = flowVelocity(wInd, hInd, heightStepNumber-1);

    end
 
    if(temp == -1) % if channel dimension combo is already unsuccessful, do not let computation with it continue
        flowTemp(wInd,hInd,heightStepNumber) = -1;
      
        break
        
    end
    
    %% Call Tucker's Function HERE, update variables (coolant side hotwall temp, Heat flux, Wall thickness) (needs updated wall and dP)
    [Q_dot, T_wallL, wall_thickness] = HeatFluxFunction(heightStepArray, width, hotWall_dP, k_w, T_target, newFluidProperties);
    if(T_wallL == -1)
        flowTemp(wInd,hInd,heightStepNumber) = -1;
        break
    end
    %wall_thickness = 0.000254;
    wall_thicknesses(wInd, hInd, heightStepNumber) = wall_thickness;
    
    %% Reynold's Number
    % Hydraulic Diameter (m)
    hyd_diam = (2*width*height)/(height+width);
    

    % % Density (kg/m^3)
    % if(heightStepNumber == 1)
    %     density = rho_start;
    % else
    %     density = 287.67129 * .53365016^(-(1+(1-temp/574.262)^.628866));
    % end 
    density = py.CoolProp.CoolProp.PropsSI('D','P',pressure,'T',temp,'Dodecane');

    % Dynamic viscosity [Pa·s]
    
    %dyn_visc = 94.544*exp(-.014*temp);
    dyn_visc = py.CoolProp.CoolProp.PropsSI('VISCOSITY','P',pressure,'T',temp,'Dodecane');%dynamic viscosity is normal viscosity, for kinematic viscosity, divide by density


    % Calculate Reynolds number
    Re = (density * velocity * hyd_diam) / dyn_visc;
    


    %% Calculate Prandtl Number
    % Thermal conductivity [W/(m·K)]
    %kf = .000005*temp + .105;
    kf = py.CoolProp.CoolProp.PropsSI('CONDUCTIVITY','P',pressure,'T',temp,'Dodecane');
    
    %^Use empirical data/curves found from papers in regen channel

    % Specific heat capacity [J/(kg·K)]
    %cp = 32.068*exp(.0023*temp);
    cp = py.CoolProp.CoolProp.PropsSI('C','P',pressure,'T',temp,'Dodecane');
    %^Use empirical data/curves found from papers in regen channel

    % Prandtl number [-]
    Pr = cp * dyn_visc / kf;



    %% Sieder Tate Nusselt's Number
    % Get viscosity at wall temperature for Sieder-Tate correction
    %mu_wall = 88.748 *exp(-.013*flowTemp(wInd, hInd, heightStepNumber-1));
    %NEEDS CHANGING TO WALL LIQUID SIDE TEMP
    %mu_wall = 94.544*exp(-.014*T_wallL); %relation using liquid side temp
    mu_wall = py.CoolProp.CoolProp.PropsSI('VISCOSITY','P',pressure,'T',T_wallL,'Dodecane');


    % Calculate Nusselt number using Sieder-Tate correlation
    Nu = 0.027 * Re^(4/5) * Pr^(1/3) * (dyn_visc/mu_wall)^0.14;
    
    % Calculate convective heat transfer coefficient [W/(m^2·K)]
    h_l = Nu * kf / hyd_diam;
    


    %% Calculate Required flow temperature
    % Calculate angle of channel slice
    angle_channel = (width/(pi*(chamberDiameter+2*(wall_thickness))))*360;
    
    % Area of Fin (m^2)
    A_fin = currentHeightStep*2*height;

    % Area of Wall on Coolant side (m^2)
    A_wallL = currentHeightStep*width;

    % Calculate fin width at start of channels (m)
    fin_width = ((pi*(chamberDiameter+2*(wall_thickness))) - (channel_number*width))/channel_number;
    
    % Area of Wall on Hotwall side (m^2)
    A_wallG = pi*currentHeightStep *chamberDiameter/channel_number;
    %currentHeightStep *(width+fin_width);

    % Finning Parameter- measure of convection from fins
    fin_param = sqrt((2*h_l*(fin_width+currentHeightStep))/(k_w*currentHeightStep*fin_width));

    % Fin Efficiency- fin convection efficiency
    fin_efficiency = tanh(fin_param*height)/(fin_param*height);

    % Calculate mass of coolant in height step, m_coolant (kg)
    m_coolant = density*currentHeightStep*((pi*((chamberDiameter/2)+wall_thickness+ height)^2)-(pi*(((chamberDiameter/2)+(wall_thickness))^2)))*(angle_channel/360);

    % Calculate required specific heat transfer rate, qdotL_total (J/kg*s)
    % W/m^2 *m^2 / kg
    qdotL_total = Q_dot*((pi*(chamberDiameter)*(angle_channel/360))*currentHeightStep)/m_coolant;

    % Calculate required coolant temp, T_L_req (K) 
    T_L_req =  -((qdotL_total*A_wallG)/(h_l*(((fin_efficiency*A_fin)+A_wallL))))+T_wallL;



    %% Calculate Coolant Temp increase, delta_T (K)
    delta_T = qdotL_total(heightStepNumber)*currentHeightStep*(width+(2*height*fin_efficiency))/(mass_flow*cp);

    flowTemp(wInd, hInd, heightStepNumber) = temp+delta_T; % Previous step flow temp plus deltaT
    
    
    
    
    if (double(flowTemp(wInd, hInd, heightStepNumber)) > T_L_req)
        flowTemp(wInd, hInd, heightStepNumber) = -1; % if coolant temp is too high, nullify

    end


    %% Calculate Coolant Pressure Drop
    if(Re < 2100) % Calculate Friction Factor, cf
        cf = 16/Re;
    elseif(5000 < Re && 200000 > Re)
        cf = 0.046/(Re^0.2);
    else
        cf = 0.0014 + (0.125/(Re^0.32));

    end

    delta_P = 2 * cf * currentHeightStep * density * (velocity^2)/ hyd_diam; %Frictional static pressure drop across the channel
    flowPressure(wInd, hInd, heightStepNumber) = pressure - delta_P; %Flow pressure
    
    %% Calculate Coolant Velocity Increase via Bernoulli's
    flowVelocity(wInd, hInd, heightStepNumber) = sqrt((delta_P/density)+(gravity*(currentHeightStep))+velocity^2); % Flow velocity
    %% 
    
   
end