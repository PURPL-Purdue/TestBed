function [flowTemp, flowVel] = calculateWallTemp(numChannels, heightStepNumber, heightStepArray, heatTransferArray, wallTempArray, hotwallSurfaceAreaArray, flowTempMatrix, flowVelocityMatrix, flowPressureMatrix, heightArray, widthArray)
    
    %% Inlet Condition Values
    T_start= 298; % K
    P_start = 3447000; % Pa
    v_start = 0; %m/s
    rho_start = 0; %kg/m^3
    m_flow_total = 0.26716591; %kg/s
    channel_number = numChannels;
    mass_flow = m_flow_total/channel_number; % Precalcuated mass flow based on # of channels in Malestrom
    

    wall_thickness = 0; % thickness of hotwall (m)
    chamberDiameter = 0.0762; % diameter of chamber (m)
    k_w = 0; % thermal conductivity of the wall (W/m*K)


    % initialize/re-define matrices & arrays
    flowTemp = flowTempMatrix;
    flowVel = flowVelocityMatrix;
    flowPressure = flowPressureMatrix;
    height_steps = heightStepArray;
    Q_dot = heatTransferArray;

    wInd = 0;
    hInd = 0;
    htCoeff = % heat transfer coefficient (calculated using bartz equation)
    for width = widthArray
        wInd = wInd + 1;
        for height = heightArray
            hInd = hInd + 1;
            
            %% Reynold's Number
            % Hydraulic Diameter (m)
            hyd_diam = (2*width*height)/(height+width);
            
            % Velocity (m/s)
            if(heightStepNumber == 1)
                velocity = v_start;
            else
                velocity = flowVel(wInd, hIn, heightStepNumber-1);
            end

            % Density (kg/m^3)
            if(heightStepNumber == 1)
                density = rho_start
            else
                density = CoolProp.PropsSI('D', 'T', flowTemp(wInd, hIn, heightStepNumber-1), 'P', flowPressure(wInd, hIn, heightStepNumber-1), coolant);
            end

            % Dynamic viscosity [Pa·s]
            if(heightStepNumber == 1)
                dyn_visc = CoolProp.PropsSI('V', 'T', T_start, 'P', P_start, coolant);
            else
                dyn_visc = CoolProp.PropsSI('V', 'T', flowTemp(wInd, hIn, heightStepNumber-1), 'P', flowPressure(wInd, hIn, heightStepNumber-1), coolant);
            end

            % Calculate Reynolds number
            Re = (density * velocity * hyd_diam) / dyn_visc;
            


            %% Calculate Prandtl Number
            % Thermal conductivity [W/(m·K)]
            kf = CoolProp.PropsSI('L', 'T', T_bulk, 'P', P_coolant, coolant);
            %^Use empirical data/curves found from papers in regen channel
            
            % Specific heat capacity [J/(kg·K)]
            cp = CoolProp.PropsSI('C', 'T', T_bulk, 'P', P_coolant, coolant);
            %^Use empirical data/curves found from papers in regen channel

            % Prandtl number [-]
            Pr = cp * dyn_visc / kf;



            %% Sieder Tate Nusselt's Number
            % Get viscosity at wall temperature for Sieder-Tate correction
            mu_wall = CoolProp.PropsSI('V', 'T', T_wall, 'P', P_coolant, coolant);

            % Calculate Nusselt number using Sieder-Tate correlation
            Nu = 0.027 * Re^(4/5) * Pr^(1/3) * (dyn_visc/mu_wall)^0.14;
            
            % Calculate convective heat transfer coefficient [W/(m^2·K)]
            h_l = Nu * kf / D_h;
            


            %% Calculate Required flow temperature
            % Area of Fin (m^2)
            A_fin = height_steps(heightStepNumber)*2*height;

            % Area of Wall on Coolant side (m^2)
            A_wallL = height_steps(heightStepNumber)*width;

            % Area of Wall on Hotwall side (m^2)
            A_wallG = height_steps(heightStepNumber) *(width+fin_width);

            % Calculate fin width at start of channels (m)
            fin_width = ((pi*(chamberDiameter+2*(wall_thickness))) - (channel_number*width))/channel_number;
            
            % Finning Parameter- measure of convection from fins
            fin_param = sqrt((2*h_l*(fin_width+height_steps(heightStepNumber)))/(k_w*height_steps(heightStepNumber)*fin_width));

            % Fin Efficiency- fin convection efficiency
            fin_efficiency = tanh(fin_param*height)/(fin_param*height);

            % Calculate angle of channel slice
            angle_channel = (width/(pi*(chamberDiameter+2*(wall_thickness))))*360;

            % Calculate mass of coolant in height step, m_coolant (kg)
            m_coolant = density*height_steps(heightStepNumber)*((pi*((chamberDiameter+wall_thickness*2)/2)^2)-(pi*((chamberDiameter/2)^2)))*(angle_channel/360);

            % Calculate required heat transfer rate, q_l_total
            

            % Calculate required coolant temp, T_L_req (K)



            flowVel(wInd, hIn, heightStepNumber) = 1; % Flow velocity
            flowTemp(wInd, hIn, heightStepNumber) = 1; %Flow temp
            flowPressure(wInd, hIn, heightStepNumber) = 1; %Flow pressure
            coolentMass; %%Calculated from volume and denisty
            fluidConactArea;
            convtCoeff;
        end
    end

end