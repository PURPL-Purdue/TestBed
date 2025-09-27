function [flowTemp, flowVel, flowPressure] = calculateWallTemp(numChannels, heightStepNumber, heightStepArray, hotwallTempArray, fluidInfoArray, flowTempMatrix, flowVelocityMatrix, flowPressureMatrix, heightArray, widthArray)
    %% Inlet Condition Values
    T_start= 298; % K
    P_start = 3447000; % Pa
    v_start = 0; %m/s
    rho_start = 1000; %kg/m^3, changed coolant density to water at standard temp
    m_flow_total = 0.26716591; %kg/s
    channel_number = numChannels;
    mass_flow = m_flow_total/channel_number; % Precalcuated mass flow based on # of channels in Malestrom
    
    chamberDiameter = 0.0762; % diameter of chamber (m)
    chamberPressure = 3447378.6466; % Chamber Pressure (Pa)
    k_w = 0; % thermal conductivity of the wall (W/m*K)
    gravity = 9.83; %m/s^2
    

    % initialize/re-define matrices & arrays
    flowTemp = flowTempMatrix;
    flowVel = flowVelocityMatrix;
    flowPressure = flowPressureMatrix;
    height_steps = heightStepArray;
    currentHeightStep = height_steps(heightStepNumber);
    T_target = hotwallTempArray ; % target gas-side hotwall temp
    fluidInformation = fluidInfoArray;

    wInd = 0;
    hInd = 0;

    for width = widthArray
        wInd = wInd + 1;

        for height = heightArray
            hInd = hInd + 1;
            
            if (heightStepNumber==1)

                hotWall_dP = P_start - chamberPressure; %calculate dP for structures (Pa)

            else

                hotWall_dP = flowPressure(wInd, hInd, heightStepNumber-1) - chamberPressure; %calculate dP for structures (Pa)

            end
            %% Call Tucker's Function HERE, update variables (coolant side hotwall temp, Heat flux, Wall thickness) (needs updated wall and dP)
            [Q_dot, T_wallL, wall_thicknesses] = HeatFlux(width, hotwall_dP, currentHeightStep, k_w, T_target, fluidInformation);
            
            wall_thickness = wall_thicknesses(heightStepNumber);

            if(flowTemp(wInd,hInd,heightStepNumber-1) == 999999999999999999999999999999999999999999) % if channel dimension combo is already unsuccessful, do not let computation with it continue
                flowTemp(wInd,hInd,heightStepNumber) = 999999999999999999999999999999999999999999;
                break
            end

            %% Reynold's Number
            % Hydraulic Diameter (m)
            hyd_diam = (2*width*height)/(height+width);
            
            % Velocity (m/s)
            if(heightStepNumber == 1)
                velocity = v_start;
            else
                velocity = flowVel(wInd, hInd, heightStepNumber-1);
            end

            % Density (kg/m^3)
            if(heightStepNumber == 1)
                density = rho_start;
            else
                density = 287.67129 * .53365016^(-(1+(1-flowTemp(wInd, hInd, heightStepNumber-1)/574.262)^.628866));
            end 

            % Dynamic viscosity [Pa·s]
            if(heightStepNumber == 1)
                  dyn_visc = 94.544*exp(-.014*T_start);
            else
                dyn_visc = 94.544*exp(-.014*flowTemp(wInd, hInd, heightStepNumber-1));
            end
    
    
            % Calculate Reynolds number
            Re = (density * velocity * hyd_diam) / dyn_visc;
            


            %% Calculate Prandtl Number
           % Thermal conductivity [W/(m·K)]
            kf = .000005*flowTemp(wInd, hInd, heightStepNumber-1) + .105;
            %^Use empirical data/curves found from papers in regen channel
        
            % Specific heat capacity [J/(kg·K)]
            cp = 32.068*exp(.0023*flowTemp(wInd, hInd, heightStepNumber-1));
            %^Use empirical data/curves found from papers in regen channel

            % Prandtl number [-]
            Pr = cp * dyn_visc / kf;



            %% Sieder Tate Nusselt's Number
            % Get viscosity at wall temperature for Sieder-Tate correction
            %mu_wall = 88.748 *exp(-.013*flowTemp(wInd, hInd, heightStepNumber-1));
            %NEEDS CHANGING TO WALL LIQUID SIDE TEMP



            % Calculate Nusselt number using Sieder-Tate correlation
            Nu = 0.027 * Re^(4/5) * Pr^(1/3) * (dyn_visc/mu_wall)^0.14;
            
            % Calculate convective heat transfer coefficient [W/(m^2·K)]
            h_l = Nu * kf / D_h;
            


            %% Calculate Required flow temperature
            % Area of Fin (m^2)
            A_fin = currentHeightStep*2*height;

            % Area of Wall on Coolant side (m^2)
            A_wallL = currentHeightStep*width;

            % Area of Wall on Hotwall side (m^2)
            A_wallG = currentHeightStep *(width+fin_width);

            % Calculate fin width at start of channels (m)
            fin_width = ((pi*(chamberDiameter+2*(wall_thickness))) - (channel_number*width))/channel_number;
            
            % Finning Parameter- measure of convection from fins
            fin_param = sqrt((2*h_l*(fin_width+currentHeightStep))/(k_w*currentHeightStep*fin_width));

            % Fin Efficiency- fin convection efficiency
            fin_efficiency = tanh(fin_param*height)/(fin_param*height);

            % Calculate angle of channel slice
            angle_channel = (width/(pi*(chamberDiameter+2*(wall_thickness))))*360;

            % Calculate mass of coolant in height step, m_coolant (kg)
            m_coolant = density*currentHeightStep*((pi*((chamberDiameter+wall_thickness*2)/2)^2)-(pi*((chamberDiameter/2)^2)))*(angle_channel/360);

            % Calculate required specific heat transfer rate, qdotL_total (J/kg*s)
            qdotL_total = Q_dot(wInd,hInm, heightStepNumber)*((pi*(chamberDiameter)*(angle_channel/360))*currentHeightStep)/m_coolant;

            % Calculate required coolant temp, T_L_req (K)
            T_L_req =  -((qdotL_total*A_wallG)/(h_l*((fin_efficiency*A_fin)+A_wallL)))+T_wallL;



            %% Calculate Coolant Temp increase, delta_T (K)
            delta_T = qdotL_total*currentHeightStep*(width+(2*height*fin_efficiency))/(mass_flow*cp);

            if (heightStepNumber == 1)
                flowTemp(wInd, hInd, heightStepNumber) = T_start + delta_T; % initial flow temp plus deltaT

            else
                flowTemp(wInd, hInd, heightStepNumber) = flowTemp(wInd, hInd, heightStepNumber-1)+delta_T; % Previous step flow temp plus deltaT

            end

            if (flowTemp(wInd, hInd, heightStepNumber) > T_L_req)
                flowTemp(wInd, hInd, heightStepNumber) = 999999999999999999999999999999999999999999; % if coolant temp is too high, nullify

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
            
            flowPressure(wInd, hInd, heightStepNumber) = flowPressure(wInd, hInd, heightStepNumber-1) - delta_P; %Flow pressure
            
            %% Calculate Coolant Velocity Increase via Bernoulli's
            flowVel(wInd, hInd, heightStepNumber) = sqrt((delta_P/density)+(gravity*(currentHeightStep))+velocity^2); % Flow velocity
            
        end
    end

end