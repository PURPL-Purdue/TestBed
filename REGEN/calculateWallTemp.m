function [flowTemp, flowVel] = calculateWallTemp(heightStepNumber, heatFluxArray, wallTempArray, hotwallSurfaceAreaArray, flowTempMatrix, heightArray, widthArray)
    
    %% Inlet Condition Values
    T_start= 298; % K
    P_start = 3447000; % Pa
    v_start = 0;
    rho_start = 0;




    flowTemp = zeros([length(widthArray), length(heightArray)]); %Preallocate array sizes
    flowVel = flowTemp; %Preallocate array sizes
    flowPressure = flowTemp;


    MASS_FLOW; % Precalcuated mass flow based on # of channels in Malestrom
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
            Re = (rho * velocity * hyd_diam) / dyn_visc;
            


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
            


            %%

            flowVel(wInd, hIn, heightStepNumber) = 1; % Flow velocity
            flowTemp(wInd, hIn, heightStepNumber) = 1; %Flow temp
            flowPressure(wInd, hIn, heightStepNumber) = 1; %Flow pressure
            coolentMass; %%Calculated from volume and denisty
            fluidConactArea;
            convtCoeff;
        end
    end

end