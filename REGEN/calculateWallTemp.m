function [flowTemp, flowVel] = calculateWallTemp(heightStepNumber, heatFlux, wallTempArray, hotwallSurfaceAreaArray, flowTempMatrix, heightArray, widthArray)
       
    D_h = params.D_h;                % Hydraulic diameter [m]
    velocity = params.velocity(heightStepNumber); % Coolant velocity [m/s]
    T_bulk = params.T_bulk(heightStepNumber);     % Bulk coolant temperature [K]
    T_wall = params.T_wall(heightStepNumber);     % Wall temperature [K]
    coolant = params.coolant;         % Coolant name
    P_coolant = params.P_coolant(heightStepNumber); % Coolant pressure [Pa]

        % Dynamic viscosity [Pa·s]
    mu_bulk = CoolProp.PropsSI('V', 'T', T_bulk, 'P', P_coolant, coolant);
    
    % Kinematic viscosity [m^2/s]
    nu = mu_bulk / rho;
    
    % Thermal conductivity [W/(m·K)]
    kf = CoolProp.PropsSI('L', 'T', T_bulk, 'P', P_coolant, coolant);
    
    % Specific heat capacity [J/(kg·K)]
    cp = CoolProp.PropsSI('C', 'T', T_bulk, 'P', P_coolant, coolant);
    
    % Prandtl number [-]
    Pr = cp * mu_bulk / kf;
    
    % Get viscosity at wall temperature for Sieder-Tate correction
    mu_wall = CoolProp.PropsSI('V', 'T', T_wall, 'P', P_coolant, coolant);
    
    % Calculate Reynolds number
    Re = (rho * velocity * D_h) / mu_bulk;
    
    % Calculate Nusselt number using Sieder-Tate correlation
    Nu = 0.027 * Re^(4/5) * Pr^(1/3) * (mu_bulk/mu_wall)^0.14;
    
    % Calculate convective heat transfer coefficient [W/(m^2·K)]
    h = Nu * kf / D_h;

    flowTemp = zeros([length(widthArray), length(heightArray)]); %Preallocate array sizes
    flowVel = flowTemp; %Preallocate array sizes
    MASS_FLOW; % Precalcuated mass flow based on # of channels in Malestrom
    wInd = 0;
    hInd = 0;
    htCoeff = % heat transfer coefficient (calculated using bartz equation)
    for width = widthArray
        wInd = wInd + 1;
        for height = heightArray
            hInd = hInd + 1;
            flowVel(wInd, hIn) = 1; % Flow velocity
            flowTemp(wInd, hIn) = 1; %Flow temp
            coolentMass; %%Calculated from volume and denisty
            fluidConactArea;
            convtCoeff;
        end
    end

end