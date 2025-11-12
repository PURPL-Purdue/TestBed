function [sigma_vm] = vonMisesStress(wInd, hInd, heightStepNumber, chamberDiameter, newFluidProperties, flowPressure, flowVelocity, frictionFactor, density)

    exhaustMach = newFluidProperties(heightStepNumber,4);
    exhaustViscosity = newFluidProperties(heightStepNumber,7);
    exhaustDensity = newFluidProperties(heightStepNumber,11);
    exhaustSoS = newFluidProperties(heightStepNumber,12);

% Longitudinal stress calculation
    exhaustPressure = (newFluidProperties(heightStepNumber,9))/100000; % Pa
    
    dP = flowPressure(wInd, hInd, heightStepNumber) - exhaustPressure;

    wallThick = 0.0015556992; % m

      % sigma_long = longitudinalStressFromPressure(p, r, t)
    % Thin-wall longitudinal stress: sigma_long = p * r / (2*t)
    % Inputs:
    %   p - internal pressure (Pa)
    %   r - inner radius (m)
    %   t - wall thickness (m)

    sigma_long = (-dP * (chamberDiameter/2)) / (2 * wallThick); % longitudinal stress value

% Coolant side tangential stress calculation
    % coolant side wall shear stress
    velocity = flowVelocity(wInd, hInd, heightStepNumber);

    tau_coolant = (frictionFactor/4)*(1/2)*(density)*(velocity^2); % coolant side tangential stress value

% Exhaust side tangential stress calculation
    eps_rough = 8e-7;

    exhaustVelocity = exhaustMach * exhaustSoS;

    LStar = 45; % THIS IS NOT THE REAL NUMBER!!!!

    exhaustReynolds = (exhaustDensity * exhaustVelocity * LStar) / exhaustViscosity;

    exhaustFriction = 1.325/(log(((eps_rough)/(3.7 *chamberDiameter))+(5.74)/((exhaustReynolds)^0.9)));

    tau_exhaust = (exhaustFriction/4)*(1/2)*(exhaustDensity)*(exhaustVelocity^2); % exhaust side tangential stress value

% von Mises stress calculation
    sigma_vm = sqrt(((tau_exhaust + tau_coolant)^2 + (sigma_long)^2)/2);
    
end
