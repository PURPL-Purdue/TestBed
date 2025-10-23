function sigma_long = longitudinalStressFromPressure(flowPressure, chamberDiameterArray)
    
    vm_Data = readmatrix("CEAOutFz_10-22-25.xlsx");

    
    exhaustPressure = (vm_Data(:,9))/100000; % Pa
    
    dP = flowPressure - exhaustPressure;

    chamberDiameter = chamberDiameterArray(heightStepNumber); %?? unit?

    wallThick = 0.0015556992; % m

      % sigma_long = longitudinalStressFromPressure(p, r, t)
    % Thin-wall longitudinal stress: sigma_long = p * r / (2*t)
    % Inputs:
    %   p - internal pressure (Pa)
    %   r - inner radius (m)
    %   t - wall thickness (m)

    sigma_long = (-dP * (chamberDiameter/2)) / (2 * wallThick);
end