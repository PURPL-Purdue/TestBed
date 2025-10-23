%% Local helper: compute wall shear from Darcy-Weisbach
function tau_exhaust  = exhaustWallShear(chamberDiameterArray, heightStepNumber)

    eps_rough = 8e-7;

    vm_Data = readmatrix("CEAOutFz_10-22-25.xlsx");

    exhaustDensity = vm_Data(:,11);
    exhaustMach = vm_Data(:,4);
    exhaustSoS = vm_Data(:,12);
    exhaustArea = vm_Data(:,2);
    exhaustViscosity = vm_Data(:,7);

    chamberDiameter = chamberDiameterArray(heightStepNumber);

    exhaustVelocity = exhaustMach .* exhaustSoS;

    LStar = 73.4776667586; % THIS IS NOT THE REAL NUMBER

    exhaustReynolds = (exhaustDensity .* exhaustVelocity .* LStar) ./ exhaustViscosity;

    exhaustFriction = 1.325/(log(((eps_rough)/(3.7 .*chamberDiameter))+(5.74)./((exhaustReynolds).^0.9)));

    tau_exhaust = (exhaustFriction/4).*(1/2).*(exhaustDensity).*(exhaustVelocity.^2);

end