%% Local helper: compute wall shear from Darcy-Weisbach
function tau_exhaust  = exhaustWallShear(chamberDiameterArray, heightStepNumber, newFluidProperties)

    eps_rough = 8e-7;

    exhaustDensity = newFluidProperties(:,11);
    exhaustMach = newFluidProperties(:,4);
    exhaustSoS = newFluidProperties(:,12);
    exhaustArea = newFluidProperties(:,2);
    exhaustViscosity = newFluidProperties(:,7);

    chamberDiameter = chamberDiameterArray(heightStepNumber);

    exhaustVelocity = exhaustMach .* exhaustSoS;

    LStar = 45; % THIS IS NOT THE REAL NUMBER

    exhaustReynolds = (exhaustDensity .* exhaustVelocity .* LStar) ./ exhaustViscosity;

    exhaustFriction = 1.325/(log(((eps_rough)/(3.7 .*chamberDiameter))+(5.74)./((exhaustReynolds).^0.9)));

    tau_exhaust = (exhaustFriction/4).*(1/2).*(exhaustDensity).*(exhaustVelocity.^2);

end