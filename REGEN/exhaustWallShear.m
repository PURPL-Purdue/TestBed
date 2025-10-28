%% Local helper: compute wall shear from Darcy-Weisbach
function tau_exhaust  = exhaustWallShear(chamberDiameterArray, heightStepNumber, newFluidProperties)

    eps_rough = 8e-7;

    exhaustDensity = newFluidProperties(heightStepNumber,11);
    exhaustMach = newFluidProperties(heightStepNumber,4);
    exhaustSoS = newFluidProperties(heightStepNumber,12);
    exhaustArea = newFluidProperties(heightStepNumber,2);
    exhaustViscosity = newFluidProperties(heightStepNumber,7);

    chamberDiameter = chamberDiameterArray(heightStepNumber);

    exhaustVelocity = exhaustMach * exhaustSoS;

    LStar = 45; % THIS IS NOT THE REAL NUMBER

    exhaustReynolds = (exhaustDensity * exhaustVelocity * LStar) / exhaustViscosity;

    exhaustFriction = 1.325/(log(((eps_rough)/(3.7 *chamberDiameter))+(5.74)/((exhaustReynolds)^0.9)));

    tau_exhaust = (exhaustFriction/4)*(1/2)*(exhaustDensity)*(exhaustVelocity^2);

    display(tau_exhaust)

end