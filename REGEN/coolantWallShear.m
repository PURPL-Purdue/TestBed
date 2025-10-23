%% Local helper: compute wall shear from Darcy-Weisbach
function tau_coolant = coolantWallShear(density, flowVelocity, frictionFactor)
    % coolant side wall shear stress
    tau_coolant = (frictionFactor./4).*(1/2).*(density).*(flowVelocity.^2);

end