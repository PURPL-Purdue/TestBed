function [sigma_vm] = vonMisesStress(tau_exhaust, tau_coolant, sigma_long)
    sigma_vm = sqrt(((tau_exhaust + tau_coolant)^2 + (sigma_long)^2)/2);
    display(sigma_vm)

end
