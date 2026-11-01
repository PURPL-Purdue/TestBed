% Manifold Area Calculations - Constant Velocity
% Authors: Chanak Gautam, Neil Babu
% First Created: 12/14/2025
% Last Updated: 10/6/2026

% Calculations done in SI, final areas in in^2

% Conversions
in_to_m = 0.0254;
lb_to_kg = 0.4536;
m2_to_in2 = 1 / in_to_m^2;
cm2_to_m2 = 1e-4;

% Flow rate
Mf = 3.0736 * lb_to_kg; % Ethanol mass flow rate [kg/s]

% Density
eth_rho = 789; % kg/m^3

% Number of orifices
eth_N = 6;

% Inlet radii
eth_inlet_rad = 0.305 * in_to_m;

% Minumum manifold area
min_CSA = 0.1 * cm2_to_m2;

% Ethanol manifold
eth_A_in = pi * eth_inlet_rad^2;
eth_v = Mf / (eth_rho * eth_A_in);
eth_mdot_orf = Mf / eth_N;

% Only calculate one half of the ring
eth_k = 0:eth_N/2;
eth_mdot_k = Mf/2 - eth_k * eth_mdot_orf;

% Minimum manifold area
eth_A_k = eth_mdot_k / (eth_rho * eth_v);
eth_A_k = max(eth_A_k, min_CSA) * m2_to_in2;

% Displaying result
disp(eth_A_k)

% Mirror results for the other half
eth_deg = [-fliplr(eth_k(2:end)) eth_k] * 360 / eth_N;
eth_area = [fliplr(eth_A_k(2:end)) eth_A_k];

% Ethanol plot
figure(1)
plot(eth_deg, eth_area, '-o', 'LineWidth', 2)
xlabel('Degrees Around Manifold [deg]')
ylabel('Manifold Cross-Sectional Area [in^2]')
title('Ethanol Manifold Area vs Angle')
grid on