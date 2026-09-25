%% Stress Calculation for Variable Area Pintle Injector
% Authors: Chanak Gautam
% First Created: 12/14/2025
% Last Updated: 4/26/2026
% Calculations done in Imperial

% pressure differentials
min_psi = 76;
max_psi = 380;
lox_psi = 456;   

% radial dimensions (inches)
outer_radius = 0.843; % outer radius of pintle wall (in)
inner_radius = 0.643; % inner radius of pintle bore (in)
wall_thickness = 0.1; 
sleeve_rad = outer_radius + wall_thickness;

% differential pressure driving stress on components
dP1 = lox_psi - min_psi; % internal minus external pressure @ min chamber
dP2 = lox_psi - max_psi;

%% Pintle Shaft

% hoop stress
pintle_hoop_stress1 = dP1 * (outer_radius) / (wall_thickness);  % psi
pintle_hoop_stress2 = dP2 * (outer_radius) / (wall_thickness);

% end-cap (radial)
endcap_stress1 = (inner_radius * dP1) / (2 * wall_thickness); % psi
endcap_stress2 = (inner_radius * dP2) / (2 * wall_thickness); % psi

% longitudinal stress approximation 
longitudinal_stress = 0.34 * pintle_hoop_stress1 + endcap_stress1; % psi
longitudinal_stress2 = 0.34 * pintle_hoop_stress2 + endcap_stress2; % psi

% von mises equivalent stress 
von_mises_min = sqrt(pintle_hoop_stress1 ^ 2 + longitudinal_stress ^2 - pintle_hoop_stress1 * longitudinal_stress);
von_mises_max = sqrt(pintle_hoop_stress2 ^ 2 + longitudinal_stress2 ^2 - pintle_hoop_stress2 * longitudinal_stress2);

%% Sleeve

sleeve_hoop_stress1 = dP1 * (sleeve_rad) / (wall_thickness);
sleeve_hoop_stress2 = dP2 * (sleeve_rad) / (wall_thickness); % psi

% lol
von_mises_sleeve_min = sqrt(sleeve_hoop_stress1 ^ 2); % psi
von_mises_sleeve_max = sqrt(sleeve_hoop_stress2 ^ 2); % psi

%% Faceplate
%loading case from roarks (11.2, 2e)
% see drive


fprintf("Pintle Shaft:\n");
fprintf("Pintle von mises (min pc): %.2d\n", von_mises_min)
fprintf("FOS: %.2d\n", 20000 / von_mises_min);
fprintf("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
fprintf("Pintle von mises (max pc): %.2d\n", von_mises_max)
fprintf("FOS: %.2d\n", 20000 / von_mises_max);
fprintf("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
fprintf("Note: Yield strength of Copper ≈ 20 ksi\n\n");
fprintf("Sleeve:\n");
fprintf("Sleeve von mises (min pc): %.2d\n", von_mises_sleeve_min);
fprintf("FOS: %.2d\n", 20000 / von_mises_sleeve_min);
fprintf("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
fprintf("Sleeve von mises (max pc): %.2d\n", von_mises_sleeve_max);
fprintf("FOS: %.2d\n", 20000 / von_mises_sleeve_max);
fprintf("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
fprintf("Note: Yield strength of SS ≈ 30 ksi\n");



