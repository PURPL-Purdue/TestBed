%% PURPL
%% Coaxial Swirl Injector Sizing

%% Mass Flows & Elements
m_dot_LOX_total = 0.565272;
m_dot_f_total = 0.353295;
elements = 4;

%% Inner Element Sizing____________________________________________________

%% CONSTANTS
g = 32.3;                                                                % gravitational acceleration (ft/s^2)
LOX_rho = 62.4;                                                          % water density (lb/ft^3)
visc_LOX = 1.0764 * 10^-5;                                               % kinematic viscosity of water (ft^2/s)

%% INPUTS
m_dot_LOX = m_dot_LOX_total / elements;                                  % water mass flow rate (lbm/s)
inner_inlets = 3;                                                        % no. of tangential inner_inlets 
spray_angle = 30;                                                        % desired spray angle (degrees)
delta_p_LOX = 17280;                                                     % pressure drop over injector (lbf/ft^2)

%% ASSUMPTIONS
K_initial = 3;                                                           % estimate geometric characteristic parameter
wall_thickness = .005104;                                                % (ft)
R_nozzle = 2;                                                            % coefficient of nozzle opening

%% DIMENSION CALCS
condition_inner = true;
iterations_inner = 1;
while condition_inner
    % fill efficiency
    E = fzero(@(E) K_initial - (((1-E) * sqrt(2)) / (E * sqrt(E))), .4);
    
    % dischsarge coefficient
    mu = E * sqrt(E / (2 - E)); 

    % outlet diameter
    d_0 = sqrt(4 * m_dot_LOX / (pi * mu * sqrt(2 * LOX_rho * delta_p_LOX * g)));
    
    % outlet radius
    outlet_rad = d_0 / 2;
    
    % swirl arm distance
    R_inner = (R_nozzle * d_0) / 2;

    % inlet diameter
    inlet_diameter_inner = sqrt((2 * R_inner * d_0) / (inner_inlets * K_initial));

    % inlet radius
    inlet_rad = inlet_diameter_inner / 2;
    
    % LOX Reynolds number
    Re = (4 * m_dot_LOX)/ (pi * LOX_rho * visc_LOX * sqrt(inner_inlets) * inlet_diameter_inner);

    % friction coefficient
    lambda = exp((25.8/log(Re)^2.58) - 2); 

    % recalculating K accounting for frictional losses
    K_lambda = (R_inner * (d_0 / 2)) / ((inner_inlets * (inlet_rad)^2) + (lambda / 2) * R_inner * (R_inner - (d_0 / 2))); 

    % recalculating filling efficiency
    E_eq = fzero(@(E_eq) K_lambda - (((1-E_eq) * sqrt(2)) / (E_eq * sqrt(E_eq))), .4);

    % recalcaulting discharge coefficient
    mu_eq = E_eq * sqrt(E_eq / (2 - E_eq)); 

    % recalculating dishcarge coefficient taking into account angular momentum losses
    mu_i = mu_eq / (sqrt(1 + mu_eq^2 * (K_initial^2 / R_nozzle^2)));

    % recalculating outlet diameter
    d_0_i = sqrt(4 * m_dot_LOX / (pi * mu_i * sqrt(2 * LOX_rho * delta_p_LOX * g)));

    % reclculating swirl arm distance
    R_inner = (R_nozzle * d_0_i) / 2;

    % recalculating K 
    K_i = (2 * R_inner * d_0_i) / (inner_inlets * inlet_diameter_inner^2); 

    % calculating ratio of gas core diameter to swirl outer diameter
    f = @(S) ((sqrt(1 - mu^2 * K_initial^2)) ...
    - (S * sqrt(S^2 - mu^2 * K_initial^2)) ...
    - (mu^2 * K_initial^2 * log((1 ...
    + sqrt(1 - mu^2 * K_initial^2)) ...
    / (S + sqrt(S^2 - mu^2 ...
    * K_initial^2)))) - mu); 

    S = fzero(f, .8);                                                    % fzero method borrowed from Jacob Bell's code

    % actual spray angle
    spray_angle_new = atand((2 * mu * K_i) / sqrt((1 + S)^2 -(4 * mu_i^2 * K_initial^2)));

    % checking if K is within a certain range
    if abs((K_initial - K_i) / K_i < .03)
        K_final = K_i;
        condition_inner = false;
    else
        iterations_inner = iterations_inner + 1;
        K_initial = K_i;
    end
end

chamber_diameter_inner = (2 * R_inner) + inlet_diameter_inner;
chamber_length_inner = 3 * R_inner;
nozzle_diameter_inner = d_0_i;
nozzle_length_inner = 2 * d_0;
inlet_length_inner = inlet_rad * 3;

%% Outer Element Sizing____________________________________________________

%% CONSTANTS
g = 32.3;                                                                  % gravitational constant (ft/s^2)
f_rho = 62.4;                                                              % water density (lb/ft^3)
visc_f = 1.0764 * 10^-5;                                                   % kinematic viscosity of water (ft^2/s)

%% INPUTS
m_dot_f = m_dot_f_total / elements;                                        % water mass flow rate (lbm/s)
delta_p_f = 17280;                                                         % fuel pressure drop over injector (lbf/ft^2) 
outer_inlets = 4;                                                          % no.of tangential inlets in outer element

%% ASSUMPTIONS
R_nozzle_outer = 3.5;                                                      % coefficient of outer nozzle opening 
film_thickness = .00098425;                                                

%% DIMENSIONS CALCS 
permitted_vortex_rad = (nozzle_diameter_inner / 2) + .00098425;                  
d_outer = 2 * permitted_vortex_rad + 2 * film_thickness;                   % first estimation of outlet diameter
a_eff = m_dot_f / sqrt(2 * f_rho * delta_p_f * g);                         % effective flow area
outlet_area = (pi / 4) * d_outer^2;                                        
mu_outer = a_eff / outlet_area;                                            % discharge coefficient area

iterations_outer = 1;
condition_outer = true;
while condition_outer
    % filling efficiency
    f1 = @(E_outer) (E_outer * sqrt(E_outer ...
        / (2 - E_outer)) - mu_outer); 

    E_outer = fzero(f1, .7);

    % geometric characteristic of outer element
    K_outer = (((1 - E_outer) * sqrt(2)) ...
        / (E_outer * sqrt(E_outer)));

    % discharge coefficient
    mu_outer = E_outer * sqrt(E_outer / (2 - E_outer));

    % ratio of the gas core to the outer outlet diameter
    f2 = @(S_outer) ((sqrt(abs(1 - mu_outer^2 * K_outer^2))) ...
    - (S_outer * sqrt(abs(S_outer^2 - mu_outer^2 * K_outer^2))) ...
    - (mu_outer^2 * K_outer^2 * log((1 ...
    + sqrt(abs(1 - mu_outer^2 * K_outer^2))) ...
    / (S_outer + sqrt(abs(S^2 - mu_outer^2 ...
    * K_outer^2))))) - mu_outer);

    S_outer = fzero(f2, .7);

    % recalculating film thickness
    film_thickness = .5 * (d_outer - (S_outer * d_outer));

    % recaulculating outlet diameter
    d_outer_new = (2 * permitted_vortex_rad) + (2 * film_thickness);

     if abs((d_outer_new - d_outer) / d_outer) < .005
         condition_outer = false;
     else
         d_outer = d_outer_new;
     end
     iterations_outer = iterations_outer + 1;
end

R_outer = (d_outer * R_nozzle_outer) / 2;
inlet_diameter_outer = sqrt((2 * R_outer * d_outer) / (outer_inlets * K_outer));
inlet_length_outer = (inlet_diameter_outer / 2) * 3;
chamber_diameter_outer = (2 * R_outer) + inlet_diameter_outer;
chamber_length_outer = R_outer;
nozzle_diameter_outer = d_outer_new;
nozzle_length_outer = 2 * d_outer; 
outer_check = (2 * R_outer * d_outer) / (outer_inlets * inlet_diameter_outer^2);

%% OUTPUTS_________________________________________________________________

fprintf("\nVERSION: WATERFLOW")
fprintf("\nELEMENTS: %d\n", elements)
fprintf("\n---------- INNER ELEMENT DIMENSIONS ---------\n")
fprintf("Tangential inlets: %d\n", inner_inlets)
fprintf("Swirl Arm (in): %f\n", R_inner * 12)
fprintf("Inlet Length (in): %f\n", inlet_length_inner * 12)
fprintf("Inlet Diameter (in): %f\n", inlet_diameter_inner * 12)
fprintf("Chamber Diameter (in): %f\n", chamber_diameter_inner * 12)
fprintf("Chamber Length(in): %f\n", chamber_length_inner * 12)
fprintf("Nozzle diameter (in): %f\n", nozzle_diameter_inner * 12)
fprintf("Nozzle Length (in): %f\n", nozzle_length_inner * 12)
fprintf("Desired Spray Angle (degrees): %f\n", spray_angle)
fprintf("Actual Spray Angle (degrees): %f\n", spray_angle_new)
fprintf("Wall Thickness (in): %f\n", wall_thickness * 12)
fprintf("----------------------------------------------\n")
fprintf("Iterations: %d\n", iterations_inner)
fprintf("\n---------- OUTER ELEMENT DIMENSIONS ----------\n")
fprintf("Tangential inlets: %d\n", outer_inlets)
fprintf("Inlet Length (in): %f\n", inlet_length_outer * 12)
fprintf("Inlet Diameter (in): %f\n", inlet_diameter_outer * 12)
fprintf("Chamber Diameter (in): %f\n", chamber_diameter_outer * 12)
fprintf("Chamber Length (in): %f\n", chamber_length_outer * 12)
fprintf("Nozzle Diameter (in): %f\n", nozzle_diameter_outer * 12)
fprintf("Nozzle Length (in): %f\n", nozzle_length_outer * 12)
fprintf("Wall Thickness (in): %f\n", wall_thickness * 12)
fprintf("-----------------------------------------------\n")
fprintf("Iterations: %d\n", iterations_outer)
