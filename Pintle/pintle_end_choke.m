in_to_m = 0.0254;
lb_to_kg = 0.4536;

% Input Parameters
flow_disturbance = 0.05; % Maximal Flow Disturbance Ratio
tolerance = 4.0 * (in_to_m / 1000.0); % +- in Thou
stiffness = 0.2; % Ratio of Pressure Drop to Chamber Pressure
shaft_ratio = 0.6; % Dominik's Vibes
head_ratio = 0.95; % Pintle Head : Sleeve Radius
sleeve_thickness = 125 * (in_to_m / 1000.0); % Thickness of Sleeve Wall

chamber_radius = (4.5 / 2) * in_to_m;
chamber_length = (4.5 / 2) * in_to_m;

of_ratio = 0.5;
max_mass_flow = 22.8 * lb_to_kg;

% Constants
psi_to_Pa = 6894.757;

mdot_f = max_mass_flow;
u_f = 1136;
rho_f = 768.09;

mdot_lox = of_ratio * mdot_f;
u_lox = 1014;
rho_lox = 1205;

P_c = 600 * psi_to_Pa;
P_mf = 900 * psi_to_Pa;
P_mlox = 720 * psi_to_Pa;


C_d = 1.0;
youngs = 117 * 10^9;

% Calculations
fprintf("\n\n\nDesign Dimensions:\n");
fprintf("--------------------------------------\n");
v_lox = C_d * sqrt(2 * stiffness * P_c / rho_lox);
v_f = C_d * sqrt(2 * stiffness * P_c / rho_f);

beta = atan(2 * chamber_radius / chamber_length);

theta_pt = acos((mdot_lox * v_lox * (1/cos(beta) - 1)) / (mdot_f * v_f * sqrt((1/cos(beta) - 1)^2 + 1))) - atan(1/cos(beta) - 1);

fprintf("Half-Spray Angle (degrees): %.2f\n", beta * 180 / pi);
fprintf("Lox Velocity (m/s): %.2f\n", v_lox);
fprintf("RP1 Velocity (m/s): %.2f\n", v_f);
fprintf("Pintle Head Angle (degrees): %.2f\n", theta_pt * 180 / pi);
fprintf("--------------------------------------\n");

A_pg = mdot_f / (C_d * sqrt(2 * stiffness * P_c * rho_f));

R_cg = (2*shaft_ratio^2*sleeve_thickness + sqrt(4*shaft_ratio^2*sleeve_thickness^2 + 4 * (1 - shaft_ratio^2) * A_pg / pi)) / 2 / (1 - shaft_ratio^2);

R_sv = sleeve_thickness + R_cg;
R_pr = R_sv * shaft_ratio;
R_pt = R_sv * head_ratio;

fprintf("Fuel Area (thou^2): %.2e\n", A_pg * (1000.0 / in_to_m)^2);
fprintf("Central Gap Radius (thou): %.2f\n", R_cg * 1000.0 / in_to_m);
fprintf("Sleeve / Pintle Head Radius (thou): %.2f\n", R_sv * 1000.0 / in_to_m);
fprintf("Pintle Shaft Radius (thou): %.2f\n", R_pr * 1000.0 / in_to_m);
fprintf("--------------------------------------\n");

A_ann = mdot_lox / (C_d * sqrt(2 * rho_lox * stiffness * P_c));

R_tot = sqrt(A_ann / pi + R_sv^2);
delta_ann = R_tot - R_sv;

fprintf("LOx Area (thou^2): %.2e\n", A_ann * (1000.0 / in_to_m)^2);
fprintf("Total Radius (thou): %.2f\n", R_tot * 1000.0 / in_to_m);
fprintf("Annular Gap (thou): %.2f\n", delta_ann * 1000.0 / in_to_m);
fprintf("--------------------------------------\n");

L_a = (R_cg - sqrt(R_cg^2 - A_pg * sin(theta_pt) / pi)) / sin(theta_pt);
L_open = L_a / cos(theta_pt);

fprintf("Actual Opening Gap (thou): %.2f\n", L_a * 1000.0 / in_to_m);
fprintf("Lopen (thou): %.2f\n", L_open * 1000.0 / in_to_m);
fprintf("--------------------------------------\n");

P_f = P_mf + rho_f * v_f^2 / 2;

head_y_strain = (R_pt^2 - R_pr^2) / (youngs * tan(theta_pt)) / ((P_f + P_c)*R_pt^2 - P_f*R_pr^2);
head_x_strain = pi * (R_pt^2 - R_pr^2) / youngs / tan(theta_pt);

fprintf("Head Height Deformation (%%%%strain, thou): %.2f%%%%\t%.2e\n", 10000.0 * head_y_strain, head_y_strain * (R_pt-R_pr) * tan(theta_pt) * 1000.0 / in_to_m);
fprintf("Head Radius Deformation (%%%%strain, thou): %.2f%%%%\t%.2e\n", 10000.0 * head_x_strain, head_x_strain * R_pt * 1000.0 / in_to_m);

