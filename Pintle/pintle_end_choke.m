in_to_m = 0.0254;
lb_to_kg = 0.4536;
N_to_lbf = 0.2248;

% Input Parameters
flow_disturbance = 0.05; % Maximal Flow Disturbance Ratio
tolerance = 4.0 * (in_to_m / 1000.0); % +- in Thou
stiffness = 0.2; % Ratio of Pressure Drop to Chamber Pressure
shaft_ratio = 0.6; % Pintle Shaft to Sleeve Radius (Dom's Vibes)
sleeve_thickness = 125 * (in_to_m / 1000.0); % Thickness of Sleeve Wall

chamber_radius = (4.5 / 2) * in_to_m;
chamber_length = (4.5 / 2) * in_to_m;

of_ratio = 1.50087; % Architecture Number
max_mass_flow = 8.82281 * lb_to_kg;

throttle = 0.5; % Minimum Throttle

servo_angle = 180.0 * (pi / 180.0); % Maximum Servo Angle

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

static_friction = 0.800; % Steel-on-Steel

% Calculations
fprintf("\n\n\nDesign Dimensions:\n");
fprintf("--------------------------------------\n");
v_lox = C_d * sqrt(2 * stiffness * P_c / rho_lox);
v_f = C_d * sqrt(2 * stiffness * P_c / rho_f);

beta = atan((4/3) * chamber_radius / chamber_length);

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
R_pt = R_cg;

fprintf("Fuel Area (thou^2): %.2e\n", A_pg * (1000.0 / in_to_m)^2);
fprintf("Central Gap Radius (thou): %.2f\n", R_cg * 1000.0 / in_to_m);
fprintf("Sleeve Radius (thou): %.2f\n", R_sv * 1000.0 / in_to_m);
fprintf("Pintle Head Radius (thou): %.2f\n", R_pt * 1000.0 / in_to_m);
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

L_min = (R_cg - sqrt(R_cg^2 - A_pg * throttle * sin(theta_pt) / pi)) / sin(theta_pt);
L_open_min = L_min / cos(theta_pt);
L_open_diff = L_open - L_open_min;
fprintf("Actual Opening Gap (thou): %.2f\n", L_a * 1000.0 / in_to_m);
fprintf("Lopen for Max Throttle (100%%) (thou): %.2f\n", L_open * 1000.0 / in_to_m);
fprintf("Lopen for Min Throttle (%02.0f%%) (thou): %.2f\n", throttle * 100.0, L_open_min * 1000.0 / in_to_m);
fprintf("Pintle Transverse Distance (thou): %.2f\n", L_open_diff * 1000.0 / in_to_m);
fprintf("--------------------------------------\n");

P_f = P_mf + rho_f * v_f^2 / 2;

F_D = pi * (P_f - P_c) * R_pt^2 - P_f * pi * R_pr^2;
thread_angle = atan(L_open_diff / (servo_angle * R_pr));
torque_max = F_D * (static_friction * cos(thread_angle) + sin(thread_angle)) / (cos(thread_angle) - static_friction * sin(thread_angle)) * R_pr;

fprintf("Axial Load (lbf, + = Into Chamber): %+.2f\n", F_D * N_to_lbf);
fprintf("Thread Angle (degrees): %.2f\n", thread_angle * 180.0 / pi);
fprintf("Maximum Thread Pitch (TPI): %.2f\n", 1.0 / (2 * pi * R_pr / in_to_m * tan(thread_angle))); 
fprintf("Min Thread Pitch (TPI): %.2f\n", 1.0 / (2 * pi * R_pr / in_to_m * static_friction));
fprintf("Max Torque on Servo (lb * in): %.2f\n", torque_max * N_to_lbf / in_to_m);
fprintf("--------------------------------------\n");

head_y_comp = ((P_f + P_c) * R_pt^2 - P_f * R_pr^2)/ R_pt / R_pr / youngs;
shaft_strain = F_D / (pi * R_pr^2) / youngs;
fprintf("Head Height Compression (%%%%, thou): %.2f\t%.2e\n", 10000.0 * head_y_comp, head_y_comp * (R_pt-R_pr) * tan(theta_pt) * 1000.0 / in_to_m);
fprintf("Shaft Strain (%%%%): %.2f\n", shaft_strain * 10000.0);

