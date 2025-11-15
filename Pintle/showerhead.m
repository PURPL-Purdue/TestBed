in_to_m = 0.0254;
lb_to_kg = 0.4536;
N_to_lbf = 0.2248;

% Input Parameters
flow_disturbance = 0.05; % Maximal Flow Disturbance Ratio
tolerance = 4.0 * (in_to_m / 1000.0); % +- in Thou
stiffness = 0.2; % Ratio of Pressure Drop to Chamber Pressure
sleeve_thickness = 100 * (in_to_m / 1000.0); % Thickness of Sleeve Wall
pintle_wall = 50 * (in_to_m / 1000.0);
min_material = 1.0 - .85; % Fraction of Material Needed to be Left in Pintle
thread_backing = 50.0 * (in_to_m / 1000.0); % How Much Material Needs to be Behind a Thread

chamber_radius = (4.5 / 2) * in_to_m;
chamber_length = (4.5 / 2) * in_to_m;

of_ratio = 1.50087; % Architecture Number
max_mass_flow = 8.82281 * lb_to_kg;

throttle = 0.2; % Minimum Throttle

servo_angle = 180.0 * (pi / 180.0); % Maximum Servo Angle

safety = 2.0;

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
fprintf("Sleeve Thickness (thou): %.0f\n", sleeve_thickness * 1000.0 / in_to_m);
fprintf("Pintle Wall Thickness (thou): %.0f\n", pintle_wall * 1000.0 / in_to_m);
fprintf("Minimum Material Left After Slotting (%%): %.0f\n", min_material * 100.0);
fprintf("--------------------------------------\n");

A_pg = mdot_f / (C_d * sqrt(2 * stiffness * P_c * rho_f));
A_ann = mdot_lox / (C_d * sqrt(2 * rho_lox * stiffness * P_c));

R_od = sqrt(A_pg / pi) + pintle_wall;
R_tot = sqrt(A_ann / pi + R_od^2 + 2 * R_od * sleeve_thickness + sleeve_thickness^2);
R_sv = R_od + sleeve_thickness;
delta_ann = R_tot - R_sv;

R_finlet = sqrt(A_pg / pi);
R_oxinlet = sqrt(A_ann / pi);

fprintf("Mass Flow LOx (lb/s): %.2f\n", mdot_lox / lb_to_kg);
fprintf("Mass Flow RP1 (lb/s): %.2f\n", mdot_f / lb_to_kg);
fprintf("Fuel Area (thou^2): %.2e\n", A_pg * (1000.0 / in_to_m)^2);
fprintf("LOx Area (thou^2): %.2e\n", A_ann * (1000.0 / in_to_m)^2);
fprintf("Fuel Inlet Diameter Min. (thou): %.0f\n", 2 * R_finlet * 1000.0 / in_to_m);
fprintf("LOx Inlet Diameter Min. (thou): %.0f\n", 2 * R_oxinlet * 1000.0 / in_to_m);
fprintf("--------------------------------------\n");
fprintf("Total Diameter (thou): %.0f\n", 2 * R_tot * 1000.0 / in_to_m);
fprintf("Pintle Diameter (thou): %.0f\n", 2 * R_od * 1000.0 / in_to_m);
fprintf("Pintle Bore Diameter (thou): %.0f\n", 2 * (R_od - pintle_wall) * 1000.0 / in_to_m);
fprintf("Sleeve Tip Diameter (thou): %.0f\n", 2 * R_sv * 1000.0 / in_to_m);
fprintf("Annular Gap (thou): %.0f\n", delta_ann * 1000.0 / in_to_m);
fprintf("--------------------------------------\n");

fuel_manifold_height = safety * A_pg / (2 * pi * (R_od - pintle_wall));
lox_manifold_height = safety * A_ann / (2 * pi * R_sv);

fprintf("Fuel Manifold Height Min. (thou): %.0f\n", fuel_manifold_height * 1000.0 / in_to_m);
fprintf("LOx Manifold Height Min. (thou): %.0f\n", lox_manifold_height * 1000.0 / in_to_m);
fprintf("--------------------------------------\n");

slot_height = A_pg * 1 / (2 * pi * (R_od - pintle_wall) * (1 - min_material));
sleeve_angle = asin(1.0 / (slot_height * (1-throttle)) * (sqrt(R_tot^2 - A_ann*throttle/pi)-R_sv));

fprintf("Slot Height @ Max (thou): %.0f\n", slot_height * 1000.0 / in_to_m);
fprintf("Slot Height @ %.0f%% (thou): %.0f\n", throttle * 100.0, slot_height * throttle * 1000.0 / in_to_m);
fprintf("Sleeve Transverse Distance: (thou): %.0f\n", slot_height * (1 - throttle) * 1000.0 / in_to_m);
fprintf("Amount of Material Left on Circumference (thou): %.0f\n", 2 * pi * R_od * min_material * 1000.0 / in_to_m);
fprintf("Sleeve Angle from Vertical (degrees): %.0f\n", sleeve_angle * 180.0 / pi);
fprintf("--------------------------------------\n");

F_U = P_c * pi * (R_sv^2 - R_od^2) + P_mlox * pi * tan(sleeve_angle) * (R_tot^2 - R_sv^2);

torque_fn = @(g)torque(g, servo_angle, static_friction, F_U, R_od + thread_backing, throttle, slot_height, R_sv + sin(sleeve_angle) * slot_height * (1-throttle));
[gear_ratio, servo_torque] = fminbnd(torque_fn, 0.05, 20.0);
thread_angle = atan(slot_height * (1 - throttle) / (servo_angle * gear_ratio * R_tot));

fprintf("Force On Sleeve Upwards, (lbf): %.0f\n", F_U * N_to_lbf);
fprintf("Gear Ratio 1:%.1f (%.2f)\n", 1 / gear_ratio, gear_ratio);
fprintf("Minimum Thread Angle (degrees): %.1f\n", thread_angle * 180.0 / pi);
fprintf("Maximum Servo Torque (lb*in): %.1f\n", servo_torque * N_to_lbf / in_to_m);

function t = torque(gear_ratio, servo_angle, friction, force, thread_radius, throttle, slot_height, gear_radius)
	thread_angle = atan(slot_height * (1 - throttle) / (servo_angle * gear_ratio * thread_radius));
	t = gear_ratio * thread_radius^2 / gear_radius * force * (friction * cos(thread_angle) + sin(thread_angle)) ...
		/ (cos(thread_angle) - friction * sin(thread_angle));
end

