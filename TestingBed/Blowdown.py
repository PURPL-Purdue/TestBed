# Hi! This is a Tank Blowdown/Pressure Ladder Calculator!
import numpy as np
P_0 = 200; #psi, Tank initial pressure
P_atm = 14; #psi, Tank initial pressure
T_0 = 70; #F,    Tank initial temperature
timestep = 0.1; # analysis time steps
t_total = 10; # tank blowdown total time
V = 1000; #ft^3, Tank Volume
A1 = 1; #in^2,   Orifice 1 Area
Cd1 = 0.9; #NA,  Orifice 1 Cd
A2 = 1; #in^2,   Orifice 2 Area
Cd2 = 0.9; #NA,  Orifice 2 Cd
A3 = 1; #in^2,   Orifice 3 Area
Cd3 = 0.9; #NA,  Orifice 3 Cd
A4 = 999; #in^2,   Orifice 4 Area
Cd4 = 0.9; #NA,  Orifice 4 Cd
A5 = 999; #in^2,   Orifice 5 Area
Cd5 = 0.9; #NA,  Orifice 5 Cd
A6 = 999; #in^2,   Orifice 6 Area
Cd6 = 0.9; #NA,  Orifice 6 Cd
Solver = 'Compressible' # 'Compressible' solves Isothermal & Isentropic Tank Blowdown with Pressure ladder, 'Incompressible' Solves Constant Pressure Tank Pressure Ladder

mdot_guess = 0;
Ladder = [A1,A2,A3,A4,A5,A6];
Ladder = Ladder[Ladder != 99]
Tank_P = np.zeros(t_total/timestep); # Array of Pressures over time
Tank_T = np.zeros(t_total/timestep); # Array of Temperatures over time
Tank_P(1) = P_0; # Set first index to starting tank pressure 
Tank_T(1) = T_0; # Set first index to starting tank temperature

if Solver == 'Compressible':
    time = np.linspace(0.1, 10, 100); # Create time steps
    for t in time: # Begin loop of time steps

        dP = []; # Set array of pressure drops across orifice
        
        while abs(sum(dP) - P_o) < 0.001:
            mdot_guess = mdot_guess + 0.001; # set initial guess of last pressure
            for A in Ladder.reverse():
                d
                