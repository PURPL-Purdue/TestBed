# Hi! This is a Tank Blowdown/Pressure Ladder Calculator!
import numpy as np
from scipy.optimize import brentq


P_0 = 200; #psi, Tank initial pressure
P_atm = 14; #psi, Tank initial pressure
T_0 = 70; #F,    Tank initial temperature
gamma = 1.4; # Specific heat ratio (for compressible)
R = 287.05; #J/kg/K

timestep = 0.1; # analysis time steps
t_total = 10; # tank blowdown total time
V = 1000; #ft^3, Tank Volume

#Set un-used areas to 99
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

CdLadder = [Cd1,Cd2,Cd3,Cd4,Cd5,Cd6];
CdLadder = CdLadder[CdLadder !=99]

#unit conversions
P_0 = P_0*6894.76; #Pa
P_atm = P_atm*6894.76; #Pa
T_0 = ((T_0-32)*5/9)+273.15; #F
Ladder = Ladder *0.00064516; #m
r_crit = (2 / (gamma + 1))**(gamma / (gamma - 1)); # crit pressure ratio

Ladder_P = np.zeros(len(Ladder),t_total/timestep); # Array of Pressures over time
Ladder_T = np.zeros(len(Ladder),t_total/timestep); # Array of Temperatures over time
Ladder_P(1,1) = P_0; # Set first index to starting tank pressure 
Ladder_T(1,1) = T_0; # Set first index to starting tank temperature

Ladder_P(len(Ladder_P),1) = P_atm; # Set final pressure to atmospheric

if Solver == 'Compressible':
    time = np.linspace(0.1, 10, 100); # Create time steps

    for t in time: # Begin loop of time steps

        dP = np.zeros(len(Ladder)); # Set array of pressure drops across orifice
        dP(1) = P_atm; # Set first index to equal atm pressure

        while abs(sum(dP) - P_o) < 0.001:
            
            mdot_guess = mdot_guess + 0.001; # set initial guess of last pressure
            
            for A in Ladder.reverse():


                # Step 1: Assume choked flow and solve for P1
                choked_factor = np.sqrt(gamma / (R * Tank_T(t/timestep))) * (2 / (gamma + 1))**((gamma + 1) / (2 * (gamma - 1)));

                P1_choked = mdot_guess / (CdLadder(t/timestep) * choked_factor)

                # Step 2: Check pressure ratio
                ratio_choked = Ladder_P(len(Ladder), t/timestep) / P1_choked

                if ratio_choked <= r_crit:
                    # Flow is choked
                    P1 = P1_choked
                    flow_regime = "Choked"

                else:
                    # Step 3: Solve unchoked mass flow equation
                    def massflow_residual(P1):
                        r = P2 / P1

                        mdot_calc = CdA * P1 * np.sqrt(
                            (2 * gamma / (R * T1 * (gamma - 1))) *
                            (r**(2/gamma) -
                            r**((gamma + 1)/gamma))
                        )

                        return mdot_calc - mdot

                    # Unchoked bounds: P2 < P1 < P2/r_crit
                    P1 = brentq(
                        massflow_residual,
                        P2,
                        P2 / r_crit
                    )
                    flow_regime = "Unchoked"
