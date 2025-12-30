# Torch Sizing Code
# Authors: Dominik Sloup
# First Created: 06/06/2025
# Last Updated: 06/06/2025
# Calculations done in SI units

import numpy as np
import math
from rocketcea.cea_obj import CEA_Obj
from pyfluids import Fluid, FluidsList, Input

# ──────────────────────────────────────────────────────────────
#  CONSTANTS
# ──────────────────────────────────────────────────────────────

g0 = 32.174 # Gravitational constant (ft/s^2)
gamma_o2 = 1.4 # Oxygen Specific heat ratio
gamma_h2 = 1.41 # Hydrogen Specific heat ratio
gamma_ch4 = 1.3 # Methane Specific heat ratio
T = 293.15 # Ambient temperature (K)


# ──────────────────────────────────────────────────────────────
#  UNIT CONVERSIONS
# ──────────────────────────────────────────────────────────────

lbf_to_N = 4.4482
psi_to_pa = 6894.76
ft_to_m = 0.3048
lb_to_kg = 0.453592
m_to_in = 39.3701
class Torch_Sizer():
    def __init__(self):

        # ──────────────────────────────────────────────────────────────
        #  INPUTS
        # ──────────────────────────────────────────────────────────────

        self.F = 10 # Nominal thrust of the torch (lbf)
        self.p_c = 140 # Nominal chamber pressure (psi)
        self.OF = 7 # Nominal OF Rati
        self.fuel_choice = 'kerosene'

        self.p_ox = 600 # Oxidizer feed pressure (psi)
        self.ox_orifice_num = 1 # Number of oxidizer orifices (N/A)
        self.Cd_ox = 0.8 # Oxidizer orifice anticipated C_d (N/A)

        self.p_fu = 740 # Fuel feed pressure (psi)
        self.fu_orifice_num = 1 # Number of fuel orifices (N/A)
        self.Cd_fu = 0.8 # Fuel orifice anticipated C_d (N/A)

        self.cstar_eff = 0.7 # C star efficiency factor (N/A)

        self.p_ox_pa = self.p_ox * psi_to_pa # Oxidizer feed pressure (Pa)
        self.p_fu_pa = self.p_fu * psi_to_pa # Fuel feed pressure (Pa)
        self.p_c_pa = self.p_c * psi_to_pa # Chamber pressure (Pa)
        self.F_N = self.F * lbf_to_N # Desired thrust (N)


    

    # ──────────────────────────────────────────────────────────────
    #  FUNCTION DEFINITIONS
    # ──────────────────────────────────────────────────────────────


    '''
    Helper function to run a CLI to accept setpoints and update 
    parameters to reflect chosen setpoints.
    '''
    def CLI_inputs(self):
        
        print("Please input nominal thrust (lbf): ")
        self.F = float(input()) # Nominal thrust of the torch (lbf)
        print("Please input nominal chamber pressure (psi): ")
        self.p_c = float(input()) # Nominal chamber pressure (psi)
        print("Please input nominal O/F Mass Ratio: ")
        self.OF = float(input()) # Nominal OF Ratio

        print("Please input your fuel choice (kerosene, methane, hydrogen)")
        self.fuel_choice = input()
        while self.fuel_choice != "kerosene" and self.fuel_choice != "methane" and self.fuel_choice != "hydrogen":
            print("Invalid fuel choice, please input (kerosene, methane, hydrogen)")
            self.fuel_choice = input()

        print("Please input your ox inlet pressure (psi)")
        self.p_ox = float(input()) # Oxidizer feed pressure (psi)

        print("Please input the number of ox inlets")
        self.ox_orifice_num = float(input()) # Number of oxidizer orifices (N/A)

        print("Please input your fuel inlet pressure (psi)")
        self.p_fu = float(input()) # Fuel feed pressure (psi)
        print("Please input the number of fuel inlets")
        self.fu_orifice_num = float(input()) # Number of fuel orifices (N/A)

        # Metric Conversions
        self.p_ox_pa = self.p_ox * psi_to_pa # Oxidizer feed pressure (Pa)
        self.p_fu_pa = self.p_fu * psi_to_pa # Fuel feed pressure (Pa)
        self.p_c_pa = self.p_c * psi_to_pa # Chamber pressure (Pa)
        self.F_N = self.F * lbf_to_N # Desired thrust (N)

    '''
    Main function to generate values based on setpoint parameters
    Returns: 
    - cstar_real, the characteristic velocity of the rocket problem
    - m_dot_fu, the mass flow rate of the fuel (kg/s)
    - m_dot_ox, the mass flow rate of the oxygen (kg/s)
    - A_fu, the area of the fuel orifice(s) (m^2)
    - d_fu, the diameter of the fuel orifice(s) (in)
    - A_ox, the area of the ox orifice(s) (m^2)
    - d_ox, the diameter of the ox orifice(s) (in)
    - d_t, the diameter of the throat (in)
    '''
    def main_solver(self):

        # Oxygen density at manifold pressure and ambient temperature (kg/m^3)
        rho_ox = Fluid(FluidsList.Oxygen).with_state(Input.pressure(self.p_ox_pa), Input.temperature(T-273.15)).density 

        if self.fuel_choice == 'kerosene':
            engine = CEA_Obj(oxName='GOX', fuelName='JetA')
            rho_fu = 810 # Kerosene density (kg/m^3) @ STP
            m_dot_fu, m_dot_ox, cstar_real = self.get_mdots(engine,self.F,self.p_c,self.OF,self.cstar_eff)
            A_fu = self.get_liquid_inj_area(m_dot_fu, self.p_fu_pa,self.p_c_pa, rho_fu, self.Cd_fu)

        elif self.fuel_choice == 'methane':
            engine = CEA_Obj(oxName='GOX', fuelName='CH4(g)')
            rho_fu = Fluid(FluidsList.Methane).with_state(Input.pressure(self.p_fu_pa), Input.temperature(T-273.15)).density
            m_dot_fu, m_dot_ox, cstar_real = self.get_mdots(engine,self.F,self.p_c,self.OF,self.cstar_eff)
            A_fu = self.get_gas_inj_area(m_dot_fu, self.p_fu_pa, self.p_c_pa, rho_fu, self.Cd_fu, gamma_ch4)

        elif self.fuel_choice == 'hydrogen':
            engine = CEA_Obj(oxName='GOX', fuelName='H2(g)')
            rho_fu = Fluid(FluidsList.Hydrogen).with_state(Input.pressure(self.p_fu_pa), Input.temperature(T-273.15)).density
            m_dot_fu, m_dot_ox, cstar_real = self.get_mdots(engine,self.F,self.p_c,self.OF,self.cstar_eff)
            A_fu = self.get_gas_inj_area(m_dot_fu, self.p_fu_pa, self.p_c_pa, rho_fu, self.Cd_fu, gamma_h2)

        A_ox = self.get_gas_inj_area(m_dot_ox, self.p_ox_pa, self.p_c_pa, rho_ox, self.Cd_ox, gamma_o2) 
        A_t = (m_dot_fu + m_dot_ox) * (cstar_real * ft_to_m) / self.p_c_pa

        d_fu = 2 * np.sqrt((A_fu / self.fu_orifice_num) / np.pi) # Fuel orifice diameter (m)
        d_ox = 2 * np.sqrt((A_ox / self.ox_orifice_num) / np.pi) # Ox orifice diameter (m)
        d_t = 2 * np.sqrt(A_t / np.pi) # Throat diameter (m)

        return cstar_real, m_dot_fu, m_dot_ox, A_fu, d_fu, A_ox, d_ox, d_t

    def get_mdots(self, engine,F,p_c,OF,cstar_eff):

        cstar = engine.get_Cstar(Pc = p_c, MR = OF)
        cstar_real = cstar * cstar_eff
        Cf, IspVac, IspSL = engine.get_PambCf(Pamb=14.7, Pc=p_c, MR=OF, eps=1)
        Isp = (cstar_real * Cf / g0)
        m_dot_lb = F / Isp # Total combined mass flow (lb/s)
        m_dot = m_dot_lb * lb_to_kg # Total combined mass flow (kg/s)
        m_dot_fu = m_dot / (1 + OF) # Fuel mass flow (kg/s)
        m_dot_ox = m_dot - m_dot_fu # Oxidizer mass flow (kg/s)

        return m_dot_fu, m_dot_ox, cstar_real

    def get_gas_inj_area(self, mdot, p_feed, p_c, rho, Cd, gamma):

        p_cr = (2 / (gamma + 1)) ** (gamma / (gamma -1))

        if p_c / p_feed < p_cr: # Choked condition
            area = mdot / (Cd * np.sqrt(gamma * rho * p_feed * (2 / (gamma + 1)) ** ((gamma + 1)/(gamma - 1))))

        else: # Unchoked condition
            area = mdot / (Cd * np.sqrt(2 * rho * p_feed * (gamma / (gamma-1)) * (((p_c/p_feed) ** (2/gamma)) - ((p_c/p_feed) ** ((gamma + 1)/gamma)))))

        return area

    def get_liquid_inj_area(self, mdot, p_feed, p_c, rho, Cd):

        return mdot / (Cd * np.sqrt(2 * rho * (p_feed - p_c))) # Total fuel injection area (m^2)

if __name__ == '__main__':
    sizer = Torch_Sizer()
    sizer.CLI_inputs()

    # ──────────────────────────────────────────────────────────────
    #  RESULTS
    # ──────────────────────────────────────────────────────────────

    cstar_real, m_dot_fu, m_dot_ox, A_fu, d_fu, A_ox, d_ox, d_t = sizer.main_solver()

    print(f"Expected characteristic velocity: {math.floor(cstar_real * ft_to_m)} m/s")
    print(f"Total mass flow: {m_dot_fu + m_dot_ox:.4f} kg/s")
    print(f"Fuel mass flow: {m_dot_fu:.4f} kg/s")
    print(f"Ox mass flow: {m_dot_ox:.4f} kg/s")

    print(f"\nThroat diameter: {d_t * m_to_in:.3f} in")

    print(f"\nTotal fuel injection area: {A_fu:.6f} m^2")
    print(f"Total fuel orifice area (in^2): {(d_fu * m_to_in) ** 2 * np.pi:.3f} in^2")
    print(f"Fuel orifice diameter: {d_fu * m_to_in:.3f} in")

    print(f"\nTotal oxidizer injection area: {A_ox:.6f} m^2")
    print(f"Total oxidizer orifice area (in^2): {(d_ox * m_to_in) ** 2 * np.pi:.3f} in^2")
    print(f"Oxidizer orifice diameter: {d_ox * m_to_in:.3f} in")