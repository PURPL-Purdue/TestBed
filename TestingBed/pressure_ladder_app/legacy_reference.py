import numpy as np
from scipy.optimize import brentq

# ============================================================
# USER INPUTS
# ============================================================
P_0 = 200.0       # psia, constant tank pressure
P_atm = 14.0      # psia, outlet pressure
P_vap = 10.0      # psia, liquid vapor pressure
T_0 = 70.0        # deg F, tank stagnation temperature
gamma = 1.4       # ideal gas specific heat ratio
R = 287.05        # J/(kg K), ideal gas constant
density = 62.428  # lbm/ft^3, liquid only
mu = 0.001        # Pa*s, liquid dynamic viscosity
mu_gas = 1.81e-5  # Pa*s, gas dynamic viscosity (air near room T)

timestep = 0.1    # reserved for future blowdown solver
t_total = 10.0
V = 1000.0        # ft^3, reserved for blowdown

# Set unused ORIFICE AREAS to 999.
A1, Cd1, beta1 = 1, 0.9, 0.5
A2, Cd2, beta2 = 1, 0.9, 0.5
A3, Cd3, beta3 = 1, 0.9, 0.5
A4, Cd4, beta4 = 999, 999, 999
A5, Cd5, beta5 = 999, 999, 999
A6, Cd6, beta6 = 999, 999, 999

# Empirical liquid cavitation parameters (illustrative only).
FL1, FL2, FL3 = 0.85, 0.85, 0.85
FL4, FL5, FL6 = 0.85, 0.85, 0.85
Fstart1, Fstart2, Fstart3 = 0.65, 0.65, 0.65
Fstart4, Fstart5, Fstart6 = 0.65, 0.65, 0.65

# Pipe diameter (in), length AFTER each orifice (ft), and roughness (in).
# Common pipe diameter required by the loss model below.
# D=2.25676 in corresponds to beta ~0.5 for a 1 in^2 bore.
D1, D2, D3, D4, D5, D6 = 2.25676, 2.25676, 2.25676, 999, 999, 999
L1, L2, L3, L4, L5, L6 = 2.0, 2.0, 2.0, 999, 999, 999
eps1, eps2, eps3, eps4, eps5, eps6 = 0.0018, 0.0018, 0.0018, 999, 999, 999

FlowSolver = 'Compressible'  # 'Incompressible' or 'Compressible'
PressureSolver = 'Steady State'

# ============================================================
# BUILD ARRAYS USING ONE ACTIVE MASK
# ============================================================
Ladder = np.array([A1, A2, A3, A4, A5, A6], dtype=float)
active = Ladder != 999
Ladder = Ladder[active]
CdLadder = np.array([Cd1, Cd2, Cd3, Cd4, Cd5, Cd6], dtype=float)[active]
betaLadder = np.array([beta1, beta2, beta3, beta4, beta5, beta6], dtype=float)[active]
FLLadder = np.array([FL1, FL2, FL3, FL4, FL5, FL6], dtype=float)[active]
FstartLadder = np.array([Fstart1, Fstart2, Fstart3, Fstart4, Fstart5, Fstart6], dtype=float)[active]
PipeD = np.array([D1, D2, D3, D4, D5, D6], dtype=float)[active]
PipeL = np.array([L1, L2, L3, L4, L5, L6], dtype=float)[active]
PipeRough = np.array([eps1, eps2, eps3, eps4, eps5, eps6], dtype=float)[active]
N = len(Ladder)
if N == 0:
    raise ValueError('No active orifices')
if np.any(Ladder <= 0) or np.any(CdLadder <= 0):
    raise ValueError('Active area and Cd must be positive')

# ============================================================
# CONVERT TO SI
# ============================================================
PSI_TO_PA = 6894.757293
P_0 *= PSI_TO_PA
P_atm *= PSI_TO_PA
P_vap *= PSI_TO_PA
T_0 = (T_0 - 32) * 5/9 + 273.15
density *= 16.018463
Ladder *= 0.00064516
PipeD *= 0.0254
PipeL *= 0.3048
PipeRough *= 0.0254

if P_0 <= P_atm:
    raise ValueError('Tank pressure must exceed outlet pressure')
if PressureSolver != 'Steady State':
    raise NotImplementedError('Only steady state is implemented')
if FlowSolver not in ('Incompressible', 'Compressible'):
    raise ValueError('Unsupported FlowSolver')

# ============================================================
# INITIALIZE SHARED RESULTS
# ============================================================
Ladder_P = np.zeros(N + 1)  # recovered node pressures (incompressible)
Ladder_T = np.full(N + 1, T_0)  # assumed interstage stagnation T (compressible)
dP = np.zeros(N)
cav = np.zeros(N, dtype=bool)
choked = np.zeros(N, dtype=bool)
mdot = 0.0

# ============================================================
# INCOMPRESSIBLE: ORIFICE PERMANENT LOSS + PIPE FRICTION
# ============================================================
if FlowSolver == 'Incompressible':
    if P_0 <= P_vap:
        raise ValueError('Tank pressure must exceed vapor pressure')
    if mu <= 0 or np.any(PipeD <= 0) or np.any(PipeL < 0) or np.any(PipeRough < 0):
        raise ValueError('Invalid pipe geometry, roughness, or liquid viscosity')
    if not np.allclose(PipeD, PipeD[0], rtol=1e-6):
        raise ValueError('This version requires one common pipe ID')
    if np.any((FLLadder <= 0) | (FLLadder > 1)):
        raise ValueError('FL must be in (0, 1]')
    if np.any((FstartLadder <= 0) | (FstartLadder > FLLadder)):
        raise ValueError('Fstart must be in (0, FL]')

    PipeArea = np.pi * PipeD**2 / 4
    OrificeD = np.sqrt(4 * Ladder / np.pi)
    betaGeom = OrificeD / PipeD
    if np.any((betaGeom <= 0) | (betaGeom >= 1)):
        raise ValueError('Orifice diameter must be smaller than pipe diameter')
    if not np.allclose(betaGeom, betaLadder, rtol=0.02):
        raise ValueError('Input beta must agree with orifice and pipe diameters')
    betaLadder = betaGeom

    K = CdLadder * Ladder * np.sqrt(2 * density / (1 - betaLadder**4))

    def OrificeLossFactor(i):
        b, Cd = betaLadder[i], CdLadder[i]
        x = np.sqrt(1 - b**4 * (1 - Cd**2))
        return (x - Cd * b**2) / (x + Cd * b**2)

    def DarcyFrictionFactor(Re, diameter, roughness):
        if Re <= 0:
            return 0.0
        if Re < 2300:
            return 64 / Re
        f_turb = (-1.8 * np.log10((roughness/(3.7*diameter))**1.11 + 6.9/Re))**(-2)
        if Re >= 4000:
            return f_turb
        w = (Re - 2300) / (4000 - 2300)
        return (1-w) * (64/Re) + w * f_turb

    def CalculateLosses(mdot_guess, i):
        dP_meter = (mdot_guess / K[i])**2
        dP_orifice = OrificeLossFactor(i) * dP_meter
        V_pipe = mdot_guess / (density * PipeArea[i])
        Re = density * V_pipe * PipeD[i] / mu
        f_D = DarcyFrictionFactor(Re, PipeD[i], PipeRough[i])
        dP_pipe = f_D * (PipeL[i] / PipeD[i]) * (density * V_pipe**2 / 2)
        return dP_meter, dP_orifice, dP_pipe, Re, f_D

    def CalculateLadder(mdot_guess, return_details=False):
        P = np.zeros(N + 1)
        P[-1] = P_atm
        meter_arr = np.zeros(N)
        orifice_arr = np.zeros(N)
        pipe_arr = np.zeros(N)
        Re_arr = np.zeros(N)
        f_arr = np.zeros(N)
        choke_governs = np.zeros(N, dtype=bool)
        for i in range(N - 1, -1, -1):
            dm, do, dp, Re, f = CalculateLosses(mdot_guess, i)
            P_after_orifice = P[i+1] + dp
            P_normal = P_after_orifice + do
            P_choke = P_vap + dm / FLLadder[i]**2
            P[i] = max(P_normal, P_choke)
            choke_governs[i] = P_choke > P_normal + 1e-7
            meter_arr[i], orifice_arr[i], pipe_arr[i] = dm, do, dp
            Re_arr[i], f_arr[i] = Re, f
        if return_details:
            return P, meter_arr, orifice_arr, pipe_arr, Re_arr, f_arr, choke_governs
        return P

    def PressureResidual(mdot_guess):
        return CalculateLadder(mdot_guess)[0] - P_0

    mdot_upper = np.min(K * FLLadder * np.sqrt(P_0 - P_vap))
    mdot = brentq(PressureResidual, 0, mdot_upper, xtol=1e-12)
    Ladder_P, dP_meter, dP_orifice, dP_pipe, Reynolds, Friction, choked = CalculateLadder(mdot, True)
    dP = Ladder_P[:-1] - Ladder_P[1:]
    for i in range(N):
        dP_start = FstartLadder[i]**2 * (Ladder_P[i] - P_vap)
        cav[i] = (dP_meter[i] >= dP_start) or choked[i]

# ============================================================
# COMPRESSIBLE: RECOVERED ORIFICE LOSS + COMPRESSIBLE PIPE LOSS
# ============================================================
else:
    if gamma <= 1 or R <= 0 or T_0 <= 0 or mu_gas <= 0:
        raise ValueError('Invalid gas properties or viscosity')
    if np.any(PipeD <= 0) or np.any(PipeL < 0) or np.any(PipeRough < 0):
        raise ValueError('Invalid pipe geometry')
    if not np.allclose(PipeD, PipeD[0], rtol=1e-6):
        raise ValueError('This approximation requires a common pipe diameter')

    PipeArea = np.pi * PipeD**2 / 4
    betaGeom = np.sqrt(4 * Ladder / np.pi) / PipeD
    if np.any((betaGeom <= 0) | (betaGeom >= 1)):
        raise ValueError('Orifice diameter must be below pipe diameter')
    if not np.allclose(betaGeom, betaLadder, rtol=0.02):
        raise ValueError('Input beta disagrees with area and pipe diameter')
    betaLadder = betaGeom

    # Empirical approximation: fraction of orifice-plate metering
    # differential which is permanently lost after pressure recovery.
    # This incompressible ISO-type relationship is NOT a validated
    # compressible recovery model for high pressure-ratio flows.
    x = np.sqrt(1 - betaLadder**4 * (1 - CdLadder**2))
    loss_fraction = (x - CdLadder * betaLadder**2) / (x + CdLadder * betaLadder**2)
    if np.any((loss_fraction <= 0) | (loss_fraction > 1)):
        raise ValueError('Invalid recovery fractions')

    r_crit = (2 / (gamma + 1))**(gamma / (gamma - 1))
    C_star = np.sqrt(gamma / (R*T_0)) * (2/(gamma+1))**(
        (gamma+1)/(2*(gamma-1)))

    def DarcyFrictionFactor(Re, diameter, roughness):
        if Re <= 0:
            return 0.0
        if Re < 2300:
            return 64/Re
        turbulent = (-1.8*np.log10((roughness/(3.7*diameter))**1.11
                                    + 6.9/Re))**(-2)
        if Re >= 4000:
            return turbulent
        weight = (Re - 2300)/1700
        return (1-weight)*64/Re + weight*turbulent

    # ONE orifice model: use a virtual lower pressure at the
    # orifice metering/jet section, calculated from the net loss.
    # DO NOT add a second orifice loss after this calculation.
    def CompressibleMdot(P_up, P_recovered, i):
        if P_up <= P_recovered or P_up <= 0:
            return 0.0
        delta_perm = P_up - P_recovered
        delta_meter = delta_perm / loss_fraction[i]
        P_virtual = max(P_up - delta_meter, 0.0)
        r_eff = max(P_virtual/P_up, r_crit)
        term = r_eff**(2/gamma) - r_eff**((gamma+1)/gamma)
        return (CdLadder[i] * Ladder[i] * P_up *
                np.sqrt(2*gamma*term/(R*T_0*(gamma-1))))

    def SolveOrifice(mdot_guess, P_recovered, i):
        if mdot_guess <= 0:
            return P_recovered
        # A choked-capacity lower bound for the upstream pressure.
        P_choke = mdot_guess/(CdLadder[i]*Ladder[i]*C_star)
        lower = P_recovered
        upper = max(lower*1.1, lower + 100.0)
        while CompressibleMdot(upper, P_recovered, i) < mdot_guess:
            upper *= 2
            if upper > 1e12:
                raise ValueError('Could not bracket gas orifice pressure')
        return brentq(lambda Pu: CompressibleMdot(Pu, P_recovered, i)
                     - mdot_guess, lower, upper, xtol=1e-8)

    # Gas pipe isothermal, low-Mach pressure-squared approximation:
    # Pin^2 - Pout^2 = f_D (L/D) G^2 R T
    # G = mdot / Apipe. This uses local ideal-gas density rather
    # than the tank density. It excludes acceleration effects;
    # near sonic conditions, use Fanno flow instead.
    def SolvePipe(mdot_guess, P_down, i):
        G = mdot_guess / PipeArea[i]
        Re = G * PipeD[i] / mu_gas
        f_D = DarcyFrictionFactor(Re, PipeD[i], PipeRough[i])
        pressure_in = np.sqrt(P_down**2 + f_D * PipeL[i] / PipeD[i]
                              * G**2 * R * T_0)
        return pressure_in, Re, f_D

    def CalculateLadder(mdot_guess, return_details=False):
        P = np.zeros(N+1)
        P[-1] = P_atm
        pipe_inlets = np.zeros(N)
        orifice_losses = np.zeros(N)
        pipe_losses = np.zeros(N)
        meter_differentials = np.zeros(N)
        Re_arr = np.zeros(N)
        f_arr = np.zeros(N)
        choke_arr = np.zeros(N, dtype=bool)
        mach_out_arr = np.zeros(N)

        for i in range(N-1, -1, -1):
            # Recover required pressure upstream of downstream pipe.
            P_after, Re, f = SolvePipe(mdot_guess, P[i+1], i)
            # Solve upstream pressure across restriction, with recovery.
            P[i] = SolveOrifice(mdot_guess, P_after, i)
            dp_perm = P[i] - P_after
            dp_meter = dp_perm / loss_fraction[i]
            P_virtual = max(P[i] - dp_meter, 0.0)
            choke_arr[i] = P_virtual/P[i] <= r_crit*(1+1e-9)
            pipe_inlets[i] = P_after
            orifice_losses[i] = dp_perm
            pipe_losses[i] = P_after - P[i+1]
            meter_differentials[i] = dp_meter
            Re_arr[i] = Re
            f_arr[i] = f
            # Isothermal sound speed for this pipe model: sqrt(R*T)
            mach_out_arr[i] = (mdot_guess / PipeArea[i]) * np.sqrt(R*T_0) / P[i+1]

        if return_details:
            return (P, pipe_inlets, meter_differentials,
                    orifice_losses, pipe_losses, Re_arr,
                    f_arr, choke_arr, mach_out_arr)
        return P

    def PressureResidual(mdot_guess):
        return CalculateLadder(mdot_guess)[0] - P_0

    # Upper bound: smallest theoretical choked capacity at full tank P0.
    mdot_upper = np.min(CdLadder * Ladder * P_0 * C_star)
    mdot = brentq(PressureResidual, 0.0, mdot_upper, xtol=1e-12)
    (Ladder_P, PipeInletP, dP_meter, dP_orifice, dP_pipe,
     Reynolds, Friction, choked, MachPipeOut) = CalculateLadder(mdot, True)
    dP = Ladder_P[:-1] - Ladder_P[1:]

# ============================================================
# OUTPUT
# ============================================================
print('\n========== FLOW RESULTS ==========')
print(f'Solver: {FlowSolver}')
print(f'Mass flow: {mdot:.6f} kg/s')
print(f'Mass flow: {mdot*2.20462262:.6f} lbm/s')
if FlowSolver == 'Incompressible':
    print(f'Liquid volume flow: {mdot/density*15850.323:.3f} GPM')
else:
    rho_tank = P_0 / (R*T_0)
    print(f'Tank gas density: {rho_tank:.4f} kg/m^3')
    print(f'Tank volume flow: {mdot/rho_tank:.6f} m^3/s')
    print(f'Critical pressure ratio: {r_crit:.5f}')

print('\n========== PRESSURE LADDER ==========')
for i, p in enumerate(Ladder_P):
    print(f'Node {i}: {p/PSI_TO_PA:.3f} psia')

print('\n========== STAGE RESULTS ==========')
for i in range(N):
    print(f'\nStage {i+1}:')
    print(f'  Total node pressure drop: {dP[i]/PSI_TO_PA:.3f} psi')
    if FlowSolver == 'Incompressible':
        print(f'  Orifice meter differential: {dP_meter[i]/PSI_TO_PA:.3f} psi')
        print(f'  Nominal orifice permanent loss: {dP_orifice[i]/PSI_TO_PA:.3f} psi')
        print(f'  Pipe friction loss: {dP_pipe[i]/PSI_TO_PA:.3f} psi')
        print(f'  Reynolds number: {Reynolds[i]:.0f}')
        print(f'  Darcy friction factor: {Friction[i]:.5f}')
        print(f'  Cavitation onset (screening): {cav[i]}')
        print(f'  Cavitation choking constraint governs: {choked[i]}')
        if choked[i]:
            print('  NOTE: Interstage pressure in cavitating regime is approximate.')
    else:
        print(f'  Permanent orifice loss: {dP_orifice[i]/PSI_TO_PA:.3f} psi')
        print(f'  Pipe pressure loss: {dP_pipe[i]/PSI_TO_PA:.3f} psi')
        print(f'  Orifice metering differential (estimated): {dP_meter[i]/PSI_TO_PA:.3f} psi')
        print(f'  Pressure after orifice recovery: {PipeInletP[i]/PSI_TO_PA:.3f} psia')
        print(f'  Pipe Reynolds number: {Reynolds[i]:.0f}')
        print(f'  Darcy friction factor: {Friction[i]:.5f}')
        print(f'  Estimated pipe outlet Mach (isothermal): {MachPipeOut[i]:.3f}')
        print(f'  Choked orifice (approximate): {choked[i]}')

if FlowSolver == 'Incompressible' and np.any(choked):
    print('\nWARNING: Choked liquid-stage recovered pressures are not uniquely predicted')
    print('by the simplified max(P_normal, P_choke) approximation.')
if FlowSolver == 'Compressible':
    if np.any(MachPipeOut > 0.3):
        print('\nWARNING: Pipe Mach > 0.3. The low-Mach isothermal pipe')
        print('approximation is not reliable; use a Fanno-flow or')
        print('coupled compressible momentum/energy model for this case.')
    if np.any(choked):
        print('\nWARNING: Gas pressure recovery at choked stages uses a')
        print('low-pressure-ratio empirical extrapolation. Validate')
        print('against compressible orifice data (e.g., xT/expansibility).')
