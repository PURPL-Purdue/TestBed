"""Steady, serial pressure-ladder solver. SI units internally.
Gas: ideal, calorically perfect, adiabatic constant-area pipe with Fanno friction.
Each stage is an orifice followed by a straight pipe. No junctions/branches.
"""
from dataclasses import dataclass, asdict
from math import pi, sqrt, log, isfinite
from scipy.optimize import brentq


@dataclass
class Stage:
    orifice_d_mm: float = 12.7
    cd: float = 0.8
    pipe_d_mm: float = 25.4
    pipe_l_m: float = 0.5
    roughness_mm: float = 0.045
    fl: float = 0.85
    fstart: float = 0.65


@dataclass
class Config:
    fluid: str = 'Gas'  # Gas or Liquid
    inlet_psia: float = 200.0
    outlet_psia: float = 14.7
    temperature_k: float = 293.15
    gamma: float = 1.4
    gas_constant: float = 287.05
    gas_viscosity: float = 1.81e-5
    liquid_density: float = 998.0
    liquid_viscosity: float = 0.001
    vapor_psia: float = 0.34

PSI = 6894.757293168


def friction_factor(re, d, roughness):
    if re <= 0: return 0.0
    if re < 2300: return 64.0/re
    ft = (-1.8*log10_safe((roughness/(3.7*d))**1.11+6.9/re))**-2
    if re >= 4000: return ft
    w = (re-2300)/1700
    return (1-w)*64/re+w*ft


def log10_safe(x):
    from math import log10
    return log10(x)


def fanno_f(m, g):
    """Fanno parameter (Darcy f_D L/D remaining to sonic)."""
    return (1-m*m)/(g*m*m)+(g+1)/(2*g)*log(((g+1)*m*m)/(2+(g-1)*m*m))


def p0_over_p(m, g):
    return (1+(g-1)*m*m/2)**(g/(g-1))


def static_pressure_from_mdot(mdot, area, t0, R, gamma, p):
    """Return subsonic Mach for static pressure p, or 1 if pipe exit chokes."""
    G=mdot/area
    # G = p M sqrt(gamma/(R*T0)) sqrt(1 + (gamma-1)M²/2)
    k=(gamma-1)/2
    q=G/(p*sqrt(gamma/(R*t0)))
    if q>=sqrt(1+k): return 1.0, True
    return sqrt(2*q*q/(1+sqrt(1+4*k*q*q))), False


def pipe_reverse(mdot, p_requested, s, c):
    """Gas pipe: back-calculate inlet from outlet using Fanno; use sonic exit floor.
    Report effective outlet static pressure if requested backpressure is too low.
    """
    d=s.pipe_d_mm*1e-3
    L=s.pipe_l_m
    area=pi*d*d/4
    g=c.gamma
    R=c.gas_constant
    t0=c.temperature_k
    G=mdot/area
    re=G*d/c.gas_viscosity
    f=friction_factor(re,d,s.roughness_mm*1e-3)
    mach_out, choke=static_pressure_from_mdot(mdot,area,t0,R,g,p_requested)
    p_exit=G/(mach_out*sqrt(g/(R*t0))*sqrt(1+(g-1)*mach_out**2/2))
    loss=f*L/d
    target=fanno_f(mach_out,g)+loss
    mach_in=brentq(lambda m:fanno_f(m,g)-target,max(1e-16,mach_out*1e-6),mach_out) if loss>1e-13 else mach_out
    # For constant T0, mass flux implies p proportional to 1/(M sqrt(1+aM²)).
    p_in=G/(mach_in*sqrt(g/(R*t0))*sqrt(1+(g-1)*mach_in**2/2))
    return {'p_in':p_in,'p_out':p_exit,'mach_in':mach_in,'mach_out':mach_out,
            'fanno_capacity':fanno_f(mach_in,g),'fanno_used':loss,'pipe_choked':choke,
            're':re,'friction':f,'pipe_dp':p_in-p_exit,
            'p0_in':p_in*p0_over_p(mach_in,g),'p0_out':p_exit*p0_over_p(mach_out,g)}


def recovery_fraction(beta, cd):
    x=sqrt(1-beta**4*(1-cd*cd))
    return (x-cd*beta*beta)/(x+cd*beta*beta)


def gas_orifice_reverse(mdot,p_recovered,s,c):
    """Approximate inverse orifice with equivalent jet pressure and recovered static p.
    Ideal-gas formula, empirical recovery fraction; not a validated high-Mach junction.
    """
    area=pi*(s.orifice_d_mm*1e-3)**2/4
    beta=s.orifice_d_mm/s.pipe_d_mm
    alpha=recovery_fraction(beta,s.cd)
    g=c.gamma; R=c.gas_constant; T=c.temperature_k
    rcrit=(2/(g+1))**(g/(g-1))
    critical_coeff=sqrt(g/(R*T))*(2/(g+1))**((g+1)/(2*(g-1)))
    def rate(up):
        if up<=p_recovered: return 0.
        pvirtual=max(0.,up-(up-p_recovered)/alpha)
        r=max(rcrit,pvirtual/up)
        return s.cd*area*up*sqrt(max(0.,2*g/(R*T*(g-1))*(r**(2/g)-r**((g+1)/g))))
    upper=max(p_recovered*1.05,mdot/(s.cd*area*critical_coeff)*1.1,1.)
    for _ in range(100):
        if rate(upper)>=mdot: break
        upper*=2
    else: raise ValueError('Could not bracket orifice pressure')
    up=brentq(lambda p:rate(p)-mdot,p_recovered,upper,xtol=1e-7)
    virtual=max(0.,up-(up-p_recovered)/alpha)
    return {'p_in':up,'p_out':p_recovered,'orifice_dp':up-p_recovered,
            'meter_dp':(up-p_recovered)/alpha,'orifice_choked':virtual/up<=rcrit*(1+1e-8),
            'recovery_fraction':alpha}


def liquid_stage_reverse(mdot,p_down,s,c):
    d=s.pipe_d_mm*1e-3; a=pi*d*d/4
    rho=c.liquid_density; mu=c.liquid_viscosity
    v=mdot/(rho*a)
    re=rho*v*d/mu
    f=friction_factor(re,d,s.roughness_mm*1e-3)
    pipe_dp=f*s.pipe_l_m/d*rho*v*v/2
    pipe_in=p_down+pipe_dp
    ao=pi*(s.orifice_d_mm*1e-3)**2/4
    beta=s.orifice_d_mm/s.pipe_d_mm
    K=s.cd*ao*sqrt(2*rho/(1-beta**4))
    meter=(mdot/K)**2
    perm=recovery_fraction(beta,s.cd)*meter
    normal=pipe_in+perm
    choke_required=c.vapor_psia*PSI+meter/s.fl**2
    up=max(normal,choke_required)
    onset=meter>=s.fstart**2*(up-c.vapor_psia*PSI)
    return {'p_in':up,'p_out':p_down,'pipe_in':pipe_in,'pipe_dp':pipe_dp,
            'orifice_dp':up-pipe_in,'meter_dp':meter,'orifice_choked':choke_required>normal+1e-7,
            'pipe_choked':False,'cavitation_onset':onset,'re':re,'friction':f}


def validate(c, stages):
    if not stages: raise ValueError('Add at least one stage.')
    if len(stages)>30: raise ValueError('A maximum of 30 serial stages is supported.')
    if c.fluid not in ('Gas','Liquid'): raise ValueError('Fluid must be Gas or Liquid.')
    if not (c.inlet_psia>c.outlet_psia>0): raise ValueError('Inlet pressure must exceed positive outlet pressure (absolute).')
    if c.temperature_k<=0: raise ValueError('Absolute temperature must be positive.')
    if c.fluid=='Gas' and not (c.gamma>1 and c.gas_constant>0 and c.gas_viscosity>0):
        raise ValueError('Invalid gas properties.')
    if c.fluid=='Liquid' and not (c.liquid_density>0 and c.liquid_viscosity>0 and 0<=c.vapor_psia<c.inlet_psia):
        raise ValueError('Invalid liquid properties.')
    for i,s in enumerate(stages,1):
        if not (0<s.orifice_d_mm<s.pipe_d_mm and 0<s.cd<=1 and s.pipe_l_m>=0 and 0<=s.roughness_mm< s.pipe_d_mm/2):
            raise ValueError(f'Invalid geometry or Cd in stage {i}.')
        if c.fluid=='Liquid' and not (0<s.fstart<=s.fl<=1):
            raise ValueError(f'Invalid cavitation factors in stage {i}.')


def evaluate(mdot,c,stages,details=False):
    p=c.outlet_psia*PSI
    rev=[]
    for i in range(len(stages)-1,-1,-1):
        s=stages[i]
        if c.fluid=='Gas':
            pipe=pipe_reverse(mdot,p,s,c)
            ori=gas_orifice_reverse(mdot,pipe['p_in'],s,c)
            stage={'stage':i+1,**pipe,**{k:v for k,v in ori.items() if k not in ('p_out',)},
                   'pipe_in':pipe['p_in'],'pipe_out':pipe['p_out']}
            p=ori['p_in']
        else:
            stage={'stage':i+1,**liquid_stage_reverse(mdot,p,s,c)}
            p=stage['p_in']
        rev.append(stage)
    return (p,list(reversed(rev))) if details else p


def solve(c,stages):
    validate(c,stages)
    target=c.inlet_psia*PSI
    # Upper bound ensures a finite root without assuming restriction controls.
    hi=1e-8
    for _ in range(95):
        if evaluate(hi,c,stages)>=target: break
        hi*=2
    else: raise ValueError('Unable to bracket mass flow; check geometry and properties.')
    mdot=brentq(lambda m:evaluate(m,c,stages)-target,1e-13,hi,xtol=1e-12)
    p_in, rows=evaluate(mdot,c,stages,True)
    node=[r['p_in']/PSI for r in rows]+[(rows[-1]['pipe_out'] if c.fluid=='Gas' else c.outlet_psia*PSI)/PSI]
    for r in rows:
        for key in ('p_in','p_out','pipe_in','pipe_out','orifice_dp','pipe_dp','meter_dp'):
            if key in r: r[key+'_psia_or_psi']=r[key]/PSI
    alerts=[]
    for r in rows:
        if r['pipe_choked']: alerts.append(f"Stage {r['stage']}: Fanno sonic limit at pipe outlet (Mach 1). Downstream pressure cannot influence upstream through that segment under this approximation.")
        if r['orifice_choked']: alerts.append(f"Stage {r['stage']}: orifice sonic/cavitation flow limitation.")
    if c.fluid=='Gas':
        alerts.append('Engineering screening only: orifice recovery uses an empirical incompressible-style fraction; junction total-pressure changes, shocks, heat transfer, and coupled multi-choke physics are not resolved.')
    else:
        alerts.append('Cavitation model is a screening approximation, not a two-phase flow solver.')
    return {'mass_flow_kg_s':mdot,'mass_flow_lbm_s':mdot*2.20462262185,
            'nodes_psia':node,'stages':rows,'alerts':alerts,
            'pipe_choke_points':[r['stage'] for r in rows if r['pipe_choked']],
            'orifice_choke_points':[r['stage'] for r in rows if r['orifice_choked']],
            'inlet_psia_calculated':p_in/PSI,'specified_backpressure_psia':c.outlet_psia}
