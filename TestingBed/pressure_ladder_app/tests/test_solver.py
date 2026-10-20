import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import pytest
from math import isclose
from solver import Config, Stage, solve, pipe_reverse, fanno_f, friction_factor, validate


def test_fanno_identity():
    assert fanno_f(1,1.4)==pytest.approx(0)
    assert fanno_f(0.2,1.4)>fanno_f(0.8,1.4)>0


def test_fanno_length_relation():
    c=Config()
    s=Stage(pipe_l_m=10)
    r=pipe_reverse(0.015,100000,s,c)
    assert fanno_f(r['mach_in'],c.gamma)-fanno_f(r['mach_out'],c.gamma)==pytest.approx(r['friction']*s.pipe_l_m/(s.pipe_d_mm/1000),rel=1e-7)
    assert r['p_in']>=r['p_out']


def test_gas_outputs_and_pressure():
    c=Config()
    result=solve(c,[Stage(),Stage()])
    assert result['mass_flow_kg_s']>0
    assert result['nodes_psia'][0]==pytest.approx(c.inlet_psia,rel=1e-6)
    assert all(r['mach_out']<=1 for r in result['stages'])


def test_geometry_affects_flow():
    c=Config()
    short=solve(c,[Stage(pipe_l_m=0.1)])['mass_flow_kg_s']
    long=solve(c,[Stage(pipe_l_m=30)])['mass_flow_kg_s']
    assert long<short
    small=solve(c,[Stage(pipe_d_mm=18,pipe_l_m=15)])['mass_flow_kg_s']
    big=solve(c,[Stage(pipe_d_mm=30,pipe_l_m=15)])['mass_flow_kg_s']
    assert small<big


def test_pipe_sonic_detection():
    r=pipe_reverse(0.02,1000,Stage(pipe_l_m=20),Config())
    assert r['pipe_choked']
    assert r['mach_out']==pytest.approx(1)
    assert r['mach_in']<1


def test_liquid_and_cavitation():
    c=Config(fluid='Liquid')
    r=solve(c,[Stage()])
    assert r['mass_flow_kg_s']>0
    assert r['nodes_psia'][0]==pytest.approx(c.inlet_psia,rel=1e-6)
    assert 'cavitation_onset' in r['stages'][0]


def test_reject_invalid():
    with pytest.raises(ValueError): solve(Config(),[Stage(pipe_d_mm=5)])
    with pytest.raises(ValueError): solve(Config(inlet_psia=3,outlet_psia=10),[Stage()])


def test_choked_exit_pressure_floor():
    r=solve(Config(),[Stage()])
    assert r['nodes_psia'][-1]>=r['specified_backpressure_psia']
