"""Streamlit frontend: run with `streamlit run app.py`."""
import json
import pandas as pd
import streamlit as st
from dataclasses import asdict
from solver import Config, Stage, solve

st.set_page_config(page_title='Pressure Ladder Studio', page_icon='🚀', layout='wide')
st.title('Pressure Ladder Studio')
st.caption('Serial flow system designer · Fanno gas pipes · liquid cavitation screening · SI-based calculation')

if 'stages' not in st.session_state:
    st.session_state.stages=[asdict(Stage()) for _ in range(3)]
if 'config' not in st.session_state:
    st.session_state.config=asdict(Config())

with st.sidebar:
    st.header('1 · Fluid & boundary conditions')
    c=st.session_state.config
    fluid=st.radio('Flow regime', ['Gas','Liquid'], index=['Gas','Liquid'].index(c['fluid']))
    c['fluid']=fluid
    c['inlet_psia']=st.number_input('Supply pressure · psia',min_value=0.01,max_value=100000.0,value=float(c['inlet_psia']),step=10.0)
    c['outlet_psia']=st.number_input('Backpressure · psia',min_value=0.01,max_value=100000.0,value=float(c['outlet_psia']),step=1.0)
    c['temperature_k']=st.number_input('Stagnation temperature · K',min_value=1.0,max_value=5000.0,value=float(c['temperature_k']),step=5.0)
    if fluid=='Gas':
        c['gamma']=st.slider('Specific heat ratio γ',1.01,1.8,float(c['gamma']),0.01)
        c['gas_constant']=st.number_input('Specific gas constant · J/(kg·K)', min_value=0.01,max_value=10000.0,value=float(c['gas_constant']))
        c['gas_viscosity']=st.number_input('Gas dynamic viscosity · Pa·s',min_value=1e-9,max_value=1.0,value=float(c['gas_viscosity']),format='%.8f')
    else:
        c['liquid_density']=st.number_input('Liquid density · kg/m³',min_value=0.01,max_value=30000.0,value=float(c['liquid_density']))
        c['liquid_viscosity']=st.number_input('Liquid dynamic viscosity · Pa·s',min_value=1e-9,max_value=100.0,value=float(c['liquid_viscosity']),format='%.6f')
        c['vapor_psia']=st.number_input('Liquid vapor pressure · psia',min_value=0.0,max_value=100000.0,value=float(c['vapor_psia']))

st.subheader('2 · Build your serial flow path')
st.write('Every row is **an orifice followed by a straight pipe**. Add stages in physical flow order. Double-click a cell to enter precise values.')
left, right=st.columns([1,4])
with left:
    if st.button('＋ Add stage',use_container_width=True):
        st.session_state.stages.append(asdict(Stage()))
        st.session_state.pop('geometry_editor',None)
        st.session_state.pop('result',None)
        st.rerun()
    if st.button('－ Remove last stage',disabled=len(st.session_state.stages)<=1,use_container_width=True):
        st.session_state.stages.pop()
        st.session_state.pop('geometry_editor',None)
        st.session_state.pop('result',None)
        st.rerun()
with right:
    st.caption('Lengths, diameters, roughness, and discharge coefficients are all independent between stages.')

cols={
 'orifice_d_mm':st.column_config.NumberColumn('Orifice Ø (mm)',min_value=0.001,format='%.3f'),
 'cd':st.column_config.NumberColumn('Cd',min_value=0.001,max_value=1.0,format='%.3f'),
 'pipe_d_mm':st.column_config.NumberColumn('Pipe ID (mm)',min_value=0.001,format='%.3f'),
 'pipe_l_m':st.column_config.NumberColumn('Pipe length (m)',min_value=0.0,format='%.3f'),
 'roughness_mm':st.column_config.NumberColumn('Roughness (mm)',min_value=0.0,format='%.5f'),
 'fl':st.column_config.NumberColumn('Liquid FL',min_value=0.001,max_value=1.0,format='%.3f'),
 'fstart':st.column_config.NumberColumn('Liquid Fstart',min_value=0.001,max_value=1.0,format='%.3f')}
keys=['orifice_d_mm','cd','pipe_d_mm','pipe_l_m','roughness_mm']+(['fl','fstart'] if fluid=='Liquid' else [])
df=pd.DataFrame(st.session_state.stages)[keys]
edited=st.data_editor(df,hide_index=False,num_rows='fixed',use_container_width=True,column_config={k:cols[k] for k in keys},key='geometry_editor')
for i,record in enumerate(edited.to_dict(orient='records')):
    st.session_state.stages[i].update(record)

with st.expander('Alternative quick adjustments · sliders'):
    stage_idx=st.selectbox('Choose stage', list(range(len(st.session_state.stages))),format_func=lambda i:f'Stage {i+1}')
    current=st.session_state.stages[stage_idx]
    # Buttons make updates deterministic relative to the editable geometry grid.
    diameter=st.slider('Set pipe ID · mm',0.5,200.0,float(min(max(current['pipe_d_mm'],0.5),200.0)),0.5)
    length=st.slider('Set pipe length · m',0.0,100.0,float(min(max(current['pipe_l_m'],0),100)),0.1)
    discharge=st.slider('Set orifice Cd',0.1,1.0,float(min(max(current['cd'],0.1),1.0)),0.01)
    if st.button('Apply sliders to selected stage'):
        current.update(pipe_d_mm=diameter,pipe_l_m=length,cd=discharge)
        st.session_state.pop('geometry_editor',None)
        st.session_state.pop('result',None)
        st.rerun()

st.subheader('3 · Solve and inspect')
if st.button('▶ Run pressure ladder',type='primary',use_container_width=True):
    try:
        result=solve(Config(**c),[Stage(**s) for s in st.session_state.stages])
        st.session_state.result=result
    except (ValueError,OverflowError,ZeroDivisionError) as exc:
        st.session_state.result=None
        st.error(str(exc))

result=st.session_state.get('result')
if result:
    a,b,d=st.columns(3)
    a.metric('Mass flow',f"{result['mass_flow_kg_s']:.6g} kg/s")
    b.metric('Orifice choke candidates',str(len(result['orifice_choke_points'])))
    d.metric('Fanno pipe choke candidates',str(len(result['pipe_choke_points'])))
    st.write('**Flow-path static pressures (psia)**')
    if result['nodes_psia'][-1]>result['specified_backpressure_psia']+0.001:
        st.info(f"Choked outlet: pipe exit static pressure {result['nodes_psia'][-1]:.3f} psia exceeds specified backpressure {result['specified_backpressure_psia']:.3f} psia.")
    st.line_chart(pd.DataFrame({'Pressure (psia)':result['nodes_psia']}).rename_axis('Node index'))
    table=pd.DataFrame(result['stages'])
    labels={'stage':'Stage','p_in_psia_or_psi':'Upstream psia', 'pipe_in_psia_or_psi':'Pipe inlet psia',
            'pipe_out_psia_or_psi':'Pipe outlet psia', 'orifice_dp_psia_or_psi':'Orifice Δp psi',
            'pipe_dp_psia_or_psi':'Pipe Δp psi','mach_in':'Pipe inlet Mach','mach_out':'Pipe outlet Mach',
            'fanno_used':'fD L/D','re':'Re','friction':'Darcy f','pipe_choked':'Pipe choked',
            'orifice_choked':'Orifice choked','cavitation_onset':'Cavitation onset'}
    present=[k for k in labels if k in table]
    st.dataframe(table[present].rename(columns=labels),hide_index=True,use_container_width=True)
    for alert in result['alerts']:
        st.warning(alert)
    st.download_button('Download results · JSON',json.dumps(result,indent=2).encode('utf-8'),file_name='pressure_ladder_results.json',mime='application/json')

st.divider()
config_file={'config':c,'stages':st.session_state.stages}
st.download_button('Save system configuration · JSON',json.dumps(config_file,indent=2),file_name='pressure_ladder_configuration.json',mime='application/json')
config_upload=st.file_uploader('Load saved configuration (JSON)',type=['json'])
if config_upload and st.button('Import configuration'):
    try:
        data=json.load(config_upload)
        config=Config(**data['config'])
        stages=[Stage(**s) for s in data['stages']]
        from solver import validate
        validate(config,stages)
        st.session_state.config=asdict(config)
        st.session_state.stages=[asdict(s) for s in stages]
        st.session_state.pop('result',None)
        st.session_state.pop('geometry_editor',None)
        st.rerun()
    except (ValueError,KeyError,TypeError) as e:
        st.error(f'Invalid configuration: {e}')

st.caption('Assumptions: steady state; serial stages; gas ideal and adiabatic; constant-area Fanno pipe with Darcy friction; recovery correlations are approximate. No branched networks, shocks, transient blowdown, or validated two-phase cavitation.')
