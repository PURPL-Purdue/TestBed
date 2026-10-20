# Pressure Ladder Studio

Browser-based, steady-state serial fluid-flow pressure ladder simulator. Supports an unlimited-in-the-UI number of stages up to a safety limit of 30, each with an orifice and downstream pipe.

## Start locally

```bash
python -m venv .venv
# Windows PowerShell:
.venv\Scripts\Activate.ps1
# macOS/Linux: source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

Streamlit displays a localhost address. To publish under a shared URL, push this folder to GitHub and select `app.py` on Streamlit Community Cloud. No secrets needed.

## Input / output
- Supply and backpressure: **absolute psia** (not psig).
- Temperature: K (gas stagnation temperature assumed constant).
- Pipe IDs, orifice diameters, and roughness: mm. Pipe lengths: m.
- Gas: gamma, specific R (J/kg/K), dynamic viscosity (Pa s).
- Liquid: density (kg/m³), viscosity (Pa s), vapor pressure (psia), empirical FL and Fstart.
- Each row: restriction, then pipe; no tees, manifolds or branches.
- Output includes node pressures, mass flow, per-stage Darcy factors, Mach, Fanno f_D L/D, choking indicators, and JSON export.

## Physics and limitations
The gas pipe uses the **subsonic Fanno equation**, with Darcy friction factor:

```
F(M) = (1-M^2)/(gamma M^2) + (gamma+1)/(2 gamma)
       * ln[(gamma+1)M^2/(2+(gamma-1)M^2)]
F(M_in) - F(M_out) = f_D L/D
```

It computes local gas Reynolds number from `mdot*D/(A*mu)` and a laminar / Haaland turbulent Darcy factor. For a specified trial mass flow, downstream static pressure determines outlet Mach, clipped at M=1 when backpressure is below the sonic exit requirement; Fanno is then integrated upstream. A pipe can therefore be flagged as an **additional choke candidate** independently of the orifices. Friction is assumed constant along each pipe (constant mass flux, viscosity and geometry). Each gas restriction uses an equivalent-jet isentropic mass-flow model plus empirical permanent-loss fraction inherited from the original implementation.

**Important engineering limitations:** This is a preliminary steady-state screening tool, not a rigorous network solver or certified design calculator. At choking, the simplified backward-pressure coupling (especially if there are several separated chokes) is not a full coupled compressible network solution. It does not resolve shocks, non-isentropic mixing, nozzle/vena-contracta geometry, full stagnation-pressure bookkeeping through junctions, heat transfer, real-fluid EOS, liquid flashing, downstream pressure recovery interactions, branching or time-dependent tank blowdown. In particular the gas pressure-recovery factor is adapted from an incompressible relationship and is **not validated for high pressure-ratio gas restrictions**. A pipe with nonzero length can have its sonic limit at the exit; the tool does not model normal shocks or supersonic pipe flow. Sonic detection at a pipe exit can be saturated at backpressures below critical. Do not rely on it alone for flight hardware, high-pressure safety decisions, or cavitation-critical systems. Validate against test data or higher-fidelity CFD.

## Files
- `app.py`: Streamlit GUI
- `solver.py`: importable physics model (no Streamlit dependency)
- `tests/test_solver.py`: unit tests including Fanno inverse identity, geometry effects, and choke flags
- `legacy_reference.py`: archived original simulator from Oct 14, 2026
- `.streamlit/config.toml`: visual theme
- `sample_system.json`: example UI configuration

Run `pytest -q` for numerical tests, and `streamlit run app.py` to inspect UI.
