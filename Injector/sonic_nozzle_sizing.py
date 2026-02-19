# Maelstrom Sonic Nozzle Sizing Code
# Authors: Chanak Gautam
# First Created: 2/17/2026
# Last Updated: 2/17/2026
# Calculations done in SI units

import numpy as np
import matplotlib.pyplot as plt
from pyfluids import Fluid, FluidsList, Input
from pathlib import Path
from ruamel.yaml import YAML
import yaml

def find_yaml(filename="Maelstrom.yaml", start_dir=None):
    start_dir = Path(start_dir or Path.cwd())
    for path in start_dir.rglob(filename):
        return path
    raise FileNotFoundError(f"{filename} not found starting from {start_dir}")

yaml_path = find_yaml()

    # Load YAML with formatting preserved
yaml = YAML()
yaml.preserve_quotes = True

with open(yaml_path, "r") as f:
    data = yaml.load(f)

# ------------------------------------------------------------
# UNIT CONSTANTS
# ------------------------------------------------------------
psi_to_pa = 6894.76
in_to_m = 0.0254
lb_to_kg = 0.453592

# ------------------------------------------------------------
# FLUID + GEOMETRY PARAMETERS
# ------------------------------------------------------------

mdot = 0.8 * data["gox_design_mdot"] * lb_to_kg # 80% GOx mass flow rate (kg/s)
p_1 = data["gox_feed_pressure"] * psi_to_pa # Dowstream pressure relative to nozzle (Pa)
t_0 = 300 # Ambient temperature (K)

