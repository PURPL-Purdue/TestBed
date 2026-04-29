# Maelstrom K-Bottle Runtime Calcs
# Authors: Chanak Gautam
# First Created: 4/27/2026
# Last Updated: 4/27/2026
# Calculations done in Imperial units

lbmol_to_mol = 453.59237

Z_ox = 0.9 # compressibility factor @ ~3000 psi
ox_molar_mass = 32 # lb/lb-mol
ambient_temp = 529.67 # Rankine (R)
R_oxygen = 10.73 # psi*ft^3/lb
k_bottle_vol = 1.76 # ft^3
k_bottle_num = 5
k_bottle_bulk = 2400 # psi
system_dp = 350 # psi
gox_mdot = 1.222  # lb/s, pulling from yaml crashes my computer lol
gox_feed = 650 # psi

minimum_bottle = (650 + 350) # psi
usable_dp = k_bottle_bulk - minimum_bottle # psi

# PV = nRT
# n = PV / RT
# mass = n * Molar Mass

gox_mol = ((usable_dp * k_bottle_vol) / (Z_ox * R_oxygen * ambient_temp)) # mol
gox_mass = gox_mol * ox_molar_mass # lb

runtime = gox_mass / gox_mdot # s

print("GOX SUPPLY CALCULATIONS")
print(f"USABLE GOX MASS PER BOTTLE (lb): {gox_mass:.2f}")
print(f"USABLE TOTAL GOX MASS (lb): {gox_mass * k_bottle_num:.2f}")
print(f"USEFUL RUNTIME  PER BOTTLE (s): {runtime:.2f}")
print(f"USEFUL TOTAL RUNTIME (s): {runtime * k_bottle_num:.2f}")





