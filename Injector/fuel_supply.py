# Maelstrom Fuel Supply Runtime Calcs
# Authors: Chanak Gautam
# First Created: 4/27/2026
# Last Updated: 4/27/2026
# Calculations done in Imperial units

ethanol_density = 6.586 # lb/gal
ipa_density = 6.55 # lb/gal
kero_density = 6.75 # lb/gal

# conservative max mass flow from test matrix
ethanol_mdot = 1.61 # lb/s
ipa_modt = 1.61 # lb/s
kero_mdot = 1.56 # lb/s

tank_capacity = 1.8 # gal

ethanol_mass = ethanol_density * tank_capacity # lb
ipa_mass = ipa_density * tank_capacity # lb
kero_mass = kero_density * tank_capacity # lb

ethanol_runtime = kero_mass / kero_mdot # s
ipa_runtime = ipa_mass / ipa_modt # s
kero_runtime = ethanol_mass / ethanol_mdot # s

print("FUEL SUPPLY CALCULATIONS (CONSERVATIVE MAX MASS FLOW)")
print(f"Ethanol Mass (lb): {ethanol_mass:.2}")
print(f"IPA Mass (lb): {ipa_mass:.2}")
print(f"Kerosense Mass (lb): {kero_mass:.2}")
print(f"Ethanol Runtime (s): {ethanol_runtime:.2}")
print(f"IPA Runtime (s): {ipa_runtime:.2}")
print(f"Kerosense Runtime (s): {kero_runtime:.2}")





