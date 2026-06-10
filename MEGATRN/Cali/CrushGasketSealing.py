import math

# =====================================================
# === INPUT PARAMETERS ================================
# =====================================================
# --- Pressure Zone Configurations ---
P_core = 200.0                  # Internal core pressure (psi)
P_coolant = 300.0               # Inter-gasket coolant jacket cavity pressure (psi)

# --- Inner Gasket Properties (Solid Flat Copper) ---
G_mean_inner = 4.26             # Mean diameter of the inner gasket (inches)
N_over_2_inner = 0.03           # Seating width parameter N/2 for inner gasket (inches)
m_inner = 4.75                  
y_inner = 13000.0               

# --- Outer Gasket Properties (Solid Flat Copper) ---
G_mean_outer = 4.86             # Mean diameter of the outer gasket (inches)
N_over_2_outer = 0.03           # Seating width parameter N/2 for outer gasket (inches)
m_outer = 4.75                  
y_outer = 13000.0               

# --- Bolt & Flange Geometry ---
bolt_diam = 0.5                 # Nominal rod diameter (inches)
bolt_yield_psi = 125000         # Rod yield strength (psi)
flange_thick_mm = 30.0          # Flange axial depth/thickness (mm)
estimated_bolt_circle_diam = 5.90551 # Adjusted Bolt Circle Diameter to clear outer gasket (inches)


# =====================================================
# === MULTI-ZONE PRESSURE FLANGE CALCULATIONS ========
# =====================================================
flange_thick_in = flange_thick_mm / 25.4

# 1. Width Profiles - Inner Gasket
b0_in = N_over_2_inner
b_in = b0_in if b0_in <= 0.25 else 0.5 * math.sqrt(b0_in)
G_in = G_mean_inner if b0_in <= 0.25 else (G_mean_inner + b0_in) - (2.0 * b_in)

# 2. Width Profiles - Outer Gasket
b0_out = N_over_2_outer
b_out = b0_out if b0_out <= 0.25 else 0.5 * math.sqrt(b0_out)
G_out = G_mean_outer if b0_out <= 0.25 else (G_mean_outer + b0_out) - (2.0 * b_out)

# 3. Summed Hydrostatic End Forces (H)
# Inner zone pressure acts on the area inside the inner gasket
H_core = (math.pi * (G_in ** 2) / 4.0) * P_core
# Inter-gasket pressure acts on the annular area between the two reaction diameters
H_annulus = (math.pi * ((G_out ** 2) - (G_in ** 2)) / 4.0) * P_coolant
H_total = H_core + H_annulus

# 4. Summed Joint Compression Loads (Hp)
# Operating compression for both gaskets reacting to their surrounding fluid environments
Hp_inner = (2.0 * b_in * math.pi * G_in) * m_inner * abs(P_coolant - P_core)
Hp_outer = (2.0 * b_out * math.pi * G_out) * m_outer * P_coolant
Hp_total = Hp_inner + Hp_outer

# 5. Combined Required Bolt Operating Load (Wm1)
Wm1 = H_total + Hp_total

# 6. Combined Required Initial Gasket Seating Load (Wm2)
# Both copper rings must be crushed simultaneously during dry assembly
Wm2_inner = math.pi * b_in * G_in * y_inner
Wm2_outer = math.pi * b_out * G_out * y_outer
Wm2_total = Wm2_inner + Wm2_outer

# Controlling structural threshold
W_controlling = max(Wm1, Wm2_total)

# 7. Bolt Sizing & Count Assessment
bolt_allowable_stress = bolt_yield_psi / 4.0
required_bolt_area_Am = W_controlling / bolt_allowable_stress
single_bolt_area = (math.pi * (bolt_diam ** 2)) / 4.0
min_bolts_by_strength = math.ceil(required_bolt_area_Am / single_bolt_area)

# 8. Maximum Circumferential Bolt Spacing Criteria
# Controlled by the gasket material parameter (using outermost gasket value for spacing conservative check)
Bs_max = (2.0 * bolt_diam) + ((6.0 * flange_thick_in) / (m_outer + 0.5))

# Perimeter distribution mapping
bolt_circle_perimeter = math.pi * estimated_bolt_circle_diam
min_bolts_by_spacing = math.ceil(bolt_circle_perimeter / Bs_max)

# Force combined criteria constraints
final_bolt_count = max(min_bolts_by_strength, min_bolts_by_spacing)
if final_bolt_count < 4:
    final_bolt_count = 4
if final_bolt_count % 2 != 0:
    final_bolt_count += 1

actual_bolt_spacing = bolt_circle_perimeter / final_bolt_count


# =====================================================
# === OUTPUTS =========================================
# =====================================================
print("==================================================")
print("     ASME VIII DUAL-PRESSURE FLANGE REPORT       ")
print("==================================================")
print(f"Core Fluid Pressure           : {P_core:.1f} psi")
print(f"Coolant Cavity Pressure       : {P_coolant:.1f} psi")
print(f"Flange Thickness              : {flange_thick_mm:.1f} mm ({flange_thick_in:.4f} in)")
print("--------------------------------------------------")
print(f"Hydrostatic End Force (Core)  : {H_core:.2f} lbf")
print(f"Hydrostatic End Force (Cavity): {H_annulus:.2f} lbf")
print(f"Total Hydrostatic Force (H)   : {H_total:.2f} lbf")
print("--------------------------------------------------")
print(f"Operating Bolt Load (Wm1)     : {Wm1:.2f} lbf")
print(f"Gasket Seating Load (Wm2)     : {Wm2_total:.2f} lbf")
print(f"Controlling Structural Load   : {W_controlling:.2f} lbf")
print("--------------------------------------------------")
print(f"Required Total Bolt Area (Am) : {required_bolt_area_Am:.4f} sq.in")
print(f"Min Bolts Required by Strength: {min_bolts_by_strength}")
print("--------------------------------------------------")
print(f"MAX ALLOWED BOLT SPACING (Bs) : {Bs_max:.4f} inches")
print(f"Design Bolt Circle Diam (C)   : {estimated_bolt_circle_diam:.3f} inches")
print("--------------------------------------------------")
print(f"FINAL REQUIRED BOLT COUNT     : {final_bolt_count} bolts")
print(f"Actual Spacing with this Count: {actual_bolt_spacing:.4f} inches")
print("==================================================")