import numpy as np
import matplotlib.pyplot as plt

# Constants
gamma = 1.4
Cd = 0.8
mdot = 0.1       # kg/s
p_c = 100        # chamber pressure
rho = 1.0        # replace with real density if needed

# Pressure range
p_feed_values = np.arange(20, 500, 1)

# Critical ratio
p_cr = (2 / (gamma + 1)) ** (gamma / (gamma - 1))

# Storage
choked_area  = []
unchoked_area = []
correct_area = []

for p_feed in p_feed_values:

    # --- CHOKED formula ---
    A_choked = mdot / (
        Cd * np.sqrt(
            gamma * rho * p_feed *
            (2 / (gamma + 1)) ** ((gamma + 1)/(gamma - 1))
        )
    )

    # --- UNCHOKED formula ---
    A_unchoked = mdot / (
        Cd * np.sqrt(
            2 * rho * p_feed *
            (gamma/(gamma-1)) *
            ((p_c/p_feed)**(2/gamma) - (p_c/p_feed)**((gamma + 1)/gamma))
        )
    )

    choked_area.append(A_choked)
    unchoked_area.append(A_unchoked)

    # --- Physically correct selection ---
    if p_c / p_feed < p_cr:      # CHOKED region
        correct_area.append(A_choked)
    else:                        # UNCHOKED region
        correct_area.append(A_unchoked)

# Convert to arrays
choked_area  = np.array(choked_area)
unchoked_area = np.array(unchoked_area)
correct_area = np.array(correct_area)

# -------------------- SIDE-BY-SIDE PLOTS --------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14,5))

# ---------- LEFT PLOT: full comparison ----------
ax1.plot(p_feed_values, choked_area,   label="Choked Formula",   color="red")
ax1.plot(p_feed_values, unchoked_area, label="Unchoked Formula", color="blue")
ax1.axvline(x=p_c / p_cr, color='black', linestyle='--', label="p_feed = p_c / p_cr")
ax1.set_xlabel("Feed Pressure p_feed")
ax1.set_ylabel("Injection Area [m²]")
ax1.set_title("Choked vs Unchoked (Full Range)")
ax1.grid(True)
ax1.legend()

# ---------- RIGHT PLOT: correct combined curve ----------
ax2.plot(p_feed_values, correct_area, color="green", label="Correct Flow Regime")
ax2.axvline(x=p_c / p_cr, color='black', linestyle='--', label="Choke Transition")
ax2.set_xlabel("Feed Pressure p_feed")
ax2.set_ylabel("Injection Area [m²]")
ax2.set_title("Correct Injection Area")
ax2.grid(True)
ax2.legend()

plt.tight_layout()
plt.show()
