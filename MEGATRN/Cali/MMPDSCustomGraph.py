import matplotlib.pyplot as plt

# Create figure and axis
fig, ax = plt.subplots(figsize=(8, 5))

# Set axis limits
ax.set_xlim(50, 700)
ax.set_ylim(0, 110)

# Axis labels
ax.set_xlabel('Temperature (°F)', fontsize=15)
ax.set_ylabel('Percentage of Room Temperature $F_{ty}$', fontsize=15)


ax.set_title('Percentage of Room Temperature vs Temperature on Tensile Yield', fontsize=15)

# Grid formatting
ax.grid(True, linestyle='--', alpha=0.6)

# Optional tick spacing
ax.set_xticks(range(50, 701, 50))
ax.set_yticks(range(0, 111, 10))

plt.tight_layout()
plt.show()