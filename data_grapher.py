import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import glob
import os

# ────────────────────────────────────────────────
# Load all matching CSVs in the current folder
# ────────────────────────────────────────────────
file_paths = sorted(glob.glob("torch_hot_fire_*_high_*.csv"))
dataframes = []
test_names = []

for path in file_paths:
    df = pd.read_csv(path)
    df['Elapsed (ms)'] = df.index
    test_name = os.path.splitext(os.path.basename(path))[0].replace("_", " ").title()
    dataframes.append(df)
    test_names.append(test_name)

# ────────────────────────────────────────────────
# Define event sequence
# ────────────────────────────────────────────────
events_default = {
    'SN-O2-01': [1000, 4500],
    'Spark Plug': [2000, 3500],
    'SN-H2-01': [3000, 4500],
    'SN-N2-01': [4600, 9000]
}
events_first = {
    'SN-O2-01': [1000, 4000],
    'Spark Plug': [2000, 3500],
    'SN-H2-01': [3000, 4000],
    'SN-N2-01': [4000, 9000]
}

# ────────────────────────────────────────────────
# Prep function
# ────────────────────────────────────────────────
def prepare_df(df, threshold=10):
    active = df[df['PT-TI-01 Pressure'] > threshold]
    if active.empty:
        return df.copy(), 0, len(df)
    t_min = max(active['Elapsed (ms)'].min() - 500, 0)
    t_max = min(active['Elapsed (ms)'].max() + 500, df['Elapsed (ms)'].max())
    return df[(df['Elapsed (ms)'] >= t_min) & (df['Elapsed (ms)'] <= t_max)].copy(), t_min, t_max

def scale_events(event_dict, test_duration_ms, shift_ms=100):
    scale = test_duration_ms / 10000
    return {
        dev: [max(t * scale - shift_ms, 0) for t in times]
        for dev, times in event_dict.items()
    }

# ────────────────────────────────────────────────
# Process all datasets
# ────────────────────────────────────────────────
processed = []
height_ratios = []

for i, df in enumerate(dataframes):
    dfc, tmin, tmax = prepare_df(df)
    duration = tmax - tmin
    events = scale_events(events_first if i == 0 else events_default, duration)
    height_ratios.append(duration)
    processed.append((dfc, tmin, tmax, events, test_names[i]))

# Shared formatter
formatter = ticker.FuncFormatter(lambda x, pos: f"{int(x * 10)}")

# ────────────────────────────────────────────────
# Figure 1 – PT-TI-01 + events
# ────────────────────────────────────────────────
fig1, axs1 = plt.subplots(len(processed), 1, figsize=(14, 3 * len(processed)),
                          gridspec_kw={'height_ratios': height_ratios})

if len(processed) == 1:
    axs1 = [axs1]

for ax, (dfc, tmin, tmax, seq, title) in zip(axs1, processed):
    dfc['t_plot'] = dfc['Elapsed (ms)'] - tmin
    ax.plot(dfc['t_plot'], dfc['PT-TI-01 Pressure'], color='steelblue')
    ymax = dfc['PT-TI-01 Pressure'].max()
    ybase = ymax + 0.02 * ymax
    ystep = -0.03 * ymax
    for device, times in seq.items():
        for idx, t in enumerate(times):
            ax.axvline(t, linestyle='--', color='gray')
            ax.text(t + 5, ybase + idx * ystep, device, rotation=90,
                    verticalalignment='bottom', fontsize=8)
    ax.set_xlim(0, tmax - tmin)
    ax.xaxis.set_major_formatter(formatter)
    ax.set_title(title)
    ax.set_xlabel("Elapsed Time (ms) ×10")
    ax.set_ylabel("PT-TI-01 (psi)")
    ax.grid(True)

fig1.suptitle("Figure 1 – PT-TI-01 Chamber Pressure with Event Markers", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.96])

# ────────────────────────────────────────────────
# Figure 2 – Chamber + Feed Pressures
# ────────────────────────────────────────────────
fig2, axs2 = plt.subplots(len(processed), 1, figsize=(14, 3 * len(processed)),
                          gridspec_kw={'height_ratios': height_ratios})

if len(processed) == 1:
    axs2 = [axs2]

for ax, (dfc, tmin, tmax, _, title) in zip(axs2, processed):
    dfc['t_plot'] = dfc['Elapsed (ms)'] - tmin
    ax.plot(dfc['t_plot'], dfc['PT-TI-01 Pressure'], label='Chamber (PT-TI-01)', color='green')
    ax.plot(dfc['t_plot'], dfc['PT-O2-05 Pressure'], label='Ox Feed (PT-O2-05)', color='steelblue')
    ax.plot(dfc['t_plot'], dfc['PT-H2-03 Pressure'], label='Fuel Feed (PT-H2-03)', color='orangered')
    ax.set_xlim(0, tmax - tmin)
    ax.xaxis.set_major_formatter(formatter)
    ax.set_title(title)
    ax.set_xlabel("Elapsed Time (ms) ×10")
    ax.set_ylabel("Pressure (psi)")
    ax.grid(True)

axs2[-1].legend(loc='upper right', fontsize=9)
fig2.suptitle("Figure 2 – Chamber, Oxidizer, and Fuel Feed Pressures", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.96])

plt.show()