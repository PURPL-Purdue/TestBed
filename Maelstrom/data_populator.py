from pathlib import Path
import numpy as np
import pandas as pd
from data_regression import data_regression

# ---------------------------
# CONFIG
# ---------------------------
INPUT_CSV  = "batch_inputs.csv"
OUTPUT_CSV = "batch_outputs.csv"

# Outputs you care about (mirrors your console summary)
OUTPUT_KEYS = [
    # FLOWS
    "m_dot_fu_mean",
    "m_dot_ox_mean",
    "OF_mean",
    # PRESSURES
    "p_fu_mean_psig",
    "p_ox_mean_psig",
    "pc_mean_psig",
    "overall_pressure_dev_psig",
    # DENSITIES
    "rho_fu_mean",
    "rho_ox_mean",
    # C*
    "cstar_expected_mean_mps",
    "cstar_actual_mean_mps",
    "cstar_eff_mean",
]

# ---------------------------
# READ INPUT (auto-detect delimiter ; vs , and skip Excel "sep=" line)
# ---------------------------
raw = Path(INPUT_CSV).read_text(encoding="utf-8", errors="ignore").splitlines()

skip = 0
if raw and raw[0].strip().lower().startswith("sep="):
    skip = 1

# find header line (after optional sep=)
header_line = raw[skip] if len(raw) > skip else ""
# choose delimiter based on header content
sep = ";" if header_line.count(";") >= header_line.count(",") else ","

df_in = pd.read_csv(
    INPUT_CSV,
    sep=sep,
    skiprows=skip,
    dtype={
        "file_name": str,
        "fuel_choice": str,
        "fuel_CdA": float,
        "ox_CdA": float,
        "A_star": float,
        "L_star": float,  # passthrough only
    },
)

# drop any Excel ghost columns
df_in = df_in.loc[:, ~df_in.columns.astype(str).str.contains(r"^Unnamed")]

# ---------------------------
# RUN BATCH
# ---------------------------
rows_out = []

for idx, row in df_in.iterrows():
    file_name = str(row.get("file_name", "")).strip()

    # Skip blank / NaN filenames safely
    if file_name == "" or file_name.lower() == "nan":
        print(f"Skipping row {idx}: blank file_name")
        rows_out.append({k: np.nan for k in OUTPUT_KEYS})
        continue

    try:
        summary, _ = data_regression(
            file_name=file_name,
            fuel_choice=str(row["fuel_choice"]).strip(),
            fuel_CdA=float(row["fuel_CdA"]),
            ox_CdA=float(row["ox_CdA"]),
            A_star=float(row["A_star"]),
            print_summary=False,
            clear_terminal=False,
        )
        rows_out.append({k: summary.get(k, np.nan) for k in OUTPUT_KEYS})

    except Exception as e:
        print(f"Row {idx} FAILED ({file_name}): {type(e).__name__}: {e}")
        rows_out.append({k: np.nan for k in OUTPUT_KEYS})

# ---------------------------
# MERGE + CLEAN (prevent duplicate column names)
# ---------------------------
df_out = pd.concat([df_in.reset_index(drop=True), pd.DataFrame(rows_out)], axis=1)

# Drop duplicate column names (keep computed results)
df_out = df_out.loc[:, ~df_out.columns.duplicated(keep="last")]

# ---------------------------
# FORMATTING
# ---------------------------
# mdot: 3 decimals
for c in ["m_dot_fu_mean", "m_dot_ox_mean"]:
    if c in df_out.columns:
        df_out[c] = pd.to_numeric(df_out[c], errors="coerce").round(3)

# OF + cstar eff: 3 decimals
for c in ["OF_mean", "cstar_eff_mean"]:
    if c in df_out.columns:
        df_out[c] = pd.to_numeric(df_out[c], errors="coerce").round(3)

# pressures: integers
for c in ["p_fu_mean_psig", "p_ox_mean_psig", "pc_mean_psig", "overall_pressure_dev_psig"]:
    if c in df_out.columns:
        df_out[c] = pd.to_numeric(df_out[c], errors="coerce").round(0).astype("Int64")

# densities: integers
for c in ["rho_fu_mean", "rho_ox_mean"]:
    if c in df_out.columns:
        df_out[c] = pd.to_numeric(df_out[c], errors="coerce").round(0).astype("Int64")

# cstar: integers
for c in ["cstar_expected_mean_mps", "cstar_actual_mean_mps"]:
    if c in df_out.columns:
        df_out[c] = pd.to_numeric(df_out[c], errors="coerce").round(0).astype("Int64")

# ---------------------------
# COLUMN ORDER (inputs first, then outputs)
# ---------------------------
final_cols = [
    "file_name",
    "fuel_choice",
    "fuel_CdA",
    "ox_CdA",
    "A_star",
    "L_star",  # untouched passthrough
] + OUTPUT_KEYS

df_out = df_out[[c for c in final_cols if c in df_out.columns]]

# ---------------------------
# WRITE OUTPUT
# IMPORTANT: write with ; for EU Excel so columns split correctly
# ---------------------------
df_out.to_csv(OUTPUT_CSV, index=False, sep=",")
print(f"Wrote {OUTPUT_CSV} (sep='{';'}')")