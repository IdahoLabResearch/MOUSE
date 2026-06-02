import re
import glob
from pathlib import Path

import openmc
import pandas as pd


def natural_sort_key(s: str):
    """Natural sort key: n0, n1, ..., n10 instead of n0, n1, n10, n2..."""
    return [int(text) if text.isdigit() else text for text in re.split(r'(\d+)', s)]


def compute_pin_peaking_factors(current_dir="."):
    """
    Compute peaking factors for all OpenMC depletion statepoints.

    Supports both:
      1) distribcell-based tally (LTMR-style)
      2) mesh-based tally (GCMR-style)

    For each statepoint:
      - Reads tally 'pin_power_kappa'
      - Sums kappa-fission power by pin/cell or mesh voxel
      - Computes PF = P_i / mean(P)

    Returns:
      summary       : DataFrame with columns [Step, Max_PF, Region_ID_Max]
      per_step_data : dict[step] -> DataFrame [Region_ID, Peaking_Factor, Step]
    """

    base = Path(current_dir)

    sp_files = glob.glob(str(base / "openmc_simulation_n*.h5"))
    sp_files = sorted(sp_files, key=natural_sort_key)

    if not sp_files:
        print("\n[PF] No depletion statepoint files found in:", base)
        print("[PF] Expected files like 'openmc_simulation_n0.h5', 'openmc_simulation_n1.h5', ...\n")
        return pd.DataFrame(), {}

    tally_name = "pin_power_kappa"
    results = []
    per_step_data = {}

    print("\n================ PEAKING FACTOR RESULTS ================\n")

    for sp_file in sp_files:
        sp_path = Path(sp_file)
        basename = sp_path.name

        m = re.search(r"n(\d+)\.h5", basename)
        if m:
            step_raw = int(m.group(1))
            step = step_raw + 1
        else:
            step = basename

        sp = openmc.StatePoint(str(sp_path))
        t = sp.get_tally(name=tally_name)
        df = t.get_pandas_dataframe(paths=False)

        # OpenMC returns MultiIndex (tuple) column names for multi-filter tallies
        # (e.g. MeshFilter + MaterialFilter in GCMR).  Flatten to plain strings so
        # that string operations below work regardless of OpenMC version.
        if any(isinstance(c, tuple) for c in df.columns):
            df.columns = [
                " ".join(str(x) for x in c if str(x).strip()).strip()
                if isinstance(c, tuple) else c
                for c in df.columns
            ]

        # Detect tally format: distribcell (LTMR) or mesh-based (GCMR).
        # Mesh column names depend on the OpenMC version and how many meshes have been
        # created in the session (the mesh ID counter increments globally, so when
        # build_openmc_model_GCMR is called twice — e.g. for Isothermal Temperature
        # Coefficients — the second mesh gets ID 2, not 1).  We detect the mesh
        # columns dynamically by regex rather than hard-coding "mesh 1 x" etc.
        x_cols = [c for c in df.columns if re.match(r"mesh \d+ x$", c)]
        y_cols = [c for c in df.columns if re.match(r"mesh \d+ y$", c)]
        z_cols = [c for c in df.columns if re.match(r"mesh \d+ z$", c)]
        flat_mesh_cols = [c for c in df.columns if re.match(r"mesh \d+$", c)]

        if "distribcell" in df.columns:
            tally_type = "distribcell"
            per_region = df.groupby("distribcell")["mean"].sum()
        elif x_cols and y_cols and z_cols:
            # New OpenMC format: separate x / y / z columns
            tally_type = "mesh_xyz"
            per_region = df.groupby([x_cols[0], y_cols[0], z_cols[0]])["mean"].sum()
        elif flat_mesh_cols:
            # Old OpenMC format: single "mesh N" column (flat voxel index)
            tally_type = "mesh_flat"
            per_region = df.groupby(flat_mesh_cols[0])["mean"].sum()
        else:
            raise ValueError(
                "[PF] Unsupported tally format. Expected a distribcell or mesh tally. "
                f"Found columns: {list(df.columns)}"
            )

        # Remove zero-power bins/cells
        per_region = per_region[per_region > 0]

        if len(per_region) == 0:
            print(f"[PF] WARNING: no positive-power regions found in step {step}. Skipping.")
            continue

        pf = per_region / per_region.mean()

        # Build region IDs aligned with the zero-power-filtered index
        if tally_type == "distribcell":
            region_ids = pf.index.tolist()
        elif tally_type == "mesh_xyz":
            region_ids = [f"({i},{j},{k})" for i, j, k in pf.index.tolist()]
        else:  # mesh_flat
            region_ids = [str(idx) for idx in pf.index.tolist()]

        out = pd.DataFrame({
            "Region_ID": region_ids,
            "Peaking_Factor": pf.values,
            "Step": step
        })

        per_step_data[step] = out

        print(f"--- Peaking factors for depletion step {step} ---")
        print(out[["Region_ID", "Peaking_Factor"]].to_string(index=False))
        print()

        results.append({
            "Step": step,
            "Max_PF": float(pf.max()),
            "Region_ID_Max": out.loc[out["Peaking_Factor"].idxmax(), "Region_ID"]
        })

    summary = pd.DataFrame(results).sort_values("Step")

    print("========== Peaking Factor Summary ==========")
    print(summary.to_string(index=False))
    print("============================================\n")

    return summary, per_step_data