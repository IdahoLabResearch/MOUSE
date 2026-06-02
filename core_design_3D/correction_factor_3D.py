# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

import openmc
import openmc.deplete
import openmc.mgxs
import matplotlib.pyplot as plt
import numpy as np
import glob
import csv
import re


def natural_sort_key(s):
    """Sort keys in a natural order (e.g., n0, n1, ..., n11)."""
    return [int(text) if text.isdigit() else text for text in re.split(r'(\d+)', s)]

def keff_3d(depletion_2d_results_file, total_height):

    # Find all state point files generated during depletion
    statepoint_files = sorted(glob.glob('openmc_simulation_n*.h5'), key=natural_sort_key)
    time_steps = []
    keff_3d_values = []
    keff_3d_values_uncertainty = []

    depletion_results = openmc.deplete.Results("depletion_results.h5")
    time, _ = depletion_results.get_keff()
    time_days = [t / 86400 for t in time]

    with open('depletion_output3.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['keff_3D', 'keff_3D_Uncertainty'])

        for idx, sp_file in enumerate(statepoint_files):
            sp = openmc.StatePoint(sp_file)
            keff_3d_val = sp.keff.nominal_value
            keff_3d_uncertainty = sp.keff.std_dev
            time_steps.append(time_days[idx])
            keff_3d_values.append(keff_3d_val)
            keff_3d_values_uncertainty.append(keff_3d_uncertainty)

            print(f"Time Step: {idx + 1}")
            print(f"keff_3D: {keff_3d_val:.5f}+/-{keff_3d_uncertainty:.5f}")
            writer.writerow([f"{keff_3d_val:.5f}", f"{keff_3d_uncertainty:.5f}"])

    plt.figure()
    plt.errorbar(time_steps, keff_3d_values, yerr=keff_3d_values_uncertainty,
                 marker='o', linestyle='-', color='r', label='keff_3D')
    plt.xlabel('Time [days]')
    plt.ylabel('k-effective')
    plt.title('keff_3D vs. Time')
    plt.grid(True)
    plt.legend()
    plt.savefig('keff3D_comparison_vs_Time.png')
    plt.show()

    cycle_length = None
    for i in range(1, len(keff_3d_values)):
        k1, k2 = keff_3d_values[i - 1], keff_3d_values[i]
        t1, t2 = time_steps[i - 1], time_steps[i]
        if (k1 < 1.0 <= k2) or (k2 < 1.0 <= k1):
            slope = (k2 - k1) / (t2 - t1)
            cycle_length = t1 + (1.0 - k1) / slope
            break

    if cycle_length is not None:
        round_cycle_length = round(cycle_length, 0)
        print(f"Estimated fuel cycle length: {round_cycle_length} days")
    else:
        print("k = 1.0 not reached within the given time steps.")
        raise ValueError("Cannot compute fuel cycle length: k=1.0 was never reached.")

    return round_cycle_length, time_steps, keff_3d_values