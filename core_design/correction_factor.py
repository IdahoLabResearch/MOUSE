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

def corrected_keff_2d(depletion_2d_results_file, total_height): 

    geometry = openmc.Geometry.from_xml()
    root_universe = geometry.root_universe
    root = root_universe
    group_edges = np.array([1e-5, 6.7e-2, 3.2e-1, 1, 4, 9.88, 4.81e1, 4.54e2, 4.9e4, 1.83e5, 8.21e5, 4e7])   # Three energy groups
    groups = openmc.mgxs.EnergyGroups(group_edges)

    mgxs_lib = openmc.mgxs.Library(geometry)
    mgxs_lib.energy_groups = groups
    mgxs_lib.mgxs_types = ['absorption', 'diffusion-coefficient', 'transport', 'scatter matrix', 'total', 'scatter']
    mgxs_lib.domain_type = 'universe'
    mgxs_lib.domains = [root_universe]
    mgxs_lib.build_library()

    # Find all state point files generated during depletion
    statepoint_files = sorted(glob.glob('openmc_simulation_n*.h5'), key=natural_sort_key)
    # Initialize lists to store time steps and keff_3D values
    time_steps = []
    keff_2d_corrected_values = []
    keff_2d_values = []

    # Read the depletion results file to extract time steps
    depletion_results = openmc.deplete.Results("depletion_results.h5")
    time, _ = depletion_results.get_keff()
    time_days = [t / 86400 for t in time]  # Convert time to days

    # Open CSV file for writing results
    with open('depletion_output3.csv', 'w', newline='') as csvfile:

        writer = csv.writer(csvfile)
        writer.writerow(['keff_2D', 'P_NL', 'keff_3D', 'keff_3D_Uncertainty'])
    
        # Iterate over all state point files
        for idx, sp_file in enumerate(statepoint_files):
            # Load the state point
            sp = openmc.StatePoint(sp_file)
            try:
                mgxs_lib.load_from_statepoint(sp)
            except LookupError as e: # If the tallies are not retreived from one of the statepoint files
                print(f"Error loading MGXS from statepoint: {e}")
                continue

            keff_2d = sp.keff.nominal_value
            keff_2d_uncertainty = sp.keff.std_dev
            abs_xs_mg = mgxs_lib.get_mgxs(root_universe, 'absorption')
            trans_xs_mg = mgxs_lib.get_mgxs(root_universe, 'transport')
            total_xs_mg = mgxs_lib.get_mgxs(root_universe, 'total')
            scatter_xs_mg = mgxs_lib.get_mgxs(root_universe, 'scatter')

            #diffcoeff_mg = mgxs_lib.get_mgxs(root_universe, 'diffusion-coefficient')

            abs_xs_array = abs_xs_mg.get_xs(nuclide='total', mgxs_type='absorption', collapse=True)
            #diffcoeff_array = diffcoeff_mg.get_xs(nuclide='total', mgxs_type='diffusion-coefficient', collapse=True)
            trans_xs_array = trans_xs_mg.get_xs(nuclide='total', mgxs_type='transport', collapse=True)
            total_xs_array = total_xs_mg.get_xs(nuclide='total', mgxs_type='total', collapse=True)
            scatter_xs_array = scatter_xs_mg.get_xs(nuclide='total', mgxs_type='scatter', collapse=True)


            abs_xs_1g = float(np.mean(abs_xs_array))
            #diffcoeff_1g = float(np.mean(diffcoeff_array))
            trans_xs_1g = float(np.mean(trans_xs_array))
            total_xs_1g = float(np.mean(total_xs_array))
            scatter_xs_1g = float(np.mean(scatter_xs_array))

            diffcoeff_1g = 1/(3*trans_xs_1g)
            L_sqrt = diffcoeff_1g / abs_xs_1g
            
            extrapolated_height = total_height + (2*diffcoeff_1g)

            Bg_sqrt = (np.pi / extrapolated_height) ** 2 
            P_nl = 1 / (1 + (L_sqrt * Bg_sqrt))

            keff_2d_corrected = P_nl * keff_2d
            keff_2d_corrected_uncertainty = P_nl * keff_2d_uncertainty
        
            # Store the time step and keff_3D value
            time_steps.append(time_days[idx])  # Use the actual time in days
            keff_2d_corrected_values.append(keff_2d_corrected)
            keff_2d_values.append(keff_2d)

            print(f"Time Step: {idx + 1}")
            print(f"keff_2D: {keff_2d:.5f}+/-{keff_2d_uncertainty:.5f}")
            #print(f"abs_xs:{abs_xs_1g:.5f}")
            #print(f"diff_coeff:{diffcoeff_1g:.5f}")
            print(f"P_NL:{P_nl:.5f}")
            print(f"keff_2D_corrected: {keff_2d_corrected:.5f}+/-{keff_2d_corrected_uncertainty:.5f}")
            #print(f"trans_xs:{trans_xs_1g:.5f}")
            #print(f"total_xs:{total_xs_1g:.5f}")
            #print(f"scatter_xs:{scatter_xs_1g:.5f}")
        
            writer.writerow([f"{keff_2d:.5f}", f"{P_nl:.5f}", f"{keff_2d_corrected:.5f}", f"{keff_2d_corrected_uncertainty:.5f}"])
    

    # Plot keff_2D and keff_3D vs. Actual Time Steps
    plt.figure()
    plt.plot(time_steps, keff_2d_values, marker='o', linestyle='-', color='r', label='keff_2D')
    plt.plot(time_steps, keff_2d_corrected_values, marker='o', linestyle='-', color='g', label='corrected_keff_2D')
    plt.xlabel('Time [days]')
    plt.ylabel('k-effective')
    plt.title('Comparison of keff_2D and corrected_keff_2D vs. Time')
    plt.grid(True)
    plt.legend()
    plt.savefig('keff_comparison_vs_Time.png')
    plt.show()
    
    # Initialize variables for cycle length calculation
    cycle_length = None

    # Iterate through the corrected keff values to find the cycle length
    for i in range(1, len(keff_2d_corrected_values)):
        k1 = keff_2d_corrected_values[i - 1]
        k2 = keff_2d_corrected_values[i]
        t1 = time_steps[i - 1]
        t2 = time_steps[i]

        # Check if k1 and k2 bracket the value of 1.0
        if (k1 < 1.0 <= k2) or (k2 < 1.0 <= k1):
            # Perform linear interpolation to find the time when k = 1.0
            slope = (k2 - k1) / (t2 - t1)
            cycle_length = t1 + (1.0 - k1) / slope
            break  # Break once the cycle length is found

    if cycle_length is not None:
        round_cycle_length = round(cycle_length, 0)
        print(f"Estimated fuel cycle length: {round_cycle_length} days")
    else:
        print("k = 1.0 not reached within the given time steps.")

    return round_cycle_length         

