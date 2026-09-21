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


def _build_mgxs_library_from_xml():
    """Build the MGXS library used by the 2D axial leakage correction."""
    geometry = openmc.Geometry.from_xml()
    root_universe = geometry.root_universe

    group_edges = np.array([
        1e-5, 6.7e-2, 3.2e-1, 1, 4, 9.88,
        4.81e1, 4.54e2, 4.9e4, 1.83e5, 8.21e5, 4e7
    ])
    groups = openmc.mgxs.EnergyGroups(group_edges)

    mgxs_lib = openmc.mgxs.Library(geometry)
    mgxs_lib.energy_groups = groups
    mgxs_lib.mgxs_types = [
        'absorption',
        'diffusion-coefficient',
        'transport',
        'scatter matrix',
        'total',
        'scatter'
    ]
    mgxs_lib.domain_type = 'universe'
    mgxs_lib.domains = [root_universe]
    mgxs_lib.build_library()

    return root_universe, mgxs_lib


def _correct_statepoint_keff(
    statepoint_file,
    total_height,
    core_radius=None,
    root_universe=None,
    mgxs_lib=None
):
    """Return raw and axial leakage corrected keff for one statepoint."""
    if root_universe is None or mgxs_lib is None:
        root_universe, mgxs_lib = _build_mgxs_library_from_xml()

    sp = openmc.StatePoint(statepoint_file)
    try:
        mgxs_lib.load_from_statepoint(sp)

        keff_2d = float(sp.keff.nominal_value)
        keff_2d_uncertainty = float(sp.keff.std_dev)

        abs_xs_mg = mgxs_lib.get_mgxs(root_universe, 'absorption')
        trans_xs_mg = mgxs_lib.get_mgxs(root_universe, 'transport')

        abs_xs_array = abs_xs_mg.get_xs(
            nuclide='total',
            mgxs_type='absorption',
            collapse=True
        )
        trans_xs_array = trans_xs_mg.get_xs(
            nuclide='total',
            mgxs_type='transport',
            collapse=True
        )

        abs_xs_1g = float(np.mean(abs_xs_array))
        trans_xs_1g = float(np.mean(trans_xs_array))

        diffcoeff_1g = 1.0 / (3.0 * trans_xs_1g)
        diffusion_length_squared = diffcoeff_1g / abs_xs_1g

        extrapolated_height = total_height + (2.0 * diffcoeff_1g)
        buckling_axial = (np.pi / extrapolated_height) ** 2
        p_nl_axial = 1.0 / (
            1.0 + diffusion_length_squared * buckling_axial
        )

        if core_radius is not None and core_radius > 0.0:
            extrapolated_radius = core_radius + (2.0 * diffcoeff_1g)
            buckling_radial = (2.405 / extrapolated_radius) ** 2
            buckling_total = buckling_axial + buckling_radial
            p_nl_total = 1.0 / (
                1.0 + diffusion_length_squared * buckling_total
            )
        else:
            p_nl_total = np.nan

        keff_corrected = p_nl_axial * keff_2d
        keff_corrected_uncertainty = (
            p_nl_axial * keff_2d_uncertainty
        )

        return {
            'keff_2d': keff_2d,
            'keff_2d_uncertainty': keff_2d_uncertainty,
            'keff_corrected': keff_corrected,
            'keff_corrected_uncertainty': keff_corrected_uncertainty,
            'axial_non_leakage_probability': p_nl_axial,
            'total_non_leakage_probability': p_nl_total,
        }
    finally:
        sp.close()


def corrected_keff_static(statepoint_file, total_height, core_radius=None):
    """
    Apply the 2D axial leakage correction to one static OpenMC statepoint.

    Unlike :func:`corrected_keff_2d`, this function does not estimate a fuel
    cycle length and therefore does not require the keff curve to cross 1.0.
    It is intended for lifecycle snapshot calculations such as shutdown
    margin and isothermal temperature coefficient evaluations.
    """
    return _correct_statepoint_keff(
        statepoint_file,
        total_height,
        core_radius=core_radius
    )


def corrected_keff_2d(
    depletion_2d_results_file,
    total_height,
    core_radius=None,
    plotting=True
):
    """
    Apply a leakage correction to 2D depletion keff values using a simple
    non-leakage probability approximation.

    Parameters
    ----------
    depletion_2d_results_file : openmc.deplete.Results
        OpenMC depletion results file object.
    total_height : float
        Total axial height used for the leakage correction [cm].
        Typically: Active Height + 2 * Axial Reflector Thickness
    core_radius : float or None
        Effective radial core radius [cm] for total leakage estimation.
        If None, only axial leakage is computed and total leakage metrics are
        returned as np.nan.

    Returns
    -------
    round_cycle_length : float
        Estimated fuel cycle length [days], based on the corrected keff curve crossing 1.0.
    time_steps : list[float]
        Cumulative depletion time points [days].
    keff_2d_values : list[float]
        Uncorrected 2D keff values.
    keff_2d_corrected_values : list[float]
        Axial-leakage-corrected keff values.
    bol_axial_non_leakage_probability : float
        Beginning-of-life axial non-leakage probability.
    estimated_axial_leakage_bol_pct : float
        Estimated beginning-of-life axial leakage [%].
    bol_total_non_leakage_probability : float
        Beginning-of-life total non-leakage probability including axial and radial
        buckling. Returned as np.nan when core_radius is not provided.
    estimated_total_leakage_bol_pct : float
        Estimated beginning-of-life total leakage [%]. Returned as np.nan when
        core_radius is not provided.
    cycle_length_status : str
        Indicates whether the cycle length was interpolated, extrapolated,
        identified as subcritical at BOL, or limited to the simulated endpoint.
    cycle_length_extrapolated : bool
        True only when the reported cycle length was extrapolated beyond the
        depletion schedule from the final two corrected keff points.
    """

    geometry = openmc.Geometry.from_xml()
    root_universe = geometry.root_universe

    group_edges = np.array([
        1e-5, 6.7e-2, 3.2e-1, 1, 4, 9.88,
        4.81e1, 4.54e2, 4.9e4, 1.83e5, 8.21e5, 4e7
    ])
    groups = openmc.mgxs.EnergyGroups(group_edges)

    mgxs_lib = openmc.mgxs.Library(geometry)
    mgxs_lib.energy_groups = groups
    mgxs_lib.mgxs_types = [
        'absorption',
        'diffusion-coefficient',
        'transport',
        'scatter matrix',
        'total',
        'scatter'
    ]
    mgxs_lib.domain_type = 'universe'
    mgxs_lib.domains = [root_universe]
    mgxs_lib.build_library()

    statepoint_files = sorted(glob.glob('openmc_simulation_n*.h5'), key=natural_sort_key)

    time_steps = []
    keff_2d_corrected_values = []
    keff_2d_values = []

    bol_axial_non_leakage_probability = np.nan
    bol_total_non_leakage_probability = np.nan
    estimated_axial_leakage_bol_pct = np.nan
    estimated_total_leakage_bol_pct = np.nan
    bol_metrics_set = False

    time, _ = depletion_2d_results_file.get_keff()
    time_days = [t / 86400 for t in time]

    with open('depletion_output3.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            'keff_2D',
            'P_nl_axial',
            'P_nl_total',
            'keff_3D',
            'keff_3D_Uncertainty',
            'Estimated_Axial_Leakage_pct_BOL',
            'Estimated_Total_Leakage_pct_BOL'
        ])

        for idx, sp_file in enumerate(statepoint_files):
            sp = openmc.StatePoint(sp_file)

            try:
                mgxs_lib.load_from_statepoint(sp)
            except LookupError as e:
                print(f"Error loading MGXS from statepoint: {e}")
                continue

            keff_2d = sp.keff.nominal_value
            keff_2d_uncertainty = sp.keff.std_dev

            abs_xs_mg = mgxs_lib.get_mgxs(root_universe, 'absorption')
            trans_xs_mg = mgxs_lib.get_mgxs(root_universe, 'transport')
            total_xs_mg = mgxs_lib.get_mgxs(root_universe, 'total')
            scatter_xs_mg = mgxs_lib.get_mgxs(root_universe, 'scatter')

            abs_xs_array = abs_xs_mg.get_xs(
                nuclide='total',
                mgxs_type='absorption',
                collapse=True
            )
            trans_xs_array = trans_xs_mg.get_xs(
                nuclide='total',
                mgxs_type='transport',
                collapse=True
            )
            total_xs_array = total_xs_mg.get_xs(
                nuclide='total',
                mgxs_type='total',
                collapse=True
            )
            scatter_xs_array = scatter_xs_mg.get_xs(
                nuclide='total',
                mgxs_type='scatter',
                collapse=True
            )

            abs_xs_1g = float(np.mean(abs_xs_array))
            trans_xs_1g = float(np.mean(trans_xs_array))
            total_xs_1g = float(np.mean(total_xs_array))
            scatter_xs_1g = float(np.mean(scatter_xs_array))

            diffcoeff_1g = 1 / (3 * trans_xs_1g)
            diffusion_length_squared = diffcoeff_1g / abs_xs_1g

            extrapolated_height = total_height + (2 * diffcoeff_1g)
            buckling_axial = (np.pi / extrapolated_height) ** 2
            p_nl_axial = 1 / (1 + diffusion_length_squared * buckling_axial)

            if core_radius is not None and core_radius > 0.0:
                extrapolated_radius = core_radius + (2 * diffcoeff_1g)
                buckling_radial = (2.405 / extrapolated_radius) ** 2
                buckling_total = buckling_axial + buckling_radial
                p_nl_total = 1 / (1 + diffusion_length_squared * buckling_total)
            else:
                p_nl_total = np.nan

            keff_2d_corrected = p_nl_axial * keff_2d
            keff_2d_corrected_uncertainty = p_nl_axial * keff_2d_uncertainty

            if not bol_metrics_set:
                bol_axial_non_leakage_probability = p_nl_axial
                estimated_axial_leakage_bol_pct = (1.0 - p_nl_axial) * 100.0

                if not np.isnan(p_nl_total):
                    bol_total_non_leakage_probability = p_nl_total
                    estimated_total_leakage_bol_pct = (1.0 - p_nl_total) * 100.0
                else:
                    bol_total_non_leakage_probability = np.nan
                    estimated_total_leakage_bol_pct = np.nan

                bol_metrics_set = True

            time_steps.append(round(float(time_days[idx]), 4))
            keff_2d_corrected_values.append(keff_2d_corrected)
            keff_2d_values.append(keff_2d)

            print(f"Time Step: {idx + 1}")
            print(f"keff_2D: {keff_2d:.5f}+/-{keff_2d_uncertainty:.5f}")
            print(f"P_nl_axial: {p_nl_axial:.5f}")
            if not np.isnan(p_nl_total):
                print(f"P_nl_total: {p_nl_total:.5f}")
            print(f"keff_2D_corrected: {keff_2d_corrected:.5f}+/-{keff_2d_corrected_uncertainty:.5f}")

            if idx == 0:
                print(f"Estimated Axial Leakage (BOL): {estimated_axial_leakage_bol_pct:.3f} %")
                if not np.isnan(estimated_total_leakage_bol_pct):
                    print(f"Estimated Total Leakage (BOL): {estimated_total_leakage_bol_pct:.3f} %")
                else:
                    print("Estimated Total Leakage (BOL): not computed (core_radius not provided)")

            writer.writerow([
                f"{keff_2d:.5f}",
                f"{p_nl_axial:.5f}",
                f"{p_nl_total:.5f}" if not np.isnan(p_nl_total) else "",
                f"{keff_2d_corrected:.5f}",
                f"{keff_2d_corrected_uncertainty:.5f}",
                f"{estimated_axial_leakage_bol_pct:.5f}" if idx == 0 else "",
                f"{estimated_total_leakage_bol_pct:.5f}" if idx == 0 and not np.isnan(estimated_total_leakage_bol_pct) else ""
            ])

    # plt.figure()
    # plt.plot(time_steps, keff_2d_values, marker='o', linestyle='-', color='r', label='keff_2D')
    # plt.plot(time_steps, keff_2d_corrected_values, marker='o', linestyle='-', color='g', label='corrected_keff_2D')
    # plt.xlabel('Time [days]')
    # plt.ylabel('k-effective')
    # plt.title('Comparison of keff_2D and corrected_keff_2D vs. Time')
    # plt.grid(True)
    # plt.legend()
    # plt.savefig('keff_comparison_vs_Time.png')
    # plt.show()


    if str(plotting).upper() == 'Y' or plotting is True:
        # Plot only the operating-period portion of the depletion history.
        # This does not change the depletion calculation or cycle-length result.
        plot_limit_days = 2000.0
        plot_indices = [
            i for i, t in enumerate(time_steps)
            if t <= plot_limit_days
        ]

        # Include the first point after the limit so the trend is visible.
        if plot_indices:
            last_index = min(plot_indices[-1] + 2, len(time_steps))
        else:
            last_index = min(2, len(time_steps))

        plt.figure()
        plt.plot(
            time_steps[:last_index],
            keff_2d_values[:last_index],
            marker='o',
            linestyle='-',
            color='r',
            label='keff_2D'
        )
        plt.plot(
            time_steps[:last_index],
            keff_2d_corrected_values[:last_index],
            marker='o',
            linestyle='-',
            color='g',
            label='corrected_keff_2D'
        )
        plt.axhline(y=1.0, color='k', linestyle='--', label='k = 1')
        plt.xlabel('Time [days]')
        plt.ylabel('k-effective')
        plt.title('Comparison of keff_2D and corrected_keff_2D vs. Time')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig('keff_comparison_vs_Time.png', dpi=300)
        plt.close()

    if not keff_2d_corrected_values:
        raise ValueError("No corrected depletion keff values were produced.")

    cycle_length = None
    cycle_length_status = None
    cycle_length_extrapolated = False

    # A design that is already subcritical cannot have a positive operating
    # cycle. Keep the workflow running so it can be reported and screened out.
    if keff_2d_corrected_values[0] <= 1.0:
        cycle_length = 0.0
        cycle_length_status = 'BOL subcritical'
        print(
            "WARNING: corrected BOL keff is at or below 1.0; the reported "
            "operating cycle length is 0 days."
        )
    else:
        # A fuel-cycle endpoint is a downward crossing of k = 1.
        for i in range(1, len(keff_2d_corrected_values)):
            k1 = keff_2d_corrected_values[i - 1]
            k2 = keff_2d_corrected_values[i]
            t1 = time_steps[i - 1]
            t2 = time_steps[i]

            if k1 >= 1.0 and k2 <= 1.0:
                if k2 == k1:
                    cycle_length = t2
                else:
                    slope = (k2 - k1) / (t2 - t1)
                    cycle_length = t1 + (1.0 - k1) / slope
                cycle_length_status = 'Interpolated k=1 crossing'
                break

    # If the schedule ends while the core is still supercritical, extrapolate
    # only when the final trend is decreasing. This is a scoping estimate; the
    # status is retained so downstream studies can distinguish it from a
    # bracketed crossing.
    if cycle_length is None:
        if len(time_steps) >= 2:
            t1, t2 = time_steps[-2], time_steps[-1]
            k1 = keff_2d_corrected_values[-2]
            k2 = keff_2d_corrected_values[-1]
            slope = (k2 - k1) / (t2 - t1)
        else:
            slope = np.nan

        if np.isfinite(slope) and slope < 0.0:
            extrapolated_cycle_length = time_steps[-1] + (
                1.0 - keff_2d_corrected_values[-1]
            ) / slope
            if extrapolated_cycle_length >= time_steps[-1]:
                cycle_length = extrapolated_cycle_length
                cycle_length_status = 'Extrapolated beyond depletion schedule'
                cycle_length_extrapolated = True
                print(
                    "WARNING: corrected keff remained above 1.0 through the "
                    "depletion schedule. Cycle length was extrapolated from "
                    "the final two corrected keff points."
                )

        if cycle_length is None:
            # A flat/rising terminal trend cannot support a defensible k=1
            # extrapolation. Use the last simulated time as a lower bound so
            # the cost workflow still completes, and flag the result clearly.
            cycle_length = float(time_steps[-1])
            cycle_length_status = 'Simulated endpoint lower bound'
            print(
                "WARNING: corrected keff remained above 1.0 and the final "
                "trend was not decreasing. Fuel lifetime is reported as the "
                "last simulated time (a lower bound); no extrapolation was "
                "performed."
            )

    round_cycle_length = round(float(cycle_length), 0)
    print(
        f"Estimated fuel cycle length: {round_cycle_length} days "
        f"[{cycle_length_status}]"
    )

    return (
        round_cycle_length,
        time_steps,
        keff_2d_values,
        keff_2d_corrected_values,
        bol_axial_non_leakage_probability,
        estimated_axial_leakage_bol_pct,
        bol_total_non_leakage_probability,
        estimated_total_leakage_bol_pct,
        cycle_length_status,
        cycle_length_extrapolated
    )
