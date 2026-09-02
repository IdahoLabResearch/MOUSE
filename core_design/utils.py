# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
import numpy as np
import openmc
import openmc.deplete
import watts
import traceback  # print full stack traces for OpenMC failures
import glob
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from core_design.correction_factor import corrected_keff_2d, corrected_keff_static
#from core_design.correction_factor import corrected_keff_steady_state       #Use this line and comment the previuos line if you run steady state 
from core_design.peaking_factor import compute_pin_peaking_factors

import pandas


# Dedicated particle count for each lifecycle temperature-coefficient snapshot.
# Depletion and cold-shutdown calculations retain their normal particle count.
TEMPERATURE_COEFFICIENT_PARTICLES = 100000

# Current MOUSE material densities at the operating reference temperature.
NAK_REFERENCE_DENSITY_G_CM3 = 0.85
ZRH_REFERENCE_DENSITY_G_CM3 = 5.6

# Mean linear expansion coefficient reported for unalloyed epsilon-ZrH1.83.
# The model uses H/Zr = 1.85. Geometry is held fixed in this density-only
# treatment.
ZRH_LINEAR_EXPANSION_PER_K = 9.15e-6


def _nak_eutectic_density_correlation_g_cm3(temperature_k):
    """Return NaK-78 density using Sodium-NaK Handbook Eqs. 1.5, 1.8, 1.9."""
    temperature_c = float(temperature_k) - 273.15
    potassium_weight_fraction = 0.778
    atomic_weight_k = 39.102
    atomic_weight_na = 22.9898
    potassium_atom_fraction = (
        potassium_weight_fraction / atomic_weight_k
        / (
            potassium_weight_fraction / atomic_weight_k
            + (1.0 - potassium_weight_fraction) / atomic_weight_na
        )
    )

    specific_volume_k = 1.0 / (
        0.8415
        - 2.172e-4 * temperature_c
        - 2.7e-8 * temperature_c ** 2
        + 4.77e-12 * temperature_c ** 3
    )
    specific_volume_na = 1.0 / (
        0.9453 - 2.2473e-4 * temperature_c
    )
    density_kg_m3 = 1000.0 / (
        potassium_atom_fraction * specific_volume_k
        + (1.0 - potassium_atom_fraction) * specific_volume_na
    )
    return density_kg_m3 / 1000.0


def _temperature_density_overrides(reference_temperature_k, temperature_k):
    """Return NaK and ZrH densities anchored to the MOUSE reference values."""
    nak_reference_correlation = _nak_eutectic_density_correlation_g_cm3(
        reference_temperature_k
    )
    nak_target_correlation = _nak_eutectic_density_correlation_g_cm3(
        temperature_k
    )
    nak_density = (
        NAK_REFERENCE_DENSITY_G_CM3
        * nak_target_correlation
        / nak_reference_correlation
    )

    linear_scale = 1.0 + ZRH_LINEAR_EXPANSION_PER_K * (
        float(temperature_k) - float(reference_temperature_k)
    )
    zrh_density = ZRH_REFERENCE_DENSITY_G_CM3 / linear_scale ** 3

    return {
        '_NaK Density Override': float(nak_density),
        '_ZrH Density Override': float(zrh_density),
    }


def circle_area(r):
    return (np.pi) * r ** 2


def cylinder_volume(r, h):
    return circle_area(r) * h


def sphere_volume(r):
    return 4 / 3 * np.pi * r ** 3


def circle_perimeter(r):
    return 2 * (np.pi) * r


def sphere_area(radius):
    area = 4 * np.pi * (radius ** 2)
    return area


def cylinder_radial_shell(r, h):
    # Calculates the lateral surface area of a cylinder
    return circle_perimeter(r) * h


def calculate_lattice_radius(params):
    """
    Backward-compatible helper.

    For LTMR, this returns the hex apothem used historically under the
    name 'Lattice Radius'. Prefer calculate_hex_apothem() in new code.
    """
    return calculate_hex_apothem(params)


def calculate_hex_edge_length(params):
    """
    Hex edge length used for the LTMR assembly boundary and drum layout.
    This matches the geometry logic used in the OpenMC LTMR template.
    """
    pin_pitch = 2 * params['Fuel Pin Radii'][-1] + params['Pin Gap Distance']
    return pin_pitch * (params['Number of Rings per Assembly'] - 1) + pin_pitch * 0.6


def calculate_hex_apothem(params):
    """
    Distance from the hex center to the middle of a flat face.
    This is the most useful single size measure for the LTMR hex lattice.
    """
    hex_edge_length = calculate_hex_edge_length(params)
    return np.sin(np.pi / 3) * hex_edge_length


def calculate_core_radius_from_hex(params):
    """
    Circular outer radius used by the 2D OpenMC radial model and by the
    leakage approximation.

    It is defined from the fuel-lattice hex apothem plus the radial reflector
    thickness.
    """
    return calculate_hex_apothem(params) + params['Radial Reflector Thickness']


def calculate_heat_flux(params):
    fuel_number = params['Fuel Pin Count']
    heat_transfer_surface = cylinder_radial_shell(
        params['Fuel Pin Radii'][-1],
        params['Active Height']
    ) * fuel_number * 1e-4  # convert from cm2 to m2

    return params['Power MWt'] / heat_transfer_surface  # MW/m^2


def calculate_pins_in_assembly(params, pin_type):
    # Get the ring configuration from the parameters
    rings = params['Pins Arrangement']
    # Keep only the last 'Number of Rings per Assembly' rings as specified in the parameters
    rings = rings[-params['Number of Rings per Assembly']:]
    return sum(row.count(pin_type) for row in rings)


def create_cells(regions: dict, materials: list) -> dict:
    return {
        key: openmc.Cell(name=key, fill=mat, region=value)
        for (key, value), mat in zip(regions.items(), materials)
    }


def calculate_number_of_rings(rings_over_one_edge):
    # Total number of positions given the number of rings along one edge
    return 2 * rings_over_one_edge * (rings_over_one_edge - 1) + \
        2 * sum(range(1, rings_over_one_edge - 1)) + \
        2 * rings_over_one_edge - 1


def calculate_number_fuel_elements_hpmr(rings_over_one_edge):
    total_number_of_rings = calculate_number_of_rings(rings_over_one_edge)
    number_of_heatpipe_pins = calculate_number_of_rings(int(np.ceil(rings_over_one_edge / 2)))
    return total_number_of_rings - number_of_heatpipe_pins


def number_of_heatpipes_hmpr(params):
    tot_rings_per_assembly = calculate_number_of_rings(params['Number of Rings per Assembly'])
    params['Number of Heatpipes per Assembly'] = tot_rings_per_assembly - params['Fuel Pin Count per Assembly']
    params['Number of Heatpipes'] = params['Number of Heatpipes per Assembly'] * params['Fuel Assemblies Count']


def calculate_total_number_of_TRISO_particles(params):
    compact_fuel_vol = cylinder_volume(params['Compact Fuel Radius'], params['Active Height'])
    one_particle_volume = sphere_volume(params['Fuel Pin Radii'][-1])
    number_of_particles_per_compact_fuel_vol = np.floor(
        params['Packing Fraction'] * compact_fuel_vol / one_particle_volume
    )
    params['Number Of TRISO Particles Per Compact Fuel'] = number_of_particles_per_compact_fuel_vol
    total_number_of_particles = number_of_particles_per_compact_fuel_vol * \
        calculate_number_of_rings(params['Assembly Rings'] - 1) * \
        calculate_number_of_rings(params['Core Rings'])
    params['Total Number of TRISO Particles'] = total_number_of_particles
    return total_number_of_particles


def calculate_heat_flux_TRISO(params):
    number_of_triso_particles = calculate_total_number_of_TRISO_particles(params)
    total_area_triso = number_of_triso_particles * sphere_area(params['Fuel Pin Radii'][0]) * 1e-4  # cm^2 to m^2
    heat_flux = params['Power MWt'] / total_area_triso
    return heat_flux


def create_universe_plot(materials_database, universe, plot_width, num_pixels, font_size, title, fig_size, output_file_name):
    import matplotlib.colors as mcolors

    potential_colors = { 
        'UZrH_alloy': 'red',
        'ZrH': 'yellow',
        'UO2': 'green',
        'UC': 'purple',
        'UCO': 'orange',
        'UN': 'cyan',
        'YHx': 'magenta',
        'NaK': 'blue',
        'Helium': 'grey',
        'Be': 'brown',
        'BeO': 'pink',
        'Zr': 'lime',
        'SS304': 'black',
        'B4C_natural': 'olive',
        'B4C_enriched': 'deepskyblue',
        'SiC': 'teal',
        'Graphite': 'coral',
        'buffer_graphite': 'gold',
        'PyC': 'salmon',
        'homog_TRISO': 'maroon',
        'heatpipe': 'seashell',
        'monolith_graphite': 'navy',
        'UZr': 'darkred',
        'ZrC': 'slategray',
        'MgO': 'lightyellow',
        'WB': 'darkgray',
        'W2B': 'dimgray',
        'WB4': 'lightgray',
        'WC': 'silver',
    }

    used_colors = set(mcolors.to_hex(c) for c in potential_colors.values())
    color_pool = [
        name for name, hex_val in mcolors.CSS4_COLORS.items()
        if mcolors.to_hex(hex_val) not in used_colors
    ]

    for mat_name in materials_database:
        if mat_name not in potential_colors:
            if not color_pool:
                raise ValueError(
                    f"Could not auto-assign a color for material '{mat_name}': "
                    f"no unique colors remaining in the CSS4 pool. "
                    f"Please manually add a color for this material in potential_colors."
                )
            auto_color = color_pool.pop(0)
            potential_colors[mat_name] = auto_color
            used_colors.add(mcolors.to_hex(auto_color))
            print(
                f"\033[93m--- WARNING: Material '{mat_name}' does not have a color specified "
                f"in potential_colors. Automatically assigned color: '{auto_color}'. "
                f"Please add a permanent entry for this material in the potential_colors "
                f"dictionary in create_universe_plot (utils.py) to suppress this warning.\033[0m"
            )

    colors = {
        materials_database[mat_name]: color
        for mat_name, color in potential_colors.items()
        if mat_name in materials_database
    }

    universe_plot = universe.plot(
        width=(plot_width, plot_width),
        pixels=(num_pixels, num_pixels),
        color_by='material',
        colors=colors
    )
    # Use a slightly smaller font for tick labels so 5 ticks fit cleanly
    label_font = font_size
    tick_font  = max(8, int(font_size * 0.75))

    universe_plot.set_xlabel('x [cm]', fontsize=label_font)
    universe_plot.set_ylabel('y [cm]', fontsize=label_font)
    universe_plot.set_title(title, fontsize=label_font)

    # For plots whose half-width is at least 1 cm, show 5 integer-valued
    # ticks per axis: -half, -quarter, 0, +quarter, +half (with half_int
    # rounded UP so the data is always fully enclosed). For sub-cm plots
    # (e.g. the zoomed fuel assembly or the TRISO particle, where half
    # is on the order of 0.05 cm), integer ticks are nonsensical — keep
    # matplotlib's default automatic ticking.
    half = plot_width / 2.0
    if half >= 1.0:
        half_int    = max(1, int(np.ceil(half)))
        quarter_int = max(0, int(round(half_int / 2.0)))
        ticks = [-half_int, -quarter_int, 0, quarter_int, half_int]
        universe_plot.set_xticks(ticks)
        universe_plot.set_yticks(ticks)
        universe_plot.set_xlim(-half_int, half_int)
        universe_plot.set_ylim(-half_int, half_int)
    else:
        universe_plot.set_xlim(-half, half)
        universe_plot.set_ylim(-half, half)
    universe_plot.tick_params(axis='x', labelsize=tick_font)
    universe_plot.tick_params(axis='y', labelsize=tick_font)

    fig = universe_plot.figure
    fig.set_size_inches(fig_size, fig_size)

    universe_materials = [cell.fill for cell in universe.get_all_cells().values()]
    used_materials = set(universe_materials)

    legend_patches = [
        mpatches.Patch(color=color, label=mat_name)
        for mat_name, color in potential_colors.items()
        if mat_name in materials_database and materials_database[mat_name] in used_materials
    ]
    universe_plot.legend(
        handles=legend_patches,
        fontsize=font_size,
        loc='center left',
        bbox_to_anchor=(1, 0.5)
    )
    fig.savefig(output_file_name, bbox_inches='tight')


def openmc_depletion(params, lattice_geometry, settings):

    openmc.config['cross_sections'] = params['cross_sections_xml_location']

    operator = openmc.deplete.CoupledOperator(
        openmc.Model(geometry=lattice_geometry, settings=settings),
        chain_file=params['simplified_chain_thermal_xml']
    )

    if 'Burnup Steps' in params:
        burnup_steps_list_MWd_per_Kg = params['Burnup Steps']
        burnup_step = np.array(burnup_steps_list_MWd_per_Kg)
        burnup = np.diff(burnup_step, prepend=0.0)  # step-wise burnup increments

        integrator = openmc.deplete.PredictorIntegrator(
            operator,
            burnup,
            1000000 * params['Power MWt'],
            timestep_units='MWd/kg'
        )
    elif 'Time Steps' in params:
        time_steps_list = params['Time Steps']
        power_list = [params['Power MWt'] * 1e6] * len(time_steps_list)
        integrator = openmc.deplete.CECMIntegrator(operator, time_steps_list, power_list)

    print("Starting depletion")
    integrator.integrate()
    print("Depletion complete")

    depletion_2d_results_file = openmc.deplete.Results("./depletion_results.h5")

    # corrected_keff_2d returns:
    # 1) fuel lifetime in days
    # 2) cumulative depletion time points in days
    # 3) raw 2D keff values
    # 4) 3D-corrected keff values
    # 5) beginning-of-life axial non-leakage probability
    # 6) beginning-of-life estimated axial leakage percentage
    # 7) beginning-of-life total non-leakage probability (NaN if core radius is unavailable)
    # 8) beginning-of-life estimated total leakage percentage (NaN if core radius is unavailable)
    (
        fuel_lifetime_days,
        time_steps,
        keff_2d_values,
        keff_2d_values_corrected,
        bol_axial_non_leakage_probability,
        bol_axial_leakage_percent,
        bol_total_non_leakage_probability,
        bol_total_leakage_percent,
    ) = corrected_keff_2d(
        depletion_2d_results_file,
        params['Active Height'] + 2 * params['Axial Reflector Thickness'],
        core_radius=params.get('Core Radius', np.nan)
    )

    try:
        pf_summary, pf_per_step = compute_pin_peaking_factors(".")
        if pf_summary.empty:
            raise ValueError("No peaking factor results were produced.")

        # Peaking-factor Step is one based, while time_steps is zero based.
        # Exclude statepoints after the corrected operating keff reaches 1.0.
        pf_summary = pf_summary.copy()
        pf_summary['Time_days'] = pf_summary['Step'].apply(
            lambda step: float(time_steps[int(step) - 1])
            if 1 <= int(step) <= len(time_steps) else np.nan
        )
        invalid_steps = pf_summary['Time_days'].isna()
        if invalid_steps.any():
            invalid_step_values = pf_summary.loc[invalid_steps, 'Step'].tolist()
            raise ValueError(
                "Peaking factor steps do not map to depletion times: "
                f"{invalid_step_values}"
            )

        pf_summary = pf_summary.loc[
            pf_summary['Time_days'] <= fuel_lifetime_days
        ].copy()
        if pf_summary.empty:
            raise ValueError(
                "No peaking factor statepoint occurs at or before the "
                "calculated end of operating life."
            )

        idx_max = pf_summary['Max_PF'].idxmax()
        params['Max Peaking Factor'] = float(
            pf_summary.loc[idx_max, 'Max_PF']
        )
        params['Step with Max Peaking Factor'] = int(
            pf_summary.loc[idx_max, 'Step']
        )
        params['Region ID with Max Peaking Factor'] = (
            pf_summary.loc[idx_max, 'Region_ID_Max']
        )
        params['Max Peaking Factors per Step'] = [
            float(value) for value in pf_summary['Max_PF']
        ]
        params['PF Summary'] = pf_summary.to_dict(orient='list')

        print(
            "[PF] Limited peaking factor evaluation to "
            f"{len(pf_summary)} saved statepoints at or before EOL "
            f"({fuel_lifetime_days:.4f} days)."
        )

    except Exception as e:
        print("[PF] WARNING: compute_pin_peaking_factors failed:", e)
        pf_summary = None
        pf_per_step = None

    orig_material = depletion_2d_results_file.export_to_materials(0)

    def uranium_mass(material):
        try:
            return material.get_mass('U235') + material.get_mass('U238')
        except Exception:
            return 0.0

    operating_fuel_material = max(orig_material, key=uranium_mass)
    if uranium_mass(operating_fuel_material) <= 0.0:
        raise ValueError(
            "Could not identify the depleted fuel material from the operating "
            "depletion results."
        )
    params['_Operating Fuel Material ID'] = int(operating_fuel_material.id)
    mass_U235 = operating_fuel_material.get_mass('U235')
    mass_U238 = operating_fuel_material.get_mass('U238')

    params['keff 2D'] = [float(k) for k in keff_2d_values]
    params['keff 3D (2D corrected)'] = [float(k) for k in keff_2d_values_corrected]
    params['Depletion Time Steps'] = [float(t) for t in time_steps]

    _, keff_with_uncertainty = depletion_2d_results_file.get_keff()
    keff_2d_uncertainties = [
        float(value)
        for value in keff_with_uncertainty[:len(keff_2d_values), 1]
    ]
    if len(keff_2d_uncertainties) != len(keff_2d_values):
        raise ValueError(
            "The number of operating keff uncertainties does not match the "
            "number of corrected depletion points."
        )

    keff_corrected_uncertainties = [
        raw_uncertainty * corrected_keff / raw_keff
        for raw_uncertainty, raw_keff, corrected_keff in zip(
            keff_2d_uncertainties,
            keff_2d_values,
            keff_2d_values_corrected
        )
    ]
    params['keff 2D Uncertainty'] = keff_2d_uncertainties
    params['keff 3D (2D corrected) Uncertainty'] = [
        float(value) for value in keff_corrected_uncertainties
    ]

    # Beginning-of-life axial leakage metrics from the buckling correction model
    params['BOL Axial Non-Leakage Probability'] = bol_axial_non_leakage_probability
    params['Estimated Axial Leakage (%)'] = bol_axial_leakage_percent

    # Total leakage uses both axial and radial buckling.
    # If Core Radius is not available, these values are returned as NaN.
    params['BOL Total Non-Leakage Probability'] = bol_total_non_leakage_probability
    params['Estimated Total Leakage (%)'] = bol_total_leakage_percent

    return fuel_lifetime_days, mass_U235, mass_U238, pf_summary


def run_depletion_analysis(params):
    openmc.run()
    lattice_geometry = openmc.Geometry.from_xml()
    settings = openmc.Settings.from_xml()
    fuel_lifetime_days, mass_U235, mass_U238, pf_summary = \
        openmc_depletion(params, lattice_geometry, settings)

    params['Fuel Lifetime'] = fuel_lifetime_days
    params['Mass U235'] = mass_U235
    params['Mass U238'] = mass_U238
    params['Uranium Mass'] = (mass_U235 + mass_U238) / 1000

# Use this function and comment the previous one if you run steady state 
# def run_steady_state_analysis(params):
#     import glob

#     openmc.run()

#     statepoint_file = sorted(glob.glob("statepoint.*.h5"))[-1]

#     keff_2d, keff_3d_corrected, p_nl_axial = corrected_keff_steady_state(
#         statepoint_file,
#         params['Active Height'] + 2 * params['Axial Reflector Thickness'],
#         core_radius=params.get('Core Radius', np.nan)
#     )

#     params['keff 2D'] = [float(keff_2d)]
#     params['keff 3D (2D corrected)'] = [float(keff_3d_corrected)]
#     params['Depletion Time Steps'] = [0.0]

#     params['BOL Axial Non-Leakage Probability'] = p_nl_axial
#     params['Estimated Axial Leakage (%)'] = (1.0 - p_nl_axial) * 100.0

#     params['Fuel Lifetime'] = np.nan
#     params['Mass U235'] = np.nan
#     params['Mass U238'] = np.nan
#     params['Uranium Mass'] = np.nan


def _sum_nuclide_mass(materials, nuclide):
    total_mass = 0.0
    for material in materials:
        try:
            total_mass += material.get_mass(nuclide)
        except Exception:
            pass
    return total_mass


def _estimate_keff_crossing_time_days(time_days, keff_values):
    if len(time_days) != len(keff_values):
        raise ValueError("time_days and keff_values must have the same length.")
    if len(time_days) < 2:
        raise ValueError("At least two depletion points are required to estimate fuel lifetime.")

    for idx in range(1, len(keff_values)):
        k1 = keff_values[idx - 1]
        k2 = keff_values[idx]
        if (k1 - 1.0) * (k2 - 1.0) <= 0.0:
            if k2 == k1:
                return time_days[idx]
            return time_days[idx - 1] + (1.0 - k1) * (time_days[idx] - time_days[idx - 1]) / (k2 - k1)

    k1 = keff_values[-2]
    k2 = keff_values[-1]
    if k2 == k1:
        raise ValueError("Cannot extrapolate fuel lifetime because the last two keff values are identical.")
    return time_days[-1] + (1.0 - k2) * (time_days[-1] - time_days[-2]) / (k2 - k1)


def openmc_depletion_3d(params, lattice_geometry, settings):
    openmc.config['cross_sections'] = params['cross_sections_xml_location']

    operator = openmc.deplete.CoupledOperator(
        openmc.Model(geometry=lattice_geometry, settings=settings),
        chain_file=params['simplified_chain_thermal_xml']
    )

    if 'Burnup Steps' in params:
        burnup_steps_list_MWd_per_Kg = params['Burnup Steps']
        burnup_step = np.array(burnup_steps_list_MWd_per_Kg)
        burnup = np.diff(burnup_step, prepend=0.0)

        integrator = openmc.deplete.PredictorIntegrator(
            operator,
            burnup,
            1000000 * params['Power MWt'],
            timestep_units='MWd/kg'
        )
    elif 'Time Steps' in params:
        time_steps_list = params['Time Steps']
        power_list = [params['Power MWt'] * 1e6] * len(time_steps_list)
        integrator = openmc.deplete.CECMIntegrator(operator, time_steps_list, power_list)
    else:
        raise ValueError("3D depletion requires either 'Burnup Steps' or 'Time Steps'.")

    print("Starting 3D depletion")
    integrator.integrate()
    print("3D depletion complete")

    depletion_3d_results_file = openmc.deplete.Results("./depletion_results.h5")
    time, keff = depletion_3d_results_file.get_keff()
    time_days = [float(t) / 86400.0 for t in time]
    keff_3d_values = [float(k) for k in keff[:, 0]]

    fuel_lifetime_days = _estimate_keff_crossing_time_days(time_days, keff_3d_values)

    original_materials = depletion_3d_results_file.export_to_materials(0)
    mass_U235 = _sum_nuclide_mass(original_materials, 'U235')
    mass_U238 = _sum_nuclide_mass(original_materials, 'U238')

    params['keff 3D'] = keff_3d_values
    params['Depletion Time Steps'] = time_days

    return fuel_lifetime_days, mass_U235, mass_U238


def run_depletion_analysis_3d(params):
    openmc.run()
    lattice_geometry = openmc.Geometry.from_xml()
    settings = openmc.Settings.from_xml()
    fuel_lifetime_days, mass_U235, mass_U238 = \
        openmc_depletion_3d(params, lattice_geometry, settings)

    params['Fuel Lifetime'] = fuel_lifetime_days
    params['Mass U235'] = mass_U235
    params['Mass U238'] = mass_U238
    params['Uranium Mass'] = (mass_U235 + mass_U238) / 1000


def monitor_heat_flux(params):
    if params['Heat Flux'] <= params['Heat Flux Criteria']:
        print("\n")
        print(f"\033[92mHeat flux: {np.round(params['Heat Flux'], 2)} MW/m^2.\033[0m")
        print("\n")
    else:
        print(f"\033[91mERROR: Heat flux is too high: {np.round(params['Heat Flux'], 2)} MW/m^2.\033[0m")


def _find_lifecycle_points(time_days, corrected_keff):
    """Locate BOL, nearest saved MOL, and the two EOL bracketing states."""
    if len(time_days) != len(corrected_keff):
        raise ValueError("Depletion times and corrected keff must have equal lengths.")

    lower_index = None
    upper_index = None
    eol_fraction = None
    eol_time = None

    for index in range(1, len(corrected_keff)):
        k_lower = corrected_keff[index - 1]
        k_upper = corrected_keff[index]
        t_lower = time_days[index - 1]
        t_upper = time_days[index]

        if k_lower == 1.0:
            lower_index = upper_index = index - 1
            eol_fraction = 0.0
            eol_time = t_lower
            break
        if k_upper == 1.0:
            lower_index = upper_index = index
            eol_fraction = 0.0
            eol_time = t_upper
            break
        if (k_lower - 1.0) * (k_upper - 1.0) < 0.0:
            lower_index = index - 1
            upper_index = index
            eol_fraction = (1.0 - k_lower) / (k_upper - k_lower)
            eol_time = t_lower + eol_fraction * (t_upper - t_lower)
            break

    if lower_index is None:
        raise ValueError(
            "Cannot define lifecycle evaluation points because corrected "
            "operating keff did not reach 1.0. Extend the burnup schedule."
        )

    mol_target_time = 0.5 * eol_time
    mol_index = int(np.argmin(np.abs(np.asarray(time_days) - mol_target_time)))

    return {
        'bol_index': 0,
        'mol_index': mol_index,
        'eol_lower_index': lower_index,
        'eol_upper_index': upper_index,
        'eol_fraction': float(eol_fraction),
        'eol_time': float(eol_time),
        'mol_target_time': float(mol_target_time),
    }


def _interpolate_value(lower, upper, fraction):
    return (1.0 - fraction) * lower + fraction * upper


def _interpolate_uncertainty(lower, upper, fraction):
    """Interpolate independent one standard deviation endpoint estimates."""
    return np.sqrt(
        ((1.0 - fraction) * lower) ** 2
        + (fraction * upper) ** 2
    )


def _lifecycle_values(values_by_index, lifecycle):
    lower_index = lifecycle['eol_lower_index']
    upper_index = lifecycle['eol_upper_index']
    fraction = lifecycle['eol_fraction']
    if lower_index == upper_index:
        eol_value = values_by_index[lower_index]
    else:
        eol_value = _interpolate_value(
            values_by_index[lower_index],
            values_by_index[upper_index],
            fraction
        )
    return [
        float(values_by_index[lifecycle['bol_index']]),
        float(values_by_index[lifecycle['mol_index']]),
        float(eol_value),
    ]


def _lifecycle_uncertainties(values_by_index, lifecycle):
    lower_index = lifecycle['eol_lower_index']
    upper_index = lifecycle['eol_upper_index']
    fraction = lifecycle['eol_fraction']
    if lower_index == upper_index:
        eol_uncertainty = values_by_index[lower_index]
    else:
        eol_uncertainty = _interpolate_uncertainty(
            values_by_index[lower_index],
            values_by_index[upper_index],
            fraction
        )
    return [
        float(values_by_index[lifecycle['bol_index']]),
        float(values_by_index[lifecycle['mol_index']]),
        float(eol_uncertainty),
    ]


def _latest_statepoint_file():
    statepoint_files = glob.glob('statepoint.*.h5')
    if not statepoint_files:
        raise FileNotFoundError("Static OpenMC run did not produce a statepoint file.")
    return max(statepoint_files, key=os.path.getmtime)


def _run_static_snapshot(
    build_openmc_model,
    params,
    materials_xml,
    temperature,
    shutdown_state,
    seed=None,
    particles=None,
    material_density_overrides=None
):
    """Run one nondepleting lifecycle snapshot and return corrected keff data."""
    original_temperature = params['Common Temperature']
    original_shutdown_state = params['Shutdown Margin Calc']
    plotting_was_present = 'plotting' in params
    original_plotting = params.get('plotting', 'N')
    seed_was_present = '_OpenMC Seed' in params
    original_seed = params.get('_OpenMC Seed')
    particles_were_present = 'Particles' in params
    original_particles = params.get('Particles')
    density_override_keys = (
        '_NaK Density Override',
        '_ZrH Density Override',
    )
    original_density_overrides = {
        key: (key in params, params.get(key))
        for key in density_override_keys
    }

    try:
        params['_Depleted Fuel Materials XML'] = materials_xml
        params['Common Temperature'] = temperature
        params['Shutdown Margin Calc'] = shutdown_state
        params['plotting'] = 'N'
        if seed is not None:
            params['_OpenMC Seed'] = int(seed)
        if particles is not None:
            params['Particles'] = int(particles)
        if material_density_overrides is not None:
            for key in density_override_keys:
                if key not in material_density_overrides:
                    raise KeyError(
                        f"Missing required density override '{key}'."
                    )
                params[key] = float(material_density_overrides[key])

        build_openmc_model(params)
        openmc.run()

        return corrected_keff_static(
            _latest_statepoint_file(),
            params['Active Height'] + 2 * params['Axial Reflector Thickness'],
            core_radius=params.get('Core Radius', np.nan)
        )
    finally:
        params.pop('_Depleted Fuel Materials XML', None)
        params['Common Temperature'] = original_temperature
        params['Shutdown Margin Calc'] = original_shutdown_state
        if plotting_was_present:
            params['plotting'] = original_plotting
        else:
            params.pop('plotting', None)
        if seed_was_present:
            params['_OpenMC Seed'] = original_seed
        else:
            params.pop('_OpenMC Seed', None)
        if particles_were_present:
            params['Particles'] = original_particles
        else:
            params.pop('Particles', None)
        for key, (was_present, original_value) in (
            original_density_overrides.items()
        ):
            if was_present:
                params[key] = original_value
            else:
                params.pop(key, None)


def _temperature_coefficient(k_base, sigma_base, k_high, sigma_high, delta_t):
    coefficient = (1.0 / k_base - 1.0 / k_high) * 1e5 / delta_t
    uncertainty = 1e5 / delta_t * np.sqrt(
        (sigma_base / k_base ** 2) ** 2
        + (sigma_high / k_high ** 2) ** 2
    )
    return float(coefficient), float(uncertainty)


def _shutdown_margin(k_shutdown, sigma_shutdown):
    margin = (1.0 / k_shutdown - 1.0) * 1e5
    uncertainty = 1e5 * sigma_shutdown / k_shutdown ** 2
    return float(margin), float(uncertainty)


def _run_lifecycle_snapshot_calculations(build_openmc_model, params):
    """Evaluate temperature coefficient and shutdown margin at BOL, MOL, EOL."""
    time_days = [float(value) for value in params['Depletion Time Steps']]
    operating_corrected_keff = [
        float(value) for value in params['keff 3D (2D corrected)']
    ]
    lifecycle = _find_lifecycle_points(time_days, operating_corrected_keff)
    selected_indices = sorted({
        lifecycle['bol_index'],
        lifecycle['mol_index'],
        lifecycle['eol_lower_index'],
        lifecycle['eol_upper_index'],
    })

    depletion_results = openmc.deplete.Results('./depletion_results.h5')
    snapshot_material_files = {}
    for index in selected_indices:
        snapshot_file = f'depleted_materials_lifecycle_{index}.xml'
        depletion_results.export_to_materials(index).export_to_xml(
            path=snapshot_file
        )
        snapshot_material_files[index] = snapshot_file

    params['Lifecycle Evaluation Times'] = [
        float(time_days[lifecycle['bol_index']]),
        float(time_days[lifecycle['mol_index']]),
        float(lifecycle['eol_time']),
    ]
    params['Lifecycle Target Times'] = [
        float(time_days[lifecycle['bol_index']]),
        float(lifecycle['mol_target_time']),
        float(lifecycle['eol_time']),
    ]

    print("\nLifecycle safety evaluation points:")
    print(f"  BOL: {params['Lifecycle Evaluation Times'][0]:.4f} days")
    print(
        f"  MOL: {params['Lifecycle Evaluation Times'][1]:.4f} days "
        f"(target {lifecycle['mol_target_time']:.4f} days)"
    )
    print(
        f"  EOL: {lifecycle['eol_time']:.4f} days, interpolated between "
        f"depletion indices {lifecycle['eol_lower_index']} and "
        f"{lifecycle['eol_upper_index']}"
    )

    base_temperature_results = {}
    high_temperature_results = {}
    shutdown_results = {}
    operating_temperature = float(params['Common Temperature'])
    temperature_particles_by_index = {}
    base_temperature_seeds = {}
    high_temperature_seeds = {}
    base_density_overrides = None
    high_density_overrides = None

    if params['Isothermal Temperature Coefficients']:
        elevated_temperature = (
            operating_temperature + float(params['Temperature Perturbation'])
        )
        base_density_overrides = _temperature_density_overrides(
            operating_temperature,
            operating_temperature
        )
        high_density_overrides = _temperature_density_overrides(
            operating_temperature,
            elevated_temperature
        )
        params['Temperature Coefficient Density Aware'] = True
        params['Temperature Coefficient Density Temperatures'] = [
            operating_temperature,
            elevated_temperature,
        ]
        params['Temperature Coefficient NaK Densities'] = [
            base_density_overrides['_NaK Density Override'],
            high_density_overrides['_NaK Density Override'],
        ]
        params['Temperature Coefficient ZrH Densities'] = [
            base_density_overrides['_ZrH Density Override'],
            high_density_overrides['_ZrH Density Override'],
        ]
        print("\nDensity-aware lifecycle temperature coefficient model:")
        print(
            f"  NaK: {base_density_overrides['_NaK Density Override']:.8f} "
            f"g/cm3 at {operating_temperature:.1f} K -> "
            f"{high_density_overrides['_NaK Density Override']:.8f} g/cm3 "
            f"at {elevated_temperature:.1f} K"
        )
        print(
            f"  ZrH: {base_density_overrides['_ZrH Density Override']:.8f} "
            f"g/cm3 at {operating_temperature:.1f} K -> "
            f"{high_density_overrides['_ZrH Density Override']:.8f} g/cm3 "
            f"at {elevated_temperature:.1f} K"
        )

    for case_number, index in enumerate(selected_indices):
        if params['Isothermal Temperature Coefficients']:
            base_seed = 104729 + 2000003 * case_number
            high_seed = 15485863 + 2000033 * case_number
            case_particles = TEMPERATURE_COEFFICIENT_PARTICLES
            base_temperature_seeds[index] = base_seed
            high_temperature_seeds[index] = high_seed
            temperature_particles_by_index[index] = case_particles

            print(
                f"\n[Lifecycle] Base-temperature ARO static case at "
                f"depletion index {index}, {operating_temperature:.1f} K, "
                f"seed {base_seed}, {case_particles} particles."
            )
            base_temperature_results[index] = _run_static_snapshot(
                build_openmc_model,
                params,
                snapshot_material_files[index],
                operating_temperature,
                shutdown_state=False,
                seed=base_seed,
                particles=case_particles,
                material_density_overrides=base_density_overrides
            )

            print(
                f"\n[Lifecycle] Elevated-temperature ARO static case at "
                f"depletion index {index}, "
                f"{operating_temperature + params['Temperature Perturbation']:.1f} K, "
                f"seed {high_seed}, {case_particles} particles."
            )
            high_temperature_results[index] = _run_static_snapshot(
                build_openmc_model,
                params,
                snapshot_material_files[index],
                operating_temperature + params['Temperature Perturbation'],
                shutdown_state=False,
                seed=high_seed,
                particles=case_particles,
                material_density_overrides=high_density_overrides
            )

        if params['Shutdown Margin Calc']:
            print(
                f"\n[Lifecycle] Cold ARI static case at depletion index "
                f"{index}."
            )
            shutdown_results[index] = _run_static_snapshot(
                build_openmc_model,
                params,
                snapshot_material_files[index],
                params['Cold Shutdown Temperature'],
                shutdown_state=True
            )

    operating_raw = [float(value) for value in params['keff 2D']]
    operating_raw_uncertainty = [
        float(value) for value in params['keff 2D Uncertainty']
    ]
    operating_corrected_uncertainty = [
        float(value)
        for value in params['keff 3D (2D corrected) Uncertainty']
    ]

    params['keff 2D ARO'] = operating_raw
    params['keff 2D ARO Uncertainty'] = operating_raw_uncertainty
    params['keff 3D (2D corrected) ARO'] = operating_corrected_keff
    params['keff 3D (2D corrected) ARO Uncertainty'] = (
        operating_corrected_uncertainty
    )

    labels = ['BOL', 'MOL', 'EOL']

    if params['Isothermal Temperature Coefficients']:
        delta_t = float(params['Temperature Perturbation'])
        params['Temperature Coefficient Particles'] = [
            temperature_particles_by_index[index]
            for index in selected_indices
        ]
        params['Temperature Coefficient Static Depletion Indices'] = list(
            selected_indices
        )
        params['Temperature Coefficient Base Seeds'] = [
            base_temperature_seeds[index] for index in selected_indices
        ]
        params['Temperature Coefficient High Seeds'] = [
            high_temperature_seeds[index] for index in selected_indices
        ]
        temp_coeff_raw = {}
        temp_coeff_raw_uncertainty = {}
        temp_coeff_corrected = {}
        temp_coeff_corrected_uncertainty = {}

        for index in selected_indices:
            base_result = base_temperature_results[index]
            high_result = high_temperature_results[index]
            raw_value, raw_uncertainty = _temperature_coefficient(
                base_result['keff_2d'],
                base_result['keff_2d_uncertainty'],
                high_result['keff_2d'],
                high_result['keff_2d_uncertainty'],
                delta_t
            )
            corrected_value, corrected_uncertainty = _temperature_coefficient(
                base_result['keff_corrected'],
                base_result['keff_corrected_uncertainty'],
                high_result['keff_corrected'],
                high_result['keff_corrected_uncertainty'],
                delta_t
            )
            temp_coeff_raw[index] = raw_value
            temp_coeff_raw_uncertainty[index] = raw_uncertainty
            temp_coeff_corrected[index] = corrected_value
            temp_coeff_corrected_uncertainty[index] = corrected_uncertainty

        raw_lifecycle = _lifecycle_values(temp_coeff_raw, lifecycle)
        raw_uncertainty_lifecycle = _lifecycle_uncertainties(
            temp_coeff_raw_uncertainty,
            lifecycle
        )
        corrected_lifecycle = _lifecycle_values(temp_coeff_corrected, lifecycle)
        corrected_uncertainty_lifecycle = _lifecycle_uncertainties(
            temp_coeff_corrected_uncertainty,
            lifecycle
        )

        base_raw = {
            index: result['keff_2d']
            for index, result in base_temperature_results.items()
        }
        base_raw_uncertainty = {
            index: result['keff_2d_uncertainty']
            for index, result in base_temperature_results.items()
        }
        base_corrected = {
            index: result['keff_corrected']
            for index, result in base_temperature_results.items()
        }
        base_corrected_uncertainty = {
            index: result['keff_corrected_uncertainty']
            for index, result in base_temperature_results.items()
        }
        high_raw = {
            index: result['keff_2d']
            for index, result in high_temperature_results.items()
        }
        high_raw_uncertainty = {
            index: result['keff_2d_uncertainty']
            for index, result in high_temperature_results.items()
        }
        high_corrected = {
            index: result['keff_corrected']
            for index, result in high_temperature_results.items()
        }
        high_corrected_uncertainty = {
            index: result['keff_corrected_uncertainty']
            for index, result in high_temperature_results.items()
        }

        params['keff 2D base temp'] = _lifecycle_values(base_raw, lifecycle)
        params['keff 2D base temp Uncertainty'] = _lifecycle_uncertainties(
            base_raw_uncertainty,
            lifecycle
        )
        params['keff 3D (2D corrected) base temp'] = _lifecycle_values(
            base_corrected,
            lifecycle
        )
        params['keff 3D (2D corrected) base temp Uncertainty'] = (
            _lifecycle_uncertainties(base_corrected_uncertainty, lifecycle)
        )
        params['keff 2D high temp'] = _lifecycle_values(high_raw, lifecycle)
        params['keff 2D high temp Uncertainty'] = _lifecycle_uncertainties(
            high_raw_uncertainty,
            lifecycle
        )
        params['keff 3D (2D corrected) high temp'] = _lifecycle_values(
            high_corrected,
            lifecycle
        )
        params['keff 3D (2D corrected) high temp Uncertainty'] = (
            _lifecycle_uncertainties(high_corrected_uncertainty, lifecycle)
        )
        params['Temp Coeff 2D Lifecycle'] = raw_lifecycle
        params['Temp Coeff 2D Lifecycle Uncertainty'] = (
            raw_uncertainty_lifecycle
        )
        params['Temp Coeff 3D (2D corrected) Lifecycle'] = corrected_lifecycle
        params['Temp Coeff 3D (2D corrected) Lifecycle Uncertainty'] = (
            corrected_uncertainty_lifecycle
        )

        base_p_nl = {
            index: result['axial_non_leakage_probability']
            for index, result in base_temperature_results.items()
        }
        high_p_nl = {
            index: result['axial_non_leakage_probability']
            for index, result in high_temperature_results.items()
        }
        params['Axial Non Leakage Probability base temp'] = (
            _lifecycle_values(base_p_nl, lifecycle)
        )
        params['Axial Non Leakage Probability high temp'] = (
            _lifecycle_values(high_p_nl, lifecycle)
        )
        params['Axial Non Leakage Probability temperature change'] = [
            float(high_value - base_value)
            for base_value, high_value in zip(
                params['Axial Non Leakage Probability base temp'],
                params['Axial Non Leakage Probability high temp']
            )
        ]

        raw_limiting_index = int(np.argmax(raw_lifecycle))
        corrected_limiting_index = int(np.argmax(corrected_lifecycle))
        params['Temp Coeff 2D'] = raw_lifecycle[raw_limiting_index]
        params['Temp Coeff 2D Uncertainty'] = (
            raw_uncertainty_lifecycle[raw_limiting_index]
        )
        params['Temp Coeff 2D Limiting Point'] = labels[raw_limiting_index]
        params['Temp Coeff 3D (2D corrected)'] = (
            corrected_lifecycle[corrected_limiting_index]
        )
        params['Temp Coeff 3D (2D corrected) Uncertainty'] = (
            corrected_uncertainty_lifecycle[corrected_limiting_index]
        )
        params['Temp Coeff 3D (2D corrected) Limiting Point'] = (
            labels[corrected_limiting_index]
        )

        print("\nIsothermal temperature coefficient lifecycle results:")
        for position, label in enumerate(labels):
            print(
                f"  {label}: 2D {raw_lifecycle[position]:.6f} +/- "
                f"{raw_uncertainty_lifecycle[position]:.6f} pcm/K; "
                f"axially corrected {corrected_lifecycle[position]:.6f} +/- "
                f"{corrected_uncertainty_lifecycle[position]:.6f} pcm/K"
            )
        print(
            "  Limiting 2D value: "
            f"{params['Temp Coeff 2D']:.6f} +/- "
            f"{params['Temp Coeff 2D Uncertainty']:.6f} pcm/K at "
            f"{params['Temp Coeff 2D Limiting Point']}"
        )
        print(
            "  Limiting axially corrected value: "
            f"{params['Temp Coeff 3D (2D corrected)']:.6f} +/- "
            f"{params['Temp Coeff 3D (2D corrected) Uncertainty']:.6f} "
            f"pcm/K at "
            f"{params['Temp Coeff 3D (2D corrected) Limiting Point']}"
        )
    else:
        params['Temp Coeff 2D'] = np.nan
        params['Temp Coeff 2D Uncertainty'] = np.nan
        params['Temp Coeff 3D (2D corrected)'] = np.nan
        params['Temp Coeff 3D (2D corrected) Uncertainty'] = np.nan

    if params['Shutdown Margin Calc']:
        sdm_raw = {}
        sdm_raw_uncertainty = {}
        sdm_corrected = {}
        sdm_corrected_uncertainty = {}

        for index in selected_indices:
            static_result = shutdown_results[index]
            raw_value, raw_uncertainty = _shutdown_margin(
                static_result['keff_2d'],
                static_result['keff_2d_uncertainty']
            )
            corrected_value, corrected_uncertainty = _shutdown_margin(
                static_result['keff_corrected'],
                static_result['keff_corrected_uncertainty']
            )
            sdm_raw[index] = raw_value
            sdm_raw_uncertainty[index] = raw_uncertainty
            sdm_corrected[index] = corrected_value
            sdm_corrected_uncertainty[index] = corrected_uncertainty

        raw_lifecycle = _lifecycle_values(sdm_raw, lifecycle)
        raw_uncertainty_lifecycle = _lifecycle_uncertainties(
            sdm_raw_uncertainty,
            lifecycle
        )
        corrected_lifecycle = _lifecycle_values(sdm_corrected, lifecycle)
        corrected_uncertainty_lifecycle = _lifecycle_uncertainties(
            sdm_corrected_uncertainty,
            lifecycle
        )

        shutdown_raw = {
            index: result['keff_2d']
            for index, result in shutdown_results.items()
        }
        shutdown_raw_uncertainty = {
            index: result['keff_2d_uncertainty']
            for index, result in shutdown_results.items()
        }
        shutdown_corrected = {
            index: result['keff_corrected']
            for index, result in shutdown_results.items()
        }
        shutdown_corrected_uncertainty = {
            index: result['keff_corrected_uncertainty']
            for index, result in shutdown_results.items()
        }

        params['keff 2D ARI'] = _lifecycle_values(shutdown_raw, lifecycle)
        params['keff 2D ARI Uncertainty'] = _lifecycle_uncertainties(
            shutdown_raw_uncertainty,
            lifecycle
        )
        params['keff 3D (2D corrected) ARI'] = _lifecycle_values(
            shutdown_corrected,
            lifecycle
        )
        params['keff 3D (2D corrected) ARI Uncertainty'] = (
            _lifecycle_uncertainties(shutdown_corrected_uncertainty, lifecycle)
        )
        params['Shutdown Margin 2D Lifecycle'] = raw_lifecycle
        params['Shutdown Margin 2D Lifecycle Uncertainty'] = (
            raw_uncertainty_lifecycle
        )
        params['Shutdown Margin 3D (2D corrected) Lifecycle'] = (
            corrected_lifecycle
        )
        params['Shutdown Margin 3D (2D corrected) Lifecycle Uncertainty'] = (
            corrected_uncertainty_lifecycle
        )

        raw_limiting_index = int(np.argmin(raw_lifecycle))
        corrected_limiting_index = int(np.argmin(corrected_lifecycle))
        params['Most Limiting Shutdown Margin 2D'] = (
            raw_lifecycle[raw_limiting_index]
        )
        params['Most Limiting Shutdown Margin 2D Uncertainty'] = (
            raw_uncertainty_lifecycle[raw_limiting_index]
        )
        params['Most Limiting Shutdown Margin 2D Point'] = (
            labels[raw_limiting_index]
        )
        params['Maximum Shutdown Margin 2D'] = float(np.max(raw_lifecycle))
        params['Most Limiting Shutdown Margin 3D (2D corrected)'] = (
            corrected_lifecycle[corrected_limiting_index]
        )
        params['Most Limiting Shutdown Margin 3D (2D corrected) Uncertainty'] = (
            corrected_uncertainty_lifecycle[corrected_limiting_index]
        )
        params['Most Limiting Shutdown Margin 3D (2D corrected) Point'] = (
            labels[corrected_limiting_index]
        )
        params['Maximum Shutdown Margin 3D (2D corrected)'] = float(
            np.max(corrected_lifecycle)
        )

        print("\nShutdown margin lifecycle results:")
        for label, value, uncertainty in zip(
            labels,
            corrected_lifecycle,
            corrected_uncertainty_lifecycle
        ):
            print(f"  {label}: {value:.3f} +/- {uncertainty:.3f} pcm")
        print(
            "  Limiting value: "
            f"{params['Most Limiting Shutdown Margin 3D (2D corrected)']:.3f} "
            f"+/- "
            f"{params['Most Limiting Shutdown Margin 3D (2D corrected) Uncertainty']:.3f} "
            f"pcm at "
            f"{params['Most Limiting Shutdown Margin 3D (2D corrected) Point']}"
        )
    else:
        params['Most Limiting Shutdown Margin 2D'] = np.nan
        params['Most Limiting Shutdown Margin 2D Uncertainty'] = np.nan
        params['Maximum Shutdown Margin 2D'] = np.nan
        params['Most Limiting Shutdown Margin 3D (2D corrected)'] = np.nan
        params['Most Limiting Shutdown Margin 3D (2D corrected) Uncertainty'] = (
            np.nan
        )
        params['Maximum Shutdown Margin 3D (2D corrected)'] = np.nan


def _run_operating_depletion_and_lifecycle(
    build_openmc_model,
    params,
    shutdown_margin_requested,
    temperature_coefficient_requested
):
    """Run the operating depletion once, then evaluate requested snapshots."""
    params['Shutdown Margin Calc'] = shutdown_margin_requested
    params['Isothermal Temperature Coefficients'] = (
        temperature_coefficient_requested
    )
    run_depletion_analysis(params)
    params['keff 2D ARO'] = list(params['keff 2D'])
    params['keff 2D ARO Uncertainty'] = list(params['keff 2D Uncertainty'])
    params['keff 3D (2D corrected) ARO'] = list(
        params['keff 3D (2D corrected)']
    )
    params['keff 3D (2D corrected) ARO Uncertainty'] = list(
        params['keff 3D (2D corrected) Uncertainty']
    )

    if (
        params['Isothermal Temperature Coefficients']
        or params['Shutdown Margin Calc']
    ):
        _run_lifecycle_snapshot_calculations(build_openmc_model, params)
    else:
        params['Temp Coeff 2D'] = np.nan
        params['Temp Coeff 2D Uncertainty'] = np.nan
        params['Temp Coeff 3D (2D corrected)'] = np.nan
        params['Temp Coeff 3D (2D corrected) Uncertainty'] = np.nan
        params['Most Limiting Shutdown Margin 2D'] = np.nan
        params['Most Limiting Shutdown Margin 2D Uncertainty'] = np.nan
        params['Maximum Shutdown Margin 2D'] = np.nan
        params['Most Limiting Shutdown Margin 3D (2D corrected)'] = np.nan
        params['Most Limiting Shutdown Margin 3D (2D corrected) Uncertainty'] = (
            np.nan
        )
        params['Maximum Shutdown Margin 3D (2D corrected)'] = np.nan


def run_openmc(build_openmc_model, heat_flux_monitor, params):
    params.setdefault('Shutdown Margin Calc', False)
    params.setdefault('Isothermal Temperature Coefficients', False)
    params.setdefault('Cold Shutdown Temperature', 300)

    original_shutdown_margin_calc = params['Shutdown Margin Calc']
    original_isothermal_temperature_coefficients = (
        params['Isothermal Temperature Coefficients']
    )
    original_common_temperature = params['Common Temperature']

    if params['Isothermal Temperature Coefficients']:
        if 'Temperature Perturbation' not in params:
            raise ValueError(
                "\n\nINPUT ERROR\n"
                "'Temperature Perturbation' is required when "
                "'Isothermal Temperature Coefficients' is True.\n"
            )
        if params['Temperature Perturbation'] <= 0.0:
            raise ValueError("'Temperature Perturbation' must be greater than zero.")

        print(
            f"Using {TEMPERATURE_COEFFICIENT_PARTICLES} particles per batch "
            "for every lifecycle temperature-coefficient snapshot. The "
            "operating depletion and cold-shutdown calculations retain their "
            "normal particle settings."
        )

    try:
        print(f"\n\nThe results/plots are saved at: {watts.Database().path}\n\n")

        params['Shutdown Margin Calc'] = False
        params['Common Temperature'] = original_common_temperature
        openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True)
        openmc_plugin(
            params,
            function=lambda: _run_operating_depletion_and_lifecycle(
                build_openmc_model,
                params,
                original_shutdown_margin_calc,
                original_isothermal_temperature_coefficients
            )
        )

    except Exception:
        print("\n\n\033[91mAn error occurred while running the OpenMC simulation:\033[0m\n\n")
        traceback.print_exc()
        raise

    finally:
        params.pop('_Depleted Fuel Materials XML', None)
        params.pop('_Operating Fuel Material ID', None)
        params['Shutdown Margin Calc'] = original_shutdown_margin_calc
        params['Isothermal Temperature Coefficients'] = (
            original_isothermal_temperature_coefficients
        )
        params['Common Temperature'] = original_common_temperature


def run_openmc_3d(build_openmc_model, heat_flux_monitor, params):
    """Run the existing explicit 3D MOUSE depletion workflow."""
    params.setdefault('Shutdown Margin Calc', False)
    params.setdefault('Isothermal Temperature Coefficients', False)
    params.setdefault('Cold Shutdown Temperature', 300)

    original_shutdown_margin_calc = params['Shutdown Margin Calc']
    original_common_temperature = params['Common Temperature']

    try:
        print(f"\n\nThe results/plots are saved at: {watts.Database().path}\n\n")

        if params['Isothermal Temperature Coefficients']:
            print(
                "[3D] WARNING: Isothermal Temperature Coefficients are not "
                "implemented for run_openmc_3d yet."
            )
        params['Temp Coeff 3D'] = np.nan

        if params['Shutdown Margin Calc']:
            params['Common Temperature'] = params['Cold Shutdown Temperature']
            params['Shutdown Margin Calc'] = True
            openmc_plugin = watts.PluginOpenMC(
                build_openmc_model,
                show_stderr=True
            )
            openmc_plugin(
                params,
                function=lambda: run_depletion_analysis_3d(params)
            )
            params['keff 3D ARI'] = params['keff 3D']

            params['Common Temperature'] = original_common_temperature
            params['Shutdown Margin Calc'] = False
            openmc_plugin = watts.PluginOpenMC(
                build_openmc_model,
                show_stderr=True
            )
            openmc_plugin(
                params,
                function=lambda: run_depletion_analysis_3d(params)
            )
            params['keff 3D ARO'] = params['keff 3D']

            sdm_3d_per_step = [
                ((1.0 - k_shutdown) / k_shutdown) * 1e5
                for k_shutdown in params['keff 3D ARI']
            ]
            params['Most Limiting Shutdown Margin 3D'] = np.min(
                sdm_3d_per_step
            )
            params['Maximum Shutdown Margin 3D'] = np.max(sdm_3d_per_step)
        else:
            params['Most Limiting Shutdown Margin 3D'] = np.nan
            params['Maximum Shutdown Margin 3D'] = np.nan

            openmc_plugin = watts.PluginOpenMC(
                build_openmc_model,
                show_stderr=True
            )
            openmc_plugin(
                params,
                function=lambda: run_depletion_analysis_3d(params)
            )
            params['keff 3D ARO'] = params['keff 3D']

    except Exception:
        print(
            "\n\n\033[91mAn error occurred while running the 3D OpenMC "
            "simulation:\033[0m\n\n"
        )
        traceback.print_exc()
        raise

    finally:
        params['Shutdown Margin Calc'] = original_shutdown_margin_calc
        params['Common Temperature'] = original_common_temperature


def cyclic_rotation(input_array, k):
    return input_array[-k:] + input_array[:-k]


def flatten_list(nested_list):
    return [item for sublist in nested_list for item in sublist]
