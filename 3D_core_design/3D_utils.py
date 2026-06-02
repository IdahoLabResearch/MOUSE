# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
import numpy as np
import openmc
import openmc.deplete
import watts
import traceback  # print full stack traces for OpenMC failures
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from core_design.correction_factor import corrected_keff_2d, keff_3d
from core_design.peaking_factor import compute_pin_peaking_factors

import pandas
import copy


def _get_mpi_rank():
    try:
        from mpi4py import MPI
        return MPI.COMM_WORLD.Get_rank()
    except ImportError:
        return 0

def _mpi_barrier():
    try:
        from mpi4py import MPI
        MPI.COMM_WORLD.Barrier()
    except ImportError:
        pass

def _mpi_bcast(obj, root=0):
    try:
        from mpi4py import MPI
        return MPI.COMM_WORLD.bcast(obj, root=root)
    except ImportError:
        return obj
    
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
        chain_file=params['simplified_chain_thermal_xml'],
        diff_burnable_mats=True
    )

    if 'Burnup Steps' in params:
        burnup_steps_list_MWd_per_Kg = params['Burnup Steps']
        burnup_step = np.array(burnup_steps_list_MWd_per_Kg)
        burnup = np.diff(burnup_step, prepend=0.0)
        integrator = openmc.deplete.PredictorIntegrator(
            operator, burnup, 1000000 * params['Power MWt'], timestep_units='MWd/kg'
        )
    elif 'Time Steps' in params:
        time_steps_list = params['Time Steps']
        power_list = [params['Power MWt'] * 1e6] * len(time_steps_list)
        integrator = openmc.deplete.CECMIntegrator(operator, time_steps_list, power_list)

    print("Start Depletion")
    integrator.integrate()
    print("End Depletion")

    rank = _get_mpi_rank()

    if rank == 0:

        depletion_2d_results_file = openmc.deplete.Results("./depletion_results.h5")

        if params.get('3D', False):
            fuel_lifetime_days, time_steps, keff_2d_values_corrected = keff_3d(
                depletion_2d_results_file,
                params['Active Height'] + 2 * params['Axial Reflector Thickness']
            )
            keff_2d_values                    = keff_2d_values_corrected
            bol_axial_non_leakage_probability = np.nan
            bol_axial_leakage_percent         = np.nan
            bol_total_non_leakage_probability = np.nan
            bol_total_leakage_percent         = np.nan
        else:
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
            idx_max = pf_summary['Max_PF'].idxmax()
            params['Max Peaking Factor']                = pf_summary.loc[idx_max, 'Max_PF']
            params['Step with Max Peaking Factor']      = pf_summary.loc[idx_max, 'Step']
            params['Region ID with Max Peaking Factor'] = pf_summary.loc[idx_max, 'Region_ID_Max']
            params['Max Peaking Factors per Step']      = pf_summary['Max_PF'].tolist()
            params['PF Summary']                        = pf_summary.to_dict(orient='list')
        except Exception as e:
            print("[PF] WARNING: compute_pin_peaking_factors failed:", e)
            pf_summary = None
            pf_per_step = None

        orig_material = depletion_2d_results_file.export_to_materials(0)
        fissile_area = circle_area(params['Fuel Pin Radii'][params['Fuel Pin Materials'].index(params['Fuel'])]) \
                     - circle_area(params['Fuel Pin Radii'][params['Fuel Pin Materials'].index(params['Fuel']) - 1])
        volume_per_pin = fissile_area * params['Active Height']
        for mat in orig_material:
            if params['Fuel'] in mat.name:
                mat.volume = volume_per_pin
        mass_U235 = sum(mat.get_mass('U235') for mat in orig_material if params['Fuel'] in mat.name)
        mass_U238 = sum(mat.get_mass('U238') for mat in orig_material if params['Fuel'] in mat.name)

        data_k = pandas.DataFrame()
        data_k['keff 2D'] = keff_2d_values
        data_k['keff 3D (2D corrected)'] = keff_2d_values_corrected
        data_k.to_csv('./objectives_keff.csv', index_label='time')

        results = (
            fuel_lifetime_days, mass_U235, mass_U238, pf_summary,
            keff_2d_values, keff_2d_values_corrected, time_steps,
            bol_axial_non_leakage_probability, bol_axial_leakage_percent,
            bol_total_non_leakage_probability, bol_total_leakage_percent,
        )
    else:
        results = None

    _mpi_barrier()
    results = _mpi_bcast(results, root=0)
    (
        fuel_lifetime_days, mass_U235, mass_U238, pf_summary,
        keff_2d_values, keff_2d_values_corrected, time_steps,
        bol_axial_non_leakage_probability, bol_axial_leakage_percent,
        bol_total_non_leakage_probability, bol_total_leakage_percent,
    ) = results

    params['keff 2D']                            = [float(k) for k in keff_2d_values]
    params['keff 3D (2D corrected)']             = [float(k) for k in keff_2d_values_corrected]
    params['Depletion Time Steps']               = [float(t) for t in time_steps]
    params['BOL Axial Non-Leakage Probability']  = bol_axial_non_leakage_probability
    params['Estimated Axial Leakage (%)']        = bol_axial_leakage_percent
    params['BOL Total Non-Leakage Probability']  = bol_total_non_leakage_probability
    params['Estimated Total Leakage (%)']        = bol_total_leakage_percent

    return fuel_lifetime_days, mass_U235, mass_U238, pf_summary


def run_depletion_analysis(params):
    _mpi_barrier()
    lattice_geometry = openmc.Geometry.from_xml()
    settings = openmc.Settings.from_xml()
    fuel_lifetime_days, mass_U235, mass_U238, pf_summary = \
        openmc_depletion(params, lattice_geometry, settings)

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


def _run_isothermal_temperature_coefficients(build_openmc_model, params):
    """
    Run the two OpenMC cases needed for the isothermal temperature coefficient:
    1) ARO at Common Temperature + Temperature Perturbation
    2) ARO at Common Temperature

    Stores the following results in params:
      - keff 2D high temp
      - keff 3D (2D corrected) high temp
      - keff 2D ARO
      - keff 3D (2D corrected) ARO
      - Temp Coeff 2D
      - Temp Coeff 3D (2D corrected)
    """
    temp_T = copy.deepcopy(params['Common Temperature'])
    params['Common Temperature'] = temp_T + params['Temperature Perturbation']

    openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True)
    openmc_plugin(params, function=lambda: run_depletion_analysis(params))
    params['keff 2D high temp'] = params['keff 2D']
    params['keff 3D (2D corrected) high temp'] = params['keff 3D (2D corrected)']

    params['Common Temperature'] = temp_T

    openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True)
    openmc_plugin(params, function=lambda: run_depletion_analysis(params))
    params['keff 2D ARO'] = params['keff 2D']
    params['keff 3D (2D corrected) ARO'] = params['keff 3D (2D corrected)']

    params['Temp Coeff 2D'] = np.max([
        (y - x) / (y * x) / params['Temperature Perturbation'] * 1e5
        for x, y in zip(params['keff 2D ARO'], params['keff 2D high temp'])
    ])
    params['Temp Coeff 3D (2D corrected)'] = np.max([
        (y - x) / (y * x) / params['Temperature Perturbation'] * 1e5
        for x, y in zip(
            params['keff 3D (2D corrected) ARO'],
            params['keff 3D (2D corrected) high temp']
        )
    ])


def run_openmc(build_openmc_model, heat_flux_monitor, params):
 
    params.setdefault('Shutdown Margin Calc', False)
    params.setdefault('Isothermal Temperature Coefficients', False)
    params.setdefault('Cold Shutdown Temperature', 300)
    params.setdefault('3D', False)
 
    original_shutdown_margin_calc                  = params['Shutdown Margin Calc']
    original_isothermal_temperature_coefficients   = params['Isothermal Temperature Coefficients']
    original_common_temperature                    = params['Common Temperature']
 
    if params['Isothermal Temperature Coefficients']:
        if 'Temperature Perturbation' not in params:
            raise ValueError(
                "\n\n--- INPUT ERROR ---\n"
                "'Temperature Perturbation' is not defined in params.\n"
                "This parameter is required when 'Isothermal Temperature Coefficients' is True.\n"
                "Please add it to your params (e.g. 'Temperature Perturbation': 100  # Kelvin)\n"
                "Typical range: 50-300K depending on your Monte Carlo statistical noise level.\n"
            )
 
    try:
        print(f"\n\nThe results/plots are saved at: {watts.Database().path}\n\n")
 
        if params['Shutdown Margin Calc']:
 
            if params['Isothermal Temperature Coefficients']:
                params['Shutdown Margin Calc'] = False
                params['Common Temperature']   = original_common_temperature
                _run_isothermal_temperature_coefficients(build_openmc_model, params)
                params['Shutdown Margin Calc'] = True
            else:
                params['Temp Coeff 2D']                = np.nan
                params['Temp Coeff 3D (2D corrected)'] = np.nan
 
            params['Common Temperature'] = params['Cold Shutdown Temperature']
            openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True)
            openmc_plugin(params, function=lambda: run_depletion_analysis(params))
            params['keff 2D ARI']                = params['keff 2D']
            params['keff 3D (2D corrected) ARI'] = params['keff 3D (2D corrected)']
 
            params['Shutdown Margin Calc'] = False
            params['Common Temperature']   = original_common_temperature
            openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True)
            openmc_plugin(params, function=lambda: run_depletion_analysis(params))
            params['keff 2D ARO']                = params['keff 2D']
            params['keff 3D (2D corrected) ARO'] = params['keff 3D (2D corrected)']
 
            sdm_2d_per_step = [
                ((1.0 - k_s) / k_s) * 1e5
                for k_s in params['keff 2D ARI']
            ]
            sdm_3d_per_step = [
                ((1.0 - k_s) / k_s) * 1e5
                for k_s in params['keff 3D (2D corrected) ARI']
            ]
            params['Most Limiting Shutdown Margin 2D']              = np.min(sdm_2d_per_step)
            params['Maximum Shutdown Margin 2D']                    = np.max(sdm_2d_per_step)
            params['Most Limiting Shutdown Margin 3D (2D corrected)'] = np.min(sdm_3d_per_step)
            params['Maximum Shutdown Margin 3D (2D corrected)']     = np.max(sdm_3d_per_step)
 
        else:
            params['Most Limiting Shutdown Margin 2D']                = np.nan
            params['Maximum Shutdown Margin 2D']                      = np.nan
            params['Most Limiting Shutdown Margin 3D (2D corrected)'] = np.nan
            params['Maximum Shutdown Margin 3D (2D corrected)']       = np.nan
 
            if params['Isothermal Temperature Coefficients']:
                _run_isothermal_temperature_coefficients(build_openmc_model, params)
            else:
                params['Temp Coeff 2D']                = np.nan
                params['Temp Coeff 3D (2D corrected)'] = np.nan
 
                openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True)
                openmc_plugin(params, function=lambda: run_depletion_analysis(params))
                params['keff 2D ARO']                = params['keff 2D']
                params['keff 3D (2D corrected) ARO'] = params['keff 3D (2D corrected)']
 
    except Exception as e:
        print("\n\n\033[91mAn error occurred while running the OpenMC simulation:\033[0m\n\n")
        traceback.print_exc()
        raise
 
    finally:
        params['Shutdown Margin Calc']                  = original_shutdown_margin_calc
        params['Isothermal Temperature Coefficients']   = original_isothermal_temperature_coefficients
        params['Common Temperature']                    = original_common_temperature


def cyclic_rotation(input_array, k):
    return input_array[-k:] + input_array[:-k]


def flatten_list(nested_list):
    return [item for sublist in nested_list for item in sublist]