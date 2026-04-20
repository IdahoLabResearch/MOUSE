# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
import numpy as np
import openmc
import openmc.deplete
import watts
import traceback # tracing errors
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from core_design.correction_factor import corrected_keff_2d
from core_design.peaking_factor import compute_pin_peaking_factors

import pandas,copy

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
    return (np.pi) * r **2

def cylinder_volume(r, h):
    return circle_area(r) * h

def sphere_volume(r):
    return 4/3 * np.pi * r **3

def circle_perimeter(r):
    return 2*(np.pi)*r

def sphere_area(radius):
    area = 4 * np.pi * (radius ** 2)
    return area


def cylinder_radial_shell(r, h):
    # calculating the outer area of a cylinder
    return circle_perimeter(r) * h

def calculate_lattice_radius(params):
    pin_diameter = 2 * params['Fuel Pin Radii'][-1]
    lattice_radius = pin_diameter * params['Number of Rings per Assembly']  +\
     params["Pin Gap Distance"] * (params['Number of Rings per Assembly'] - 1)              
    return lattice_radius

def calculate_heat_flux(params):
    fuel_number =  params['Fuel Pin Count']  
    heat_transfer_surface = cylinder_radial_shell(params['Fuel Pin Radii'][-1], params['Active Height'] ) * fuel_number  * 1e-4 # convert from cm2 to m2
    
    return params['Power MWt']/heat_transfer_surface # MW/m^2

def calculate_pins_in_assembly(params, pin_type):
     # Get the rings configuration from the parameters
    rings = params['Pins Arrangement']
    # Only keep the last 'Number of Rings per Assembly' number of rings as specified in the parameters
    rings = rings[-params['Number of Rings per Assembly']:]
    return sum(row.count(pin_type) for row in rings )

def create_cells(regions:dict, materials:list)->dict:
    return {key:openmc.Cell(name=key, fill=mat, region=value) for (key,value), mat in zip(regions.items(), materials)}


def calculate_number_of_rings(rings_over_one_edge):
    # total number of rings given the rings over one edge
    return 2 * rings_over_one_edge * (rings_over_one_edge -1) +\
        2 * sum(range(1, rings_over_one_edge -1)) +\
            2*rings_over_one_edge-1
 
def calculate_number_fuel_elements_hpmr(rings_over_one_edge):
    total_number_of_rings = calculate_number_of_rings(rings_over_one_edge)
    number_of_heatpipe_pins = calculate_number_of_rings(int(np.ceil(rings_over_one_edge/2)))
    return total_number_of_rings - number_of_heatpipe_pins

def number_of_heatpipes_hmpr(params):
    tot_rings_per_assembly = calculate_number_of_rings(params['Number of Rings per Assembly'])
    params['Number of Heatpipes per Assembly'] = tot_rings_per_assembly - params['Fuel Pin Count per Assembly']
    params['Number of Heatpipes'] = params['Number of Heatpipes per Assembly'] * params['Fuel Assemblies Count'] 

def calculate_total_number_of_TRISO_particles(params):
    compact_fuel_vol = cylinder_volume(params['Compact Fuel Radius'], params['Active Height'])
    one_particle_volume = sphere_volume(params['Fuel Pin Radii'][-1])
    number_of_particles_per_compact_fuel_vol = np.floor(params['Packing Fraction'] *compact_fuel_vol / one_particle_volume) 
    params['Number Of TRISO Particles Per Compact Fuel'] =number_of_particles_per_compact_fuel_vol
    total_number_of_particles = number_of_particles_per_compact_fuel_vol * calculate_number_of_rings(params['Assembly Rings']) *\
     calculate_number_of_rings(params['Core Rings'])
    params['Total Number of TRISO Particles'] = total_number_of_particles
    return total_number_of_particles

def calculate_heat_flux_TRISO(params):
    number_of_triso_particles = calculate_total_number_of_TRISO_particles(params)
    total_area_triso = number_of_triso_particles * sphere_area(params['Fuel Pin Radii'][0]) * 1e-4 #  # cm^2 to m^2
    heat_flux = params['Power MWt'] / total_area_triso
    return heat_flux


def create_universe_plot(materials_database, universe, plot_width, num_pixels, font_size, title, fig_size, output_file_name):
    import matplotlib.colors as mcolors

    # -----------------------------------------------------------------------------------------
    # Known material colors — each material must have a UNIQUE color.
    # If you add a new material to the materials database, add a corresponding
    # unique color here. If you forget, the code will automatically assign one
    # and warn you (see auto-assignment logic below).
    #
    # IMPORTANT: Do not use color aliases that resolve to the same hex value.
    # For example, 'cyan' and 'aqua' are identical in matplotlib (#00FFFF).
    # Always verify new colors are visually distinct from existing ones.
    # -----------------------------------------------------------------------------------------
    potential_colors = {
        'TRIGA_fuel':       'red',
        'ZrH':              'yellow',
        'UO2':              'green',
        'UC':               'purple',
        'UCO':              'orange',
        'UN':               'cyan',          # #00FFFF
        'YHx':              'magenta',
        'NaK':              'blue',
        'Helium':           'grey',          # #808080
        'Be':               'brown',
        'BeO':              'pink',
        'Zr':               'lime',
        'SS304':            'black',
        'B4C_natural':      'olive',
        'B4C_enriched':     'deepskyblue',   # was 'aqua' — FIXED: 'aqua'=='cyan', now unique
        'SiC':              'teal',
        'Graphite':         'coral',
        'buffer_graphite':  'gold',
        'PyC':              'salmon',
        'homog_TRISO':      'maroon',
        'heatpipe':         'seashell',
        'monolith_graphite':'navy',
        'UZr':              'darkred',
        'ZrC':              'slategray',
        'MgO':              'lightyellow',
        'WB':               'darkgray',
        'W2B':              'dimgray',
        'WB4':              'lightgray',
        'WC':               'silver',        # was 'gray' — FIXED: 'gray'=='grey', now unique
    }

    # -----------------------------------------------------------------------------------------
    # Auto-color assignment for materials not in potential_colors.
    # If a material exists in the database but has no assigned color, the code
    # automatically picks a unique color from a pool of distinct CSS4 colors,
    # avoiding all colors already in use. A warning is printed so the developer
    # knows to add a permanent color entry above.
    # -----------------------------------------------------------------------------------------

    # Build a pool of candidate colors — all CSS4 named colors not already used
    used_colors = set(mcolors.to_hex(c) for c in potential_colors.values())
    color_pool = [
        name for name, hex_val in mcolors.CSS4_COLORS.items()
        if mcolors.to_hex(hex_val) not in used_colors
    ]

    # Check for any materials in the database that are missing from potential_colors
    for mat_name in materials_database:
        if mat_name not in potential_colors:
            if not color_pool:
                raise ValueError(
                    f"Could not auto-assign a color for material '{mat_name}': "
                    f"no unique colors remaining in the CSS4 pool. "
                    f"Please manually add a color for this material in potential_colors."
                )
            # Assign the first available unique color from the pool
            auto_color = color_pool.pop(0)
            potential_colors[mat_name] = auto_color
            # Update used_colors so the next auto-assignment is also unique
            used_colors.add(mcolors.to_hex(auto_color))
            print(
                f"\033[93m--- WARNING: Material '{mat_name}' does not have a color specified "
                f"in potential_colors. Automatically assigned color: '{auto_color}'. "
                f"Please add a permanent entry for this material in the potential_colors "
                f"dictionary in create_universe_plot (utils.py) to suppress this warning.\033[0m"
            )

    # Create the plot_colors dictionary only with existing materials
    colors = {materials_database[mat_name]: color for mat_name, color in potential_colors.items() if mat_name in materials_database}

    # Create the plot
    universe_plot = universe.plot(width=(plot_width, plot_width),
                                  pixels=(num_pixels, num_pixels), color_by='material', colors=colors)
    universe_plot.set_xlabel('x [cm]', fontsize=font_size)
    universe_plot.set_ylabel('y [cm]', fontsize=font_size)
    universe_plot.set_title(title, fontsize=font_size)

    universe_plot.tick_params(axis='x', labelsize=font_size)
    universe_plot.tick_params(axis='y', labelsize=font_size)
   
    # Retrieve the figure from the Axes object
    fig = universe_plot.figure
    fig.set_size_inches(fig_size, fig_size)

    # Extract the materials present in the universe
    universe_materials = [cell.fill for cell in universe.get_all_cells().values()]
    used_materials = set(universe_materials)
    
    # Create legend patches for only the used materials
    legend_patches = [mpatches.Patch(color=color, label=mat_name) 
                      for mat_name, color in potential_colors.items() 
                      if mat_name in materials_database and materials_database[mat_name] in used_materials]
    # Add the legend to the plot, positioning it outside the plot area
    universe_plot.legend(handles=legend_patches, fontsize=font_size, loc='center left', bbox_to_anchor=(1, 0.5))
    # Save the figure to a file
    fig.savefig(output_file_name, bbox_inches='tight')




    
def openmc_depletion(params, lattice_geometry, settings):
    
    openmc.config['cross_sections'] = params['cross_sections_xml_location'] 
    # depletion operator, performing transport simulations, is created using the geometry and settings xml files
    operator = openmc.deplete.CoupledOperator(openmc.Model(geometry=lattice_geometry, 
            settings=settings),
            chain_file= params['simplified_chain_thermal_xml'])
    if 'Burnup Steps' in params:
        burnup_steps_list_MWd_per_Kg = params['Burnup Steps']
    
        #MWd/kg (MW-day of energy deposited per kilogram of initial heavy metal)
        burnup_step = np.array(burnup_steps_list_MWd_per_Kg)      
        burnup = np.diff (burnup_step, prepend =0.0 )
        
        # Deplete using a first-order predictor algorithm.
        integrator = openmc.deplete.PredictorIntegrator(operator, burnup,
                                                    1000000 * params['Power MWt'] , timestep_units='MWd/kg')
    elif 'Time Steps' in params:                                               
        time_steps_list = params['Time Steps'] 
        power_list = [params['Power MWt'] * 1e6] * len(time_steps_list)
        integrator = openmc.deplete.CECMIntegrator(operator, time_steps_list, power_list)

    print("Start Depletion")
    integrator.integrate()
    print("End Depletion")

    # Only rank 0 handles post-processing, then broadcasts results
    rank = _get_mpi_rank()

    if rank == 0:
        depletion_2d_results_file = openmc.deplete.Results("./depletion_results.h5")
        fuel_lifetime_days, keff_2d_values, keff_2d_values_corrected = corrected_keff_2d(
            depletion_2d_results_file, 
            params['Active Height'] + 2 * params['Axial Reflector Thickness']
        )

        try:
            pf_summary, pf_per_step = compute_pin_peaking_factors(".")
            idx_max = pf_summary['Max_PF'].idxmax()
            params['Max Peaking Factor']             = pf_summary.loc[idx_max, 'Max_PF']
            params['Step with Max Peaking Factor']   = pf_summary.loc[idx_max, 'Step']
            params['Rod ID with Max Peaking Factor'] = pf_summary.loc[idx_max, 'Rod_ID_Max']
            params['Max Peaking Factors per Step']   = pf_summary['Max_PF'].tolist()
            params['PF Summary']                     = pf_summary.to_dict(orient='list')
        except Exception as e:
            print("[PF] WARNING: compute_pin_peaking_factors failed:", e)
            pf_summary = None
            pf_per_step = None

        orig_material = depletion_2d_results_file.export_to_materials(0)
        mass_U235 = orig_material[0].get_mass('U235')
        mass_U238 = orig_material[0].get_mass('U238')

        data_k = pandas.DataFrame()
        data_k['keff 2D'] = keff_2d_values
        data_k['keff 3D (2D corrected)'] = keff_2d_values_corrected
        data_k.to_csv('./objectives_keff.csv', index_label='time')

        results = (fuel_lifetime_days, mass_U235, mass_U238, pf_summary, 
                keff_2d_values, keff_2d_values_corrected)
    else:
        results = None

    # Broadcast results to all ranks
    _mpi_barrier()
    results = _mpi_bcast(results, root=0)
    fuel_lifetime_days, mass_U235, mass_U238, pf_summary, keff_2d_values, keff_2d_values_corrected = results

    params['keff 2D'] = keff_2d_values
    params['keff 3D (2D corrected)'] = keff_2d_values_corrected

    return fuel_lifetime_days, mass_U235, mass_U238, pf_summary


def run_depletion_analysis(params):
    # openmc.run()
    lattice_geometry = openmc.Geometry.from_xml()
    settings = openmc.Settings.from_xml()
    fuel_lifetime_days, mass_U235, mass_U238, pf_summary = \
        openmc_depletion(params, lattice_geometry, settings)

    params['Fuel Lifetime'] = fuel_lifetime_days  # days
    params['Mass U235'] = mass_U235  # grams
    params['Mass U238'] = mass_U238  # grams
    params['Uranium Mass'] = (mass_U235 + mass_U238) / 1000  # kg


def monitor_heat_flux(params):
    if params['Heat Flux'] <= params['Heat Flux Criteria']:
        print("\n")
        print(f"\033[92mHEAT FLUX is: {np.round(params['Heat Flux'],2)} MW/m^2.\033[0m")
        print("\n")

    else:
        print(f"\033[91mERROR: HIGH HEAT FLUX IS TOO HIGH: {np.round(params['Heat Flux'],2)} MW/m^2.\033[0m")   
        return "High Heat Flux"

def run_openmc(build_openmc_model, heat_flux_monitor, params):

    params.setdefault('SD Margin Calc', False)
    params.setdefault('Isothermal Temperature Coefficients', False)
    original_sd_margin_calc = params['SD Margin Calc']
    original_itc = params['Isothermal Temperature Coefficients']

    if params['Isothermal Temperature Coefficients']:
        if 'Temperature Perturbation' not in params.keys():
            raise ValueError(
                "\n\n--- INPUT ERROR ---\n"
                "'Temperature Perturbation' is not defined in params.\n"
                "This parameter is required when 'Isothermal Temperature Coefficients' is True.\n"
                "Please add it to your params (e.g. 'Temperature Perturbation': 100  # Kelvin)\n"
                "Typical range: 50-300K depending on your Monte Carlo statistical noise level.\n"
            )

    if heat_flux_monitor == "High Heat Flux":
        print("ERROR: HIGH HEAT FLUX")
    else:
        try:
            print(f"\n\nThe results/plots are saved at: {watts.Database().path}\n\n")
            if params['SD Margin Calc']:
                if params['Isothermal Temperature Coefficients']:
                    params['SD Margin Calc'] = False
                    temp_T = copy.deepcopy(params['Common Temperature'])
                    params['Common Temperature'] = params['Common Temperature'] + params['Temperature Perturbation']
                    openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True)  
                    openmc_plugin(params, function=lambda: run_depletion_analysis(params))

                    params['keff 2D high temp'] = params['keff 2D']
                    params['keff 3D (2D corrected) high temp'] = params['keff 3D (2D corrected)']

                    params['Common Temperature'] = temp_T

                    openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True)  
                    openmc_plugin(params, function=lambda: run_depletion_analysis(params))
                    params['Temp Coeff 2D'] = np.max([(y - x) / (y*x) / (params['Temperature Perturbation'])*1e5 for x,y in zip(params['keff 2D'],params['keff 2D high temp'])])
                    params['Temp Coeff 3D (2D corrected)'] = np.max([(y - x) / (y*x) / (params['Temperature Perturbation'])*1e5 for x,y in zip(params['keff 3D (2D corrected)'],params['keff 3D (2D corrected) high temp'])])
                    params['SD Margin Calc'] = True
                else:
                    params['Temp Coeff 2D'] = np.nan
                    params['Temp Coeff 3D (2D corrected)'] = np.nan

                openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True)  
                openmc_plugin(params, function=lambda: run_depletion_analysis(params))
                params['keff 2D ARI'] = params['keff 2D']
                params['keff 3D (2D corrected) ARI'] = params['keff 3D (2D corrected)']
                params['SD Margin Calc'] = False
                openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True)  
                openmc_plugin(params, function=lambda: run_depletion_analysis(params))
                params['SDM 2D'] = np.max([(y - x)*1e5 for x,y in zip(params['keff 2D'],params['keff 2D ARI'])])
                params['SDM 3D (2D corrected)'] = np.max([(y - x)*1e5 for x,y in zip(params['keff 3D (2D corrected)'],params['keff 3D (2D corrected) ARI'])])
            else:
                params['SDM 2D'] = np.nan
                params['SDM 3D (2D corrected)'] = np.nan
                if params['Isothermal Temperature Coefficients']:
                    temp_T = copy.deepcopy(params['Common Temperature'])
                    params['Common Temperature'] = params['Common Temperature'] + params['Temperature Perturbation']
                    openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True)  
                    openmc_plugin(params, function=lambda: run_depletion_analysis(params))

                    params['keff 2D high temp'] = params['keff 2D']
                    params['keff 3D (2D corrected) high temp'] = params['keff 3D (2D corrected)']

                    params['Common Temperature'] = temp_T

                    openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True)  
                    openmc_plugin(params, function=lambda: run_depletion_analysis(params))
                    params['Temp Coeff 2D'] = np.max([(y - x) / (y*x) / (params['Temperature Perturbation'])*1e5 for x,y in zip(params['keff 2D'],params['keff 2D high temp'])])
                    params['Temp Coeff 3D (2D corrected)'] = np.max([(y - x) / (y*x) / (params['Temperature Perturbation'])*1e5 for x,y in zip(params['keff 3D (2D corrected)'],params['keff 3D (2D corrected) high temp'])])
                else:
                    params['Temp Coeff 2D'] = np.nan
                    params['Temp Coeff 3D (2D corrected)'] = np.nan
                    openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True) 
                    openmc_plugin(params, function=lambda: run_depletion_analysis(params))

        except Exception as e:
            print("\n\n\033[91mAn error occurred while running the OpenMC simulation:\033[0m\n\n")
            traceback.print_exc()
            raise  # fix: re-raise so the outer try/except in the main script can catch it

        finally:
            params['SD Margin Calc'] = original_sd_margin_calc
            params['Isothermal Temperature Coefficients'] = original_itc

def cyclic_rotation(input_array, k):
    return input_array[-k:] + input_array[:-k]


def flatten_list(nested_list):
    return [item for sublist in nested_list for item in sublist]  