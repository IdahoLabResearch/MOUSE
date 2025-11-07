# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
import numpy as np
import openmc
import openmc.deplete
import watts
import traceback # tracing errors
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from core_design.correction_factor import corrected_keff_2d

import pandas,copy

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
    # Define potential colors for materials
    potential_colors = {
        'TRIGA_fuel': 'red',
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
        'B4C_enriched': 'aqua',
        'SiC': 'teal',
        'Graphite': 'coral',
        'buffer_graphite': 'gold',
        'PyC': 'salmon',
        'homog_TRISO': 'maroon',
        'heatpipe': 'seashell',
        'monolith_graphite': 'navy'
    }
    
    # Check for materials in the database that do not have a color specified
    undefined_colors = [mat_name for mat_name in materials_database if mat_name not in potential_colors]
    if undefined_colors:
        raise ValueError(f"\033[91m\n\nError: The following materials do not have colors specified in the potential_colors dictionary:\n {', '.join(undefined_colors)}\n\n\033[0m")
    
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

    depletion_2d_results_file = openmc.deplete.Results("./depletion_results.h5")  # Example file path
 
    fuel_lifetime_days,keff_2d_values,keff_2d_values_corrected = corrected_keff_2d(depletion_2d_results_file, params['Active Height'] + 2 * params['Axial Reflector Thickness'])
    orig_material = depletion_2d_results_file.export_to_materials(0)
    mass_U235 = orig_material[0].get_mass('U235')
    mass_U238 = orig_material[0].get_mass('U238')
    data_k = pandas.DataFrame()
    data_k['keff 2D'] = keff_2d_values
    data_k['keff 3D (2D corrected)'] = keff_2d_values_corrected
    data_k.to_csv('./objectives_keff.csv',index_label='time')

    params['keff 2D'] = keff_2d_values
    params['keff 3D (2D corrected)'] = keff_2d_values_corrected
    
    return fuel_lifetime_days, mass_U235, mass_U238


def run_depletion_analysis(params):
    openmc.run()
    lattice_geometry = openmc.Geometry.from_xml()
    settings = openmc.Settings.from_xml()
    depletion_results = openmc_depletion(params, lattice_geometry, settings)
    params['Fuel Lifetime'] = depletion_results[0]  # days
    params['Mass U235'] = depletion_results[1]  # grams
    params['Mass U238'] = depletion_results[2]  # grams
    params['Uranium Mass'] = (params['Mass U235']  + params['Mass U238']) / 1000 # Kg


def monitor_heat_flux(params):
    if params['Heat Flux'] <= params['Heat Flux Criteria']:
        print("\n")
        print(f"\033[92mHEAT FLUX is: {np.round(params['Heat Flux'],2)} MW/m^2.\033[0m")
        print("\n")

    else:
        print(f"\033[91mERROR: HIGH HEAT FLUX IS TOO HIGH: {np.round(params['Heat Flux'],2)} MW/m^2.\033[0m")   
        return "High Heat Flux"

def run_openmc(build_openmc_model, heat_flux_monitor, params):
    if heat_flux_monitor == "High Heat Flux":
        print("ERROR: HIGH HEAT FLUX")
    else:    
        try:
            print(f"\n\nThe results/plots are saved at: {watts.Database().path}\n\n")
            #openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True)  # running the LTMR Model
            #openmc_plugin(params, function=lambda: run_depletion_analysis(params)) 
            if params['SD Margin Calc']:
                if params['Isothermal Temperature Coefficients']:
                    params['SD Margin Calc'] = False
                    temp_T = copy.deepcopy(params['Common Temperature'])
                    params['Common Temperature'] = params['Common Temperature'] + 300
                    openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True)  # running the LTMR Model
                    openmc_plugin(params, function=lambda: run_depletion_analysis(params)) 

                    params['keff 2D high temp'] = params['keff 2D']
                    params['keff 3D (2D corrected) high temp'] = params['keff 3D (2D corrected)']

                    params['Common Temperature'] = temp_T

                    openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True)  # running the LTMR Model
                    openmc_plugin(params, function=lambda: run_depletion_analysis(params)) 
                    params['ITC 2D'] = [(y - x) / (y*x) / (300)*1e5 for x,y in zip(params['keff 2D'],params['keff 2D high temp'])]
                    params['ITC 3D (2D corrected)'] =  [(y - x) / (y*x) / (300)*1e5 for x,y in zip(params['keff 3D (2D corrected)'],params['keff 3D (2D corrected) high temp'])]#(params['keff 3D (2D corrected)'] - params['keff 3D (2D corrected) ARI'])*1e5
                    params['SD Margin Calc'] = True
                    
                openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True)  # running the LTMR Model
                openmc_plugin(params, function=lambda: run_depletion_analysis(params)) 
                params['keff 2D ARI'] = params['keff 2D']
                params['keff 3D (2D corrected) ARI'] = params['keff 3D (2D corrected)']
                params['SD Margin Calc'] = False
                openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True)  # running the LTMR Model
                openmc_plugin(params, function=lambda: run_depletion_analysis(params)) 
                params['SDM 2D'] = [(y - x)*1e5 for x,y in zip(params['keff 2D'],params['keff 2D ARI'])]
                params['SDM 3D (2D corrected)'] =  [(y - x)*1e5 for x,y in zip(params['keff 3D (2D corrected)'],params['keff 3D (2D corrected) ARI'])]#(params['keff 3D (2D corrected)'] - params['keff 3D (2D corrected) ARI'])*1e5
            else:
                if params['Isothermal Temperature Coefficients']:
                    temp_T = copy.deepcopy(params['Common Temperature'])
                    params['Common Temperature'] = params['Common Temperature'] + 300
                    openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True)  # running the LTMR Model
                    openmc_plugin(params, function=lambda: run_depletion_analysis(params)) 

                    params['keff 2D high temp'] = params['keff 2D']
                    params['keff 3D (2D corrected) high temp'] = params['keff 3D (2D corrected)']

                    params['Common Temperature'] = temp_T

                    openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True)  # running the LTMR Model
                    openmc_plugin(params, function=lambda: run_depletion_analysis(params)) 
                    params['ITC 2D'] = [(y - x) / (y*x) / (300)*1e5 for x,y in zip(params['keff 2D'],params['keff 2D high temp'])]
                    params['ITC 3D (2D corrected)'] =  [(y - x) / (y*x) / (300)*1e5 for x,y in zip(params['keff 3D (2D corrected)'],params['keff 3D (2D corrected) high temp'])]#(params['keff 3D (2D corrected)'] - params['keff 3D (2D corrected) ARI'])*1e5
                else:
                    openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True)  # running the LTMR Model
                    openmc_plugin(params, function=lambda: run_depletion_analysis(params)) 
        except Exception as e:
            print("\n\n\033[91mAn error occurred while running the OpenMC simulation:\033[0m\n\n")
            traceback.print_exc()


def cyclic_rotation(input_array, k):
    return input_array[-k:] + input_array[:-k]


def flatten_list(nested_list):
    return [item for sublist in nested_list for item in sublist]  