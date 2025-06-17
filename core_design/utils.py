import numpy as np
import openmc
import openmc.deplete
import watts
import traceback # tracing errors





def circle_area(r):
    return (np.pi) * r **2

def cylinder_volume(r, h):
    return circle_area(r) * h

def sphere_volume(r):
    return 4/3 * np.pi * r **3

def circle_perimeter(r):
    return 2*(np.pi)*r


def cylinder_radial_shell(r, h):
    # calculating the outer area of a cylinder
    return circle_perimeter(r) * h

# def hexagon_area(side):
#     return (3 * np.sqrt(3))/2 * side**2

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

# def hex_area(ftf):
#     apothem = ftf/2
#     side = apothem / (np.sqrt(3)/2)
#     perimeter = 6*side
#     area = apothem * perimeter /2
#     return area


def calculate_reflector_mass_LTMR(params):
    hex_area =  2.598 * params['Lattice Radius'] * params['Lattice Radius']
    core_radius = params['Core Radius']
    area_of_all_drums = params['All Drums Area'] 
    drum_height = params['Drum Height']
    # I assume for now that the drums are always fully inside the reflector
    
    area_reflector = 3.14 * core_radius * core_radius - hex_area  - area_of_all_drums # cm2
    vol_reflector = area_reflector * drum_height # cm^3
    mass_reflector = vol_reflector * 3.02/1000 # mass in Kg
    params['Reflector Mass'] = mass_reflector 


def calculate_number_of_core_rings(core_rings_over_one_edge):
    return 2 * core_rings_over_one_edge * (core_rings_over_one_edge -1) +\
        2 * sum(range(1, core_rings_over_one_edge -1)) +\
            2*core_rings_over_one_edge-1













# # def create_moderator_pin_regions(params):
# #         ## Reflector
# #     moderator_radii = {'moderator': params['moderator_pin_radii'][0],
# #                     'cladding': params['moderator_pin_radii'][1]
# #                     }
# #     shells = [openmc.ZCylinder(r=r) for r in moderator_radii.values()]
    
# #     region = {'moderator': -shells[0],
# #             'cladding': +shells[0] & -shells[1],
# #             'coolant': +shells[1]
# #     }
# #     return region























def create_universe_plot(pin_universe, pin_plot_width, num_pixels, font_size,\
    title, fig_size, output_file_name):
    
    pin_plot = pin_universe.plot(width = ( pin_plot_width, pin_plot_width),
                                 pixels=(num_pixels, num_pixels), color_by='material')
    pin_plot.set_xlabel('x [cm]', fontsize= font_size)
    pin_plot.set_ylabel('y [cm]', fontsize= font_size)
    pin_plot.set_title(title, fontsize= font_size)

    pin_plot.tick_params(axis='x', labelsize= font_size)
    pin_plot.tick_params(axis='y', labelsize= font_size)
   
    # Retrieve the figure from the Axes object
    fig = pin_plot.figure
    fig.set_size_inches(fig_size, fig_size) 
    fig.tight_layout()
    # Save the figure to a file
    fig.savefig(output_file_name) 

    
def openmc_depletion(params, lattice_geometry, settings):
    
    openmc.config['cross_sections'] = params['cross_sections_xml_location'] 
    
    # depletion operator, performing transport simulations, is created using the geometry and settings xml files
    operator = openmc.deplete.CoupledOperator(openmc.Model(geometry=lattice_geometry, 
            settings=settings),
            chain_file= params['simplified_chain_thermal_xml'])
    burnup_steps_list_MWd_per_Kg = params['Burnup Steps']
    
    #MWd/kg (MW-day of energy deposited per kilogram of initial heavy metal)
    burnup_step = np.array(burnup_steps_list_MWd_per_Kg)     #np.array([0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 140.0]) 
    burnup = np.diff (burnup_step, prepend =0.0 )
    
    # Deplete using a first-order predictor algorithm.
    integrator = openmc.deplete.PredictorIntegrator(operator, burnup,
                                                    1000000 * params['Power MWt'] , timestep_units='MWd/kg')
    print("Start Depletion")
    integrator.integrate()
    print("End Depletion")
    results = openmc.deplete.Results("./depletion_results.h5")
    time, k = results.get_keff()

    time /= (24 * 60 * 60)  # convert back to days from second
    for j, ki in enumerate(k):
        if ki[0] < 1.0:
            i = j-1
            break
    fuel_lifetime_days = (time[i])   
    orig_material = results.export_to_materials(0)

    mass_U235 = orig_material[0].get_mass('U235')
    mass_U238 = orig_material[0].get_mass('U238')
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
            openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True)  # running the LTMR Model
            openmc_plugin(params, function=lambda: run_depletion_analysis(params)) 

        except Exception as e:
            print("\n\n\033[91mAn error occurred while running the OpenMC simulation:\033[0m\n\n")
            traceback.print_exc()



# def list_to_dict(var_list):
#     return {str(var): var for var in var_list}        











def cyclic_rotation(input_array, k):
    return input_array[-k:] + input_array[:-k]



def flatten_list(nested_list):
    return [item for sublist in nested_list for item in sublist]
    
  