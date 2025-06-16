import numpy as np
import openmc
import openmc.deplete
import watts
import traceback # tracing errors





def circle_area(r):
    return (np.pi) * r **2

# def cylinder_volume(r, h):
#     return circle_area(r) * h

# def sphere_volume(r):
#     return 4/3 * np.pi * r **3

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


# def calculate_number_of_core_rings(core_rings_over_one_edge):
#     return 2 * core_rings_over_one_edge * (core_rings_over_one_edge -1) +\
#         2 * sum(range(1, core_rings_over_one_edge -1)) +\
#             2*core_rings_over_one_edge-1







# def create_fuel_pin_regions_TRISO(params):
#     # Read what the used decided for the dimensions of the fuel pin
#     ## Fuel (these variables' names need to be reviewed!)
#     fuel_radii = {'kernel': params['fuel_pin_radii'][0],
#                 'buffer': params['fuel_pin_radii'][1],
#                 'layer_1': params['fuel_pin_radii'][2],
#                 'layer_2': params['fuel_pin_radii'][3],
#                 'layer_3': params['fuel_pin_radii'][4]
#     }

#     # # Creating surfaces
#     shells = [openmc.Sphere(r=r) for r in fuel_radii.values()]

    
#     region = {'kernel': -shells[0],
#               'buffer': +shells[0] & -shells[1],
#               'layer_1': +shells[1] & -shells[2],
#               'layer_2': +shells[2] & -shells[3],
#               'layer_3': +shells[3] & -shells[4]
#               }
 

#     return region





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







# def create_drums_universe_CGMR(params, absorber_thickness, drum_radius,
#                           control_drum_absorber_material,
#                           control_drum_reflector_material):

#     absorber_arc = np.pi/3
#     REFERENCE_ANGLE = 240 # This angle is a constant that puts the drum in the correct orientation in reference to the lattice geometry
#     rotation_angle = 0






#     cd_inner_shell = openmc.ZCylinder(r= drum_radius - absorber_thickness)
#     cd_outer_shell = openmc.ZCylinder(r= drum_radius)
#     cd_gap_shell = openmc.ZCylinder(r= params['drum_tube_radius'])

#     cutting_plane_1 = openmc.Plane(a=1, b=absorber_arc/2)
#     cutting_plane_2 = openmc.Plane(a=1, b=-absorber_arc/2)

#     drum_absorber = +cd_inner_shell & -cd_outer_shell & -cutting_plane_1 & -cutting_plane_2
#     drum_reflector = -cd_outer_shell & ~drum_absorber
#     drum_gap_hs = +cd_outer_shell & - cd_gap_shell
#     drum_outside = +cd_gap_shell
    
#     drum_absorber = openmc.Cell(name='drum_absorber', fill= control_drum_absorber_material, region=drum_absorber)
#     drum_reflector = openmc.Cell(name='drum_reflector', fill= control_drum_reflector_material, region=drum_reflector)
#     drum_gap = openmc.Cell(name='drum_gap', region=drum_gap_hs)
#     drum_exterior = openmc.Cell(name='drum_outside', fill= control_drum_reflector_material, region=drum_outside)

#     drum_reference = openmc.Universe(cells=(drum_reflector, drum_absorber, drum_gap, drum_exterior))

    
#     drum_cells = []
#     for r in range(0, 360, 60):
#         dc = openmc.Cell(name=f'drum{r}', fill=drum_reference)
#         dc.rotation = [0, 0, REFERENCE_ANGLE + -r + rotation_angle]
#         drum_cells.append(dc)

#     drums = [openmc.Universe(cells=(dc,)) for dc in drum_cells]  
#     return drums
















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
    print("OpenMC finished Running")
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




# def create_TRISO_particles_lattice_universe(params, triso_universe, materials_database):
#     active_fuel_top = params['lattice height']/2
#     active_fuel_bot = - params['lattice height']/2

#     # Generating TRISO particle assembly in cylindrical pin cell
#     active_core_maxz = openmc.ZPlane(z0=active_fuel_top, boundary_type='reflective')
#     active_core_minz = openmc.ZPlane(z0=active_fuel_bot, boundary_type='reflective')

#     compact_surf = openmc.ZCylinder(r=params['lattice radius'])

#     compact_region = -compact_surf & -active_core_maxz & +active_core_minz

#     packed_shells = openmc.model.pack_spheres(radius= params['fuel_pin_radii'][-1], region=compact_region, pf= params['packing factor'])

#     compact_triso_particles = [openmc.model.TRISO(params['fuel_pin_radii'][-1], fill=triso_universe, center=c) for c in packed_shells]
#     compact_triso_particles_number = len(compact_triso_particles)



#     compact_cell = openmc.Cell(region=compact_region)

#     lower_left, up_right = compact_cell.region.bounding_box
#     shape = (4, 4, 4)
#     pitch = (up_right - lower_left)/shape

#     triso_assembly = openmc.model.create_triso_lattice(compact_triso_particles, lower_left, pitch, shape, materials_database[params['matrix material']])
#     compact_cell.fill = triso_assembly

#     outer_fuel_region = +compact_surf & -active_core_maxz & +active_core_minz
#     outer_fuel_cell = openmc.Cell(fill= materials_database[params['Moderator']], region=outer_fuel_region)

#     fuel_universe = openmc.Universe(cells=[compact_cell, outer_fuel_cell])
#     return active_core_maxz, active_core_minz,  fuel_universe,  compact_triso_particles_number 


# def create_universe_from_core_top_and_bottom_planes(radius, active_core_maxz, active_core_minz, material_inside, material_outside):
#     surf = openmc.ZCylinder(r=radius)
#     cell = openmc.Cell(region=-surf & -active_core_maxz & +active_core_minz, fill=material_inside)
#     outside_cell = openmc.Cell(region=+surf & -active_core_maxz & +active_core_minz, fill=material_outside)
#     universe = openmc.Universe(cells=[cell, outside_cell])    
#     return universe


# def create_assembly(num_rings, lattice_pitch, inner_fill, fuel_pin , moderator_pin, outer_ring=None, simplified_output=True):
#     # Create a hexagonal lattice for the assembly
#     assembly = openmc.HexLattice()
#     # Set the center of the hexagonal lattice
#     assembly.center = (0., 0.)
#     # Set the pitch (distance between pin centers) of the lattice
#     assembly.pitch = (lattice_pitch,)
#     # Define the outer universe of the lattice: the inner fill material
#     assembly.outer = inner_fill
#     # Set the orientation of the hexagonal lattice
#     assembly.orientation = 'x'

#     # Initialize the rings with the first ring containing the fuel pin
#     rings = [[fuel_pin]]
#     # Initialize the count of fuel cells
#     fuel_cells = 1
#     # Loop to create the rings of fuel pins around the center
#     for n in range(1, num_rings-1):
#         ring_cells = 6*n
#         rings.insert(0, [fuel_pin]*ring_cells)
#         fuel_cells += ring_cells

#     if outer_ring:
#         rings.insert(0, outer_ring)
#     else:
#         # Create and insert an outer ring of moderator pins
#         rings.insert(0, [moderator_pin]*6*(num_rings-1))

#     assembly.universes = rings

#     assembly_boundary = openmc.model.hexagonal_prism(edge_length=lattice_pitch*(num_rings-1), orientation='x')

#     assembly_cell = openmc.Cell(fill=assembly, region=assembly_boundary)
#     assembly_universe = openmc.Universe(cells=[assembly_cell])

#     if simplified_output:
#         return assembly_universe
#     else:
#         return assembly_universe, fuel_cells

# def cyclic_rotation(input_array, k):
#     return input_array[-k:] + input_array[:-k]



# def flatten_list(nested_list):
#     return [item for sublist in nested_list for item in sublist]
    
  