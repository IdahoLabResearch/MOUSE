import numpy as np
import openmc
import openmc.deplete
import watts



def create_cells(regions:dict, materials:list)->dict:
    return {key:openmc.Cell(name=key, fill=mat, region=value) for (key,value), mat in zip(regions.items(), materials)}


def circle_area(r):
    return (np.pi) * r **2

def cylinder_volume(r, h):
    return circle_area(r) * h

def sphere_volume(r):
    return 4/3 * np.pi * r **3

def circle_perimeter(r):
    return 2*(np.pi)*r


def cylinder_radial_shell(r, h):
    return circle_perimeter(r) * h

def hexagon_area(side):
    return (3 * np.sqrt(3))/2 * side**2

def calculate_lattice_radius(cladding_radius, pin_gap_distance, number_of_assembly_rings):
    lattice_radius = (cladding_radius * 2 + pin_gap_distance) * number_of_assembly_rings
    return lattice_radius

def calculate_heat_flux(cladding_radius, lattice_height, rings, power_MW_th):
    fuel_number =  sum(r.count("FUEL") for r in rings)
    heat_transfer_surface = cylinder_radial_shell(cladding_radius, lattice_height ) * fuel_number  * 1e-4 # convert from cm2 to m2
    
    return power_MW_th/heat_transfer_surface # MW/m^2



def hex_area(ftf):
    apothem = ftf/2
    side = apothem / (np.sqrt(3)/2)
    perimeter = 6*side
    area = apothem * perimeter /2
    return area


def calculate_reflector_mass_LTMR(hex_area, core_radius, area_of_all_drums, drum_height):
    # I assume for now that the drums are always fully inside the reflector
    
    area_reflector = 3.14 * core_radius * core_radius - hex_area  - area_of_all_drums # cm2
    vol_reflector = area_reflector * drum_height # cm^3
    mass_reflector = vol_reflector * 3.02/1000 # mass in Kg
    return mass_reflector

def calculate_number_of_core_rings(core_rings_over_one_edge):
    return 2 * core_rings_over_one_edge * (core_rings_over_one_edge -1) +\
        2 * sum(range(1, core_rings_over_one_edge -1)) +\
            2*core_rings_over_one_edge-1


def create_fuel_pin_regions(params):
    # Read what the used decided for the dimensions of the fuel pin
    ## Fuel (these variables' names need to be reviewed!)
    fuel_radii = {'insert': params['fuel_pin_radii'][0],
                'gap1': params['fuel_pin_radii'][1],
                'fuel_meat': params['fuel_pin_radii'][2],
                'gap2': params['fuel_pin_radii'][3],
                'cladding': params['fuel_pin_radii'][4]
    }

    # # Creating surfaces
    shells = [openmc.ZCylinder(r=r) for r in fuel_radii.values()]
    
    region = {'insert': -shells[0],
            'gap1': +shells[0] & -shells[1],
            'fuel_meat': +shells[1] & -shells[2],
            'gap2': +shells[2] & -shells[3],
            'cladding': +shells[3] & -shells[4],
            'coolant': +shells[4]
    }

    return region




def create_fuel_pin_regions_TRISO(params):
    # Read what the used decided for the dimensions of the fuel pin
    ## Fuel (these variables' names need to be reviewed!)
    fuel_radii = {'kernel': params['fuel_pin_radii'][0],
                'buffer': params['fuel_pin_radii'][1],
                'layer_1': params['fuel_pin_radii'][2],
                'layer_2': params['fuel_pin_radii'][3],
                'layer_3': params['fuel_pin_radii'][4]
    }

    # # Creating surfaces
    shells = [openmc.Sphere(r=r) for r in fuel_radii.values()]

    
    region = {'kernel': -shells[0],
              'buffer': +shells[0] & -shells[1],
              'layer_1': +shells[1] & -shells[2],
              'layer_2': +shells[2] & -shells[3],
              'layer_3': +shells[3] & -shells[4]
              }
 

    return region





def create_moderator_pin_regions(params):
        ## Reflector
    moderator_radii = {'moderator': params['moderator_pin_radii'][0],
                    'cladding': params['moderator_pin_radii'][1]
                    }
    shells = [openmc.ZCylinder(r=r) for r in moderator_radii.values()]
    
    region = {'moderator': -shells[0],
            'cladding': +shells[0] & -shells[1],
            'coolant': +shells[1]
    }
    return region



def create_drums_universe(absorber_thickness, drum_radius,
                          control_drum_absorber_material,
                          control_drum_reflector_material):

    absorber_arc = np.pi/3
    REFERENCE_ANGLE = 0
    rotation_angle = 180

    cd_inner_shell = openmc.ZCylinder(r= drum_radius - absorber_thickness)
    cd_outer_shell = openmc.ZCylinder(r= drum_radius)

    cutting_plane_1 = openmc.Plane(a=1, b=absorber_arc/2)
    cutting_plane_2 = openmc.Plane(a=1, b=-absorber_arc/2)

    drum_absorber = +cd_inner_shell & -cd_outer_shell & -cutting_plane_1 & -cutting_plane_2
    drum_reflector = -cd_outer_shell & ~drum_absorber
    drum_outside = +cd_outer_shell
    drum_absorber = openmc.Cell(name='drum_absorber', fill= control_drum_absorber_material, region=drum_absorber)
    drum_reflector = openmc.Cell(name='drum_reflector', fill= control_drum_reflector_material, region=drum_reflector)
    drum_exterior = openmc.Cell(name='drum_outside', region=drum_outside)

    drum_reference = openmc.Universe(cells=(drum_reflector, drum_absorber, drum_exterior))
    
    drum_cells = []
    for r in range(0, 360, 60):
        dc = openmc.Cell(name=f'drum{r}', fill=drum_reference)
        dc.rotation = [0, 0, REFERENCE_ANGLE + r + rotation_angle]
        drum_cells.append(dc)

    drums = [openmc.Universe(cells=(dc,)) for dc in drum_cells]    
    return drums

def create_drums_universe_CGMR(params, absorber_thickness, drum_radius,
                          control_drum_absorber_material,
                          control_drum_reflector_material):

    absorber_arc = np.pi/3
    REFERENCE_ANGLE = 240 # This angle is a constant that puts the drum in the correct orientation in reference to the lattice geometry
    rotation_angle = 0






    cd_inner_shell = openmc.ZCylinder(r= drum_radius - absorber_thickness)
    cd_outer_shell = openmc.ZCylinder(r= drum_radius)
    cd_gap_shell = openmc.ZCylinder(r= params['drum_tube_radius'])

    cutting_plane_1 = openmc.Plane(a=1, b=absorber_arc/2)
    cutting_plane_2 = openmc.Plane(a=1, b=-absorber_arc/2)

    drum_absorber = +cd_inner_shell & -cd_outer_shell & -cutting_plane_1 & -cutting_plane_2
    drum_reflector = -cd_outer_shell & ~drum_absorber
    drum_gap_hs = +cd_outer_shell & - cd_gap_shell
    drum_outside = +cd_gap_shell
    
    drum_absorber = openmc.Cell(name='drum_absorber', fill= control_drum_absorber_material, region=drum_absorber)
    drum_reflector = openmc.Cell(name='drum_reflector', fill= control_drum_reflector_material, region=drum_reflector)
    drum_gap = openmc.Cell(name='drum_gap', region=drum_gap_hs)
    drum_exterior = openmc.Cell(name='drum_outside', fill= control_drum_reflector_material, region=drum_outside)

    drum_reference = openmc.Universe(cells=(drum_reflector, drum_absorber, drum_gap, drum_exterior))

    
    drum_cells = []
    for r in range(0, 360, 60):
        dc = openmc.Cell(name=f'drum{r}', fill=drum_reference)
        dc.rotation = [0, 0, REFERENCE_ANGLE + -r + rotation_angle]
        drum_cells.append(dc)

    drums = [openmc.Universe(cells=(dc,)) for dc in drum_cells]  
    return drums




def create_assembly_universe(params, fuel_pin_universe, moderator_pin_universe, pin_pitch, reflector_material, outer_universe):

    assembly = openmc.HexLattice()

    assembly.center = (0., 0.)
    assembly.pitch = (pin_pitch,)
    assembly.outer = outer_universe # coolant universe is the outer universe probably
    rings = params['rings']
 
    for i in range(len(rings)):
        for j in range(len(rings[i])):
            if rings[i][j] == 'FUEL':
                rings[i][j] = fuel_pin_universe
            elif rings[i][j] == 'MODERATOR':
                rings[i][j] = moderator_pin_universe
    
    rings = rings[-params['assembly_rings']:]
    assembly.universes = rings
    
    assembly_boundary = openmc.model.hexagonal_prism(edge_length=\
        pin_pitch*(params['assembly_rings']-1)+pin_pitch*0.6, corner_radius = (params['fuel_pin_radii'])[-1] \
            + params["pin_gap_distance"])

    fuel_assembly_cell = openmc.Cell(fill=assembly, region=assembly_boundary)
    reflector_cell = openmc.Cell(fill = reflector_material, region=~assembly_boundary)

    assembly_universe = openmc.Universe(cells=[fuel_assembly_cell, reflector_cell])
    return assembly_universe



def create_control_drums_positions( params, number_of_drums):
        
    # Placement of drums happen by tracing a line through the core apothems
    # then 2 drums are place after each apothem by deviating from this line
    # by a deviation angle
    sector = (60/180) * np.pi
    
    deviation =  np.pi* (params["deviation angle between drums"]  / 180) 
    positions = []
    for s in range(number_of_drums):
        positions.append(s*sector-deviation)
        positions.append(s*sector+deviation)

    return positions 





def create_core_geometry(params, drums, drums_positions, assembly_universe ):
    
    # The distance between the center of the control drum and the center of the hexagonal lattice
    cd_distance = 0.86602540378 * params['lattice_radius'] + params['drum_tube_radius'] 
    drum_tube_radius = params['drum_tube_radius']
    drum_universes = []
    for d in drums:
        drum_universes.append(d)
        drum_universes.append(d)

    drum_shells = []
    drum_cells = []
    for p, du in zip(drums_positions, drum_universes):
        x, y = np.cos(p)*cd_distance, np.sin(p)*cd_distance
        drum_shell = openmc.ZCylinder(x0=x, y0=y, r=drum_tube_radius)
        drum_shells.append(drum_shell)
        drum_cell = openmc.Cell(fill=du, region=-drum_shell)
        drum_cell.translation = (x, y, 0)  # translates the center of the drum universe to match the cylinder position
        drum_cells.append(drum_cell)
    
    drums_outside = +drum_shells[0]
    for d in drum_shells[1:]:
        drums_outside = drums_outside & +d

    outer_surface = openmc.ZCylinder(r=params['core_radius'] , boundary_type='vacuum')

    core_cell = openmc.Cell(fill= assembly_universe, region=-outer_surface & drums_outside)

    core_geometry = openmc.Geometry([core_cell] + drum_cells)  
    return core_geometry 


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
    burnup_steps_list_MWd_per_Kg = params['burnup_steps_MWd_per_Kg']
    
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

    drum_count = 0
    # Iterate through all cells in the geometry
    for cell in lattice_geometry.get_all_cells().values():
   
        #Check if the cell is a drum (you need to define how to identify a drum)
        # For example, if drums are identified by a specific material or cell id:
        if "drum" in cell.name.lower():  # Assuming drum cells have "drum" in their name
            drum_count += 1

    # Get the number of drums in the geometry
    num_drums = drum_count
    params['number of drums'] = num_drums 

    settings = openmc.Settings.from_xml()
    # depletion_results = openmc_depletion(params, lattice_geometry, settings)
    # params['fuel_lifetime_days'] = depletion_results[0]  # days
    # params['mass_U235'] = depletion_results[1]  # grams
    # params['mass_U238'] = depletion_results[2]  # grams

def monitor_heat_flux(params, openmc_plugin ):
    if params['heat_flux'] <= params['heat_flux_criteria']:
        print(f"\033[92mHEAT FLUX is: {params['heat_flux']} MW/m^2.\033[0m")
        openmc_plugin(params, function=lambda: run_depletion_analysis(params))
    else:
        print(f"\033[91mWARNING: HIGH HEAT FLUX: {params['heat_flux']} MW/m^2.\033[0m")   

def list_to_dict(var_list):
    return {str(var): var for var in var_list}        




def create_TRISO_particles_lattice_universe(params, triso_universe, materials_database):
    active_fuel_top = params['lattice height']/2
    active_fuel_bot = - params['lattice height']/2

    # Generating TRISO particle assembly in cylindrical pin cell
    active_core_maxz = openmc.ZPlane(z0=active_fuel_top, boundary_type='reflective')
    active_core_minz = openmc.ZPlane(z0=active_fuel_bot, boundary_type='reflective')

    compact_surf = openmc.ZCylinder(r=params['lattice radius'])

    compact_region = -compact_surf & -active_core_maxz & +active_core_minz

    packed_shells = openmc.model.pack_spheres(radius= params['fuel_pin_radii'][-1], region=compact_region, pf= params['packing factor'])

    compact_triso_particles = [openmc.model.TRISO(params['fuel_pin_radii'][-1], fill=triso_universe, center=c) for c in packed_shells]
    compact_triso_particles_number = len(compact_triso_particles)



    compact_cell = openmc.Cell(region=compact_region)

    lower_left, up_right = compact_cell.region.bounding_box
    shape = (4, 4, 4)
    pitch = (up_right - lower_left)/shape

    triso_assembly = openmc.model.create_triso_lattice(compact_triso_particles, lower_left, pitch, shape, materials_database[params['matrix material']])
    compact_cell.fill = triso_assembly

    outer_fuel_region = +compact_surf & -active_core_maxz & +active_core_minz
    outer_fuel_cell = openmc.Cell(fill= materials_database[params['moderator']], region=outer_fuel_region)

    fuel_universe = openmc.Universe(cells=[compact_cell, outer_fuel_cell])
    return active_core_maxz, active_core_minz,  fuel_universe,  compact_triso_particles_number 


def create_universe_from_core_top_and_bottom_planes(radius, active_core_maxz, active_core_minz, material_inside, material_outside):
    surf = openmc.ZCylinder(r=radius)
    cell = openmc.Cell(region=-surf & -active_core_maxz & +active_core_minz, fill=material_inside)
    outside_cell = openmc.Cell(region=+surf & -active_core_maxz & +active_core_minz, fill=material_outside)
    universe = openmc.Universe(cells=[cell, outside_cell])    
    return universe


def create_assembly(num_rings, lattice_pitch, inner_fill, fuel_pin , moderator_pin, outer_ring=None, simplified_output=True):
    # Create a hexagonal lattice for the assembly
    assembly = openmc.HexLattice()
    # Set the center of the hexagonal lattice
    assembly.center = (0., 0.)
    # Set the pitch (distance between pin centers) of the lattice
    assembly.pitch = (lattice_pitch,)
    # Define the outer universe of the lattice: the inner fill material
    assembly.outer = inner_fill
    # Set the orientation of the hexagonal lattice
    assembly.orientation = 'x'

    # Initialize the rings with the first ring containing the fuel pin
    rings = [[fuel_pin]]
    # Initialize the count of fuel cells
    fuel_cells = 1
    # Loop to create the rings of fuel pins around the center
    for n in range(1, num_rings-1):
        ring_cells = 6*n
        rings.insert(0, [fuel_pin]*ring_cells)
        fuel_cells += ring_cells

    if outer_ring:
        rings.insert(0, outer_ring)
    else:
        # Create and insert an outer ring of moderator pins
        rings.insert(0, [moderator_pin]*6*(num_rings-1))

    assembly.universes = rings

    assembly_boundary = openmc.model.hexagonal_prism(edge_length=lattice_pitch*(num_rings-1), orientation='x')

    assembly_cell = openmc.Cell(fill=assembly, region=assembly_boundary)
    assembly_universe = openmc.Universe(cells=[assembly_cell])

    if simplified_output:
        return assembly_universe
    else:
        return assembly_universe, fuel_cells

def cyclic_rotation(input_array, k):
    return input_array[-k:] + input_array[:-k]



def flatten_list(nested_list):
    return [item for sublist in nested_list for item in sublist]
    
  