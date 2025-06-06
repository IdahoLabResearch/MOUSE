# Importing libraries
import openmc
from core_design.utils import *
from core_design.openmc_materials_database import collect_materials_data

"""
An OpenMC function that accepts an instance of "parameters" 
and generates the necessary XMl files
"""

def build_openmc_model_GCMR(params):

    
    # **************************************************************************************************************************
    #                                                Sec. 1 : MATERIALS
    # **************************************************************************************************************************
    materials_database = collect_materials_data(params)
    fuel = materials_database[params['fuel']]
    reflector = materials_database[params['Reflector']]
    moderator = materials_database[params['Moderator']]
    moderator_booster = materials_database[params['moderator_booster']]

    control_drum_absorber = materials_database[params['Control Drum Absorber']]
    control_drum_reflector = materials_database[params['Control Drum Reflector']]
    coolant =  materials_database[params['Coolant']]

    # **************************************************************************************************************************
    #                                                Sec. 2 : GEOMETRY: TRISO particles
    # **************************************************************************************************************************

    # # # ## TRISO particles
    fuel_pin_region = create_fuel_pin_regions_TRISO(params)
    fuel_materials = []
    for mat in params['fuel_pin_materials']:
        if mat == None:
            fuel_materials.append(None)
        else: 
            material_1 = materials_database[mat]
            fuel_materials.append(material_1)
    # Giving the user error message if the number of materials is not the same as the number of regions
    assert len(fuel_pin_region) == len(fuel_materials), "The number of regions, {len(fuel_pin_region)} should be\
        the same as the number of introduced materials, {len(fuel_materials)}"  
    
    triso_cells = create_cells(fuel_pin_region, fuel_materials) 
    triso_universe = openmc.Universe(cells=triso_cells.values())  

    if params['plotting'] == "yes":
        # plotting
        create_universe_plot(triso_universe, 
                        pin_plot_width = 2.2 * params['fuel_pin_radii'][-1],
                        num_pixels = 500, 
                        font_size = 16,
                        title = "TRISO Particle", 
                        fig_size = 8, 
                        output_file_name = "TRISO_Particle.png")



    # ## The Fuel Universe (TRISO particles with background material in between and moderator material around the TRISO)

    active_core_maxz, active_core_minz, fuel_universe, compact_triso_particles_number  = create_TRISO_particles_lattice_universe(params, triso_universe, materials_database)
    # In[7]:
                            

    # # ## coolat channels & Booster Pins & Burnable Poison
    #small coolant channels
    small_coolant_universe = create_universe_from_core_top_and_bottom_planes(params['coolant channel radius'],\
    active_core_maxz, active_core_minz, coolant , materials_database[params['matrix material']])
    
    booster_universe = create_universe_from_core_top_and_bottom_planes(params['booster radius'],\
    active_core_maxz, active_core_minz, materials_database[params['moderator_booster']] , materials_database[params['Moderator']]) 


    # # # Construct hexagonal cells surrounded by coolant channels
    # Define the boundary of the hexagonal prism with the given edge length
    hex_boundary = openmc.model.hexagonal_prism(edge_length= params['hex lattice_radius'])
    # Create an instance of HexLattice
    fuel_lattice = openmc.HexLattice()
    # Set the center of the hexagonal lattice
    fuel_lattice.center = (0., 0.)
    # Set the pitch (distance between the centers of adjacent hexagons) of the hexagonal lattice
    fuel_lattice.pitch = (params['hex lattice_radius'],)###
    
    fuel_lattice.outer = openmc.Universe(cells=[openmc.Cell(fill= materials_database[params['Moderator']])]) # inner_fill or moderator_universe
    fuel_lattice.universes =  [[small_coolant_universe]*6, [fuel_universe]]
    fuel_lattice_hex = openmc.Universe(cells=[openmc.Cell(fill=fuel_lattice, region=hex_boundary)])

    # Booster Lattice
    booster_lattice = openmc.HexLattice()
    booster_lattice.center = (0., 0.)
    booster_lattice.pitch = (params['hex lattice_radius'],)###
    booster_lattice.outer = openmc.Universe(cells=[openmc.Cell(fill= materials_database[params['Moderator']])]) 
    booster_lattice.universes = [[small_coolant_universe]*6, [booster_universe]]
    booster_lattice_hex = openmc.Universe(cells=[openmc.Cell(fill=booster_lattice, region=hex_boundary)])

    # Coolant Lattice
    coolant_lattice = openmc.HexLattice()
    coolant_lattice.center = (0., 0.)
    coolant_lattice.pitch = (params['hex lattice_radius'],)
    coolant_lattice.outer = openmc.Universe(cells=[openmc.Cell(fill= materials_database[params['Moderator']])]) 
    coolant_lattice.universes = [[small_coolant_universe]*6, [openmc.Universe(cells=[openmc.Cell(fill= materials_database[params['Moderator']])])]]
    coolant_lattice_hex = openmc.Universe(cells=[openmc.Cell(fill=coolant_lattice, region=hex_boundary)])
                            
    # **************************************************************************************************************************
    #                                                Sec. 3 : ASSEMBLY & RINGS
    # **************************************************************************************************************************

    assembly_universe, assembly_fuel_cells = create_assembly(params['assembly_rings'] , params['lattice_pitch'],\
     openmc.Universe(cells=[openmc.Cell(fill= materials_database[params['Moderator']])]),\
     fuel_lattice_hex, booster_lattice_hex, outer_ring=None, simplified_output=False)
    
    if params['plotting'] == "yes":
    # plotting 

        create_universe_plot(assembly_universe, 
                pin_plot_width = 2.2 *params['lattice_pitch'] * params['assembly_rings']  ,
                num_pixels = 5000, 
                font_size = 16,
                title = "Fuel Assembly", 
                fig_size = 8, 
                output_file_name = "Fuel Assembly.png")



    # ## Corner and Edge Assemblies
    """
    The edge and corner assemblies are special cases of the normal assembly. 
    If nothing was done, moderator rods would appear as half-rods in the outer region of the assembled core.
    To address this we have 2 options:
    1- Keep moderator pins in the outer region
    2- Remove moderator pins in the outer region
 
    1 requires defining a graphite assembly with moderator rods in part of the outer region and use this to surround the core.
    2 requires defining edge and corner assemblies that do not have moderator pins in a part of the outer region.
 
    We will use solution 2, since it would introduce a relatively good parasitic absorber (H1) between the core and the reflector.
    The easiest way to define this special assemblies is to define a special outer ring, 
    where only part of the outer ring has moderator cells, then we perform a cyclic rotation of the list of pins, 
    which will cause the assembly to effectively rotate using the outer pins as a reference.
    The defining characteristic of these special assemblies is that edge assemblies have 4 sides with moderator pins, while corner assemblies have 3 sides with moderator pins.
     
    Due to the way OpenMC numbers the pins in a HexLattice, the easiest way is to make a reference outer ring for both the edge and the corner cases based on the numbering criteria, 
    then to use a cyclic rotation of this list to bring these reference rings into their proper position to represent assemblies at the correct rotation.
    """

    # # # Corner assembly universe
    corner_ring_ref = [coolant_lattice_hex]*((params['assembly_rings']-1)*3+1) + [booster_lattice_hex]*((params['assembly_rings']-1)*3-1)
    corner_ring_1 = cyclic_rotation(corner_ring_ref, (params['assembly_rings']-1)*3)
    corner_rings = [corner_ring_1] + [cyclic_rotation(corner_ring_1, (params['assembly_rings']-1)*i) for i in range(1,6)]
    corner_assembly_universe = [create_assembly(params['assembly_rings'], params['lattice_pitch'], openmc.Universe(cells=[openmc.Cell(fill= materials_database[params['Moderator']])]),\
     fuel_lattice_hex, booster_lattice_hex, outer_ring=cr, simplified_output=True) for cr in corner_rings]

    # # Edge assembly universe
    edge_ring_ref = [coolant_lattice_hex]*((params['assembly_rings']-1)*2+1) + [booster_lattice_hex]*((params['assembly_rings']-1)*4-1)
    edge_ring_1 = cyclic_rotation(edge_ring_ref, (params['assembly_rings']-1)*4)
    edge_rings = [edge_ring_1] + [cyclic_rotation(edge_ring_1, (params['assembly_rings']-1)*i) for i in range(1,6)]
    edge_assembly_universe = [create_assembly(params['assembly_rings'], params['lattice_pitch'], openmc.Universe(cells=[openmc.Cell(fill= materials_database[params['Moderator']])]),\
     fuel_lattice_hex, booster_lattice_hex, outer_ring=er) for er in edge_rings]

    # **************************************************************************************************************************
    #                                           Sec. 4 : User-Defined Parameters (Control Drums)
    # ************************************************************************************************************************** 

    # # # ## Drum Assembly

    drums = create_drums_universe_CGMR(params, absorber_thickness = params['drum_Absorber_thickness'],
                                  drum_radius = params['Drum_Radius'],
                          control_drum_absorber_material = control_drum_absorber,
                          control_drum_reflector_material = control_drum_reflector)
    


    # **************************************************************************************************************************
    #                                           Sec. 5 : User-Defined Parameters (Core)
    # **************************************************************************************************************************                     

    
    active_core = openmc.HexLattice()
    active_core.center = (0., 0.)
    active_core.pitch = (params['assembly_ftf'],)
    active_core.outer = openmc.Universe(cells=[openmc.Cell(fill= materials_database[params['Reflector']])])  # reflector Area

    rings = [[assembly_universe]]
    assembly_number = 1
    for n in range(1,  params['core_rings']-1):
        ring_cells = 6*n
        rings.insert(0, [assembly_universe]*ring_cells)
        assembly_number += ring_cells

    rings.insert(0, flatten_list([[ca] + [ea]*( params['core_rings']-2)\
        for (ca, ea) in zip(corner_assembly_universe, edge_assembly_universe)]))
    rings.insert(0, flatten_list([[openmc.Universe(cells =\
        [openmc.Cell(fill= materials_database[params['Reflector']])])] +\
            [cd]*( params['core_rings']-1) for cd in drums]))
    params['number of drums'] = (params['core_rings']-1) * len(drums)
    active_core.universes = rings
    outer_surface = openmc.ZCylinder(r=params['core_radius'], boundary_type='vacuum')
    active_core_cell = openmc.Cell(fill=active_core, region=-outer_surface & -active_core_maxz & +active_core_minz)
    active_core_universe = openmc.Universe(cells=[active_core_cell])

    if params['plotting'] == "yes":
            create_universe_plot(active_core_universe, 
            pin_plot_width = 2.2 *params['assembly_ftf'] *  params['core_rings'] ,
            num_pixels = 500, 
            font_size = 16,
            title = "Core", 
            fig_size = 8, 
            output_file_name = "Core.png")

    # **************************************************************************************************************************
    #                                                Sec. 6 : VOLUME INFO for Depletion
    # **************************************************************************************************************************

    core_fuel_cells = assembly_number * assembly_fuel_cells
    core_compact_volume = cylinder_volume(params['lattice radius'], params['active_height']) * core_fuel_cells
    core_triso_number = core_compact_volume / params['lattice_compact_volume'] * compact_triso_particles_number
    kernel_volume = sphere_volume(params['fuel_pin_radii'][0])
    fuel_volume = core_triso_number * kernel_volume
    fuel.volume = fuel_volume
    all_materials = fuel_materials + [ fuel, reflector,  moderator, moderator_booster, control_drum_absorber, coolant, control_drum_reflector ]
    
    # removing "None" materials
    all_materials_cleaned_list = [item for item in all_materials if item is not None]
    materials = openmc.Materials(list(set(all_materials_cleaned_list)))
   
    openmc.Materials.cross_sections = params['cross_sections_xml_location']
    materials.export_to_xml()
    full_core_geometry = openmc.Geometry([active_core_cell])
    full_core_geometry.export_to_xml()

    # **************************************************************************************************************************
    #                                                Sec. 7 : SIMULATION
    # **************************************************************************************************************************

    # OpenMC simulation parameters




    batches = 100
    inactive = 10
    particles = 500

    settings_file = openmc.Settings()
    settings_file.batches = batches
    settings_file.inactive = inactive
    settings_file.particles = particles
    settings_file.output = {'tallies': True}
    settings_file.temperature = {'default': 900.0,
                                 'method': 'interpolation',
                                 'tolerance': 50.0}

    # Define a cylindrical source distribution
    r = openmc.stats.Uniform(0, params['core_radius'])
    theta = openmc.stats.Uniform(0, 2*np.pi)
    z = openmc.stats.Uniform(- params['lattice height']/2, params['lattice height']/2)
    uniform_cyl = openmc.stats.CylindricalIndependent(r, theta, z)
    src = openmc.Source(space=uniform_cyl)
    src.only_fissionable = True

    settings_file.source = src

    settings_file.export_to_xml()




















