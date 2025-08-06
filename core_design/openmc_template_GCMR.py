# Importing libraries
import numpy as np
import openmc
from core_design.utils import create_universe_plot, create_cells, cyclic_rotation, flatten_list, cylinder_volume, sphere_volume
from core_design.openmc_materials_database import collect_materials_data

"""
An OpenMC function that accepts an instance of "parameters" 
and generates the necessary XMl files
"""

def build_openmc_model_GCMR(params):
    # **************************************************************************************************************************
    #                                                Sec. 0 : Helper Functions
    # **************************************************************************************************************************
    def create_fuel_pin_regions_TRISO(params):
        # Read what the used decided for the dimensions of the fuel pin
        ## Fuel (these variables' names need to be reviewed!)
        fuel_radii = {'kernel': params['Fuel Pin Radii'][0],
                    'buffer': params['Fuel Pin Radii'][1],
                    'layer_1': params['Fuel Pin Radii'][2],
                    'layer_2': params['Fuel Pin Radii'][3],
                    'layer_3': params['Fuel Pin Radii'][4]
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


    def create_TRISO_particles_lattice_universe(params, triso_universe, materials_database):
        active_fuel_top = 2
        active_fuel_bot = -2

        # Generating TRISO particle assembly in cylindrical pin cell
        active_core_maxz = openmc.ZPlane(z0 = active_fuel_top, boundary_type='reflective')
        active_core_minz = openmc.ZPlane(z0 = active_fuel_bot, boundary_type='reflective')

        compact_surf = openmc.ZCylinder(r=params['Compact Fuel Radius'])

        compact_region = -compact_surf & -active_core_maxz & +active_core_minz

        packed_shells = openmc.model.pack_spheres(radius= params['Fuel Pin Radii'][-1], region=compact_region, pf= params['Packing Fraction'])

        compact_triso_particles = [openmc.model.TRISO(params['Fuel Pin Radii'][-1], fill=triso_universe, center=c) for c in packed_shells]
        compact_triso_particles_number = len(compact_triso_particles)



        compact_cell = openmc.Cell(region=compact_region)

        lower_left, up_right = compact_cell.region.bounding_box
        shape = (4, 4, 4)
        pitch = (up_right - lower_left)/shape

        triso_assembly = openmc.model.create_triso_lattice(compact_triso_particles, lower_left, pitch, shape, materials_database[params['Matrix Material']])
        compact_cell.fill = triso_assembly

        outer_fuel_region = +compact_surf & -active_core_maxz & +active_core_minz
        outer_fuel_cell = openmc.Cell(fill= materials_database[params['Moderator']], region=outer_fuel_region)

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

      

    def create_drums_universe_CGMR(params, absorber_thickness, drum_radius,
                            control_drum_absorber_material,
                            control_drum_reflector_material):

        absorber_arc = np.pi/3
        REFERENCE_ANGLE = 240 # This angle is a constant that puts the drum in the correct orientation in reference to the lattice geometry
        rotation_angle = 0


        cd_inner_shell = openmc.ZCylinder(r= drum_radius - absorber_thickness)
        cd_outer_shell = openmc.ZCylinder(r= drum_radius)
        
        # The radius of the tube of the control drum
        params['Drum Tube Radius'] = params['Drum Radius'] +(params['Drum Radius']/ 45)  # cm

        cd_gap_shell = openmc.ZCylinder(r= params['Drum Tube Radius'] )

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

    # **************************************************************************************************************************
    #                                                Sec. 1 : MATERIALS
    # **************************************************************************************************************************
    materials_database = collect_materials_data(params)
    fuel = materials_database[params['Fuel']]
    reflector = materials_database[params['Reflector']]
    moderator = materials_database[params['Moderator']]
    moderator_booster = materials_database[params['Moderator Booster']]

    control_drum_absorber = materials_database[params['Control Drum Absorber']]
    control_drum_reflector = materials_database[params['Control Drum Reflector']]
    coolant =  materials_database[params['Coolant']]
    
    # **************************************************************************************************************************
    #                                                Sec. 2 : GEOMETRY: TRISO particles
    # **************************************************************************************************************************

    # # # ## TRISO particles
    fuel_pin_region = create_fuel_pin_regions_TRISO(params)
    fuel_materials = []
    for mat in params['Fuel Pin Materials']:
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

    if params['plotting'] == "Y":
        # plotting
        create_universe_plot(materials_database, triso_universe, 
                        plot_width = 2.2 * params['Fuel Pin Radii'][-1],
                        num_pixels = 500, 
                        font_size = 16,
                        title = "TRISO Particle", 
                        fig_size = 8, 
                        output_file_name = "TRISO_Particle.png")

    # The Fuel Universe (TRISO particles with background material in between and moderator material around the TRISO)
    active_core_maxz, active_core_minz, fuel_universe, compact_triso_particles_number  = create_TRISO_particles_lattice_universe(params, triso_universe, materials_database)
   
                            
    # # ## coolat channels & Booster Pins & Burnable Poison
    #small coolant channels
    small_coolant_universe = create_universe_from_core_top_and_bottom_planes(params['Coolant Channel Radius'],\
    active_core_maxz, active_core_minz, coolant , materials_database[params['Matrix Material']])
    
    booster_universe = create_universe_from_core_top_and_bottom_planes(params['Moderator Booster Radius'],\
    active_core_maxz, active_core_minz, materials_database[params['Moderator Booster']] , materials_database[params['Moderator']]) 


    # # # Construct hexagonal cells surrounded by coolant channels
    params['Hex Lattice Radius'] = params['Lattice Pitch'] /np.sqrt(3)
    # Define the boundary of the hexagonal prism with the given edge length
    hex_boundary = openmc.model.hexagonal_prism(edge_length= params['Hex Lattice Radius'])
    # Create an instance of HexLattice
    fuel_lattice = openmc.HexLattice()
    # Set the center of the hexagonal lattice
    fuel_lattice.center = (0., 0.)
    # Set the pitch (distance between the centers of adjacent hexagons) of the hexagonal lattice
    fuel_lattice.pitch = (params['Hex Lattice Radius'],)###
    
    fuel_lattice.outer = openmc.Universe(cells=[openmc.Cell(fill= materials_database[params['Moderator']])]) # inner_fill or moderator_universe
    fuel_lattice.universes =  [[small_coolant_universe]*6, [fuel_universe]]
    fuel_lattice_hex = openmc.Universe(cells=[openmc.Cell(fill=fuel_lattice, region=hex_boundary)])

    # Booster Lattice
    booster_lattice = openmc.HexLattice()
    booster_lattice.center = (0., 0.)
    booster_lattice.pitch = (params['Hex Lattice Radius'],)###
    booster_lattice.outer = openmc.Universe(cells=[openmc.Cell(fill= materials_database[params['Moderator']])]) 
    booster_lattice.universes = [[small_coolant_universe]*6, [booster_universe]]
    booster_lattice_hex = openmc.Universe(cells=[openmc.Cell(fill=booster_lattice, region=hex_boundary)])

    # Coolant Lattice
    coolant_lattice = openmc.HexLattice()
    coolant_lattice.center = (0., 0.)
    coolant_lattice.pitch = (params['Hex Lattice Radius'],)
    coolant_lattice.outer = openmc.Universe(cells=[openmc.Cell(fill= materials_database[params['Moderator']])]) 
    coolant_lattice.universes = [[small_coolant_universe]*6, [openmc.Universe(cells=[openmc.Cell(fill= materials_database[params['Moderator']])])]]
    coolant_lattice_hex = openmc.Universe(cells=[openmc.Cell(fill=coolant_lattice, region=hex_boundary)])
                            
    # **************************************************************************************************************************
    #                                                Sec. 3 : Fuel ASSEMBLY 
    # **************************************************************************************************************************

    assembly_universe, assembly_fuel_cells = create_assembly(params['Assembly Rings'] , params['Lattice Pitch'],\
     openmc.Universe(cells=[openmc.Cell(fill= materials_database[params['Moderator']])]),\
     fuel_lattice_hex, booster_lattice_hex, outer_ring=None, simplified_output=False)
    
    if params['plotting'] == "Y":
    # plotting 

        create_universe_plot(materials_database, assembly_universe, 
                plot_width =      2 *params['Lattice Pitch'] * params['Assembly Rings']  ,
                num_pixels = 5000, 
                font_size = 16,
                title = "Fuel Assembly", 
                fig_size = 8, 
                output_file_name = "Fuel Assembly.png")

        create_universe_plot(materials_database, assembly_universe, 
                plot_width =      0.3 * params['Lattice Pitch'] * params['Assembly Rings']  ,
                num_pixels = 5000, 
                font_size = 16,
                title = "Fuel Assembly", 
                fig_size = 8, 
                output_file_name = "Fuel Assembly (zoomed in).png")        



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
    corner_ring_ref = [coolant_lattice_hex]*((params['Assembly Rings']-1)*3+1) + [booster_lattice_hex]*((params['Assembly Rings']-1)*3-1)
    corner_ring_1 = cyclic_rotation(corner_ring_ref, (params['Assembly Rings']-1)*3)
    corner_rings = [corner_ring_1] + [cyclic_rotation(corner_ring_1, (params['Assembly Rings']-1)*i) for i in range(1,6)]
    corner_assembly_universe = [create_assembly(params['Assembly Rings'], params['Lattice Pitch'], openmc.Universe(cells=[openmc.Cell(fill= materials_database[params['Moderator']])]),\
     fuel_lattice_hex, booster_lattice_hex, outer_ring=cr, simplified_output=True) for cr in corner_rings]

    # # Edge assembly universe
    edge_ring_ref = [coolant_lattice_hex]*((params['Assembly Rings']-1)*2+1) + [booster_lattice_hex]*((params['Assembly Rings']-1)*4-1)
    edge_ring_1 = cyclic_rotation(edge_ring_ref, (params['Assembly Rings']-1)*4)
    edge_rings = [edge_ring_1] + [cyclic_rotation(edge_ring_1, (params['Assembly Rings']-1)*i) for i in range(1,6)]
    edge_assembly_universe = [create_assembly(params['Assembly Rings'], params['Lattice Pitch'], openmc.Universe(cells=[openmc.Cell(fill= materials_database[params['Moderator']])]),\
     fuel_lattice_hex, booster_lattice_hex, outer_ring=er) for er in edge_rings]

    # **************************************************************************************************************************
    #                                           Sec. 4 : User-Defined Parameters (Control Drums)
    # ************************************************************************************************************************** 

    # # # ## Drum Assembly

    drums = create_drums_universe_CGMR(params, absorber_thickness = params['Drum Absorber Thickness'],
                                  drum_radius = params['Drum Radius'],
                          control_drum_absorber_material = control_drum_absorber,
                          control_drum_reflector_material = control_drum_reflector)
    


    # **************************************************************************************************************************
    #                                           Sec. 5 : User-Defined Parameters (Core)
    # **************************************************************************************************************************                     


    
    active_core = openmc.HexLattice()
    active_core.center = (0., 0.)  
    # the height of the hexagonal of one fuel assembly
    active_core.pitch = (params['Assembly FTF'],)
    active_core.outer = openmc.Universe(cells=[openmc.Cell(fill= materials_database[params['Reflector']])])  # reflector Area

    rings = [[assembly_universe]]
    assembly_number = 1
    for n in range(1,  params['Core Rings']-1):
        ring_cells = 6*n
        rings.insert(0, [assembly_universe]*ring_cells)
        assembly_number += ring_cells

    rings.insert(0, flatten_list([[ca] + [ea]*( params['Core Rings']-2)\
        for (ca, ea) in zip(corner_assembly_universe, edge_assembly_universe)]))
    rings.insert(0, flatten_list([[openmc.Universe(cells =\
        [openmc.Cell(fill= materials_database[params['Reflector']])])] +\
            [cd]*( params['Core Rings']-1) for cd in drums]))
    params['number of drums'] = (params['Core Rings']-1) * len(drums)
    active_core.universes = rings
    outer_surface = openmc.ZCylinder(r=params['Core Radius'], boundary_type='vacuum')
    active_core_cell = openmc.Cell(fill=active_core, region=-outer_surface & -active_core_maxz & +active_core_minz)
    active_core_universe = openmc.Universe(cells=[active_core_cell])

    if params['plotting'] == "Y":
            create_universe_plot(materials_database, active_core_universe, 
            plot_width = 2.2 *params['Assembly FTF'] *  params['Core Rings'] ,
            num_pixels = 500, 
            font_size = 16,
            title = "Core", 
            fig_size = 8, 
            output_file_name = "Core.png")

    if params['plotting'] == "Y":
            create_universe_plot(materials_database, active_core_universe, 
            plot_width = 0.5 * params['Assembly FTF'] *  params['Core Rings'] ,
            num_pixels = 500, 
            font_size = 16,
            title = "Core", 
            fig_size = 8, 
            output_file_name = "Core (zoomed in).png")        

    # **************************************************************************************************************************
    #                                                Sec. 6 : VOLUME INFO for Depletion
    # **************************************************************************************************************************
    # The volume of a compact fuel volume defined before to have a a height of 4
    params['Lattice Compact Volume'] =  cylinder_volume(params['Compact Fuel Radius'], 4)

    core_fuel_cells = assembly_number * assembly_fuel_cells
    core_compact_volume = cylinder_volume(params['Compact Fuel Radius'], params['Active Height']) * core_fuel_cells
    core_triso_number = core_compact_volume / params['Lattice Compact Volume'] * compact_triso_particles_number
    kernel_volume = sphere_volume(params['Fuel Pin Radii'][0])
    fuel_volume = core_triso_number * kernel_volume
    fuel.volume = fuel_volume
    all_materials = fuel_materials + [ fuel, reflector,  moderator, moderator_booster, control_drum_absorber, coolant, control_drum_reflector ]
    
    # removing "None" materials
    all_materials_cleaned_list = [item for item in all_materials if item is not None]
    materials = openmc.Materials(list(set(all_materials_cleaned_list)))
   
    openmc.Materials.cross_sections = params['cross_sections_xml_location']
    materials.export_to_xml()
    core = openmc.Universe(cells=[active_core_cell])
    core_geometry = openmc.Geometry(core)
    core_geometry.export_to_xml()
    
    # # **************************************************************************************************************************
    # #                                                Sec. 1.7 : TALLIES
    # # **************************************************************************************************************************

    tallies_file = openmc.Tallies()

    group_edges = np.array([1e-5, 6.7e-2, 3.2e-1, 1, 4, 9.88, 4.81e1, 4.54e2, 4.9e4, 1.83e5, 8.21e5, 4e7])# 11 energy groups from HPMR report table no.5 in ev
    groups = openmc.mgxs.EnergyGroups(group_edges)

    mgxs_lib = openmc.mgxs.Library(core_geometry)
    mgxs_lib.energy_groups = groups
    mgxs_lib.legendre_order     = 1
    mgxs_lib.mgxs_types = ['absorption', 'diffusion-coefficient', 'transport', 'scatter matrix', 'total', 'scatter']
    mgxs_lib.domain_type = 'universe'
    mgxs_lib.domains = [core]
    mgxs_lib.build_library()
    mgxs_lib.add_to_tallies_file(tallies_file, merge=False)
    tallies_file.export_to_xml()
    # **************************************************************************************************************************
    #                                                Sec. 7 : SIMULATION
    # **************************************************************************************************************************

    # OpenMC simulation parameters




    batches = 100
    inactive = 10
    particles = 1000

    settings_file = openmc.Settings()
    settings_file.batches = batches
    settings_file.inactive = inactive
    settings_file.particles = particles
    settings_file.output = {'tallies': True}
    settings_file.temperature = {'default': params['Common Temperature'],
                                 'method': 'interpolation',
                                 'tolerance': 50.0}

    # Define a cylindrical source distribution
    r = openmc.stats.Uniform(0, params['Core Radius'])
    theta = openmc.stats.Uniform(0, 2*np.pi)
    z = openmc.stats.Uniform(- 2, 2)
    uniform_cyl = openmc.stats.CylindricalIndependent(r, theta, z)
    src = openmc.Source(space=uniform_cyl)
    src.only_fissionable = True

    settings_file.source = src

    settings_file.export_to_xml()