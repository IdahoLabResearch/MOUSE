# Importing libraries
import openmc
import numpy as np
from core_design.openmc_materials_database import collect_materials_data
from core_design.utils import create_universe_plot, circle_area



# **************************************************************************************************************************
#                                                Sec. 0 : Helper Functions
# **************************************************************************************************************************

"""
Helper functions are smaller, reusable functions ar defined to perform specific tasks,
and then used later to simplify and organize the code
"""

def create_pin_regions(params, pin_type):
    """
    Creating the pin regions
    @ In, params, dict, The parameters that are used to “fill in” input files with placeholders.
    @ In, pin_type, str, The type of pin ('moderator' or 'fuel').
    @ out, regions, dict, regions of the specified pin
    """

    if pin_type == 'moderator':
        # Extract the radii values for different regions of the moderator pin from the input parameters
        pin_radii = {'moderator': params['Moderator Pin Radii'][0],
                     'cladding': params['Moderator Pin Radii'][1]
        }
        region_keys = ['moderator', 'cladding', 'coolant']

    elif pin_type == 'fuel':
        # Extract the radii values for different regions of the fuel pin from the input parameters
        pin_radii = {'insert': params['Fuel Pin Radii'][0],
                     'gap1': params['Fuel Pin Radii'][1],
                     'fuel_meat': params['Fuel Pin Radii'][2],
                     'gap2': params['Fuel Pin Radii'][3],
                     'cladding': params['Fuel Pin Radii'][4]
        }
        region_keys = ['insert', 'gap1', 'fuel_meat', 'gap2', 'cladding', 'coolant']
    
    else:
        raise ValueError("Invalid pin type. Must be 'moderator' or 'fuel'.")

    # Creating surfaces
    # Create cylindrical surfaces for each of the specified radii
    shells = [openmc.ZCylinder(r=r) for r in pin_radii.values()]

    # Define the regions within the pin using the created cylindrical surfaces
    regions = {}
    for i, key in enumerate(region_keys[:-1]):
        if i == 0:
            regions[key] = -shells[i]
        else:
            regions[key] = +shells[i-1] & -shells[i]
    regions[region_keys[-1]] = +shells[-1]

    return regions



def create_drums_universe(params, control_drum_absorber_material, control_drum_reflector_material):
    """
    Creating the universe of control drums
    @ In, params, dict, The parameters that are used to “fill in” input files with placeholders.
    @ In, control_drum_absorber_material, openmc.material.Material, The control drum material (absorber)
    @ In, control_drum_reflector_material, openmc.material.Material, The control drum material (reflector)
    @ out,list, list of universes for control drums
    """

    
    absorber_thickness = params['Drum Absorber Thickness']
    drum_radius = params['Drum Radius'] 

    # Define the arc angle for the absorber and reference angle for rotation
    absorber_arc = np.pi/3
    REFERENCE_ANGLE = 0
    rotation_angle = 180

    # Create cylindrical surfaces for the inner and outer shells of the control drum
    cd_inner_shell = openmc.ZCylinder(r= drum_radius - absorber_thickness)
    cd_outer_shell = openmc.ZCylinder(r= drum_radius)

    # Define planes to cut the absorber arc segment
    cutting_plane_1 = openmc.Plane(a=1, b=absorber_arc/2)
    cutting_plane_2 = openmc.Plane(a=1, b=-absorber_arc/2)

    # Define the regions for the absorber and reflector
    drum_absorber_region = +cd_inner_shell & -cd_outer_shell & -cutting_plane_1 & -cutting_plane_2
    drum_reflector_region = -cd_outer_shell & ~drum_absorber_region
    drum_outside_region = +cd_outer_shell

    # Create cells for the absorber, reflector, and exterior
    drum_absorber = openmc.Cell(name='drum_absorber', fill=control_drum_absorber_material, region=drum_absorber_region)
    drum_reflector = openmc.Cell(name='drum_reflector', fill=control_drum_reflector_material, region=drum_reflector_region)
    drum_exterior = openmc.Cell(name='drum_outside', region=drum_outside_region)

    # Create the reference universe containing the drum cells
    drum_reference = openmc.Universe(cells=(drum_reflector, drum_absorber, drum_exterior))
    
    # List to store drum cells with different rotations
    drum_cells = []
    for r in range(0, 360, 60):
        # Create a new cell for each drum with a specific rotation
        dc = openmc.Cell(name=f'drum{r}', fill=drum_reference)
        dc.rotation = [0, 0, REFERENCE_ANGLE + r + rotation_angle]
        drum_cells.append(dc)

    # Create universes for each drum cell
    drums = [openmc.Universe(cells=(dc,)) for dc in drum_cells]
    
    return drums




def create_assembly_universe(params, fuel_pin_universe, moderator_pin_universe, pin_pitch, reflector_material, outer_coolant_universe):
    """
    Creating the universe of the fuel assembly
    @ In, params, dict, The parameters that are used to “fill in” input files with placeholders.
    @ In, fuel_pin_universe, openmc.universe.Universe, 
    @ In, moderator_pin_universe, openmc.universe.Universe, 
    @ In, pin_pitch, float, the center-to-center distance between adjacent fuel/moderator pins 
    @ In, reflector_material, openmc.material.Material, the material of the outer radial reflector
    @ In, outer_coolant_universe, openmc.universe.Universe, the openmc universe of the coolant in the assembly
    @ out, assembly_universe, openmc.universe.Universe, the fuel assembly universe


    'openmc.material.Material'> <class 'openmc.universe.Universe'> <class 'openmc.universe.Universe'>
    """
    
    # Create a hexagonal lattice for the fuel assembly
    assembly = openmc.HexLattice()

    # Set the center of the hexagonal lattice to the origin (0, 0)
    assembly.center = (0., 0.)
    
    # Set the pitch of the lattice cells (distance between centers of adjacent pins)
    assembly.pitch = (pin_pitch,)
    
    # Define the outer universe, which is likely the coolant region surrounding the fuel assembly
    assembly.outer = outer_coolant_universe
    
    # Get the rings configuration from the parameters
    rings = params['Pins Arrangement']

    # Only keep the last 'Number of Rings per Assembly' number of rings as specified in the parameters
    rings = rings[-params['Number of Rings per Assembly']:]
    
    # Loop through each ring and replace 'FUEL' and 'MODERATOR' strings with corresponding universes
    for i in range(len(rings)):
        for j in range(len(rings[i])):
            if rings[i][j] == 'FUEL':
                rings[i][j] = fuel_pin_universe
            elif rings[i][j] == 'MODERATOR':
                rings[i][j] = moderator_pin_universe
    

    # Assign the resulting ring configuration to the lattice universes
    assembly.universes = rings
    
    # Define the boundary of the assembly using a hexagonal prism
    assembly_boundary = openmc.model.hexagonal_prism(
        edge_length=pin_pitch * (params['Number of Rings per Assembly'] - 1) + pin_pitch * 0.6, 
        corner_radius=(params['Fuel Pin Radii'])[-1] + params["Pin Gap Distance"]
    )

    # Create a cell for the fuel assembly within the defined boundary
    fuel_assembly_cell = openmc.Cell(fill=assembly, region=assembly_boundary)
    
    # Create a cell for the reflector material surrounding the assembly
    reflector_cell = openmc.Cell(fill=reflector_material, region=~assembly_boundary)

    # Combine the fuel assembly cell and reflector cell into a universe
    assembly_universe = openmc.Universe(cells=[fuel_assembly_cell, reflector_cell])
    

    # Return the created assembly universe
    return assembly_universe


def create_control_drums_positions(number_of_drums):
    """
    Creating the positions of the control drums around the reactor
    @ In, number_of_drums, int, number of drums around the reactor
    @ out, positions, list, list of control drum positions
    """
        
    # Placement of drums happen by tracing a line through the core apothems
    # then 2 drums are place after each apothem by deviating from this line
    # by a deviation angle
    sector = (60/180) * np.pi

    deviation_angle_between_drums = 12.86 # degrees
    
    deviation =  np.pi* (deviation_angle_between_drums  / 180) 
    positions = []
    for s in range(number_of_drums):
        positions.append(s*sector-deviation)
        positions.append(s*sector+deviation)

    return positions 

def create_core_geometry(params, drums, drums_positions, assembly_universe):

    """
    Creating the geometry for the entire core
    @ In, params, dict, The parameters that are used to “fill in” input files with placeholders.
    @ In, drums, list, universes of the drums
    @ In, drums_positions, list of control drum positions
    @ In, assembly_universe, openmc.universe.Universet, fuel assembly universe
    @ out, core_geometry, openmc.geometry.Geometry, the geometry of the entire core
    """

    params['Drum Tube Radius'] = params['Drum Radius'] + params['Drum Radius'] / 90 # cm

    # The distance between the center of the control drum and the center of the hexagonal lattice
    cd_distance = 0.86602540378 * params['Lattice Radius'] + params['Drum Tube Radius']
    drum_tube_radius = params['Drum Tube Radius']
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

    outer_surface = openmc.ZCylinder(r=params['Core Radius'] , boundary_type='vacuum')

    core_cell = openmc.Cell(fill= assembly_universe, region=-outer_surface & drums_outside)

    core_geometry = openmc.Geometry([core_cell] + drum_cells) 
    return core_geometry 


def create_cells(regions:dict, materials:list)->dict:
    return {key:openmc.Cell(name=key, fill=mat, region=value) for (key,value), mat in zip(regions.items(), materials)}






# **************************************************************************************************************************
#                                                Sec. 1 : OpenMC Model
# **************************************************************************************************************************

"""
An OpenMC function that accepts an instance of "parameters" 
and generates the necessary XMl files
"""

def build_openmc_model_LTMR(params):
    """
    OpenMC Model
    @ In, params, watts.parameters.Parameters, The parameters that are used to “fill in” 
    input files with placeholders. params mostly behaves like a Python dictionary with a few extra capabilities
    """
    # **************************************************************************************************************************
    #                                                Sec. 1.1 : MATERIALS
    # **************************************************************************************************************************
   
    # reading all the materials properties
    materials_database = collect_materials_data(params)

    # reading the materials properties for the fuel, coolant, refelctor, control drum (the drum includes two materials: absoerber and reflector)   
    fuel = materials_database[params['Fuel']]
    coolant = materials_database[params['Coolant']]
    reflector = materials_database[params['Reflector']]
    control_drum_absorber = materials_database[params['Control Drum Absorber']]
    control_drum_reflector = materials_database[params['Control Drum Reflector']]
    

    # **************************************************************************************************************************
    #                                                Sec. 1.2 : GEOMETRY: Fuel Pins, Moderator Pins, Coolant
    # **************************************************************************************************************************
    
    # Creating Fuel Pins regions
    fuel_pin_regions = create_pin_regions(params, 'fuel')
    
    # Creating fuel pin materials
    fuel_materials = []
    for mat in params['Fuel Pin Materials']:
        if mat == None:
            fuel_materials.append(None)
        else: 
            material_1 = materials_database[mat]
            fuel_materials.append(material_1)
    fuel_materials.append(coolant)

     # Giving the user error message if the number of materials is not the same as the number of regions
    assert len(fuel_pin_regions) == len(fuel_materials), "The number of regions, {len(fuel_pin_regions)} should be\
        the same as the number of introduced materials, {len(fuel_materials)}"
    
    # creating the fuel pin universe
    fuel_cells = create_cells(fuel_pin_regions, fuel_materials)
    fuel_pin_universe = openmc.Universe(cells=fuel_cells.values())

    if params['plotting'] == "Y":
        # plotting
        create_universe_plot(fuel_pin_universe, 
                        pin_plot_width = 2.2 * params['Fuel Pin Radii'][-1],
                        num_pixels = 500, 
                        font_size = 16,
                        title = "Fuel Pin Universe", 
                        fig_size = 8, 
                        output_file_name = "fuel_pin_universe.png")

    
    # Creating Moderator Pin regions
    moderator_pin_regions = create_pin_regions(params, 'moderator')
    
    # Creating moderator pin materials
    moderator_materials = []
    for mat in params['Moderator Pin Materials']:
        if mat == None:
            moderator_materials.append(None)
        else: 
            material_1 = materials_database[mat]
            moderator_materials.append(material_1)
    moderator_materials.append(coolant)
    
    # Giving the user error message if the number of materials is not the same as the number of regions
    assert len(moderator_pin_regions) == len(moderator_materials), "The number of regions, {len(moderator_pin_regions)} should be\
        the same as the number of introduced materials, {len(moderator_materials)}"

    # creating them moerator pin universe
    moderator_cells = create_cells(moderator_pin_regions, moderator_materials)
    moderator_pin_universe = openmc.Universe(cells=moderator_cells.values())

    if params['plotting'] == "Y":
        # plotting
        create_universe_plot(moderator_pin_universe, 
                        pin_plot_width = 2.2 * params['Moderator Pin Radii'][-1],
                        num_pixels = 500, 
                        font_size = 16,
                        title = "Moderator Pin Universe", 
                        fig_size = 8, 
                        output_file_name = "moderator_pin_universe.png")
    
    
    # Coolant Universe
    coolant_cell = openmc.Cell(fill=coolant)
    coolant_universe = openmc.Universe(cells=(coolant_cell,))
   
    # **************************************************************************************************************************
    #                                                Sec. 1.3 : CONTROL DRUMS
    # **************************************************************************************************************************
    
    drums = create_drums_universe(params, control_drum_absorber, control_drum_reflector )
    
    # **************************************************************************************************************************
    #                                                Sec. 1.4 : Fuel Assembly
    # **************************************************************************************************************************
    
    # The center-to-center distance between adjacent fuel/moderator pins
    pin_pitch = 2 * (params['Fuel Pin Radii'][-1] ) + params["Pin Gap Distance"]

    assembly_universe = create_assembly_universe(params,
                                                 fuel_pin_universe,
                                                 moderator_pin_universe,
                                                 pin_pitch,
                                                 reflector,
                                                 coolant_universe)


    if params['plotting'] == "Y":
        # plotting
        create_universe_plot(assembly_universe, 
                        pin_plot_width = 2.2 * params['Lattice Radius'],
                        num_pixels = 2000, 
                        font_size = 16,
                        title = "Fuel Assembly Universe", 
                        fig_size = 8, 
                        output_file_name = "fuel_assembly_universe.png")
    # # **************************************************************************************************************************
    # #                                                Sec. 1.5 : VOLUME INFO for Depletion
    # # **************************************************************************************************************************
    #find where the fuel is in the fuel pin
    fuel_index = params['Fuel Pin Materials'].index(params['Fuel'])

    fissile_area = circle_area(params['Fuel Pin Radii'][fuel_index] )\
        - circle_area(params['Fuel Pin Radii'][fuel_index - 1])
    fuel.volume = fissile_area *params['Active Height'] * params['Fuel Pin Count']
   
    all_materials = fuel_materials +\
        moderator_materials + [coolant, reflector, control_drum_absorber, control_drum_reflector]
    
    # removing "None" materials
    all_materials_cleaned_list = [item for item in all_materials if item is not None]
    materials = openmc.Materials(list(set(all_materials_cleaned_list)))
   
    openmc.Materials.cross_sections = params['cross_sections_xml_location']
    materials.export_to_xml()
    

    # # **************************************************************************************************************************
    # #                                                Sec. 1.6 : CORE DRUM REPLACEMENT
    # # **************************************************************************************************************************

    control_drum_positions = create_control_drums_positions(number_of_drums = len(drums))

    core_geometry = create_core_geometry(params,
                                         drums,
                                         drums_positions = control_drum_positions,
                                         assembly_universe  = assembly_universe )


    core_geometry.export_to_xml()
    
    if params['plotting'] == "Y":
        create_universe_plot(core_geometry, 
                        pin_plot_width = 2.01 * params['Core Radius'],
                        num_pixels = 2000, 
                        font_size = 16,
                        title = "Reactor Core", 
                        fig_size = 8, 
                        output_file_name = "core.png")

        # # Plotting
        # plot = openmc.Plot.from_geometry(core_geometry)
        # plot.pixels = (2000, 2000)
        # plot.filename = 'lattice_plot'
    
        # plot.to_ipython_image()
    
    # # **************************************************************************************************************************
    # #                                                Sec. 1.7 : SIMULATION
    # # **************************************************************************************************************************
 
    point = openmc.stats.Point((0, 0, 0))
    source = openmc.Source(space=point)
    settings = openmc.Settings()
    settings.source = source
    settings.batches = 100
    settings.inactive = 50
    settings.particles = 2000

    # Set max_lost_particles to half the number of particles
    settings.max_lost_particles = 1000
    settings.export_to_xml()