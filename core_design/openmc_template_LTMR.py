# Importing libraries
import openmc
from core_design.utils import *
from core_design.openmc_materials_database import collect_materials_data

"""
An OpenMC function that accepts an instance of "parameters" 
and generates the necessary XMl files
"""

def build_openmc_model_LTMR(params):
    
    # **************************************************************************************************************************
    #                                                Sec. 1 : MATERIALS
    # **************************************************************************************************************************
    materials_database = collect_materials_data(params)
    
    fuel = materials_database[params['fuel']]
    coolant = materials_database[params['coolant']]
    reflector = materials_database[params['Reflector']]
    control_drum_absorber = materials_database[params['Control Drum Absorber']]
    control_drum_reflector = materials_database[params['Control Drum Reflector']]
    
    # **************************************************************************************************************************
    #                                                Sec. 2 : GEOMETRY: Fuel Pins, Moderator Pins, Coolant
    # **************************************************************************************************************************
    
    # Fuel Pins
    fuel_pin_region = create_fuel_pin_regions(params)
    fuel_materials = []
    for mat in params['fuel_pin_materials']:
        if mat == None:
            fuel_materials.append(None)
        else: 
            material_1 = materials_database[mat]
            fuel_materials.append(material_1)
    fuel_materials.append(materials_database[params['coolant']])

     # Giving the user error message if the number of materials is not the same as the number of regions
    assert len(fuel_pin_region) == len(fuel_materials), "The number of regions, {len(fuel_pin_region)} should be\
        the same as the number of introduced materials, {len(fuel_materials)}"
    
    fuel_cells = create_cells(fuel_pin_region, fuel_materials)
    fuel_pin_universe = openmc.Universe(cells=fuel_cells.values())

    if params['plotting'] == "yes":
        # plotting
        create_universe_plot(fuel_pin_universe, 
                        pin_plot_width = 2.2 * params['fuel_pin_radii'][-1],
                        num_pixels = 500, 
                        font_size = 16,
                        title = "Fuel Pin", 
                        fig_size = 8, 
                        output_file_name = "fuel_pin.png")

    
    # Moderator Pins
    moderator_pin_region = create_moderator_pin_regions(params)
    moderator_materials = []
    for mat in params['moderator_pin_materials']:
        if mat == None:
            moderator_materials.append(None)
        else: 
            material_1 = materials_database[mat]
            moderator_materials.append(material_1)
    moderator_materials.append(materials_database[params['coolant']])
    
    # Giving the user error message if the number of materials is not the same as the number of regions
    assert len(moderator_pin_region) == len(moderator_materials), "The number of regions, {len(moderator_pin_region)} should be\
        the same as the number of introduced materials, {len(moderator_materials)}"

    moderator_cells = create_cells(moderator_pin_region, moderator_materials)
    moderator_pin_universe = openmc.Universe(cells=moderator_cells.values())

    if params['plotting'] == "yes":
        # plotting
        create_universe_plot(moderator_pin_universe, 
                        pin_plot_width = 2.2 * params['moderator_pin_radii'][-1],
                        num_pixels = 500, 
                        font_size = 16,
                        title = "Moderator Pin", 
                        fig_size = 8, 
                        output_file_name = "moderator_pin.png")
    
    
    # Coolant
    coolant_cell = openmc.Cell(fill=coolant)
    coolant_universe = openmc.Universe(cells=(coolant_cell,))
   
    # **************************************************************************************************************************
    #                                                Sec. 3 : CONTROL DRUMS
    # **************************************************************************************************************************


    drums = create_drums_universe(absorber_thickness = params['drum_Absorber_thickness'],
                                  drum_radius = params['Drum_Radius'],
                          control_drum_absorber_material = control_drum_absorber,
                          control_drum_reflector_material = control_drum_reflector )
    
    if params['plotting'] == "yes":
        # plotting
        for ii in range (len(drums)):
            the_drum = drums[ii]
            create_universe_plot(the_drum, 
                    pin_plot_width = 2.2 * params['Drum_Radius'],
                    num_pixels = 500, 
                    font_size = 16,
                    title = "Control Drum", 
                    fig_size = 8, 
                    output_file_name = f"control_drum{ii}.png")
    
    # **************************************************************************************************************************
    #                                                Sec. 4 : ASSEMBLY & RINGS
    # **************************************************************************************************************************
    
    pin_pitch = 2 * (params['fuel_pin_radii'][-1] ) + params["pin_gap_distance"]

    assembly_universe = create_assembly_universe(params,
                                                 fuel_pin_universe,
                                                 moderator_pin_universe,
                                                 pin_pitch,
                                                 reflector_material = materials_database[params['reflector']],
                                                 outer_universe = coolant_universe)

    # **************************************************************************************************************************
    #                                                Sec. 5 : VOLUME INFO for Depletion
    # **************************************************************************************************************************
    #find where the fuel is in the fuel pin
    fuel_index = params['fuel_pin_materials'].index(params['fuel'])

    fissile_area = circle_area(params['fuel_pin_radii'][fuel_index] )\
        - circle_area(params['fuel_pin_radii'][fuel_index - 1])
    fuel.volume = fissile_area *params['active_height'] * params['fuel_pin_count']
   
    all_materials = fuel_materials +\
        moderator_materials + [coolant, reflector, control_drum_absorber, control_drum_reflector]
    
    # removing "None" materials
    all_materials_cleaned_list = [item for item in all_materials if item is not None]
    materials = openmc.Materials(list(set(all_materials_cleaned_list)))
   
    openmc.Materials.cross_sections = params['cross_sections_xml_location']
    materials.export_to_xml()
    

    # **************************************************************************************************************************
    #                                                Sec. 6 : CORE DRUM REPLACEMENT
    # **************************************************************************************************************************

    control_drum_positions = create_control_drums_positions(params,
                                                            number_of_drums = len(drums))

    core_geometry = create_core_geometry(params,
                                         drums,
                                         drums_positions = control_drum_positions,
                                         assembly_universe  = assembly_universe )


    core_geometry.export_to_xml()
    
    if params['plotting'] == "yes":
        # Plotting
        plot = openmc.Plot.from_geometry(core_geometry)
        plot.pixels = (2000, 2000)
        plot.filename = 'lattice_plot'
    
        plot.to_ipython_image()
    
    # **************************************************************************************************************************
    #                                                Sec. 7 : SIMULATION
    # **************************************************************************************************************************
 
    point = openmc.stats.Point((0, 0, 0))
    source = openmc.Source(space=point)
    settings = openmc.Settings()
    settings.source = source
    settings.batches = 100
    settings.inactive = 50
    settings.particles = 400
    settings.export_to_xml()