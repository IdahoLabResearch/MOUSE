from core_design.openmc_materials_database import *
from core_design.utils import *

def calculate_drum_volume(DRUM_RADIUS, drum_height, absorber_thickness, params):
    drum_volume = 3.14*DRUM_RADIUS * DRUM_RADIUS *drum_height
    
    drum_absorp_vol = (3.14*( DRUM_RADIUS * DRUM_RADIUS - (DRUM_RADIUS-absorber_thickness)*(DRUM_RADIUS-absorber_thickness) )*drum_height)/3
    drum_refl_vol = drum_volume - drum_absorp_vol 
    if params['reactor type'] == "LTMR":
        number_of_drums = 12 
    elif params['reactor type'] == "GCMR":
        number_of_drums = 6 * (params['core_rings']-1) 
        

    all_drums_volume = drum_volume * number_of_drums
    
    drum_absorp_vol_all = drum_absorp_vol  * number_of_drums  
    drum_refl_vol_all = drum_refl_vol  * number_of_drums  
    
    materials_database = collect_materials_data(params)
    drums_absorber_density = materials_database[params['Control Drum Absorber']].density
    drums_reflector_density =  materials_database[params['Control Drum Reflector']].density
    drum_absorp_all_mass = drum_absorp_vol_all * drums_absorber_density/1000 # (in Kg)
    drum_refl_all_mass = drum_refl_vol_all  * drums_reflector_density/1000 #  (in Kg)
    
    Control_Drums_Mass =  drum_absorp_all_mass + drum_refl_all_mass
    return all_drums_volume, drum_absorp_all_mass, drum_refl_all_mass, Control_Drums_Mass


def calculate_reflector_mass_GCMR(params):
    materials_database = collect_materials_data(params)
    tot_number_assemblies = calculate_number_of_core_rings(params['core_rings'] )
    reflector_volume = (circle_area(params['core_radius']) -\
       tot_number_assemblies * hex_area(params['assembly_ftf']) - params['tot_drum_area_all']) * params['active_height']
        
    reflector_density = materials_database[params['Reflector']].density
    reflector_mass = reflector_density * reflector_volume  / 1000 # Kg
    return reflector_mass

def calculate_moderator_mass_GCMR(params): 
    materials_database = collect_materials_data(params)
    tot_number_assemblies = calculate_number_of_core_rings(params['core_rings'] )

    # The area of one hexagonal lattice in the core
    A_hex  = hex_area(params['assembly_ftf'])
    
    # area occuplied by the fuel in one hexagonal lattice (assembly)
    num_fuel_regions_per_hex = calculate_number_of_core_rings(params['assembly_rings'] - 1 )
    area_fuel_per_hex = params['packing factor'] * circle_area(params['lattice radius']) * num_fuel_regions_per_hex
    area_coolant_per_hex = 2 * num_fuel_regions_per_hex * circle_area(params['coolant channel radius'])
    area_moderator_booster_per_hex =  0.5 * 6 * (params['assembly_rings'] - 1) * circle_area(params['booster radius'] )

    # area ocuupied by the moderators in one of one hexagonal lattices
    moderator_area = A_hex - area_fuel_per_hex - area_coolant_per_hex - area_moderator_booster_per_hex 

    # total moderator mass
    tot_moderator_mass = (tot_number_assemblies * moderator_area  * params['active_height'] * materials_database[params['Moderator']].density / 1000) # Kg
    tot_booster_mass = tot_number_assemblies * area_moderator_booster_per_hex * params['active_height'] * materials_database[params['moderator_booster']].density / 1000
    return tot_moderator_mass, tot_booster_mass


def calculate_moderator_mass(params): 
    # for the moderator pins
    materials_database = collect_materials_data(params)  
    moderator_volume = params['moderator_pin_count'] * circle_area( (params['moderator_pin_radii'])[0]) *params['active_height'] 
    moderator_mass = (1/ 1000) * moderator_volume * materials_database[(params['moderator_pin_materials'])[0]].density # Kg
    return moderator_mass
