# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
import math
from core_design.openmc_materials_database import *
from core_design.utils import *

def calculate_drums_volumes_and_masses(params):
    
    DRUM_RADIUS = params['Drum Radius']
    drum_height = params['Drum Height']
    absorber_thickness = params['Drum Absorber Thickness']
    drum_volume = 3.14*DRUM_RADIUS * DRUM_RADIUS *drum_height
    if 'coating_angle' in params.keys():
        drum_absorp_vol = (3.14*( DRUM_RADIUS * DRUM_RADIUS) - 3.14/180*params['coating_angle']*(DRUM_RADIUS-absorber_thickness)*(DRUM_RADIUS-absorber_thickness))*drum_height/3
    else:
        drum_absorp_vol = (3.14*( DRUM_RADIUS * DRUM_RADIUS - (DRUM_RADIUS-absorber_thickness)*(DRUM_RADIUS-absorber_thickness) )*drum_height)/3    
    drum_refl_vol = drum_volume - drum_absorp_vol 
    if params['reactor type'] == "LTMR":
        number_of_drums = 12 
        params['Drum Count'] = number_of_drums
    elif params['reactor type'] == "GCMR":
        if 'Drum Count' in params.keys():
            number_of_drums = params['Drum Count']
        else:
            number_of_drums = 6 * (params['Core Rings']-1) 
            params['Drum Count'] = number_of_drums
    elif params['reactor type'] == "HPMR":
            number_of_drums = 12 
            params['Drum Count'] = number_of_drums

    all_drums_volume = drum_volume * number_of_drums
    
    drum_absorp_vol_all = drum_absorp_vol  * number_of_drums  
    drum_refl_vol_all = drum_refl_vol  * number_of_drums  
    
    materials_database = collect_materials_data(params)
    drums_absorber_density = materials_database[params['Control Drum Absorber']].density
    drums_reflector_density =  materials_database[params['Control Drum Reflector']].density
    drum_absorp_all_mass = drum_absorp_vol_all * drums_absorber_density/1000 # (in Kg)
    drum_refl_all_mass = drum_refl_vol_all  * drums_reflector_density/1000 #  (in Kg)
    
    Control_Drums_Mass =  drum_absorp_all_mass + drum_refl_all_mass
    params['All Drums Volume'] = all_drums_volume
    params['Control Drum Absorber Mass'] = drum_absorp_all_mass
    params['Control Drum Reflector Mass'] = drum_refl_all_mass
    params['Control Drums Mass'] = Control_Drums_Mass # all the drums masses
    params['All Drums Area'] = params['All Drums Volume']  / params['Drum Height']



def hexagonal_area_from_ftf(ftf_distance):
    # Calculate the area directly from the flat-to-flat distance
    area = (np.sqrt(3) / 2) * ftf_distance ** 2
    return area

def calculate_reflector_mass_LTMR(params):
    hex_area =  2.598 * params['Lattice Radius'] * params['Lattice Radius']
    core_radius = params['Core Radius']
    area_of_all_drums = params['All Drums Area'] 
    drum_height = params['Drum Height']
    # I assume for now that the drums are always fully inside the reflector
    
    area_reflector = 3.14 * core_radius * core_radius - hex_area  - area_of_all_drums # cm2
    vol_reflector = area_reflector * drum_height # cm^3
    
    materials_database = collect_materials_data(params)
    reflector_density = materials_database[params['Reflector']].density
    mass_reflector = vol_reflector * reflector_density/1000 # mass in Kg
    params['Reflector Mass'] = mass_reflector 
    params['Axial Reflector Mass'] = (1/1000) * reflector_density * cylinder_volume(core_radius, params['Axial Reflector Thickness'])



def calculate_reflector_mass_GCMR(params):
    materials_database = collect_materials_data(params)
    tot_number_assemblies = calculate_number_of_rings(params['Core Rings'] )
    reflector_height = params['Active Height']
    reflector_volume = reflector_height * (circle_area(params['Core Radius']) 
                                           - tot_number_assemblies*hexagonal_area_from_ftf(params['Assembly FTF'])
                                           - params['All Drums Area']) 
        
    reflector_density = materials_database[params['Reflector']].density
    reflector_mass = reflector_density * reflector_volume  / 1000 # Kg
    params['Reflector Mass'] = reflector_mass
    params['Axial Reflector Mass'] = 2 *  (1/1000) * reflector_density * cylinder_volume(params['Core Radius'], params['Axial Reflector Thickness'])


def calculate_moderator_mass_GCMR(params): 
    materials_database = collect_materials_data(params)
    tot_number_assemblies = calculate_number_of_rings(params['Core Rings'] )

    # The area of one hexagonal lattice in the core
    A_hex  = hexagonal_area_from_ftf(params['Assembly FTF'])
    
    # area occuplied by the fuel in one hexagonal lattice (assembly)
    num_fuel_regions_per_hex = calculate_number_of_rings( params['Assembly Rings'] - 1 )

    area_fuel_per_hex = params['Packing Fraction'] * circle_area(params['Compact Fuel Radius']) * num_fuel_regions_per_hex
    area_coolant_per_hex = 2 * num_fuel_regions_per_hex * circle_area(params['Coolant Channel Radius'])
    area_moderator_booster_per_hex =  0.5 * 6 * (params['Assembly Rings'] - 1) * circle_area(params['Moderator Booster Radius'] )

    # area ocuupied by the moderators in one of one hexagonal lattices
    moderator_area = A_hex - area_fuel_per_hex - area_coolant_per_hex - area_moderator_booster_per_hex 

    # total moderator mass
    tot_moderator_mass = tot_number_assemblies * moderator_area  * params['Active Height'] * materials_database[params['Moderator']].density / 1000 # Kg
    params['Moderator Mass'] = tot_moderator_mass 
    tot_booster_mass   = tot_number_assemblies * area_moderator_booster_per_hex * params['Active Height'] * materials_database[params['Moderator Booster']].density / 1000 # Kg
    params['Moderator Booster Mass'] = tot_booster_mass

def calculate_reflector_and_moderator_mass_HPMR(params):
    materials_database = collect_materials_data(params)
    # first, determine the area of the big hexagonal of monolith surrounding the assemblies
    assembly_long_diag = 1.1547 * params['Assembly FTF']
    assembly_side_length =  params['Assembly FTF'] / (np.sqrt(3))
    big_hex_FTF = params['Number of Rings per Core'] *  assembly_long_diag  + (params['Number of Rings per Core'] - 1) * assembly_side_length
    big_hex_area = hexagonal_area_from_ftf(big_hex_FTF)
    reflector_volume = ( circle_area(params['Core Radius']) - big_hex_area ) * params['Active Height']
    reflector_density = materials_database[params['Reflector']].density
    reflector_mass    = reflector_density * reflector_volume  / 1000 # Kg
    params['Reflector Mass'] = reflector_mass
    # mass of two axial reflectors
    params['Axial Reflector Mass'] = 2 * (1/1000) * reflector_density * cylinder_volume(params['Core Radius'], params['Axial Reflector Thickness'])

    # moderator = big hex minus the fuel and heatpipes
    fuel_area     =  params['Fuel Pin Count']  * circle_area(params['Fuel Pin Radii'][-1])
    heatpipe_area =  params['Number of Heatpipes'] * circle_area(params['Heat Pipe Radii'][-1])
    params['Moderator Total Area'] = big_hex_area - fuel_area - heatpipe_area
    params['Moderator Mass'] = params['Moderator Total Area'] * params['Active Height'] * materials_database[params['Moderator']].density / 1000 #Kg
    



def calculate_reflector_and_moderator_mass_HPMR_vtb(params):
    materials_database = collect_materials_data(params)
    # first, determine the area of the big hexagonal of monolith surrounding the assemblies
    assembly_long_diag = 1.1547 * params['Assembly FTF']
    assembly_side_length =  params['Assembly FTF'] / (np.sqrt(3))
    big_hex_FTF = params['Number of Rings per Core'] *  assembly_long_diag  + (params['Number of Rings per Core'] - 1) * assembly_side_length
    big_hex_area = hexagonal_area_from_ftf(big_hex_FTF)
    reflector_volume = (circle_area(params['Core Radius']) - big_hex_area ) * params['Active Height']
    reflector_density = materials_database[params['Radial Reflector']].density
    reflector_mass    = reflector_density * reflector_volume  / 1000 # Kg
    params['Reflector Mass'] = reflector_mass
    # mass of two axial reflectors
    params['Axial Reflector Mass'] = 2 * (1/1000) * materials_database[params['Axial Reflector']].density * cylinder_volume(params['Core Radius'], params['Axial Reflector Thickness'])

    # moderator = big hex minus the fuel and heatpipes
    fuel_area     =  params['Fuel Pin Count']  * circle_area(params['Fuel Pin Radii'][-1])
    heatpipe_area =  params['Number of Heatpipes'] * circle_area(params['Heat Pipe Radii'][-1])
    moderator_booster_area = params['Number of Moderator Booster'] * circle_area(params['Moderator Booster Raddi'])
    params['Moderator Booster Mass'] = moderator_booster_area * params['Active Height'] * materials_database[params['Moderator Booster']].density / 1000 #Kg
    
    # Remove drum area from the last ring
    DRUM_RADIUS = params['Drum Radius']
    drum_area = 3.14*DRUM_RADIUS * DRUM_RADIUS
    number_of_drums = 12 
    params['Moderator Total Area'] = big_hex_area - fuel_area - heatpipe_area - moderator_booster_area - drum_area*number_of_drums
    params['Moderator Mass'] = params['Moderator Total Area'] * params['Active Height'] * materials_database[params['Moderator']].density / 1000 #Kg
    #sys.exit()

def calculate_moderator_mass(params): 
    # for the moderator pins
    materials_database = collect_materials_data(params)  
    moderator_volume = params['Moderator Pin Count'] * circle_area( (params['Moderator Pin Radii'])[0]) *params['Active Height'] 
    moderator_mass = (1/ 1000) * moderator_volume * materials_database[(params['Moderator Pin Materials'])[0]].density # Kg
    return moderator_mass
