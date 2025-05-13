from reactor_engineering_evaluation.tools import ellipsoid_shell, circle_area, materials_densities


# Vessel Calcs
def vessels_specs(params):
    
    vessel_height = params['active_height'] + 2.5* params['extra_reflector'] +\
        params['vessel_lower_plenum_height'] + params['vessel_upper_plenum_height'] +\
            params['vessel_upper_gas_gap'] # This is the first vessel
    
    vessel_volume = (ellipsoid_shell(params['vessel_radius'], params['vessel_radius'], params['vessel_bottom_depth'] )/2)\
        *params['vessel_thickness'] + (circle_area(params['vessel_radius'] + params['vessel_thickness'])\
            - circle_area(params['vessel_radius'])) * vessel_height
    vessel_mass_kg = vessel_volume  * materials_densities(params['vessel_material'])/1000



    guard_vessel_radius = params['vessel_radius'] + params['vessel_thickness'] + params['gap_between_vessel_and_guard_vessel'] 
    guard_bottom_depth = params['vessel_bottom_depth']  + params['vessel_thickness'] + params['gap_between_vessel_and_guard_vessel']
    guard_vessel_volume = (ellipsoid_shell(guard_vessel_radius, guard_vessel_radius, guard_bottom_depth)/2)*\
        params['guard_vessel_thickness'] + (circle_area(guard_vessel_radius + params['guard_vessel_thickness']) -\
            circle_area(guard_vessel_radius))*vessel_height
    guard_vessel_mass_kg = guard_vessel_volume * materials_densities(params['guard_vessel_material'])/1000

    # cooling vessel
 
    cooling_vessel_radius = guard_vessel_radius + params['gap_between_guard_vessel_and_cooling_vessel'] # cm
    cooling_bottom_depth = guard_bottom_depth + params['guard_vessel_thickness'] +\
        params['gap_between_guard_vessel_and_cooling_vessel']
    cooling_vessel_volume = (ellipsoid_shell(cooling_vessel_radius, cooling_vessel_radius, cooling_bottom_depth)/2)*\
        params['cooling_vessel_thickness'] + (circle_area(cooling_vessel_radius + params['cooling_vessel_thickness'])\
            - circle_area(cooling_vessel_radius))*vessel_height
    cooling_vessel_mass = cooling_vessel_volume *materials_densities(params['cooling_vessel_material'])/1000
    
    # gap_intake = 3

    # This is the intake vessel
    # intake_vessel_thickness = 0.3
    intake_vessel_radius = cooling_vessel_radius + params['gap_between_cooling_vessel_and_intake_vessel']
    intake_bottom_depth = cooling_bottom_depth + params['cooling_vessel_thickness'] + params['gap_between_cooling_vessel_and_intake_vessel']
    intake_vessel_volume = (ellipsoid_shell(intake_vessel_radius, intake_vessel_radius, intake_bottom_depth)/2)\
        *params['intake_vessel_thickness']  + (circle_area(intake_vessel_radius + params['intake_vessel_thickness'] ) -\
            circle_area(intake_vessel_radius))*vessel_height
        
    intake_vessel_mass = intake_vessel_volume *materials_densities(params['intake_vessel_material'] )/1000
    
    total_vessel_height = intake_bottom_depth + vessel_height
    vessels_full_radius = intake_vessel_radius + params['intake_vessel_thickness']
    
    return vessels_full_radius, vessel_height, total_vessel_height,\
        vessel_mass_kg, guard_vessel_mass_kg, cooling_vessel_mass, intake_vessel_mass




