# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

from reactor_engineering_evaluation.tools import ellipsoid_shell, circle_area, materials_densities


# Vessel Calcs
def vessels_specs(params):
    # Refers to the Inner Vessel 
    # For the GCMR: the core barrel
    vessel_height = (params['Active Height'] 
                     + 2*params['Axial Reflector Thickness'] 
                     + params['Vessel Lower Plenum Height'] 
                     + params['Vessel Upper Plenum Height'] 
                     + params['Vessel Upper Gas Gap']) # This is the first vessel
    if params['reactor type'] == "GCMR":
        # Volume based on CAD model
        # Has upper and lower head (ellipsoid)
        vessel_volume = (ellipsoid_shell(params['Vessel Radius'], params['Vessel Radius'], params['Vessel Bottom Depth']) * params['Vessel Thickness']
                        + (circle_area(params['Vessel Radius'] + params['Vessel Thickness']) - circle_area(params['Vessel Radius'])) * vessel_height)
    else:
        vessel_volume = (ellipsoid_shell(params['Vessel Radius'], params['Vessel Radius'], params['Vessel Bottom Depth'])/2)\
            * params['Vessel Thickness'] + (circle_area(params['Vessel Radius'] + params['Vessel Thickness'])\
                - circle_area(params['Vessel Radius'])) * vessel_height
    vessel_mass_kg = vessel_volume * materials_densities(params['Vessel Material'])/1000

    # Refers to the Outer Vessel
    # For the GCMR: RPV
    # For the LTMR: Guard Vessel
    guard_vessel_radius = params['Vessel Radius'] + params['Vessel Thickness'] + params['Gap Between Vessel And Guard Vessel'] 
    guard_bottom_depth = params['Vessel Bottom Depth'] + params['Vessel Thickness'] + params['Gap Between Vessel And Guard Vessel']
    if params['reactor type'] == "GCMR":
        guard_vessel_volume = (ellipsoid_shell(guard_vessel_radius, guard_vessel_radius, guard_bottom_depth) * params['Guard Vessel Thickness'] 
                              + (circle_area(guard_vessel_radius + params['Guard Vessel Thickness']) - circle_area(guard_vessel_radius)) * vessel_height)
    else:
        guard_vessel_volume = (ellipsoid_shell(guard_vessel_radius, guard_vessel_radius, guard_bottom_depth)/2)*\
            params['Guard Vessel Thickness'] + (circle_area(guard_vessel_radius + params['Guard Vessel Thickness']) -\
                circle_area(guard_vessel_radius)) * vessel_height
    guard_vessel_mass_kg = guard_vessel_volume * materials_densities(params['Guard Vessel Material'])/1000

    # Refers to the RCCS / Cooling Vessel
    cooling_vessel_radius = guard_vessel_radius + params['Gap Between Guard Vessel And Cooling Vessel'] # cm
    cooling_bottom_depth = guard_bottom_depth + params['Guard Vessel Thickness'] +\
        params['Gap Between Guard Vessel And Cooling Vessel']
    cooling_vessel_volume = (ellipsoid_shell(cooling_vessel_radius, cooling_vessel_radius, cooling_bottom_depth)/2)*\
        params['Cooling Vessel Thickness'] + (circle_area(cooling_vessel_radius + params['Cooling Vessel Thickness'])\
            - circle_area(cooling_vessel_radius)) * vessel_height
    cooling_vessel_mass = cooling_vessel_volume * materials_densities(params['Cooling Vessel Material'])/1000
    
    # Refers to the RCCS Intake Vessel
    intake_vessel_radius = cooling_vessel_radius + params['Gap Between Cooling Vessel And Intake Vessel']
    intake_bottom_depth = cooling_bottom_depth + params['Cooling Vessel Thickness'] + params['Gap Between Cooling Vessel And Intake Vessel']
    intake_vessel_volume = (ellipsoid_shell(intake_vessel_radius, intake_vessel_radius, intake_bottom_depth)/2)\
        * params['Intake Vessel Thickness'] + (circle_area(intake_vessel_radius + params['Intake Vessel Thickness']) -\
            circle_area(intake_vessel_radius)) * vessel_height
        
    intake_vessel_mass = intake_vessel_volume * materials_densities(params['Intake Vessel Material'])/1000
    
    # Includes Inner + Outer Vessel + RCCS + RCCS Intake
    total_vessel_height = intake_bottom_depth + vessel_height
    vessels_full_radius = intake_vessel_radius + params['Intake Vessel Thickness']
    total_vessels_mass = vessel_mass_kg + guard_vessel_mass_kg + cooling_vessel_mass + intake_vessel_mass
    
    # NOTE for GCMR:
    # Designs are different in that the Vessels Total Radius or Height 
    # is not the true vessel specs and includes the RCCS and the RCCS Intake.
    # Guard Vessel refers to the RPV
    # Vessel refers to the Core Barrel
    params['Vessels Total Radius'] = vessels_full_radius 
    params['Vessel Height'] = vessel_height
    params['Vessels Total Height'] = total_vessel_height
    params['Guard Vessel Radius'] = guard_vessel_radius
    params['Cooling Vessel Radius'] = cooling_vessel_radius
    params['Intake Vessel Radius'] = intake_vessel_radius
    params['Vessel Mass'] = vessel_mass_kg
    params['Guard Vessel Mass'] = guard_vessel_mass_kg
    params['Cooling Vessel Mass'] = cooling_vessel_mass
    params['Intake Vessel Mass'] = intake_vessel_mass
    params['Total Vessels Mass'] = total_vessels_mass

    if params['reactor type'] == 'GCMR':
        params['RPV Outer Radius'] = (params['Guard Vessel Radius'] + params['Guard Vessel Thickness'])
        params['RPV Outer Height'] = params['Vessel Height'] + 2*params['Gap Between Vessel And Guard Vessel'] + 2*params['Guard Vessel Thickness'] + 2*params['Vessel Bottom Depth'] + 2*params['Vessel Thickness']
