"""
This example demonstrates a detailed bottom-up cost estimate for a 15 MWe Liquid Metal Thermal Microreactor (LTMR).
OpenMC is used for core design calculations.
Back-of-the-envelope calculations were conducted for other Balance of Plant components.
Users can modify any of the parameters in the "params" dictionary below.
"""


# Import necessary libraries
import numpy as np
import time
import warnings

# Import custom modules and functions
import watts  # Import watts library for simulation workflows for one or multiple codes
import openmc  # Import OpenMC for nuclear reactor physics simulations

# Import templates and utility functions from the core design module
from core_design.openmc_template_LTMR import *
from core_design.drums import *
# from core_design.utils import *
from core_design.pins_arrangement import rings_1

# Import engineering evaluation tools and functions
from reactor_engineering_evaluation.tools import *
from reactor_engineering_evaluation.operation import *
from reactor_engineering_evaluation.fuel_calcs import *
from reactor_engineering_evaluation.vessels_calcs import *
from reactor_engineering_evaluation.BOP import *

# Import cost estimation functions
from cost.baseline_costs import *

# Suppress warnings
warnings.filterwarnings("ignore")  # Ignore all warnings

# Record the current time (in seconds)
time_start = time.time()  # Store the current time to track the script's execution duration

# Initialize user parameters using the watts library
params = watts.Parameters()


# **************************************************************************************************************************
#                                                Sec. 1 : User Configuration
# **************************************************************************************************************************

# Whether the user prefers to display the core design plots. 
# Note: Plotting is not recommended for sensitivity analysis and parametric studies.
params['plotting'] = "No"  # options are "Yes" and "No"

# If the user is using OpenMC calculations, the location of the neutron data library must be specified. 
# For example, data can be obtained from: https://openmc.org/official-data-libraries/
params['cross_sections_xml_location'] =\
    '/home/hannbn/projects/MARVEL_MRP/Github_repos/openmc_data/endfb-viii.0-hdf5/cross_sections.xml'

# If the user is using OpenMC calculations, the location of the depletion chain file must be specified. 
# For example, data can be obtained from: https://openmc.org/depletion-chains/
params['simplified_chain_thermal_xml'] =\
    '/home/hannbn/projects/MARVEL_MRP/Github_repos/openmc_data/simplified_thermal_chain11.xml'


# **************************************************************************************************************************
#                                                Sec. 2 : User-Defined Parameters (Materials)
# **************************************************************************************************************************  
params['reactor type'] = "LTMR"
# The user can change any of the following materials,
# but the new material has to be included in the materials database at "core_design/openmc_materials_database"
params['TRISO Fueled'] = "No"  # options are "Yes" and "No"
# The fuel and its properties
params['fuel'] = 'TRIGA_fuel'
params['enrichment'] = 0.1975 # The enrichment is a fraction. It has to be between 0 and 1
params["H_Zr_ratio"] = 1.6 # The proportion of hydrogen atoms relative to zirconium atoms
params['U_met_wo'] = 0.3 # The ratio between the weight of Uranium and the total weight of the fuel (less than 1)

# Coolant, reflector and control drums
params['coolant'] = 'NaK'
params['Reflector'] = 'BeO' # external reflector around the core
params['Control Drum Absorber'] = 'B4C_natural' # The absorber material in the control drums
params['Control Drum Reflector'] = 'BeO'     # The reflector material in the control drums   

# Temperature assumed for the OpenMC calculations
params['common_temperature'] = 600 # Kelvins


# **************************************************************************************************************************
#                                           Sec. 3 : User-Defined Parameters (Geometry: Fuel Pins, Moderator Pins, Coolant)
# **************************************************************************************************************************  

# The reactor core includes fuel pins and moderator pins
# The fuel pin is assumed to include 5 regions:
# 1- a filler rod surrounded by 2- a gap, which is then surrounded by 3- the fuel. This is followed by 4- another gap, 
# which is finally surrounded by 5- the cladding.

# The materials of the 5 regions of the fuel pin
params['fuel_pin_materials'] = ['Zr', None, 'TRIGA_fuel', None, 'SS304']
# The outer radii of each of the 5 regions of the fuel pin
params['fuel_pin_radii'] = [0.28575, 0.3175, 1.5113, 1.5367, 1.5875] # cm

# The moderator pin is assumed to include 2 regions: a moderator and a cladding
params['moderator_pin_materials'] = ['ZrH', 'SS304']
params['Moderator'] = params['moderator_pin_materials'][0]
params['moderator_pin_radii'] = [1.5367, 1.5875] # cm


# **************************************************************************************************************************
#                                           Sec. 4 : User-Defined Parameters (Hexagonal Lattice)
# **************************************************************************************************************************  

# The gap between pins (fuel pins or moderator pins)
params["pin_gap_distance"] = 0.1 # cm

# The number of rings (pins) along the edge of the hexagonal lattice
params['assembly_rings'] = 12

# The previous parameters are used to calculate the hexagonal lattice radius and height
params['lattice_radius'] = calculate_lattice_radius(params['fuel_pin_radii'][-1], params["pin_gap_distance"], params['assembly_rings'])
params['active_height'] = 2 * params['lattice_radius']

# The selected fuel pins/moerator pins arrangement
params['rings'] = rings_1
# The total number of fuel pins in the lattice
params['fuel_pin_count'] = sum(row.count("FUEL") for row in params['rings'])
params['moderator_pin_count'] = sum(row.count("MODERATOR") for row in params['rings'])

params['Moderator Mass'] = calculate_moderator_mass(params)

#The thickness of the reflector around the lattice
params['extra_reflector'] = 14 #cm
# The total core radius 
params['core_radius'] = params['lattice_radius'] + params['extra_reflector']


# **************************************************************************************************************************
#                                           Sec. 5 : User-Defined Parameters (Control Drums)
# ************************************************************************************************************************** 

params['drum_radius_to_lattice_radius'] = 0.22784810068 # Ratio of the control drum radius to the lattice radius
params['Drum_Radius'] = params['drum_radius_to_lattice_radius'] * params['lattice_radius'] # cm # Control drum radius

# Each of the control drums includes a reflector material (most of the drum) and an absorber material (the outer layer)
# The thickness of the absorber layer
params['drum_Absorber_thickness'] = 1 # cm


# Placement of drums happen by tracing a line through the core apothems
# then 2 drums are place around each apothem by deviating from this line
# by a deviation angle
params["deviation angle between drums"] = 12.86 # degrees

# The radius of the tube of the control drum
params['drum_tube_radius'] = params['Drum_Radius'] + params['Drum_Radius'] / 90 # cm

params['drum_height'] = params['active_height']



# **************************************************************************************************************************
#                                           Sec. 6 : User-Defined Parameters (Overall System)
# ************************************************************************************************************************** 

# The reactor power (thermal) MW
params['Power MWt'] = 20 # MWt
params['thermal_efficiency'] = 0.31
params['Power MWe'] = params['Power MWt'] * params['thermal_efficiency']

# The actual heat flux (MW/m^2)
params['heat_flux'] = calculate_heat_flux(params['fuel_pin_radii'][-1], params['active_height'], params['rings'], params['Power MWt'])
# Target Heat Flux : Approximate calculated value for a typical sodium-cooled fast reactor (SFR)
params['heat_flux_criteria'] = 0.9

# The burnup steps for depletion calculations
params['burnup_steps_MWd_per_Kg'] = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0,
                                     30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 140.0]


# **************************************************************************************************************************
#                                           Sec. 7 : Running OpenMC
# ************************************************************************************************************************** 

# A function that runs OpenMC and the depletion calculations
# If the heat flux is too high, the code stops running
# The results such as the heat flux, the fuel lifetime and the mass of the fuel are added to the "params"

try:
    openmc_plugin = watts.PluginOpenMC(build_openmc_model_LTMR, show_stderr=True)  # running the LTMR Model
    # monitor_heat_flux(params, openmc_plugin)

except Exception as e:
    # Handle any errors that occur during the simulation
    print(f"An error occurred while running the OpenMC simulation: {e}")

# TEMPORARY  ## DELETE LATER!!!!!!!!!!!!!!!!!!!!
params['fuel_lifetime_days'] = 2078 # days
params['mass_U235'] = 67711.4 # grams
params['mass_U238'] = 278650.8  # grams
params['Uranium Mass'] = (params['mass_U235'] + params['mass_U238']) / 1000 # Kg

params['all_drums_volume'], params['Control Drum Absorber Mass'], params['Control Drum Reflector Mass'], params['Control Drums Mass'] =\
    calculate_drum_volume(params['Drum_Radius'], params['drum_height'],\
        params['drum_Absorber_thickness'], params)
        
params['tot_drum_area_all'] = params['all_drums_volume'] / params['drum_height']


# **************************************************************************************************************************
#                                           Sec. 7 : Fuel Calcs
# ************************************************************************************************************************** 

params['Natural Uranium Mass'], params['fuel_tail_waste_mass_Kg'], params['SWU'] =\
        fuel_calculations(params)


# **************************************************************************************************************************
#                                           Sec. 7 : Balance of Plant
# ************************************************************************************************************************** 
params['Primary HX Mass']  =  calculate_heat_exchanger_mass(params)[0] # Kg
params['Secondary HX Mass'] = 0
params['Primary Pump'] = 'Yes' # options are "Yes" and "No"
params['Secondary Pump'] = 'No' # options are "Yes" and "No"
params['Primary Pump Mechanical Power']  = calculate_pump_mechanical_power(params)[0]
params['Pump Isentropic Efficiency'] = 0.8
# **************************************************************************************************************************
#                                           Sec. 8 : Shielding
# ************************************************************************************************************************** 

# #Shielding
params['in_vessel_shield_thickness'] = 10.16 #cm
params['in_vessel_shielding_inner_radius'] = params['core_radius'] 
params['in_vessel_shielding_outer_radius'] = params['core_radius'] + params['in_vessel_shield_thickness']
params['In Vessel Shielding Material'] = 'B4C_natural' 


params['out_of_vessel_shield_thickness'] = 39.37 #cm
params['Out Of Vessel Shielding Material'] = 'WEP' # water extended polymer (WEP)
# The out of vessel shield is not fully made of the out of vessel material (e.g. WEP) so we use an effective density factor
params['out_vessel_shield_effective_density_factor'] = 0.5


# **************************************************************************************************************************
#                                           Sec. 9 : Vessels Calculations
# ************************************************************************************************************************** 
# # Vessels parameters
params['vessel_radius'] = params['core_radius'] + params['in_vessel_shield_thickness']  
params['vessel_thickness'] = 2 # cm
params['vessel_lower_plenum_height'] = 30 # cm
params['vessel_upper_plenum_height'] = 60 # cm
params['vessel_upper_gas_gap'] = 10 # cm
params['vessel_bottom_depth'] = 35  # cm
params['vessel_material'] ='stainless_steel'

params['gap_between_vessel_and_guard_vessel'] = 2 # cm
params['guard_vessel_thickness'] = 0.5
params['guard_vessel_material'] ='stainless_steel'

params['gap_between_guard_vessel_and_cooling_vessel'] = 5 # cm
params['cooling_vessel_thickness'] = 0.5 # cm
params['cooling_vessel_material'] ='stainless_steel'

params['gap_between_cooling_vessel_and_intake_vessel'] = 3 # cm
params['intake_vessel_thickness'] = 0.5 # cm
params['intake_vessel_material'] ='stainless_steel'

# Calculating the vessels dimensions and masses (in Kg)
params['vessels_total_radius'], params['vessel_height'] , params['vessels_total_height'],\
        params['Vessel Mass'], params['Guard Vessel Mass'] ,\
            params['Cooling Vessel Mass'], params['Intake Vessel Mass'] , params['Total Vessels Mass']= vessels_specs(params)

params['In Vessel Shielding Mass'] = cylinder_annulus_mass(params['in_vessel_shielding_outer_radius'],\
    params['in_vessel_shielding_inner_radius'], params['vessel_height'], params['In Vessel Shielding Material'] )  

params['Out Of Vessel Shielding Mass'] = params['out_vessel_shield_effective_density_factor'] * cylinder_annulus_mass(params['out_of_vessel_shield_thickness']+ params['vessels_total_radius'],\
        params['out_of_vessel_shield_thickness'], params['vessels_total_height'], params['Out Of Vessel Shielding Material']) 

params ['Vessel and Guard Vessel Masses'] = params['Vessel Mass'] +  params['Guard Vessel Mass']




params['hex_area'] = 2.598 * params['lattice_radius'] * params['lattice_radius']

params['Reflector Mass'] = calculate_reflector_mass_LTMR(params['hex_area'],\
    params['core_radius'], params['tot_drum_area_all'], params['drum_height'])
# **************************************************************************************************************************
#                                           Sec. 10 : Operation
# **************************************************************************************************************************

# # operation
params['num_people_required_per_refueling'] = 5
params['num_people_required_per_startup'] = 4

params['levelization_period_years'] = 60 # in years
params['refueling_period_days'] = 15
params['number_of_unanticipated_shutdowns_per_year']= 0.8

params['duration_to_startup_after_refueling_days'] = 7
params['duration_to_startup_after_shutdown_days'] = 14
params['reactors_monitored_by_one_person'] = 5 
params['FTEs_for_security_staff'] = 5 

params['people_by_days_refueling_per_year'], params['people_by_days_startup_per_year'],\
        params['capacity_factor'] = reactor_operation(params)


# **************************************************************************************************************************
#                                           Sec. 11 : Economic Parameters
# **************************************************************************************************************************
# preconstruction cost params

# A conservative estimate for the land area 
# Ref: McDowell, B., and D. Goodman. "Advanced Nuclear Reactor Plant Parameter Envelope and
#Guidance." National Reactor Innovation Center (NRIC), NRIC-21-ENG-0001 (2021). 

params['Land Area'] = 18 # acres
params['escalation_year'] = 2023
# excavation volume needs to be detailed
params['Excavation Volume'] = 412.605 # m^3

# reactor building
reactor_building = ['Reactor building', 2.6, 2.6,  11.0, 2, 2, 2]
RB_slab_roof_volume, RB_basemat_volume, RB_walls_volume = calculate_reactor_building_structure_volume(reactor_building)


# Energy conversion building
energy_conversion_building = ['Energy conversion building', 2.6, 6.0,   5.6, 2, 2, 2]
EB_slab_roof_volume, EB_basemat_volume, EB_walls_volume = calculate_energy_conversion_building_structure_volume(energy_conversion_building)


# control building
control_building = ['Control building', 2.4, 8.4,   3.0, 2, 2, 2]
cb_slab_roof_volume, cb_basemat_volume, cb_walls_volume = calculate_control_building_structure_volume(control_building)


# Refueling building
refueling_building = ['Refueling area building', 8.0, 9.0,   5.0, 2, 2, 2]
rb_slab_roof_volume, rb_basemat_volume, rb_walls_volume = calculate_refueling_building_structure_volume(refueling_building)


# spent fuel building
spent_fuel_building = ['Spent fuel building', 6.0, 8.0,   5.6, 2, 2, 2]
sfb_slab_roof_volume, sfb_basemat_volume, sfb_walls_volume = calculate_spent_fuel_building_structure_volume(spent_fuel_building)


# emergency building
emergency_building = ['Emergency and startup power system building', 4.0, 4.0,   5.6, 2, 2, 2]
eb_slab_roof_volume, eb_basemat_volume, eb_walls_volume = calculate_emergency_building_structure_volume(emergency_building)

# storage building
storage_building = ['Storage and warehouse building', 6.0, 6.0,   5.6, 2, 2, 2]
sb_slab_roof_volume, sb_basemat_volume, sb_walls_volume = calculate_emergency_building_structure_volume(storage_building)


# radwaste building
radwaste_storage_building = ['Radwaste building', 6.0, 6.0,   5.6, 2, 2, 2]
radb_slab_rood_volume, radb_basemat_volume, radb_walls_volume = calculate_radwaste_building_structure_volume(radwaste_storage_building)


# Financing params
params['interest_rate'] = 0.06 # 

# *************************************************************************************************************************
#                                           Sec. 12 : Cost
# **************************************************************************************************************************
Cost_estimate = bottom_up_cost_estimate('cost/Cost_Database.xlsx', params) 
print(Cost_estimate.to_string(index=False))

# **************************************************************************************************************************
#                                           Sec. 13 : Post Processing
# **************************************************************************************************************************

# params.show_summary(show_metadata=True, sort_by='time') # print all the parameters

#Calculate the code execution time
elapsed_time = (time.time() - time_start) / 60
print('Execution time:', np.round(elapsed_time, 1), 'minutes')
