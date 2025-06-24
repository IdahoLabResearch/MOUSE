"""
This example demonstrates a detailed bottom-up cost estimate for a 15 MWe Gas Cooled Triso-Fueled Microreactor (CGMR).
OpenMC is used for core design calculations.
Back-of-the-envelope calculations were conducted for other Balance of Plant components.
Users can modify any of the parameters in the "params" dictionary below.
"""


# Import necessary libraries
import numpy as np
import time
import warnings
import traceback

# Import custom modules and functions
import watts  # Import watts library for simulation workflows for one or multiple codes
import openmc  # Import OpenMC for nuclear reactor physics simulations

# Import templates and utility functions from the core design module
from core_design.openmc_template_GCMR import *
from core_design.drums import *

# # Import engineering evaluation tools and functions
from reactor_engineering_evaluation.tools import *

from reactor_engineering_evaluation.fuel_calcs import *
from reactor_engineering_evaluation.vessels_calcs import *
from reactor_engineering_evaluation.BOP import *

# Import cost estimation functions
from cost.baseline_costs import *
from cost.cost_utils import *
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
params['plotting'] = "yes"  # options are "yes" and "no"

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

# The user can change any of the following materials,
# but the new material has to be included in the materials database at "core_design/openmc_materials_database"
params['reactor type'] = "GCMR"
params['TRISO Fueled'] = "Yes" 

# The fuel and its properties
params['fuel'] = 'un'
params['enrichment'] = 0.19 # The enrichment is a fraction. It has to be between 0 and 1

# Mixing UO2 and UC by atom fraction
params['UO2 atom fraction'] = 0.7
params['Reflector'] = 'BeO'
params['Moderator'] = 'Graphite'
params['moderator_booster'] = 'ZrH' 
params['Coolant'] = 'Helium'

# Temperature assumed for the OpenMC calculations
params['common_temperature'] = 850 # Kelvins

# **************************************************************************************************************************
#                                           Sec. 3 : User-Defined Parameters (Geometry: Fuel Pins, Moderator Pins, Coolant)
# **************************************************************************************************************************  

# The reactor core includes fuel pins and moderator pins
# The fuel pin is assumed to include 5 regions: kernel, buffer, 3 layers of ceramics
# The materials of the 5 regions of the fuel pin
params['fuel_pin_materials'] = ['un', 'buffer', 'PyC', 'SiC', 'PyC']

# The outer radii of each of the 5 regions of the fuel pin
params['fuel_pin_radii'] = [0.025, 0.035, 0.039, 0.0425, 0.047] # cm


# A background material that is used anywhere within the lattice but outside a TRISO particle
# The background material (params['matrix material']) is within the compact region where the TRISO particles are packed.
# The moderator (params['moderator']) is outside this compact region but within the larger cylindrical fuel region.
params['matrix material'] = 'Graphite'

# The dimensions of the cylinderical lattice containing TRISO particles
params['lattice height'] = 4 # cm

# The radius of the area that is occupied by the TRISO particles
params['lattice radius'] = 0.6225 # cm
params['lattice_compact_volume'] = cylinder_volume(params['lattice radius'],params['lattice height'])
params['packing factor'] = 0.3 

# Coolant Channel
params['coolant channel radius'] = 0.35 #cm
params['booster radius'] =  0.55 # cm

# Hexagonal Lattice
params['hex lattice_radius']  = 1.3 #cm

# the distance between the centers of adjacent hexagons within the lattice structure
params['lattice_pitch'] =  params['hex lattice_radius']*np.sqrt(3)

# **************************************************************************************************************************
#                                           Sec. 4 : User-Defined Parameters (Hexagonal Lattice)
# **************************************************************************************************************************  

# The number of rings (pins) along the edge of the hexagonal lattice
params['assembly_rings'] = 6

# the height of the hexagonal of one fuel assembly
params['assembly_ftf'] = params['lattice_pitch']*(params['assembly_rings']-1)*np.sqrt(3)

# #The thickness of the reflector around the lattice
params['core_rings'] = 5

# extra thickness for the reflector beyond the control drums
params['extra_reflector'] = params['assembly_ftf'] / 10 #cm
params['core_radius'] = params['assembly_ftf']* params['core_rings'] +  params['extra_reflector']
params['active_height'] = 2 * params['core_radius']

# **************************************************************************************************************************
#                                           Sec. 5 : User-Defined Parameters (Control Drums)
# ************************************************************************************************************************** 

# Each of the control drums includes a reflector material (most of the drum) and an absorber material (the outer layer)
params['Control Drum Reflector'] = 'BeO'   # The reflector material in the control drums  
params['Control Drum Absorber'] = 'B4C_enriched' # The absorber material in the control drums
params['Control Drum Count'] = 24
params['Drum_Radius'] = 9 #cm  

# The thickness of the absorber layer
params['drum_Absorber_thickness'] = 1 # cm

# # If there are two drums alone each edge of the hexagonal lattice, the angle between drums pairs will be 60 degrees
# params['angle_between_drums_pairs'] = 60 # degrees

# # Placement of drums happen by tracing a line through the core apothems
# # then 2 drums are place around each apothem by deviating from this line
# # by a deviation angle
# params["deviation angle between drums"] = 12.86 # degrees

# The radius of the tube of the control drum
params['drum_tube_radius'] = params['Drum_Radius'] +(params['Drum_Radius']/ 45)  # cm

params['drum_height_to_lattice_height'] = 1.24
params['drum_height'] = 1 * params['active_height']

params['all_drums_volume'], params['Control Drum Absorber Mass'], params['Control Drum Reflector Mass'], params['Control Drums Mass'] =\
    calculate_drum_volume(params['Drum_Radius'], params['drum_height'],\
        params['drum_Absorber_thickness'], params)

params['tot_drum_area_all'] = params['all_drums_volume'] / params['drum_height']

params['Reflector Mass'] = calculate_reflector_mass_GCMR(params)          
params['Moderator Mass'], params['Moderator Booster Mass'] = calculate_moderator_mass_GCMR(params) 

# **************************************************************************************************************************
#                                           Sec. 6 : User-Defined Parameters (Overall System)
# ************************************************************************************************************************** 

# # The reactor power (thermal) MW
params['Power MWt'] = 15 # MWt
params['thermal_efficiency'] = 0.4
params['Power MWe'] = params['Power MWt'] * params['thermal_efficiency']
params['Power kWe'] = params['Power MWe'] * 1e3
# TEMPRARY: NEED TO CHANGE!!!

# # The actual heat flux (MW/m^2)
params['heat_flux'] = 0.8 # calculate_heat_flux(params['fuel_pin_radii'][-1], params['lattice_height'], params['rings'], params['power_MW_th'])
# # Target Heat Flux : Approximate calculated value for a typical sodium-cooled fast reactor (SFR)
params['heat_flux_criteria'] = 0.9

# # The burnup steps for depletion calculations
params['burnup_steps_MWd_per_Kg'] = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0,
                                     30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0]



# **************************************************************************************************************************
#                                           Sec. 7 : Running OpenMC
# ************************************************************************************************************************** 

# A function that runs OpenMC and the depletion calculations
# If the heat flux is too high, the code stops running
# The results such as the heat flux, the fuel lifetime and the mass of the fuel are added to the "params"

try:
    openmc_plugin = watts.PluginOpenMC(build_openmc_model_GCMR, show_stderr=True)  # running the LTMR Model
    # monitor_heat_flux(params, openmc_plugin)


except Exception as e:
    print("An error occurred:")
    traceback.print_exc()

    # # Handle any errors that occur during the simulation
    # print("\n")  # Blank line before
    # print(f"\033[91mAn error occurred while running the OpenMC simulation: {e}\033[0m")
    # print("\n")  # Blank line after

# # ## TEMPORARY  ## DELETE LATER!!!!!!!!!!!!!!!!!!!!
params['fuel_lifetime_days'] =1305 # days
params['mass_U235'] = 61975 # grams
params['mass_U238'] = 263372.87  # grams
params['Uranium Mass'] = (params['mass_U235'] + params['mass_U238']) / 1000 # Kg


# **************************************************************************************************************************
#                                           Sec. 7 : Fuel Calcs
# ************************************************************************************************************************** 

params['Natural Uranium Mass'], params['fuel_tail_waste_mass_Kg'], params['SWU'] =\
        fuel_calculations(params)

# **************************************************************************************************************************
#                                           Sec. 7 : Primary Loop / Balance of Plant
# ************************************************************************************************************************** 

# Primary Loop Specifications
params['Primary Loop Count'] = 2
params['Primary Loop per loop load fraction'] = 0.5
params['Primary Loop Power MWt'] = params['Power MWt'] * params['Primary Loop per loop load fraction']
params['Primary Loop Purification'] = False

# Equipments
params['HX Material']               = 'SS316'
params['Primary HX Mass']           = 1000 * 0.8 * params['Primary Loop Power MWt'] # [Kg] TODO: replace with better estimate
params['Secondary HX Mass']         = 0
params['Compressor Isentropic Efficiency'] = 0.8

# Coolant Temperature Differnce between the inlet and the outlet (reactor size)
params['Primary Loop Outlet Temperature'] = 550 + 273.15 # K
params['Primary Loop Temperature Difference'] = 250 
params['Primary Loop Mass Flow Rate'] = mass_flow_rate(params['Primary Loop Power MWt'],
                                                       params['Primary Loop Temperature Difference'], 
                                                       params['Coolant'])
params['Primary Loop Pressure Drop'] = 50e3 # [Pa] TODO: implement some rough TH code for dP estimation
params['Primary Loop Compressor Power'] = calculate_circulator_mechanical_power(params)

params['BoP Count'] = 2
params['BoP per loop load fraction'] = 0.5
params['BoP Power kWe'] = params['Power kWe'] * params['BoP per loop load fraction']

# **************************************************************************************************************************
#                                           Sec. 8 : Shielding
# ************************************************************************************************************************** 
# There is no in-vessel sheilding thickness
params['in_vessel_shield_thickness'] = 0 #cm
params['in_vessel_shielding_inner_radius'] = params['core_radius'] 
params['in_vessel_shielding_outer_radius'] = params['core_radius'] + params['in_vessel_shield_thickness']
params['In Vessel Shielding Material'] = 'B4C_natural'

params['out_of_vessel_shield_thickness'] = 39.37 #cm
params['Out Of Vessel Shielding Material'] = 'WEP'
# The out of vessel shield is not fully made of the out of vessel material (e.g. WEP) so we use an effective density factor
params['out_vessel_shield_effective_density_factor'] = 0.5

# **************************************************************************************************************************
#                                           Sec. 9 : Vessels Calculations
# ************************************************************************************************************************** 
# # Vessels parameters
params['vessel_radius'] = params['core_radius'] + params['in_vessel_shield_thickness']   

params['vessel_thickness'] = 7.6 # cm
params['vessel_lower_plenum_height'] = 30 # cm
params['vessel_upper_plenum_height'] = 60 # cm
params['vessel_upper_gas_gap'] = 10 # cm
params['vessel_bottom_depth'] = 35  # cm
params['vessel_material'] ='stainless_steel'

# assuming that no guard vessel 
params['gap_between_vessel_and_guard_vessel'] = 2 # cm
params['guard_vessel_thickness'] = 5
params['guard_vessel_material'] ='stainless_steel'

params['gap_between_guard_vessel_and_cooling_vessel'] = 5 # cm
params['cooling_vessel_thickness'] = 0.5 # cm
params['cooling_vessel_material'] ='stainless_steel'

params['gap_between_cooling_vessel_and_intake_vessel'] = 0.3 # cm
params['intake_vessel_thickness'] = 0.5 # cm
params['intake_vessel_material'] ='stainless_steel'

# Calculating the vessels dimensions and masses
params['vessels_total_radius'], params['vessel_height'] , params['vessels_total_height'],\
        params['Vessel Mass'], params['Guard Vessel Mass'] ,\
            params['Cooling Vessel Mass'], params['Intake Vessel Mass'], params['Total Vessels Mass'] = vessels_specs(params)

params['In Vessel Shielding Mass'] = cylinder_annulus_mass(params['in_vessel_shielding_outer_radius'],\
    params['in_vessel_shielding_inner_radius'], params['vessel_height'], params['In Vessel Shielding Material'] )  

params['Out Of Vessel Shielding Mass'] = params['out_vessel_shield_effective_density_factor'] * cylinder_annulus_mass(params['out_of_vessel_shield_thickness']+ params['vessels_total_radius'],\
        params['out_of_vessel_shield_thickness'], params['vessels_total_height'], params['Out Of Vessel Shielding Material']) 
params ['Vessel and Guard Vessel Masses'] = params['Vessel Mass'] +  params['Guard Vessel Mass']

# **************************************************************************************************************************
#                                           Sec. 10 : Operation
# **************************************************************************************************************************

params['Operation Mode'] = "Autonomous" # "Non-Autonomous" or "Autonomous"
params['Number of Operators'] = 2

params['Levelization Period'] = 60 # [yr]
params['Refueling Period'] = 7 # [d]
params['Emergency Shutdowns Per Year']= 0.2

params['Startup Duration after Refueling'] = 2 # [d]
params['Startup Duration after Emergency Shutdown'] = 14
params['Reactors Monitored Per Operator'] = 10
params['Security Staff Per Shift'] = 1


# **************************************************************************************************************************
#                                           Sec. 11 : Economic Parameters
# **************************************************************************************************************************
# preconstruction cost params

# A conservative estimate for the land area 
# Ref: McDowell, B., and D. Goodman. "Advanced Nuclear Reactor Plant Parameter Envelope and
#Guidance." National Reactor Innovation Center (NRIC), NRIC-21-ENG-0001 (2021). 
params['Land Area'] = 18 # acres
params['escalation_year'] = 2024
# excavation volume needs to be detailed
params['Excavation Volume'] = 463.93388 # m3 

# Financing params
# Financing params
params['Interest Rate'] = 0.065 # 
params['Construction Duration'] = 12 # months 
params['Debt To Equity Ratio'] = 0.5 
params['Annual Return'] = 0.0475

params['Reactor Building Slab Roof Volume'] = 219.18168 # m^3
params['Reactor Building Basement Volume'] = 219.18168 # m^3
params['Reactor Building Exterior Walls Volume'] = 438.04376 # m^3

# Energy conversion building
params['Turbine Building Slab Roof Volume'] = 132 # m^3
params['Turbine Building Basement Volume'] = 132 # m^3
params['Turbine Building Exterior Wall'] = 192.64 # m^3

# control building
params['Control Building Slab Roof Volume'] = 8.1 # m^3
params['Control Building Basement Volume'] = 27 # m^3
params['Control Building Exterior Walls Volume'] = 19.44 # m^3

# Refueling building 
params['Refueling Building Slab Roof Volume'] = 0 # m^3
params['Refueling Building Basement Volume'] = 0 # m^3
params['Refueling Building Exterior Walls Volume'] = 0 # m^3

# spent fuel building
params['Spent Fuel Building Slab Roof Volume'] = 0 # m^3
params['Spent Fuel Building Basement Volume'] = 0
params['Spent Fuel Building Exterior Walls Volume'] = 0

params['Emergency Building Slab Roof Volume'] =  0
params['Emergency Building Basement Volume'] = 0
params['Emergency Building Exterior Walls Volume'] = 0

params['Storage Building Slab Roof Volume'] =  0
params['Storage Building Basement Volume'] =  0
params['Storage Building Exterior Walls Volume'] =  0

params['Radwaste Building Slab Roof Volume'] =  0
params['Radwaste Building Basement Volume'] =  0
params['Radwaste Building Exterior Walls Volume'] =  0

# **************************************************************************************************************************
#                                           Sec. 12 : Cost
# **************************************************************************************************************************
Cost_estimate = bottom_up_cost_estimate('cost/Cost_Database.xlsx', params) 
print(Cost_estimate.to_string(index=False))

# **************************************************************************************************************************
#                                           Sec. 13 : Post Processing
# **************************************************************************************************************************

# params.show_summary(show_metadata=True, sort_by='time') # print all the parameters

# #Calculate the code execution time
elapsed_time = (time.time() - time_start) / 60
print('Execution time:', np.round(elapsed_time, 2), 'minutes')