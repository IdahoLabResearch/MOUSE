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

# Import custom modules and functions
import watts  # Import watts library for simulation workflows for one or multiple codes
import openmc  # Import OpenMC for nuclear reactor physics simulations

# Import templates and utility functions from the core design module
from core_design.openmc_template_GCMR import *
# from core_design.utils import *
# from core_design.pins_arrangement import rings_1

# # Import engineering evaluation tools and functions
from reactor_engineering_evaluation.tools import *
from reactor_engineering_evaluation.operation import *
from reactor_engineering_evaluation.fuel_calcs import *
from reactor_engineering_evaluation.vessels_calcs import *
from reactor_engineering_evaluation.BOP import *

# Import cost estimation functions
# from cost.cost_scaling import cost_estimate  # Import cost estimation function

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
params['plotting'] = "no"  # options are "yes" and "no"

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

# The fuel and its properties
params['fuel'] = 'un'
params['enrichment'] = 0.19 # The enrichment is a fraction. It has to be between 0 and 1
# Mixing UO2 and UC by atom fraction
params['UO2 atom fraction'] = 0.7

# params['poison'] = 'b4c'
# The boron is enriched to 18.7% in the isotope B-10, specified by enrichment_target='B10' and the enrichment type wt%.
params['boron enrichment'] = 0.187

params['reflector'] = 'BeO'
params['moderator'] = 'graphite'
params['moderator_booster'] = 'YHx' 
params['coolant'] = 'helium'

# Temperature assumed for the OpenMC calculations
params['common_temperature'] = 850 # Kelvins

params['control_drum_absorber'] = 'b4c' # The absorber material in the control drums
params['control_drum_reflector'] = 'BeO'   # The reflector material in the control drums  


# **************************************************************************************************************************
#                                           Sec. 3 : User-Defined Parameters (Geometry: Fuel Pins, Moderator Pins, Coolant)
# **************************************************************************************************************************  

# The reactor core includes fuel pins and moderator pins
# The fuel pin is assumed to include 5 regions: kernel, buffer, 3 layers of ceramics

# The materials of the 5 regions of the fuel pin
params['fuel_pin_materials'] = ['un', 'buffer', 'PyC', 'SiC', 'PyC']
# The outer radii of each of the 5 regions of the fuel pin
params['fuel_pin_radii'] = [0.025, 0.035, 0.039, 0.0425, 0.047] # cm


#A background material that is used anywhere within the lattice but outside a TRISO particle
params['matrix material'] = 'graphite'
# The dimensions of the cylinderical lattice containing TRISO particles
params['lattice height'] = 4 # cm
params['lattice radius'] = 0.6225 # cm
params['lattice_compact_volume'] = cylinder_volume(params['lattice radius'],params['lattice height'])
params['packing factor'] = 0.3 

params['small coolant channel radius'] = 0.35 #cm
params['big coolant channel radius'] = 0.5 #cm
params['booster radius'] = 0.55 # cm
#Burnable Poison
params['poison radius'] = 0.6375 # cm

# Hexagonal Lattice
params['hex lattice_radius']  = 1.3 #cm
params['lattice_pitch'] =  params['hex lattice_radius']*np.sqrt(3)



# **************************************************************************************************************************
#                                           Sec. 4 : User-Defined Parameters (Hexagonal Lattice)
# **************************************************************************************************************************  

# The number of rings (pins) along the edge of the hexagonal lattice
params['assembly_rings'] = 6
params['assembly_ftf'] = params['lattice_pitch']*(params['assembly_rings']-1)*np.sqrt(3)
# # The previous parameters are used to calculate the hexagonal lattice radius and height
# params['lattice_radius'] = calculate_lattice_radius(params['fuel_pin_radii'][-1], params["pin_gap_distance"], params['assembly_rings'])
# params['lattice_height'] = 2 * params['lattice_radius']

# # The selected fuel pins/moerator pins arrangement
# params['rings'] = rings_1
# # The total number of fuel pins in the lattice
# params['fuel_pin_count'] = sum(row.count("FUEL") for row in params['rings'])

# #The thickness of the reflector around the lattice
params['core_rings'] = 5
params['extra_reflector'] = params['assembly_ftf'] / 10 #cm
params['core_radius'] =params['assembly_ftf']* params['core_rings']  + params['extra_reflector']
params['active_height'] = 2 * params['core_radius']
# # The total core radius 
# params['core_radius'] = params['lattice_radius'] + params['extra_reflector']


# **************************************************************************************************************************
#                                           Sec. 5 : User-Defined Parameters (Control Drums)
# ************************************************************************************************************************** 

params['Drum_Radius'] = 9 #cm  

# # Each of the control drums includes a reflector material (most of the drum) and an absorber material (the outer layer)
# # The thickness of the absorber layer
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
params['drum_height'] = params['drum_height_to_lattice_height'] * params['active_height']

params['all_drums_volume'], params['drum_absorp_all_mass'], params['drum_refl_all_mass'] =\
    calculate_drum_volume(params['Drum_Radius'], params['drum_height'],\
        params['drum_Absorber_thickness'])

params['tot_drum_area_all'] = params['all_drums_volume'] / params['drum_height']

# **************************************************************************************************************************
#                                           Sec. 6 : User-Defined Parameters (Overall System)
# ************************************************************************************************************************** 

# # The reactor power (thermal) MW
params['power_MW_th'] = 15 # MWt

# TEMPRARY: NEED TO CHANGE!!!

# # The actual heat flux (MW/m^2)
params['heat_flux'] = 0.8 # calculate_heat_flux(params['fuel_pin_radii'][-1], params['lattice_height'], params['rings'], params['power_MW_th'])
# # Target Heat Flux : Approximate calculated value for a typical sodium-cooled fast reactor (SFR)
params['heat_flux_criteria'] = 0.9

# # The burnup steps for depletion calculations
params['burnup_steps_MWd_per_Kg'] = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0,
                                     30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0]


# TEMPORARY  ## DELETE LATER!!!!!!!!!!!!!!!!!!!!
# params['fuel_lifetime_days'] = 1508 # days
# params['mass_U235'] = 67711.4 # grams
# params['mass_U238'] = 278650.8  # grams
# **************************************************************************************************************************
#                                           Sec. 7 : Running OpenMC
# ************************************************************************************************************************** 

# A function that runs OpenMC and the depletion calculations
# If the heat flux is too high, the code stops running
# The results such as the heat flux, the fuel lifetime and the mass of the fuel are added to the "params"

try:
    openmc_plugin = watts.PluginOpenMC(build_openmc_model_GCMR, show_stderr=True)  # running the LTMR Model
    # monitor_heat_flux(params, openmc_plugin)

    # def run_func():
    #     print("done")
    # openmc_result = openmc_plugin(params, function=run_func) 
    # monitor_heat_flux(params, openmc_plugin)
   
#     # monitor_heat_flux(params, openmc_plugin)
except Exception as e:
    # Handle any errors that occur during the simulation
    print(f"An error occurred while running the OpenMC simulation: {e}")

# ## TEMPORARY  ## DELETE LATER!!!!!!!!!!!!!!!!!!!!
params['fuel_lifetime_days'] =1305 # days
params['mass_U235'] = 61975 # grams
params['mass_U238'] = 263372.87  # grams

# **************************************************************************************************************************
#                                           Sec. 7 : Fuel Calcs
# ************************************************************************************************************************** 

params['natural_U_mass_consumption_Kg'], params['fuel_tail_waste_mass_Kg'], params['SWU_kg'] =\
        fuel_calculations(params)

# **************************************************************************************************************************
#                                           Sec. 7 : Balance of Plant
# ************************************************************************************************************************** 
params['Primary HX Mass'] , params['Intermediate HX Mass'] =  calculate_heat_exchanger_mass(params) # Kg

# **************************************************************************************************************************
#                                           Sec. 8 : Shielding
# ************************************************************************************************************************** 

# #Shielding
params['in_vessel_shield_thickness'] = 10 #cm
params['in_vessel_shielding_inner_radius'] = params['core_radius'] 
params['in_vessel_shielding_outer_radius'] = params['core_radius'] + params['in_vessel_shield_thickness']
params['in_vessel_material'] = 'boron_carbide' 


params['out_of_vessel_shield_thickness'] = 39.37 #cm
params['out_vessel_shield_material'] = 'water_extended_polymer'
# The out of vessel shield is not fully made of the out of vessel material (e.g. WEP) so we use an effective density factor
params['out_vessel_shield_effective_density_factor'] = 0.5

# **************************************************************************************************************************
#                                           Sec. 9 : Vessels Calculations
# ************************************************************************************************************************** 
# # Vessels parameters
params['vessel_radius'] = params['core_radius'] + params['in_vessel_shield_thickness']   + 20 # cm
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
params['cooling_vessel_thickness'] = 5 # cm
params['cooling_vessel_material'] ='stainless_steel'

params['gap_between_cooling_vessel_and_intake_vessel'] = 0.3 # cm
params['intake_vessel_thickness'] = 0.3 # cm
params['intake_vessel_material'] ='stainless_steel'

# Calculating the vessels dimensions and masses
params['vessels_total_radius'], params['vessel_height'] , params['vessels_total_height'],\
        params['vessel_mass_kg'], params['guard_vessel_mass_kg'] ,\
            params['cooling_vessel_mass'], params['intake_vessel_mass_kg'] = vessels_specs(params)


params['in_vessel_shielding_mass'] = cylinder_annulus_mass(params['in_vessel_shielding_outer_radius'],\
    params['in_vessel_shielding_inner_radius'], params['vessel_height'], params['in_vessel_material'] )  

params['out_of_vessel_shielding_mass'] = params['out_vessel_shield_effective_density_factor'] * cylinder_annulus_mass(params['out_of_vessel_shield_thickness']+ params['vessels_total_radius'],\
        params['out_of_vessel_shield_thickness'], params['vessels_total_height'], params['out_vessel_shield_material']) 
    






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
#                                           Sec. 11 : Post Processing
# **************************************************************************************************************************

params.show_summary(show_metadata=True, sort_by='time') # print all the parameters

#Calculate the code execution time
elapsed_time = (time.time() - time_start) / 60
print('Execution time:', np.round(elapsed_time, 1), 'minutes')



















# params['reactor type'] = "LTMR"
# params['thermal_efficiency'] = 0.31
# params['hex_area'] = 2.598 * params['lattice_radius'] * params['lattice_radius']
# params['reflector_mass'] = calculate_reflector_mass(params['hex_area'],\
#     params['core_radius'], params['tot_drum_area_all'], params['drum_height'])



# # Function to run OpenMC plugin
# def run_openmc(params):
#     openmc_plugin = watts.PluginOpenMC(build_openmc_model, show_stderr=True)
    
#     def run_func():
#         openmc.run()
#         lattice_geometry = openmc.Geometry.from_xml()
#         settings = openmc.Settings.from_xml()
#         depletion_results = openmc_depletion(params, lattice_geometry, settings)
#         params['fuel_lifetime_days'] = depletion_results[0] # days
#         params['mass_U235'] = depletion_results[1] # grams
#         params['mass_U238'] = depletion_results[2] # grams

#     if params['heat_flux'] <= params['heat_flux_criteria']:
#         print(f"\033[92mHEAT FLUX is: {params['heat_flux']} MW/m^2.\033[0m")
#         openmc_result = openmc_plugin(params, function=run_func)
        
        
#         elapsed_time = (time.time() - time_start) / 60
#         print('Execution time:', np.round(elapsed_time, 1), 'minutes')
    
#     else:
#         print(f"\033[91mHIGH HEAT FLUX: {params['heat_flux']} MW/m^2.\033[0m")



# def design_and_cost_evaluations(params):
#     run_openmc(params)
# =    

#     # return detailed_cost
    
# design_and_cost_evaluations(params)    









# # Main execution flow

#     # coa_table = design_and_cost_evaluations(params)[0]
#     # cap_cost, ann_cost, ann_cost_levelized, lcoe1 = design_and_cost_evaluations(params)[1:]
    
#     # print(cap_cost, ann_cost, ann_cost_levelized, lcoe1)
#     # html_string = styler._repr_html_()
#     # print(html_string)
#     # print(tabulate(coa_table.fillna(''), headers='keys',  showindex=False)) #, headers='keys', tablefmt='pretty'))



#     #     # Save to an Excel file
#     # with pd.ExcelWriter('styled_dataframe.xlsx', engine='openpyxl') as writer:
#     #     styler.to_excel(writer, sheet_name='Sheet1',index=False)
    



