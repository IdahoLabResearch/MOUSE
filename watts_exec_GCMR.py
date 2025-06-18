"""
This script performs a bottom-up cost estimate for a Gas Cooled Microreactor (GCMR).
OpenMC is used for core design calculations, and other Balance of Plant components are estimated.
Users can modify parameters in the "params" dictionary below.
"""

import numpy as np
import watts  # Simulation workflows for one or multiple codes
from core_design.openmc_template_GCMR import *
from core_design.pins_arrangement import LTMR_pins_arrangement
from core_design.utils import *
from core_design.drums import *
from reactor_engineering_evaluation.fuel_calcs import fuel_calculations
from reactor_engineering_evaluation.BOP import *
from reactor_engineering_evaluation.vessels_calcs import *
from reactor_engineering_evaluation.tools import *
from cost.baseline_costs import bottom_up_cost_estimate

import warnings
warnings.filterwarnings("ignore")

import time
time_start = time.time()

params = watts.Parameters()

def update_params(updates):
    params.update(updates)

# **************************************************************************************************************************
#                                                Sec. 0: Settings
# **************************************************************************************************************************
update_params({
    'plotting': "Y",  # "Y" or "N": Yes or No
    'cross_sections_xml_location': '/home/hannbn/openmc_data/endfb-viii.0-hdf5/cross_sections.xml',
    'simplified_chain_thermal_xml': '/home/hannbn/openmc_data/simplified_thermal_chain11.xml'
})

# **************************************************************************************************************************
#                                                Sec. 1: Materials
# **************************************************************************************************************************
params['Control Drum Reflector'] = 'BeO'   # The reflector material in the control drums  

update_params({
    'reactor type': "GCMR",
    'TRISO Fueled': "Yes",
    'Fuel': 'un',
    'Enrichment': 0.19,  # The enrichment is a fraction. It has to be between 0 and 1
    'UO2 atom fraction': 0.7,  # Mixing UO2 and UC by atom fraction
    'Reflector': 'BeO',
    'Matrix Material': 'Graphite', #         # matrix material is a background material  within the compact fuel element between the TRISO particles
    'Moderator': 'Graphite', # # The moderator is outside this compact fuel region 
    'Moderator Booster': 'YHx',
    'Coolant': 'Helium',
    'Common Temperature': 850,  # Kelvins
    'Control Drum Absorber': 'B4C_natural',  # The absorber material in the control drums
    'Control Drum Reflector': 'BeO'  # The reflector material in the control drums
})
# **************************************************************************************************************************
#                                           Sec. 3 : User-Defined Parameters (Geometry: Fuel Pins, Moderator Pins, Coolant)
# **************************************************************************************************************************  

# The reactor core includes fuel pins and moderator pins
# The fuel pin is assumed to include 5 regions: kernel, buffer, 3 layers of ceramics

update_params({
    # fuel pin details
    'Fuel Pin Materials': ['un', 'buffer', 'PyC', 'SiC', 'PyC'],
    'Fuel Pin Radii': [0.025, 0.035, 0.039, 0.0425, 0.047],  # cm
    'Compact Fuel Radius': 0.6225,  # cm # The radius of the area that is occupied by the TRISO particles (fuel compact/ fuel element)
    'Packing Factor': 0.3,

    # Coolant channel and booster dimensions
    'Coolant Channel Radius': 0.35,  # cm
    'Booster Radius': 0.55,  # cm

    # Hexagonal lattice dimensions
    'Hex Lattice Radius':  1.3,  # cm
   
})

params['Lattice Pitch'] = params['Hex Lattice Radius']* np.sqrt(3)
params['Lattice Compact Volume'] =  cylinder_volume(params['Compact Fuel Radius'], 4),

# **************************************************************************************************************************
#                                           Sec. 4 : User-Defined Parameters (Hexagonal Lattice)
# **************************************************************************************************************************  

# The number of rings (pins) along the edge of the hexagonal lattice
params['assembly_rings'] = 6
# the height of the hexagonal of one fuel assembly
params['assembly_ftf'] = params['Lattice Pitch']*(params['assembly_rings']-1)*np.sqrt(3)

# #The thickness of the reflector around the lattice
params['core_rings'] = 5
# extra thickness for the reflector beyond the control drums
params['extra_reflector'] = params['assembly_ftf'] / 10 #cm
params['core_radius'] = params['assembly_ftf']* params['core_rings'] +  params['extra_reflector']
params['Active Height'] = 2 * params['core_radius']



# **************************************************************************************************************************
#                                           Sec. 5 : User-Defined Parameters (Control Drums)
# ************************************************************************************************************************** 

params['Drum Radius'] = 9 #cm  

# # Each of the control drums includes a reflector material (most of the drum) and an absorber material (the outer layer)
# # The thickness of the absorber layer
params['Drum Absorber Thickness'] = 1 # cm

# # If there are two drums alone each edge of the hexagonal lattice, the angle between drums pairs will be 60 degrees
# params['angle_between_drums_pairs'] = 60 # degrees

# # Placement of drums happen by tracing a line through the core apothems
# # then 2 drums are place around each apothem by deviating from this line
# # by a deviation angle
# params["deviation angle between drums"] = 12.86 # degrees

# The radius of the tube of the control drum
params['drum_tube_radius'] = params['Drum Radius'] +(params['Drum Radius']/ 45)  # cm

params['drum_height_to_lattice_height'] = 1.24
params['Drum Height'] = 1 * params['Active Height']

calculate_drums_volumes_and_masses(params)

params['Reflector Mass'] = calculate_reflector_mass_GCMR(params)          
params['Moderator Mass'] = calculate_moderator_mass_GCMR(params) 
# **************************************************************************************************************************
#                                           Sec. 6 : User-Defined Parameters (Overall System)
# ************************************************************************************************************************** 

# # The reactor power (thermal) MW
params['Power MWt'] = 15 # MWt
params['thermal_efficiency'] = 0.4
params['Power MWe'] = params['Power MWt'] * params['thermal_efficiency']
# TEMPRARY: NEED TO CHANGE!!!

# # The actual heat flux (MW/m^2)
params['Heat Flux'] = 0.8 # calculate_heat_flux(params['fuel_pin_radii'][-1], params['lattice_height'], params['rings'], params['power_MW_th'])
# # Target Heat Flux : Approximate calculated value for a typical sodium-cooled fast reactor (SFR)
params['Heat Flux Criteria'] = 0.9

# # The burnup steps for depletion calculations
params['Burnup Steps'] = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0,
                                     30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0]



# **************************************************************************************************************************
#                                           Sec. 7 : Running OpenMC
# ************************************************************************************************************************** 

# A function that runs OpenMC and the depletion calculations
# If the heat flux is too high, the code stops running
# The results such as the heat flux, the fuel lifetime and the mass of the fuel are added to the "params"

heat_flux_monitor = monitor_heat_flux(params)
run_openmc(build_openmc_model_GCMR, heat_flux_monitor, params)

# # # ## TEMPORARY  ## DELETE LATER!!!!!!!!!!!!!!!!!!!!
# params['fuel_lifetime_days'] =1305 # days
# params['mass_U235'] = 61975 # grams
# params['mass_U238'] = 263372.87  # grams
# params['Uranium Mass'] = (params['mass_U235'] + params['mass_U238']) / 1000 # Kg

# # **************************************************************************************************************************
# #                                           Sec. 7 : Fuel Calcs
# # ************************************************************************************************************************** 

# params['Natural Uranium Mass'], params['fuel_tail_waste_mass_Kg'], params['SWU'] =\
#         fuel_calculations(params)

# # **************************************************************************************************************************
# #                                           Sec. 7 : Balance of Plant
# # ************************************************************************************************************************** 
# params['Primary HX Mass']  =  calculate_heat_exchanger_mass(params)[0] # Kg
# params['Secondary HX Mass'] = 0
# params['Compressor Pressure Ratio'] = 4 
# params['Compressor Isentropic Efficiency'] = 0.8
# # Coolant Temperature Differnce between the inlet and the outlet (reactor size)
# params['Coolant Temperature Difference'] = 360 #K or C
# params['Coolant Mass Flow Rate'] = mass_flow_rate( params['Power MWt'] , params['Coolant Temperature Difference'], params['Coolant'])
# # **************************************************************************************************************************
# #                                           Sec. 8 : Shielding
# # ************************************************************************************************************************** 

# # #Shielding
# # There is no in-vessel sheilding thickness
# params['in_vessel_shield_thickness'] = 0 #cm
# params['in_vessel_shielding_inner_radius'] = params['core_radius'] 
# params['in_vessel_shielding_outer_radius'] = params['core_radius'] + params['in_vessel_shield_thickness']
# params['In Vessel Shielding Material'] = 'B4C_natural'


# params['out_of_vessel_shield_thickness'] = 39.37 #cm
# params['Out Of Vessel Shielding Material'] = 'WEP'
# # The out of vessel shield is not fully made of the out of vessel material (e.g. WEP) so we use an effective density factor
# params['out_vessel_shield_effective_density_factor'] = 0.5

# # **************************************************************************************************************************
# #                                           Sec. 9 : Vessels Calculations
# # ************************************************************************************************************************** 
# # # Vessels parameters
# params['vessel_radius'] = params['core_radius'] + params['in_vessel_shield_thickness']   

# params['vessel_thickness'] = 7.6 # cm
# params['vessel_lower_plenum_height'] = 30 # cm
# params['vessel_upper_plenum_height'] = 60 # cm
# params['vessel_upper_gas_gap'] = 10 # cm
# params['vessel_bottom_depth'] = 35  # cm
# params['vessel_material'] ='stainless_steel'

# # assuming that no guard vessel 
# params['gap_between_vessel_and_guard_vessel'] = 0 # cm
# params['guard_vessel_thickness'] = 0
# params['guard_vessel_material'] ='stainless_steel'

# params['gap_between_guard_vessel_and_cooling_vessel'] = 5 # cm
# params['cooling_vessel_thickness'] = 0.5 # cm
# params['cooling_vessel_material'] ='stainless_steel'

# params['gap_between_cooling_vessel_and_intake_vessel'] = 0.3 # cm
# params['intake_vessel_thickness'] = 0.5 # cm
# params['intake_vessel_material'] ='stainless_steel'

# # Calculating the vessels dimensions and masses
# params['vessels_total_radius'], params['vessel_height'] , params['vessels_total_height'],\
#         params['Vessel Mass'], params['Guard Vessel Mass'] ,\
#             params['Cooling Vessel Mass'], params['Intake Vessel Mass'], params['Total Vessels Mass'] = vessels_specs(params)

# params['In Vessel Shielding Mass'] = cylinder_annulus_mass(params['in_vessel_shielding_outer_radius'],\
#     params['in_vessel_shielding_inner_radius'], params['vessel_height'], params['In Vessel Shielding Material'] )  


# params['Out Of Vessel Shielding Mass'] = params['out_vessel_shield_effective_density_factor'] * cylinder_annulus_mass(params['out_of_vessel_shield_thickness']+ params['vessels_total_radius'],\
#         params['out_of_vessel_shield_thickness'], params['vessels_total_height'], params['Out Of Vessel Shielding Material']) 
# params ['Vessel and Guard Vessel Masses'] = params['Vessel Mass'] +  params['Guard Vessel Mass']

# # **************************************************************************************************************************
# #                                           Sec. 10 : Operation
# # **************************************************************************************************************************

# params['Operation Mode'] = "Autonomous" # "Non-Autonomous" or "Autonomous"
# params['Number of Operators'] = 2

# params['Levelization Period'] = 60 # in years
# params['Refueling Period'] = 7
# params['Emergency Shutdowns Per Year']= 0.2

# params['Startup Duration after Refueling'] = 2
# params['Startup Duration after Emergency Shutdown'] = 14
# params['Reactors Monitored Per Operator'] = 10
# params['Security Staff Per Shift'] = 1


# # **************************************************************************************************************************
# #                                           Sec. 11 : Economic Parameters
# # **************************************************************************************************************************
# # preconstruction cost params

# # A conservative estimate for the land area 
# # Ref: McDowell, B., and D. Goodman. "Advanced Nuclear Reactor Plant Parameter Envelope and
# #Guidance." National Reactor Innovation Center (NRIC), NRIC-21-ENG-0001 (2021). 

# params['Land Area'] = 18 # acres
# params['escalation_year'] = 2023
# # excavation volume needs to be detailed
# params['Excavation Volume'] = 463.93388 # m3 
# # Financing params
# # Financing params
# params['Interest Rate'] = 0.065 # 
# params['Construction Duration'] = 12 # months 
# params['Debt To Equity Ratio'] = 0.5 

# params['Annual Return'] = 0.0475

# params['Reactor Building Slab Roof Volume'] = 219.18168 # m^3
# params['Reactor Building Basement Volume'] = 219.18168 # m^3
# params['Reactor Building Exterior Walls Volume'] = 438.04376 # m^3

# # Energy conversion building
# params['Turbine Building Slab Roof Volume'] = 132 # m^3
# params['Turbine Building Basement Volume'] = 132 # m^3
# params['Turbine Building Exterior Wall'] = 192.64 # m^3

# # control building
# params['Control Building Slab Roof Volume'] = 8.1 # m^3
# params['Control Building Basement Volume'] = 27 # m^3
# params['Control Building Exterior Walls Volume'] = 19.44 # m^3

# # Refueling building 
# params['Refueling Building Slab Roof Volume'] = 312 # m^3
# params['Refueling Building Basement Volume'] = 312 # m^3
# params['Refueling Building Exterior Walls Volume'] = 340 # m^3

# # spent fuel building
# params['Spent Fuel Building Slab Roof Volume'] = 384 # m^3
# params['Spent Fuel Building Basement Volume'] = 384
# params['Spent Fuel Building Exterior Walls Volume'] = 448

# params['Emergency Building Slab Roof Volume'] =  128
# params['Emergency Building Basement Volume'] = 128
# params['Emergency Building Exterior Walls Volume'] = 180


# params['Storage Building Slab Roof Volume'] =  180
# params['Storage Building Basement Volume'] =  180
# params['Storage Building Exterior Walls Volume'] =  246.4

# params['Radwaste Building Slab Roof Volume'] =  280
# params['Radwaste Building Basement Volume'] =  280
# params['Radwaste Building Exterior Walls Volume'] =  358

# # **************************************************************************************************************************
# #                                           Sec. 12 : Cost
# # **************************************************************************************************************************
# Cost_estimate = bottom_up_cost_estimate('cost/Cost_Database.xlsx', params) 
# print(Cost_estimate.to_string(index=False))

# # **************************************************************************************************************************
# #                                           Sec. 13 : Post Processing
# # **************************************************************************************************************************

# params.show_summary(show_metadata=True, sort_by='time') # print all the parameters

# # #Calculate the code execution time
# elapsed_time = (time.time() - time_start) / 60
# print('Execution time:', np.round(elapsed_time, 2), 'minutes')