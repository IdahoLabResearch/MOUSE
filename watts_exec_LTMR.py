"""
This example demonstrates a detailed bottom-up cost estimate for a 15 MWe Liquid Metal Thermal Microreactor (LTMR).
OpenMC is used for core design calculations.
Back-of-the-envelope calculations were conducted for other Balance of Plant components.
Users can modify any of the parameters in the "params" dictionary below.
"""

# # # Import necessary libraries
import numpy as np

import watts  # Import watts library for simulation workflows for one or multiple codes

# Import  core design module
from core_design.openmc_template_LTMR import *
from core_design.pins_arrangement import LTMR_pins_arrangement
from core_design.utils import *
from core_design.drums import *

from reactor_engineering_evaluation.fuel_calcs import fuel_calculations
from reactor_engineering_evaluation.BOP import *
from reactor_engineering_evaluation.vessels_calcs import * 
from reactor_engineering_evaluation.tools import *
from cost.baseline_costs import bottom_up_cost_estimate

# Suppress warnings
import warnings
warnings.filterwarnings("ignore")  # Ignore all warnings

# # Record the current time (in seconds)
import time
time_start = time.time()  # Store the current time to track the script's execution duration

# Initialize user parameters using the watts library
params = watts.Parameters()
# **************************************************************************************************************************
#                                                Sec. 1 : User Configuration
# **************************************************************************************************************************

# # Whether the user prefers to display the core design plots. 
# # Note: Plotting is not recommended for sensitivity analysis and parametric studies.
params['plotting'] = "Y"  # "Y" or "N" : Yes or No

# If the user is using OpenMC calculations, the location of the neutron data library and the depletion chain file must be specified. 
# from: https://openmc.org/official-data-libraries/
params['cross_sections_xml_location'] ='/home/hannbn/openmc_data/endfb-viii.0-hdf5/cross_sections.xml'
# from: https://openmc.org/depletion-chains/
params['simplified_chain_thermal_xml'] ='/home/hannbn/openmc_data/simplified_thermal_chain11.xml'

# **************************************************************************************************************************
#                                                Sec. 2 : User-Defined Parameters (Materials)
# **************************************************************************************************************************  
params['reactor type'] = "LTMR"
# # The user can change any of the following materials,
# # but the new material has to be included in the materials database at "core_design/openmc_materials_database"
# params['TRISO Fueled'] = "No"  # options are "Yes" and "No"

# # The fuel and its properties
params['Fuel'] = 'TRIGA_fuel'
params['Enrichment'] = 0.1975 # The enrichment is a fraction. It has to be between 0 and 1
params["H_Zr_ratio"] = 1.6 # The proportion of hydrogen atoms relative to zirconium atoms
params['U_met_wo'] = 0.3 # The ratio between the weight of Uranium and the total weight of the fuel (less than 1)

# # Coolant, reflector and control drums
params['Coolant'] = 'NaK'
params['Reflector'] = 'BeO' # external reflector around the core
params['Control Drum Absorber'] = 'B4C_natural' # The absorber material in the control drums
params['Control Drum Reflector'] = 'BeO'     # The reflector material in the control drums   

# # Temperature assumed for the OpenMC calculations
params['Common Temperature'] = 600 # Kelvins

# **************************************************************************************************************************
#                                           Sec. 3 : User-Defined Parameters (Geometry: Fuel Pins, Moderator Pins, Coolant)
# **************************************************************************************************************************  

# The reactor core includes fuel pins and moderator pins
# The fuel pin is assumed to include 5 regions:
# 1- a filler rod surrounded by 2- a gap, which is then surrounded by 3- the fuel. This is followed by 4- another gap, 
# which is finally surrounded by 5- the cladding.

# # The materials of the 5 regions of the fuel pin.      
params['Fuel Pin Materials'] = ['Zr', None, 'TRIGA_fuel', None, 'SS304']
# The outer radii of each of the 5 regions of the fuel pin
params['Fuel Pin Radii'] = [0.28575, 0.3175, 1.5113, 1.5367, 1.5875] # cm

# The moderator pin is assumed to include 2 regions: a moderator and a cladding
params['Moderator Pin Materials'] = ['ZrH', 'SS304']
params['Moderator'] = params['Moderator Pin Materials'][0]
# The inner radius of the moderator pin is specified while the outer one is the same as the fuel pin outer radius
params['Moderator Pin Inner Radius'] = 1.5367 # cm
params['Moderator Pin Radii'] = [params['Moderator Pin Inner Radius'], params['Fuel Pin Radii'][-1]]


# **************************************************************************************************************************
#                                           Sec. 4 : User-Defined Parameters (Hexagonal Lattice)
# **************************************************************************************************************************  

# The gap between pins (fuel pins or moderator pins)
params["Pin Gap Distance"] = 0.1 # cm

# # The selected fuel pins/moerator pins arrangement
params['Pins Arrangement'] = LTMR_pins_arrangement
# The number of rings (pins) along the edge of the hexagonal lattice.   
params['Number of Rings per Assembly'] = 12

# # The previous parameters are used to calculate the hexagonal lattice radius and height
params['Lattice Radius'] = calculate_lattice_radius(params)
params['Active Height'] = 2 * params['Lattice Radius'] 


# # The total number of fuel pins in the fuel assembly
params['Fuel Pin Count']      = calculate_pins_in_assembly(params, "FUEL")
params['Moderator Pin Count'] = calculate_pins_in_assembly(params, "MODERATOR")

params['Moderator Mass'] = calculate_moderator_mass(params)

# #The thickness of the reflector around the lattice
params['Reflector Thickness'] = 14 #cm
# # The total core radius 
params['Core Radius'] = params['Lattice Radius'] +params['Reflector Thickness'] 


# **************************************************************************************************************************
#                                           Sec. 5 : User-Defined Parameters (Control Drums)
# ************************************************************************************************************************** 

params['Drum Radius']= 0.23 * params['Lattice Radius'] # cm # Control drum radius

# # Each of the control drums includes a reflector material (most of the drum) and an absorber material (the outer layer)
# # The thickness of the absorber layer
params['Drum Absorber Thickness'] = 1 # cm

# # Placement of drums happen by tracing a line through the core apothems

params['Drum Height'] = params['Active Height']

# calculating the total volumes and masses of drums
calculate_drums_volumes_and_masses(params)
# calculate the relector mass
calculate_reflector_mass_LTMR(params)
# **************************************************************************************************************************
#                                           Sec. 6 : User-Defined Parameters (Overall System)
# ************************************************************************************************************************** 

# # The reactor power (thermal) MW
params['Power MWt'] = 20 # MWt
params['Thermal Efficiency'] = 0.31
params['Power MWe'] = params['Power MWt'] * params['Thermal Efficiency']

# The actual heat flux (MW/m^2)
params['Heat Flux'] = calculate_heat_flux(params)

# Target Heat Flux : Approximate calculated value for a typical sodium-cooled fast reactor (SFR)
params['Heat Flux Criteria'] =  0.9   # (MW/m^2)(MW/m^2)

# The burnup steps for depletion calculations
params['Burnup Steps'] = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0,
                                     30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 140.0] # MWd_per_Kg
# **************************************************************************************************************************
#                                           Sec. 7 : Running OpenMC
# ************************************************************************************************************************** 

# If the heat flux is too high, the code stops running
# The results such as the heat flux, the fuel lifetime and the mass of the fuel are added to the "params"

heat_flux_monitor = monitor_heat_flux(params)

# Run openmc
run_openmc(build_openmc_model_LTMR, heat_flux_monitor, params)

# # **************************************************************************************************************************
# #                                           Sec. 7 : Fuel Calcs
# # ************************************************************************************************************************** 

# calculate the masses and SWU of the fuel
fuel_calculations(params)

# # **************************************************************************************************************************
# #                                           Sec. 7 : Balance of Plant
# # ************************************************************************************************************************** 
params['Primary HX Mass']  =  calculate_heat_exchanger_mass(params)[0] # Kg
params['Secondary HX Mass'] = 0
params['Primary Pump'] = 'Yes' # options are "Yes" and "No"
params['Secondary Pump'] = 'No' # options are "Yes" and "No"
params['Primary Pump Mechanical Power']  = calculate_pump_mechanical_power(params)[0]
params['Pump Isentropic Efficiency'] = 0.8
# # **************************************************************************************************************************
# #                                           Sec. 8 : Shielding
# # ************************************************************************************************************************** 

# Shielding
params['In Vessel Shield Thickness'] = 10.16 #cm
params['In Vessel Shield Inner Radius'] = params['Core Radius'] 
params['In Vessel Shield Outer Radius Radius'] = params['Core Radius'] + params['In Vessel Shield Thickness']
params['In Vessel Shielding Material'] = 'B4C_natural' 


params['Out Of Vessel Shield Thickness'] = 39.37 #cm
params['Out Of Vessel Shield Material'] = 'WEP' # water extended polymer (WEP)
# The out of vessel shield is not fully made of the out of vessel material (e.g. WEP) so we use an effective density factor
params['Out Of Vessel Shield Effective Density Factor'] = 0.5

# **************************************************************************************************************************
#                                           Sec. 9 : Vessels Calculations
# ************************************************************************************************************************** 
# Vessels parameters
params['Vessel Radius'] = params['Core Radius'] + params['In Vessel Shield Thickness'] 

params['Vessel Thickness'] = 2 # cm
params['Vessel Lower Plenum Height'] = 30 # cm
params['Vessel Upper Plenum Height'] = 60 # cm
params['Vessel Upper Gas Gap'] = 10 # cm
params['Vessel Bottom Depth'] = 35  # cm
params['Vessel Material'] = 'stainless_steel'

params['Gap Between Vessel And Guard Vessel'] = 2 # cm
params['Guard Vessel Thickness'] = 0.5 # cm
params['Guard Vessel Material'] = 'stainless_steel'

params['Gap Between Guard Vessel And Cooling Vessel'] = 5 # cm
params['Cooling Vessel Thickness'] = 0.5 # cm
params['Cooling Vessel Material'] = 'stainless_steel'

params['Gap Between Cooling Vessel And Intake Vessel'] = 3 # cm
params['Intake Vessel Thickness'] = 0.5 # cm
params['Intake Vessel Material'] = 'stainless_steel'

# Calculating the vessels dimensions and masses (in Kg)
vessels_specs(params)
# The mass of the in-vessel and out of vessel shielding
calculate_shielding_masses(params)

# # **************************************************************************************************************************
# #                                         Sec. 10 : Operation
# # **************************************************************************************************************************

params['Operation Mode'] = "Autonomous" # "Non-Autonomous" or "Autonomous"
params['Number of Operators'] = 2

params['Levelization Period'] = 60 # in years
params['Refueling Period'] = 7
params['Emergency Shutdowns Per Year']= 0.2

params['Startup Duration after Refueling'] = 2
params['Startup Duration after Emergency Shutdown'] = 14
params['Reactors Monitored Per Operator'] = 10 # in case of autonomous control
params['Security Staff Per Shift'] = 1

# **************************************************************************************************************************
#                                           Sec. 11 : Economic Parameters
# **************************************************************************************************************************
# preconstruction cost params

# # A conservative estimate for the land area 
# # Ref: McDowell, B., and D. Goodman. "Advanced Nuclear Reactor Plant Parameter Envelope and
# # Guidance." National Reactor Innovation Center (NRIC), NRIC-21-ENG-0001 (2021). 

params['Land Area'] = 18 # acres
params['Escalation Year'] = 2023
# excavation volume needs to be detailed
params['Excavation Volume'] = 412.605 # m^3

# # reactor building
params['Reactor Building Slab Roof Volume'] = 87.12 # m^3
params['Reactor Building Basement Volume'] = 87.12 # m^3
params['Reactor Building Exterior Walls Volume'] = 228.8 # m^3

# Energy conversion building
params['Turbine Building Slab Roof Volume'] = 132 # m^3
params['Turbine Building Basement Volume'] = 132 # m^3
params['Turbine Building Exterior Wall'] = 192.64 # m^3
# control building
params['Control Building Slab Roof Volume'] = 8.1
params['Control Building Basement Volume'] = 27
params['Control Building Exterior Walls Volume'] = 19.44

# # Refueling building
params['Refueling Building Slab Roof Volume'] = 312
params['Refueling Building Basement Volume'] = 312
params['Refueling Building Exterior Walls Volume'] = 340

# spent fuel building
params['Spent Fuel Building Slab Roof Volume'] = 240
params['Spent Fuel Building Basement Volume'] = 240
params['Spent Fuel Building Exterior Walls Volume'] = 313.6

params['Emergency Building Slab Roof Volume'] =  128
params['Emergency Building Basement Volume'] = 128
params['Emergency Building Exterior Walls Volume'] = 180

params['Storage Building Slab Roof Volume'] =  200
params['Storage Building Basement Volume'] =  200
params['Storage Building Exterior Walls Volume'] =  268.8

params['Radwaste Building Slab Roof Volume'] =  200
params['Radwaste Building Basement Volume'] =  200
params['Radwaste Building Exterior Walls Volume'] =  268.8

# Financing params
params['Interest Rate'] = 0.065 # 
params['Construction Duration'] = 12 # months 
params['Debt To Equity Ratio'] = 0.5 


# Given that decommissioning costs are incurred over the life of the reactor it is important to represent
# them not as lump sum costs, but rather annual payments that are placed into a trust and earn a base return
# on the investments
params['Annual Return'] = 0.0475

# *************************************************************************************************************************
#                                           Sec. 12 : Cost
# **************************************************************************************************************************
Cost_estimate = bottom_up_cost_estimate('cost/Cost_Database.xlsx', params) 
Cost_estimate .to_excel("LTMR.xlsx", index=False)

# **************************************************************************************************************************
#                                          Sec. 13 : Post Processing
#**************************************************************************************************************************
params.show_summary(show_metadata=True, sort_by='time') # print all the parameters

#Calculate the code execution time
elapsed_time = (time.time() - time_start) / 60
print('Execution time:', np.round(elapsed_time, 1), 'minutes')