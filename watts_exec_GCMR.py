"""
This script performs a bottom-up cost estimate for a Gas Cooled Microreactor (GCMR).
OpenMC is used for core design calculations, and other Balance of Plant components are estimated.
Users can modify parameters in the "params" dictionary below.
"""

import numpy as np
import watts  # Simulation workflows for one or multiple codes
from core_design.openmc_template_GCMR import *
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
update_params({
    'reactor type': "GCMR",
    'TRISO Fueled': "Yes",
    'Fuel': 'UN',
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
#                                           Sec. 2: Geometry: Fuel Pins, Moderator Pins, Coolant, Hexagonal Lattice
# **************************************************************************************************************************  

update_params({
    # fuel pin details
    'Fuel Pin Materials': ['UN', 'buffer_graphite', 'PyC', 'SiC', 'PyC'],
    'Fuel Pin Radii': [0.025, 0.035, 0.039, 0.0425, 0.047],  # cm
    'Compact Fuel Radius': 0.6225,  # cm # The radius of the area that is occupied by the TRISO particles (fuel compact/ fuel element)
    'Packing Factor': 0.3,
    
    # Coolant channel and booster dimensions
    'Coolant Channel Radius': 0.35,  # cm
    'Moderator Booster Radius': 0.55, # cm
      'Lattice Pitch'  : 2.25,
      'Assembly Rings' : 6,
      'Core Rings' : 5,
})
params['Assembly FTF'] = params['Lattice Pitch']*(params['Assembly Rings']-1)*np.sqrt(3)
params['Reflector Thickness'] = params['Assembly FTF'] / 10  # cm
params['Core Radius'] = params['Assembly FTF']* params['Core Rings'] +  params['Reflector Thickness']
params['Active Height'] = 2 * params['Core Radius']
# **************************************************************************************************************************
#                                           Sec. 3: Control Drums
# ************************************************************************************************************************** 
update_params({
    'Drum Radius' : 9, #cm   
    'Drum Absorber Thickness': 1, # cm
    'Drum Height': params['Active Height']
    })

calculate_drums_volumes_and_masses(params)
calculate_reflector_mass_GCMR(params)          
calculate_moderator_mass_GCMR(params) 
# **************************************************************************************************************************
#                                           Sec. 4: Overall System
# ************************************************************************************************************************** 
update_params({
    'Power MWt': 15,  # MWt
    'Thermal Efficiency': 0.4,
    'Heat Flux Criteria': 0.9,  # MW/m^2 (This one needs to be reviewed)
    'Burnup Steps': [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0,
                                     30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0]  # MWd_per_Kg
    })

params['Power MWe'] = params['Power MWt'] * params['Thermal Efficiency'] 
params['Heat Flux'] =  calculate_heat_flux_TRISO(params) # MW/m^2
# **************************************************************************************************************************
#                                           Sec. 5: Running OpenMC
# ************************************************************************************************************************** 

heat_flux_monitor = monitor_heat_flux(params)
run_openmc(build_openmc_model_GCMR, heat_flux_monitor, params)

# **************************************************************************************************************************
#                                         Sec. 6: Balance of Plant
# ************************************************************************************************************************** 
params.update({
    'Primary HX Mass': calculate_heat_exchanger_mass(params)[0],  # Kg
    'Secondary HX Mass': 0,
    'Compressor Pressure Ratio': 4,
    'Compressor Isentropic Efficiency': 0.8,
    'Coolant Temperature Difference': 360  # Coolant Temperature Differnce between the inlet and the outlet (reactor side):  K or C
    })

# calculate coolant mass flow rate
mass_flow_rate(params)
# # **************************************************************************************************************************
# #                                           Sec. 8 : Shielding
# # ************************************************************************************************************************** 
update_params({
    'In Vessel Shield Thickness': 0,  # cm (no shield in vessel for GCMR)
    'In Vessel Shield Inner Radius': params['Core Radius'],
    'In Vessel Shield Material': 'B4C_natural',
    'Out Of Vessel Shield Thickness': 39.37,  # cm
    'Out Of Vessel Shield Material': 'WEP',
    'Out Of Vessel Shield Effective Density Factor': 0.5 # The out of vessel shield is not fully made of the out of vessel material (e.g. WEP) so we use an effective density factor

})
params['In Vessel Shield Outer Radius'] =  params['Core Radius'] + params['In Vessel Shield Thickness']

# **************************************************************************************************************************
#                                           Sec. 9 : Vessels Calculations
# ************************************************************************************************************************** 
update_params({
    'Vessel Radius': params['Core Radius'] +  params['In Vessel Shield Thickness'],
    'Vessel Thickness': 7.6,  # cm
    'Vessel Lower Plenum Height': 30,  # cm
    'Vessel Upper Plenum Height': 60,  # cm
    'Vessel Upper Gas Gap': 10,  # cm
    'Vessel Bottom Depth': 35,  # cm
    'Vessel Material': 'stainless_steel',
    # No guard vessel
    'Gap Between Vessel And Guard Vessel': 0,  # cm
    'Guard Vessel Thickness': 0,  # cm
    'Guard Vessel Material': 'stainless_steel',
    
    'Gap Between Guard Vessel And Cooling Vessel': 5,  # cm
    'Cooling Vessel Thickness': 0.5,  # cm
    'Cooling Vessel Material': 'stainless_steel',
    'Gap Between Cooling Vessel And Intake Vessel': 3,  # cm
    'Intake Vessel Thickness': 0.5,  # cm
    'Intake Vessel Material': 'stainless_steel'
})

vessels_specs(params)  # calculate the volumes and masses of the vessels
calculate_shielding_masses(params)  # calculate the masses of the shieldings

# # **************************************************************************************************************************
# #                                           Sec. 10 : Operation
# # **************************************************************************************************************************
update_params({
    'Operation Mode': "Autonomous", # "Non-Autonomous" or "Autonomous"
    'Number of Operators': 2,
    'Levelization Period': 60,  # years
    'Refueling Period': 7,
    'Emergency Shutdowns Per Year': 0.2,
    'Startup Duration after Refueling': 2,
    'Startup Duration after Emergency Shutdown': 14,
    'Reactors Monitored Per Operator': 10,
    'Security Staff Per Shift': 1
})

# **************************************************************************************************************************
#                                           Sec. 11 : Economic Parameters
# **************************************************************************************************************************

update_params({
    # A conservative estimate for the land area 
    # Ref: McDowell, B., and D. Goodman. "Advanced Nuclear Reactor Plant Parameter Envelope and
    #Guidance." National Reactor Innovation Center (NRIC), NRIC-21-ENG-0001 (2021). 
    'Land Area': 18,  # acres
    
    'Escalation Year': 2023,
    
    'Excavation Volume': 412.605,  # m^3
    'Reactor Building Slab Roof Volume': 87.12,  # m^3
    'Reactor Building Basement Volume': 87.12,  # m^3
    'Reactor Building Exterior Walls Volume': 228.8,  # m^3
    
    'Turbine Building Slab Roof Volume': 132,  # m^3
    'Turbine Building Basement Volume': 132,  # m^3
    'Turbine Building Exterior Walls Volume': 192.64,  # m^3
    
    'Control Building Slab Roof Volume': 8.1,  # m^3
    'Control Building Basement Volume': 27,  # m^3
    'Control Building Exterior Walls Volume': 19.44,  # m^3
    
    'Refueling Building Slab Roof Volume': 312,  # m^3
    'Refueling Building Basement Volume': 312,  # m^3
    'Refueling Building Exterior Walls Volume': 340,  # m^3
    
    'Spent Fuel Building Slab Roof Volume': 240,  # m^3
    'Spent Fuel Building Basement Volume': 240,  # m^3
    'Spent Fuel Building Exterior Walls Volume': 313.6,  # m^3
    
    'Emergency Building Slab Roof Volume': 128,  # m^3
    'Emergency Building Basement Volume': 128,  # m^3
    'Emergency Building Exterior Walls Volume': 180,  # m^3
    
    'Storage Building Slab Roof Volume': 200,  # m^3
    'Storage Building Basement Volume': 200,  # m^3
    'Storage Building Exterior Walls Volume': 268.8,  # m^3
    
    'Radwaste Building Slab Roof Volume': 200,  # m^3
    'Radwaste Building Basement Volume': 200,  # m^3
    'Radwaste Building Exterior Walls Volume': 268.8,  # m^3,
    
    'Interest Rate': 0.065,
    'Construction Duration': 12,  # months
    'Debt To Equity Ratio': 0.5,
    'Annual Return': 0.0475  # Annual return on decommissioning costs
})

# **************************************************************************************************************************
#                                           Sec. 13: Post Processing
# **************************************************************************************************************************

# Estimate costs using the cost database file and save the output to an Excel file
bottom_up_cost_estimate('cost/Cost_Database.xlsx', params, "output_GCMR.xlsx")
elapsed_time = (time.time() - time_start) / 60  # Calculate execution time
print('Execution time:', np.round(elapsed_time, 2), 'minutes')