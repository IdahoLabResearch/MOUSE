"""
This script performs a bottom-up cost estimate for a Liquid Metal Thermal Microreactor (LTMR).
OpenMC is used for core design calculations, and other Balance of Plant components are estimated.
Users can modify parameters in the "params" dictionary below.
"""
import numpy as np
import watts  # Simulation workflows for one or multiple codes
from core_design.openmc_template_LTMR import *
from core_design.pins_arrangement import LTMR_pins_arrangement
from core_design.utils import *
from core_design.drums import *
from reactor_engineering_evaluation.fuel_calcs import fuel_calculations
from reactor_engineering_evaluation.BOP import *
from reactor_engineering_evaluation.vessels_calcs import *
from reactor_engineering_evaluation.tools import *
from cost.cost_estimation import detailed_bottom_up_cost_estimate

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
    'reactor type': "LTMR", # LTMR or GCMR
    'TRISO Fueled': "No",
    'Fuel': 'TRIGA_fuel',
    'Enrichment': 0.1975,  # Fraction between 0 and 1
    "H_Zr_ratio": 1.6,  # Proportion of hydrogen to zirconium atoms
    'U_met_wo': 0.3,  # Weight ratio of Uranium to total fuel weight (less than 1)
    'Coolant': 'NaK',
    'Reflector': 'BeO',
    'Control Drum Absorber': 'B4C_natural',
    'Control Drum Reflector': 'BeO',
    'Common Temperature': 600,  # Kelvins
    'HX Material': 'SS316'
})

# **************************************************************************************************************************
#                                           Sec. 2: Geometry: Fuel Pins, Moderator Pins, Coolant, Hexagonal Lattice
# **************************************************************************************************************************  

update_params({
    'Fuel Pin Materials': ['Zr', None, 'TRIGA_fuel', None, 'SS304'],
    'Fuel Pin Radii': [0.28575, 0.3175, 1.5113, 1.5367, 1.5875],  # cm
    'Moderator Pin Materials': ['ZrH', 'SS304'],
    'Moderator': 'ZrH',  # params['Moderator Pin Materials'][0]
    'Moderator Pin Inner Radius': 1.5367,  # cm
    'Moderator Pin Radii': [1.5367, 1.5875],  # [params['Moderator Pin Inner Radius'], params['Fuel Pin Radii'][-1]]
    "Pin Gap Distance": 0.1,  # cm
    'Pins Arrangement': LTMR_pins_arrangement,
    'Number of Rings per Assembly': 12,
    'Reflector Thickness': 14,  # cm
    'Axial Reflector Thickness': 0 #cm
})

params['Lattice Radius'] = calculate_lattice_radius(params)
params['Active Height']  = 2 * params['Lattice Radius']
params['Fuel Pin Count'] = calculate_pins_in_assembly(params, "FUEL")
params['Moderator Pin Count'] =  calculate_pins_in_assembly(params, "MODERATOR")
params['Moderator Mass'] = calculate_moderator_mass(params)
params['Core Radius'] = params['Lattice Radius'] + 14  # params['Reflector Thickness']

# **************************************************************************************************************************
#                                           Sec. 3: Control Drums
# ************************************************************************************************************************** 

update_params({
    'Drum Radius': 0.23 * params['Lattice Radius'],  # cm
    'Drum Absorber Thickness': 1,  # cm
    'Drum Height': params['Active Height']
})

calculate_drums_volumes_and_masses(params)
calculate_reflector_mass_LTMR(params)

# **************************************************************************************************************************
#                                           Sec. 4: Overall System
# ************************************************************************************************************************** 

update_params({
    'Power MWt': 20,  # MWt
    'Thermal Efficiency': 0.31,
    'Heat Flux Criteria': 0.9,  # MW/m^2
    'Burnup Steps': [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 140.0]  # MWd_per_Kg
})
params['Power MWe'] = params['Power MWt'] * params['Thermal Efficiency']
params['Heat Flux'] =  calculate_heat_flux(params)
# **************************************************************************************************************************
#                                           Sec. 5: Running OpenMC
# ************************************************************************************************************************** 

heat_flux_monitor = monitor_heat_flux(params)
run_openmc(build_openmc_model_LTMR, heat_flux_monitor, params)
fuel_calculations(params)  # calculate the fuel mass and SWU

# **************************************************************************************************************************
#                                           Sec. 6: Primary Loop + Balance of Plant
# ************************************************************************************************************************** 

update_params({
    'Secondary HX Mass': 0,
    'Primary Pump': 'Yes',
    'Secondary Pump': 'No',
    'Pump Isentropic Efficiency': 0.8,
    'Primary Loop Inlet Temperature': 430 + 273.15, # K
    'Primary Loop Outlet Temperature': 520 + 273.15, # K
    'Secondary Loop Inlet Temperature': 395 + 273.15, # K
    'Secondary Loop Outlet Temperature': 495 + 273.15, # K,
})

params['Primary HX Mass'] = calculate_heat_exchanger_mass(params)  # Kg
# Update BoP Parameters
params.update({
    'BoP Count': 2, # Number of BoP present in plant
    'BoP per loop load fraction': 0.5, # based on assuming that each BoP Handles the total load evenly (1/2)
    })
params['BoP Power kWe'] = 1000 * params['Power MWe'] * params['BoP per loop load fraction']
# calculate coolant mass flow rate
mass_flow_rate(params)
calculate_primary_pump_mechanical_power(params)

# **************************************************************************************************************************
#                                           Sec. 7: Shielding
# ************************************************************************************************************************** 

update_params({
    'In Vessel Shield Thickness': 10.16,  # cm
    'In Vessel Shield Inner Radius': params['Core Radius'],
    'In Vessel Shield Material': 'B4C_natural',
    'Out Of Vessel Shield Thickness': 39.37,  # cm
    'Out Of Vessel Shield Material': 'WEP',
    'Out Of Vessel Shield Effective Density Factor': 0.5 # The out of vessel shield is not fully made of the out of vessel material (e.g. WEP) so we use an effective density factor
})

params['In Vessel Shield Outer Radius'] =  params['Core Radius'] + params['In Vessel Shield Thickness']

# **************************************************************************************************************************
#                                           Sec. 8: Vessels Calculations
# ************************************************************************************************************************** 

update_params({
    'Vessel Radius': params['Core Radius'] +  params['In Vessel Shield Thickness'],
    'Vessel Thickness': 2,  # cm
    'Vessel Lower Plenum Height': 30,  # cm
    'Vessel Upper Plenum Height': 60,  # cm
    'Vessel Upper Gas Gap': 10,  # cm
    'Vessel Bottom Depth': 35,  # cm
    'Vessel Material': 'stainless_steel',
    'Gap Between Vessel And Guard Vessel': 2,  # cm
    'Guard Vessel Thickness': 0.5,  # cm
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

# **************************************************************************************************************************
#                                           Sec. 9: Operation
# **************************************************************************************************************************

update_params({
    'Operation Mode': "Autonomous",
    'Number of Operators': 2,
    'Levelization Period': 60,  # years
    'Refueling Period': 7,
    'Emergency Shutdowns Per Year': 0.2,
    'Startup Duration after Refueling': 2,
    'Startup Duration after Emergency Shutdown': 14,
    'Reactors Monitored Per Operator': 10,
    'Security Staff Per Shift': 1
})
## Calculated based on 1 tanks
## Density=855  kg/m3, Volume=8.2402 m3 (standard tank size?)
params['Onsite Coolant Inventory'] = 1 * 855 * 8.2402 # kg
params['Annual Coolant Supply Frequency'] = 0.1 # LTMR should not require frequent refilling

total_refueling_period = params['Fuel Lifetime'] + params['Refueling Period'] + params['Startup Duration after Refueling'] # days
total_refueling_period_yr = total_refueling_period/365
params['A75: Vessel Replacement Period (cycles)']        = np.floor(30/total_refueling_period_yr)*total_refueling_period_yr
params['A75: Core Barrel Replacement Period (cycles)']   = np.floor(15/total_refueling_period_yr)
params['A75: Reflector Replacement Period (cycles)']     = 2
params['A75: Drum Replacement Period (cycles)']          = 2
# **************************************************************************************************************************
#                                           Sec. 12: Buildings & Economic Parameters
# **************************************************************************************************************************

update_params({
    'Land Area': 18,  # acres
    'Escalation Year': 2023,
    
    'Excavation Volume': 412.605,  # m^3
    'Reactor Building Slab Roof Volume': 87.12,  # m^3
    'Reactor Building Basement Volume': 87.12,  # m^3
    'Reactor Building Exterior Walls Volume': 228.8,  # m^3

    'Integrated Heat Exchanger Building Slab Roof Volume': 0,  # m^3
    'Integrated Heat Exchanger Building Basement Volume': 0,  # m^3
    'Integrated Heat Exchanger Building Exterior Walls Volume': 0,  # m^3
    'Integrated Heat Exchanger Building Superstructure Area': 0, # m^2
    
        # Manipulator Building
    'Manipulator Building Slab Roof Volume':0, # m^3
    'Manipulator Building Basement Volume': 0, # m^3
    'Manipulator Building Exterior Walls Volume': 0, # m^3

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
    'Annual Return': 0.0475,  # Annual return on decommissioning costs
    'NOAK Unit Number': 100
})

# **************************************************************************************************************************
#                                           Sec. 13: Post Processing
# **************************************************************************************************************************
params['Number of Samples'] = 10 # Accounting for cost uncertainties
# Estimate costs using the cost database file and save the output to an Excel file
detailed_bottom_up_cost_estimate('cost/Cost_Database.xlsx', params, "output_LTMR.xlsx")
elapsed_time = (time.time() - time_start) / 60  # Calculate execution time
print('Execution time:', np.round(elapsed_time, 1), 'minutes')