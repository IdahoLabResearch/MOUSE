# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

"""
This script performs a bottom-up cost estimate for a heat pipe Microreactor.
OpenMC is used for core design calculations, and other Balance of Plant components are estimated.
Users can modify parameters in the "params" dictionary below.
"""
import numpy as np
import watts  # Simulation workflows for one or multiple codes
from core_design.openmc_template_HPMR import *
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
    'cross_sections_xml_location': '/projects/MRP_MOUSE/openmc_data/endfb-viii.0-hdf5/cross_sections.xml', # on INL HPC
    'simplified_chain_thermal_xml': '/projects/MRP_MOUSE/openmc_data/simplified_thermal_chain11.xml'       # on INL HPC
})

# **************************************************************************************************************************
#                                                Sec. 1: Materials
# **************************************************************************************************************************
# These parameters are based on: https://inldigitallibrary.inl.gov/sites/sti/sti/Sort_99962.pdf
update_params({
    'reactor type': "HPMR",
    'TRISO Fueled': "Yes",
    'Fuel': 'homog_TRISO',     
    'Enrichment': 0.19985,
    'Radial Reflector': 'Graphite',
    'Axial Reflector': 'Graphite',
    'Moderator': 'monolith_graphite',
    'Coolant': 'Helium',
    'Control Drum Absorber': 'B4C_natural',
    'Control Drum Reflector': 'Graphite',
    'Cooling Device': 'heatpipe',
    'Common Temperature': 1000,  #K
    'HX Material': 'SS316'
})

# **************************************************************************************************************************
#                                           Sec. 2: Geometry: Fuel Pins, Moderator Pins, Coolant, Hexagonal Lattice
# **************************************************************************************************************************

update_params({
    'Fuel Pin Materials': ['homog_TRISO', 'Helium'],
    'Fuel Pin Radii': [1.00, 1.05], #cm
    'Heat Pipe Materials': ['heatpipe', 'Helium'],
    'Heat Pipe Radii': [1.10, 1.15],
    'Number of Rings per Assembly': 6,
    'Number of Rings per Core': 3,
    'Lattice Pitch': 3.4,
})
params['Assembly FTF'] = (params['Lattice Pitch'] * (params['Number of Rings per Assembly'] - 1) + 1.4 * params['Fuel Pin Radii'][-1]) * np.sqrt(3)
params['hexagonal Core Edge Length'] = (params['Assembly FTF'] * (params['Number of Rings per Core']-1)) + (params['Assembly FTF']/2) + 6.6
params['Radial Reflector Thickness'] = 50 #cm
params['Core Radius'] = 0.5*np.sqrt(3)*params['hexagonal Core Edge Length'] + params['Radial Reflector Thickness']
params['Active Height'] = 2 * params['Core Radius'] #cm  
params['Axial Reflector Thickness'] = params['Radial Reflector Thickness']
params['Fuel Pin Count per Assembly'] = calculate_number_fuel_elements_hpmr(params['Number of Rings per Assembly'])
params['Fuel Assemblies Count'] = (3 * params['Number of Rings per Core']**2) - (3 * params['Number of Rings per Core'])
params['Fuel Pin Count'] = params['Fuel Assemblies Count'] * params['Fuel Pin Count per Assembly']
number_of_heatpipes_hmpr(params)

# **************************************************************************************************************************
#                                           Sec. 3: Control Drums
# ************************************************************************************************************************** 

update_params({
    'Drum Count': 12,   # allowed: 6, 12, 18, 24
    'Drum Radius': 0.4 * params['Radial Reflector Thickness'], 
    'Drum Absorber Thickness': 1,  # cm
    'Drum Height': params['Active Height']
})

calculate_drums_volumes_and_masses(params)
calculate_reflector_and_moderator_mass_HPMR(params)

# **************************************************************************************************************************
#                                           Sec. 4: Overall System
# ************************************************************************************************************************** 
update_params({
    'Power MWt': 7, 
    'Thermal Efficiency': 0.36,
    'Heat Flux Criteria': 0.9,  # MW/m^2 
    'Time Steps': [t * 86400 for t in [0.01, 0.99, 3, 6, 20, 70, 100, 165, 365, 365, 365, 365, 365, 365, 365.00]]  # seconds
})
params['Power MWe'] = params['Power MWt'] * params['Thermal Efficiency']
params['Heat Flux'] = calculate_heat_flux(params)

# **************************************************************************************************************************
#                                           Sec. 5: Running OpenMC
# **************************************************************************************************************************

# --- Shutdown Margin (SDM) ---
# When True, an additional OpenMC simulation is run with all control drums rotated
# to the fully inserted (ARI - All Rods In) position. The SDM is then calculated
# as the difference in reactivity (in pcm) between the ARO and ARI configurations.
# A positive SDM means the reactor can be safely shut down with all drums inserted.
# Recommended: True for final design verification; can be set to False to save
# computation time during early design exploration.
params['Shutdown Margin Calc'] = False  # True or False

# --- Isothermal Temperature Coefficient ---
# When True, two additional OpenMC simulations are run: one at 'Common Temperature'
# and one at 'Common Temperature' + 'Temperature Perturbation'. The temperature
# coefficient is then calculated in units of pcm/K.
# A negative coefficient indicates the reactor is self-stabilizing (desired behavior).
# Recommended: True for safety analysis; can be set to False to save computation time.
# params['Isothermal Temperature Coefficients'] = True  # True or False

# --- Temperature Perturbation ---
# The temperature step (in Kelvin) used for the isothermal temperature coefficient calculation.
# Must be large enough to produce a keff difference above OpenMC Monte Carlo statistical
# noise, but small enough to stay in the linear reactivity regime.
# Typical range: 50–300 K. 100 K is chosen here as a balance between accuracy and
# avoiding nonlinear effects. 
# Units: Kelvin
# This parameter is REQUIRED only when 'Isothermal Temperature Coefficients' is True.
# params['Temperature Perturbation'] = 100  # K

heat_flux_monitor = monitor_heat_flux(params)
run_openmc(build_openmc_model_HPMR, heat_flux_monitor, params)
fuel_calculations(params)  # calculate the fuel mass and SWU

# **************************************************************************************************************************
#                                        Sec. 6: Primary Loop + Balance of Plant
# ************************************************************************************************************************** 
params.update({
    'Primary Loop Purification': True,
    'Secondary HX Mass': 0,
    'Primary Loop Count': 2,
    'Primary Loop Inlet Temperature': 650 + 273.15, # K
    'Primary Loop Outlet Temperature': 650 + 273.15, # K
    'Secondary Loop Inlet Temperature': 300 + 273.15, # K
    'Secondary Loop Outlet Temperature': 630 + 273.15, # K,
   })
params['Primary HX Mass'] = calculate_heat_exchanger_mass(params)  # Kg

params.update({
    'BoP Count': 2,
    'BoP per loop load fraction': 0.5,
    })
params['BoP Power kWe'] = 1000 * params['Power MWe'] * params['BoP per loop load fraction']

# **************************************************************************************************************************
#                                           Sec. 7 : Shielding
# ************************************************************************************************************************** 
update_params({
    'In Vessel Shield Thickness': 0,  # cm (no shield in vessel for HPMR)
    'In Vessel Shield Inner Radius': params['Core Radius'],
    'In Vessel Shield Material': 'B4C_natural',
    'Out Of Vessel Shield Thickness': 39.37,  # cm
    'Out Of Vessel Shield Material': 'WEP',
    'Out Of Vessel Shield Effective Density Factor': 0.5
})
params['In Vessel Shield Outer Radius'] = params['Core Radius'] + params['In Vessel Shield Thickness']

# **************************************************************************************************************************
#                                           Sec. 8 : Vessels Calculations
# ************************************************************************************************************************** 
update_params({
    'Vessel Radius': params['Core Radius'] + params['In Vessel Shield Thickness'],
    'Vessel Thickness': 2,  # cm — ASME Sec III Div 5 high-T creep minimum at >650°C; vessel at atmospheric (heat pipes individually sealed); refs: INL HPMR Design A, eVinci, Oklo (was 1, below Div 5 minimum)
    'Vessel Lower Plenum Height': 20,  # cm — heat-pipe header / lower core support; refs: INL HPMR Design A, MARVEL (was 2.848, unit-conv bug)
    'Vessel Upper Plenum Height': 47.152,       # cm — heat-pipe condenser interface to secondary HX
    'Vessel Upper Gas Gap': 0,
    'Vessel Bottom Depth': 32.129,
    'Vessel Material': 'stainless_steel',
    # Guard vessel intentionally removed: each heat pipe is individually sealed (Na in SS316), no bulk primary coolant inventory requiring secondary containment
    'Gap Between Vessel And Guard Vessel': 0,
    'Guard Vessel Thickness': 0,  # cm
    'Guard Vessel Material': 'low_alloy_steel',
    'Gap Between Guard Vessel And Cooling Vessel': 5,  # cm
    'Cooling Vessel Thickness': 0.5,  # cm
    'Cooling Vessel Material': 'stainless_steel',
    'Gap Between Cooling Vessel And Intake Vessel': 5,  # cm — Hejzlar & Buongiorno 2007 NED RVACS minimum (was 4)
    'Intake Vessel Thickness': 0.5,  # cm
    'Intake Vessel Material': 'stainless_steel'
})

vessels_specs(params)
calculate_shielding_masses(params)

# **************************************************************************************************************************
#                                           Sec. 9 : Operation
# **************************************************************************************************************************
update_params({
    'Operation Mode': "Remotely Monitored",
    # Number of concurrent on-site operators per shift. On-Site Staffed means these positions are staffed 24/7;
    # Remotely Monitored means they are present only for refueling and startup after refueling or emergency shutdowns.
    'Number of On-Site Operators per Shift': 2,
    'Levelization Period': 60,  # years
    'Refueling Period': 7,
    'Emergency Shutdowns Per Year': 0.2,
    'Startup Duration after Refueling': 2,
    'Startup Duration after Emergency Shutdown': 14,
    'Reactors Monitored Per Operator': 10,
    'Security Staff Per Shift': 1
})

params['Onsite Coolant Inventory'] = 0  # the helium gap is extremely thin and can be neglected
params['Replacement Coolant Inventory'] = 0
# params['Annual Coolant Supply Frequency'] = 1 if params['Primary Loop Purification'] else 6

total_refueling_period = params['Fuel Lifetime'] + params['Refueling Period'] + params['Startup Duration after Refueling'] # days
total_refueling_period_yr = total_refueling_period/365
params['A75: Outer Vessel Structure Replacement Period (years)'] = 20
params['A75: Inner Vessel Structure Replacement Period (cycles)'] = np.floor(10/total_refueling_period_yr)
params['A75: Reflector Replacement Period (cycles)']     = np.floor(10/total_refueling_period_yr)
params['A75: Reactor Control Devices Replacement Period (cycles)'] = np.floor(10/total_refueling_period_yr)
params['A75: Moderator Booster Replacement Period (cycles)'] = 1
params['A75: Moderator Replacement Period (cycles)'] = 1
params['A75: Moderator Replacement Period (cycles)'] = 1
params['Maintenance to Direct Cost Ratio']                = 0.015
params['A78: CAPEX to Decommissioning Cost Ratio'] = 0.15

# **************************************************************************************************************************
#                                           Sec. 10 : Economic Parameters
# **************************************************************************************************************************
update_params({
    'Land Area': 18,  # acres
    'Escalation Year': 2025,
    'Excavation Volume': 412.605,  # m^3
    'Reactor Building Slab Roof Volume': (9750*6502.4*1500)/1e9,  # m^3
    'Reactor Building Basement Volume': (9750*6502.4*1500)/1e9,  # m^3
    'Reactor Building Exterior Walls Volume': ((2*9750*3500*1500)+(3502.4*3500*(1500+750)))/1e9,  # m^3
    'Reactor Building Superstructure Area': ((2*3500*3500)+(2*7500*3500))/1e6, # m^2
    'Integrated Heat Exchanger Building Slab Roof Volume': 0,  # m^3
    'Integrated Heat Exchanger Building Basement Volume': 0,  # m^3
    'Integrated Heat Exchanger Building Exterior Walls Volume': 0,  # m^3
    'Integrated Heat Exchanger Building Superstructure Area': 0, # m^2
    'Turbine Building Slab Roof Volume': (12192*2438*200)/1e9,  # m^3
    'Turbine Building Basement Volume': (12192*2438*200)/1e9,  # m^3
    'Turbine Building Exterior Walls Volume': ((12192*2496*200)+(2038*2496*200))*2/1e9,  # m^3
    'Control Building Slab Roof Volume': (12192*2438*200)/1e9,  # m^3
    'Control Building Basement Volume': (12192*2438*200)/1e9,  # m^3
    'Control Building Exterior Walls Volume': ((12192*2496*200)+(2038*2496*200))*2/1e9,  # m^3
    'Manipulator Building Slab Roof Volume': (4876.8*2438.4*400)/1e9, # m^3
    'Manipulator Building Basement Volume': (4876.8*2438.4*1500)/1e9, # m^3
    'Manipulator Building Exterior Walls Volume': ((4876.8*4445*400)+(2038.4*4445*400*2))/1e9, # m^3
    'Refueling Building Slab Roof Volume': 0,  # m^3
    'Refueling Building Basement Volume': 0,  # m^3
    'Refueling Building Exterior Walls Volume': 0,  # m^3
    'Spent Fuel Building Slab Roof Volume': 0,  # m^3
    'Spent Fuel Building Basement Volume': 0,  # m^3
    'Spent Fuel Building Exterior Walls Volume': 0,  # m^3
    'Emergency Building Slab Roof Volume': 0,  # m^3
    'Emergency Building Basement Volume': 0,  # m^3
    'Emergency Building Exterior Walls Volume': 0,  # m^3
    'Storage Building Slab Roof Volume': (8400*3500*400)/1e9, # m^3
    'Storage Building Basement Volume': (8400*3500*400)/1e9, # m^3
    'Storage Building Exterior Walls Volume': ((8400*2700*400)+(3100*2700*400*2))/1e9, # m^3
    'Radwaste Building Slab Roof Volume': 0,  # m^3
    'Radwaste Building Basement Volume': 0,  # m^3
    'Radwaste Building Exterior Walls Volume': 0,  # m^3,
    'Interest Rate': 0.07,
    'Discount Rate': 0.07,
    'Construction Duration': 12,  # months
    'Debt To Equity Ratio': 1,
    'Annual Return': 0.0475,
    'NOAK Unit Number': 100,
})

# --- No Tax Credits Applied ---
# This example does not apply any ITC or PTC tax credits.
# To apply ITC, add: params['ITC credit level'] = 0.30  (see watts_exec_LTMR.py for full details)
# To apply PTC, add: params['PTC credit value'] = 15.0  (see watts_exec_GCMR.py for full details)
# When ITC/PTC is enabled, optionally cap how many units may claim the credit
# under the IRA sunset:
#     params['Number of Units Claiming ITC/PTC'] = 10
# (FOAK = unit 1; NOAK column = unit 'NOAK Unit Number'. Units past the cutoff
# fall back to the un-subsidized values.)
# Note: ITC and PTC are mutually exclusive — only one can be selected per project.

# **************************************************************************************************************************
#                                           Sec. 11: Post Processing
# **************************************************************************************************************************
params['Number of Samples'] = 1  # number of samples for cost uncertainty analysis
# Estimate costs using the cost database file and save the output to an Excel file
estimate = detailed_bottom_up_cost_estimate('cost/Cost_Database.xlsx')
elapsed_time = (time.time() - time_start) / 60  # calculate execution time
print('Execution time:', np.round(elapsed_time, 1), 'minutes')
