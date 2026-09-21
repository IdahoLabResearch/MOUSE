# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

"""
This script performs a bottom-up cost estimate for a Liquid Metal Thermal Microreactor (LTMR).
OpenMC is used for core design calculations, and other Balance of Plant components are estimated.
Users can modify parameters in the "params" dictionary below.
"""
import numpy as np
import watts  # Simulation workflows for one or multiple codes
from core_design.openmc_template_LTMR import *
from core_design.pins_arrangement import get_ltmr_pins_arrangement
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
    # Used for both reactor fuel depletion and concrete activation.
    'simplified_chain_thermal_xml': '/projects/MRP_MOUSE/openmc_data/chain_endf_b8.0.xml'                  # on INL HPC
})

# **************************************************************************************************************************
#                                                Sec. 1: Materials
# **************************************************************************************************************************
update_params({
    'reactor type': "LTMR", # LTMR or GCMR
    'TRISO Fueled': "No",
    'Fuel': 'UO2',
    'Enrichment': 0.05,  # Fraction between 0 and 1
    'Coolant': 'NaK',
    'Radial Reflector': 'Graphite',
    'Axial Reflector': 'Graphite',
    'Moderator': 'ZrH',
    'Control Drum Absorber': 'B4C_enriched',
    'Control Drum Reflector': 'Graphite',
    'Common Temperature': 600,  # Kelvins
    'HX Material': 'SS316'
})
# **************************************************************************************************************************
#                                           Sec. 2: Geometry: Fuel Pins, Moderator Pins, Coolant, Hexagonal Lattice
# **************************************************************************************************************************  

# Select one of the predefined sixfold-symmetric shutdown-rod layouts.
# Supported values are 6 and 12.
params['Number of Shutdown Rods'] = 6

update_params({
    'Fuel Pin Materials': ['Zr', None, params['Fuel'], None, 'SS304'],
    'Fuel Pin Radii': [0.28575, 0.3175, 1.5113, 1.5367, 1.5875],  # cm
    'Moderator Pin Materials': ['ZrH', 'SS304'],  
    'Moderator Pin Inner Radius': 1.5367,  # cm
    'Moderator Pin Radii': [1.5367, 1.5875],  # [params['Moderator Pin Inner Radius'], params['Fuel Pin Radii'][-1]]
    "Pin Gap Distance": 0.1,  # cm
    'Pins Arrangement': get_ltmr_pins_arrangement(
        params['Number of Shutdown Rods']
    ),
    'Number of Rings per Assembly': 14, # the number of rings can be 12 or lower as long as the heat flux criteria is not violated
    'Radial Reflector Thickness': 44.53999867260457,  # cm
    'Axial Reflector Thickness': 44.53999867260457,  # cm
})

params['Lattice Apothem'] = calculate_hex_apothem(params)
params['Lattice Radius'] = params['Lattice Apothem']
params['Assembly FTF'] = 2 * params['Lattice Apothem']
params['Active Height']  = 120
params['Shutdown Rod Height'] = params['Active Height']  # cm
params['Fuel Pin Count'] = calculate_pins_in_assembly(params, "FUEL")
params['Moderator Pin Count'] = calculate_pins_in_assembly(params, "MODERATOR")
params['Moderator Mass'] = calculate_moderator_mass(params)
params['Core Radius'] = 83.11277015716345  # cm

# **************************************************************************************************************************
#                                           Sec. 3: Control Drums
# ************************************************************************************************************************** 

update_params({
    'Number of Drums': 6,
    'Drum Radius': 22.02527406887039,  # cm
    'Drum Tube Radius': 22.26999933630228,  # cm
    'Drum Absorber Thickness': 1,  # cm
    'Drum Absorber Arc Degrees': 120,
    'Drum Height': 209.0799973452091,  # cm
    'Shutdown Rod Absorber': 'B4C_enriched',
    'Shutdown Rod Cladding': 'SS304',

    # Must fit inside the existing pin envelope
    'Shutdown Rod Absorber Radius': 1.30,  # cm
    'Shutdown Rod Clad Radius': 1.50,      # cm

})

# Explicit original geometry. The shared helper may validate and synchronize
# dependent values, but it is not allowed to change any selected dimension.
explicit_original_geometry = {
    'Core Radius': 83.11277015716345,
    'Radial Reflector Thickness': 44.53999867260457,
    'Axial Reflector Thickness': 44.53999867260457,
    'Drum Radius': 22.02527406887039,
    'Drum Tube Radius': 22.26999933630228,
    'Drum Height': 209.0799973452091,
}
update_ltmr_reflector_geometry_from_drums(params)
for parameter_name, expected_value in explicit_original_geometry.items():
    if not np.isclose(params[parameter_name], expected_value, rtol=0.0, atol=1.0e-9):
        raise RuntimeError(
            f"Explicit original geometry mismatch for {parameter_name}: "
            f"expected {expected_value}, got {params[parameter_name]}"
        )
calculate_drums_volumes_and_masses(params)
calculate_shutdown_rods_volumes_and_masses(params)
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

# Declare shielding and vessel inputs before run_openmc so the BOL leakage
# statepoint can produce the shield result before the fuel-depletion sequence.
update_params({
    'Dynamic Shielding Calculation': True,
    'Shielding Dose Limit': 0.5,  # mrem/h
    'Shielding Irradiation Years': 60.0,
    'Shielding Activation Step Days': 90.0,
    'Shielding Decay Days': 90.0,
    'Shielding Detector Distance': 30.0,
    # Volume-averaged dose in a 10 cm-radius sphere centered 30 cm above
    # the pit opening. This modestly improves photon tally statistics.
    'Shielding Detector Radius': 10.0,
    'Shielding Concrete Thickness': 200.0,
    # Three activation zones through the concrete depth on the side wall and
    # base: 0-5 cm, 5-25 cm, and 25-200 cm.
    'Shielding Concrete Activation Zone Boundaries': [5.0, 25.0],
    'Shielding Pit Clearance': 10.0,
    'Shielding Temperature': 300.0,
    'Shielding Particles': 5000,
    'Shielding Batches': 10,
    'Shielding Photon Particles': 50000,
    'Shielding Photon Batches': 10,
    'Shielding Photon Retry Multiplier': 4,
    'Shielding Photon Maximum Retries': 1,
    'Shielding Thickness Tolerance': 10.0,
    'Shielding Confirmation Increment': 5.0,
    # One-sided 95% upper estimate used only when a zero-score photon tally
    # requires conservative exponential tail extrapolation.
    'Shielding Extrapolation Confidence Multiplier': 1.645,
    'Shielding Initial Upper Thickness': 60.0,
    'Shielding Maximum Thickness': 300.0,
    'Shielding Maximum Search Iterations': 1,
    'In Vessel Shield Thickness': 10.16,  # cm
    'In Vessel Shield Inner Radius': params['Core Radius'],
    'In Vessel Shield Material': 'B4C_natural',
    'Out Of Vessel Shield Thickness': 39.37,  # initial/fallback value, cm
    'Out Of Vessel Shield Material': 'WEP',
    'Out Of Vessel Shield Effective Density Factor': 0.5,
    'Vessel Radius': params['Core Radius'] + 10.16,
    'Vessel Thickness': 2,  # cm
    'Vessel Lower Plenum Height': 50,  # cm
    'Vessel Upper Plenum Height': 47.152,  # cm
    'Vessel Upper Gas Gap': 0,
    'Vessel Bottom Depth': 32.129,
    'Vessel Material': 'stainless_steel',
    'Gap Between Vessel And Guard Vessel': 5,  # cm
    'Guard Vessel Thickness': 1,  # cm
    'Guard Vessel Material': 'stainless_steel',
    'Gap Between Guard Vessel And Cooling Vessel': 5,  # cm
    'Cooling Vessel Thickness': 0.5,  # cm
    'Cooling Vessel Material': 'stainless_steel',
    'Gap Between Cooling Vessel And Intake Vessel': 5,  # cm
    'Intake Vessel Thickness': 0.5,  # cm
    'Intake Vessel Material': 'stainless_steel',
})
params['In Vessel Shield Outer Radius'] = (
    params['Core Radius'] + params['In Vessel Shield Thickness']
)

# **************************************************************************************************************************
#                                           Sec. 5: Running OpenMC
# **************************************************************************************************************************

# --- Shutdown Margin  ---
# When True, dedicated cold static calculations evaluate shutdown margin at
# BOL, MOL, and EOL. Shutdown rods are inserted and control-drum absorbers face
# the core at 'Cold Shutdown Temperature'.
# Recommended: True for final design verification; can be set to False to save
# computation time during early design exploration.
params['Shutdown Margin Calc'] = True  # True or False

# --- Isothermal Temperature Coefficient ---
# When True, dedicated static pairs at BOL, MOL, and EOL are run at
# 'Common Temperature' and at 'Common Temperature' + 'Temperature Perturbation'.
# The temperature coefficient is calculated in units of pcm/K and includes the
# modeled NaK and ZrH density changes between the two temperatures.
# A negative coefficient indicates the reactor is self-stabilizing (desired behavior).
# Recommended: True for safety analysis; can be set to False to save computation time.
params['Isothermal Temperature Coefficients'] = True  # True or False

# --- Temperature Perturbation ---
# The temperature step (in Kelvin) used for the isothermal temperature coefficient calculation.
# Must be large enough to produce a keff difference above OpenMC Monte Carlo statistical
# noise, but small enough to stay in the linear reactivity regime.
# Typical range: 50–300 K. 100 K is chosen here as a balance between accuracy and
# avoiding nonlinear effects. 
# Units: Kelvin
# This parameter is REQUIRED only when 'Isothermal Temperature Coefficients' is True.
params['Temperature Perturbation'] = 100  # K

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
# Inputs are declared before run_openmc; the dynamic calculation updates the
# out-of-vessel thickness before depletion begins.

# **************************************************************************************************************************
#                                           Sec. 8: Vessels Calculations
# ************************************************************************************************************************** 
vessels_specs(params)  # calculate the volumes and masses of the vessels
calculate_shielding_masses(params)  # calculate the masses of the shieldings

# **************************************************************************************************************************
#                                           Sec. 9: Operation
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
## Based on https://www.edf.fr/sites/default/files/mediatheque/dp_creys_2017.pdf :
## 5,500 tonnes of sodium from the reactor vessel and secondary circuits at the Creys-Malville plant (France), which is 3,000 MWt.
# This gives a rough estimate of 1833 kg/MWt.
params['Onsite Coolant Inventory'] = 1833 * params['Power MWt']
params['Replacement Coolant Inventory'] = 0  # NaK is assumed not to require replacement
# params['Annual Coolant Supply Frequency']  # LTMR should not require frequent refilling

params['A75: Outer Vessel Structure Replacement Period (years)'] = 20
params['A75: Inner Vessel Structure Replacement Period (years)'] = 10
params['A75: Reflector Replacement Period (years)'] = 10
params['A75: Reactor Control Devices Replacement Period (years)'] = 10
params['A75: Moderator Booster Replacement Period (cycles)'] = 1
params['Maintenance to Direct Cost Ratio']                = 0.015
# A78: Annualized Decommissioning Cost
params['A78: CAPEX to Decommissioning Cost Ratio'] = 0.15

# **************************************************************************************************************************
#                                           Sec. 10: Buildings & Economic Parameters
# **************************************************************************************************************************

update_params({
    'Land Area': 18,  # acres
    'Escalation Year': 2025,

    'Excavation Volume': 412.605,  # m^3
    'Reactor Building Slab Roof Volume': (9750*6502.4*1500)/1e9,  # m^3
    'Reactor Building Basement Volume': (9750*6502.4*1500)/1e9,  # m^3
    'Reactor Building Exterior Walls Volume': ((2*9750*3500*1500)+(3502.4*3500*(1500+750)))/1e9,  # m^3
    'Reactor Building Superstructure Area': ((2*3500*3500)+(2*7500*3500))/1e6, # m^2
    
    # Connected to the Reactor Building (contains steel liner)
    'Integrated Heat Exchanger Building Slab Roof Volume': 0,  # m^3
    'Integrated Heat Exchanger Building Basement Volume': 0,  # m^3
    'Integrated Heat Exchanger Building Exterior Walls Volume': 0,  # m^3
    'Integrated Heat Exchanger Building Superstructure Area': 0, # m^2
    
    # Assumed to be High 40' CONEX Container with 20 cm wall thickness (including conex wall)
    'Turbine Building Slab Roof Volume': (12192*2438*200)/1e9,  # m^3
    'Turbine Building Basement Volume': (12192*2438*200)/1e9,  # m^3
    'Turbine Building Exterior Walls Volume': ((12192*2496*200)+(2038*2496*200))*2/1e9,  # m^3
    
    # Assumed to be High 40' CONEX Container with 20 cm wall thickness (including conex wall)
    'Control Building Slab Roof Volume': (12192*2438*200)/1e9,  # m^3
    'Control Building Basement Volume': (12192*2438*200)/1e9,  # m^3
    'Control Building Exterior Walls Volume': ((12192*2496*200)+(2038*2496*200))*2/1e9,  # m^3
    
    # Manipulator Building
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
    
    # Building to host operational spares (CO2, He, filters, etc.)
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
    'Annual Return': 0.0475,  # Annual return on decommissioning costs
    'NOAK Unit Number': 100
})

# --- ITC (Investment Tax Credit) ---
# The ITC is a one-time credit applied to the Overnight Capital Cost (OCC) of the plant.
# Under the IRA (Inflation Reduction Act), advanced nuclear facilities placed in service
# after Dec 31, 2024 may qualify for the Clean Electricity ITC (Section 48E).
# The ITC level depends on whether the project meets certain requirements:
#   - Base rate (no prevailing wage): 6% of OCC
#   - With prevailing wage + apprenticeship requirements: 30% of OCC
#   - With prevailing wage + domestic content bonus: 40% of OCC
#   - With prevailing wage + domestic content + energy community bonus: 50% of OCC
# Typical values: 0.06, 0.30, 0.40, 0.50
# Note: ITC and PTC are mutually exclusive — only one can be selected per project.
# To disable ITC, remove or comment out this parameter.
params['ITC credit level'] = 0.30  # fraction — assumes prevailing wage requirements are met

# --- IRA Sunset: Number of Units Claiming ITC/PTC ---
# Under the IRA, ITC and PTC eligibility ends at a sunset year. Once the sunset
# is reached, units placed in service after that point cannot claim the credit.
# This parameter caps how many units in the deployment sequence may avail the
# credit. A unit is eligible only if its position is <= this cutoff:
#   - FOAK column = unit 1 (always eligible if cutoff >= 1)
#   - NOAK column = unit 'NOAK Unit Number' (eligible only if NOAK Unit Number <= cutoff)
# When a unit is past the cutoff, the ITC/PTC-adjusted outputs fall back to the
# un-subsidized values, producing a step in the LCOE-vs-deployment-scale curve
# at the sunset point.
# Typical value: a fleet-size estimate consistent with deployments before the
# IRA sunset (e.g. 50, 100). Set very high to keep the original behavior of
# applying the credit to every unit.
params['Number of Units Claiming ITC/PTC'] = 10

# **************************************************************************************************************************
#                                           Sec. 11: Post Processing
# **************************************************************************************************************************
params['Number of Samples'] = 100  # number of samples for cost uncertainty analysis
# Estimate costs using the cost database file and save the output to an Excel file
estimate = detailed_bottom_up_cost_estimate('cost/Cost_Database.xlsx')
elapsed_time = (time.time() - time_start) / 60  # calculate execution time
print('Execution time:', np.round(elapsed_time, 1), 'minutes')
