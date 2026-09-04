# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

"""
This script performs a bottom-up cost estimate for a Gas Cooled Microreactor (GCMR).
Parallel screening case: E5_enr_15p0.
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
from cost.cost_estimation import detailed_bottom_up_cost_estimate
from cost.fleet_mode import (
    servicing_facility_allocation,
    servicing_facility_occ_learning_multipliers,
)

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
update_params({
    'reactor type': "GCMR",  # LTMR or GCMR
    'TRISO Fueled': "Yes",
    'Fuel': 'UCO',
    'Enrichment': 0.15,  # The enrichment is a fraction. It has to be between 0 and 1
    'UO2 atom fraction': 0.7,  # Mixing UO2 and UC by atom fraction
    'Radial Reflector': 'Graphite',
    'Axial Reflector': 'Graphite',
    'Matrix Material': 'Graphite',  # matrix material is the background material within the compact fuel element between TRISO particles
    'Moderator': 'Graphite',  # the moderator is outside the compact fuel region
    'Moderator Booster Materials': ['ZrH'],
    'Coolant': 'Helium',
    'Common Temperature': 850,  # Kelvins
    # IG-110 proxy: mean of axial/transverse CTE values in
    # ORNL/TM-2017/705, Table 2.2 (4.5 and 4.2 microstrain/K).
    'Graphite Linear Expansion Coefficient': 4.3e-6,  # 1/K
    'Control Drum Absorber': 'B4C_enriched',  # The absorber material in the control drums
    'Control Drum Reflector': 'Graphite',  # The reflector material in the control drums
    'Shutdown Rod Absorber': 'B4C_enriched',
    'Shutdown Rod Cladding': 'SS304',
    'HX Material': 'SS316', 
})

# **************************************************************************************************************************
#                                           Sec. 2: Geometry: Fuel Pins, Moderator Pins, Coolant, Hexagonal Lattice
# **************************************************************************************************************************  

update_params({
    # fuel pin details
    'Fuel Pin Materials': ['UCO', 'buffer_graphite', 'PyC', 'SiC', 'PyC'],
    'Fuel Pin Radii': [0.0250, 0.0350, 0.0390, 0.0425, 0.0465],  # cm # https://art.inl.gov/NRC%20Training%202019/04_TRISO_Fuel.pdf
    'Compact Fuel Radius': 0.6225,  # cm # The radius of the area that is occupied by the TRISO particles (fuel compact/ fuel element)
    'Packing Fraction': 0.4,
    'TRISO Packing Seed': 1,
    
    # Coolant channel and booster dimensions
    'Coolant Channel Radius': 0.35,  # cm
    'Moderator Booster Radii': [0.5],  # cm
    'Lattice Pitch': 2.25,
    'Assembly Rings': 6,
    'Core Rings': 5,

    # Central assembly
    'Central Shutdown Rod Radius': 0.85,  # cm
    'Central Shutdown Rod Clad Radius': 1.05,  # cm; 0.20 cm SS304
    'Central Shutdown Rod Ring': 2,
    'Central Shutdown Rod Count': 12,

    # Six assemblies surrounding the center
    'Surrounding Shutdown Rod Radius': 0.45,  # cm
    'Surrounding Shutdown Rod Clad Radius': 0.65,  # cm; 0.20 cm SS304
    'Surrounding Shutdown Rod Ring': 2,
    'Surrounding Shutdown Rod Count': 2,
    'Surrounding Shutdown Assembly Count': 6,

    # Explicit geometry values for this design. The geometry helper validates
    # these values and does not replace them with calculated dimensions.
    'Assembly FTF': 19.48557158514987,  # cm
    'Active Height': 200.0,  # cm
    'Radial Reflector Thickness': 9.742785792574935,  # cm
    'Axial Reflector Thickness': 9.742785792574935,  # cm
    'Core Radius': 107.17064371832429,  # cm
    'Shutdown Rod Height': 200.0,  # cm
})

# **************************************************************************************************************************
#                                           Sec. 3: Control Drums
# ************************************************************************************************************************** 
update_params({
    'Drum Count': 24,
    'Drum Radius': 9.530986101432001,  # cm
    'Drum Tube Radius': 9.742785792574935,  # cm
    'Drum Absorber Thickness': 1, # cm
    'Drum Absorber Arc Degrees': 120.0,
    'Drum Height': 219.48557158514987,  # cm
    })
calculate_drums_volumes_and_masses(params)
calculate_gcmr_shutdown_rods_volumes_and_masses(params)
calculate_reflector_mass_GCMR(params)          
calculate_moderator_mass_GCMR(params) 

# **************************************************************************************************************************
#                                           Sec. 4: Overall System
# ************************************************************************************************************************** 
update_params({
    'Power MWt': 15,  # MWt
    'Thermal Efficiency': 0.4,
    'Heat Flux Criteria': 0.9,  # MW/m^2 (needs review)
    'Burnup Steps': [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0,
                     30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0]  # MWd_per_Kg
    })

params['Power MWe'] = params['Power MWt'] * params['Thermal Efficiency'] 
params['Heat Flux'] = calculate_heat_flux_TRISO(params) # MW/m^2

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
# params['Shutdown Margin Calc'] = True  # True or False
# params['Cold Shutdown Temperature'] = 300  # K

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

# heat_flux_monitor = monitor_heat_flux(params)
# run_openmc(build_openmc_model_GCMR, heat_flux_monitor, params)

# --- Previously calculated OpenMC results ---
# To bypass OpenMC later, comment out run_openmc(...) above and uncomment these assignments.
params['Fuel Lifetime'] = 1786  # days
params['Mass U235'] = 85655.7587486539  # g
params['Mass U238'] = 484226.2801659319  # g
params['Uranium Mass'] = 569.8820389145857  # kg
fuel_calculations(params)  # calculate the fuel mass and SWU

# **************************************************************************************************************************
#                                         Sec. 6: Primary Loop + Balance of Plant
# ************************************************************************************************************************** 
params.update({
    'Primary Loop Purification': False,
    'Secondary HX Mass': 0,
    'Compressor Pressure Ratio': 4,
    'Compressor Isentropic Efficiency': 0.8,
    'Primary Loop Count': 2,  # number of primary coolant loops in the plant
    'Primary Loop per loop load fraction': 0.5,  # each loop handles an equal share of the total load
    'Primary Loop Inlet Temperature': 300 + 273.15, # K
    'Primary Loop Outlet Temperature': 550 + 273.15, # K
    'Secondary Loop Inlet Temperature': 270 + 273.15, # K — cold-end PCHE pinch 30°C with 300°C primary inlet (was 290 -> 10°C pinch, below realistic PCHE design)
    'Secondary Loop Outlet Temperature': 500 + 273.15, # K,
    'Primary Loop Pressure Drop': 50e3,  # Pa — estimated assumption
})
params['Primary HX Mass'] = calculate_heat_exchanger_mass(params)  # Kg
# calculate coolant mass flow rate
mass_flow_rate(params)
compressor_power(params)

# Update BoP Parameters
params.update({
    'BoP Count': 2,  # number of BoP systems in the plant
    'BoP per loop load fraction': 0.5,  # each BoP handles an equal share of the total load
    })
params['BoP Power kWe'] = 1000 * params['Power MWe'] * params['BoP per loop load fraction']

# Integrated Heat Transfer Vessel
params.update({
    'Integrated Heat Transfer Vessel Thickness': 0, # cm
    'Integrated Heat Transfer Vessel Material': 'SA508',
})
GCMR_integrated_heat_transfer_vessel(params)

# **************************************************************************************************************************
#                                           Sec. 7 : Shielding
# ************************************************************************************************************************** 
update_params({
    'In Vessel Shield Thickness': 0,  # cm (no shield in vessel for GCMR)
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
    'Vessel Thickness': 3,  # cm — ASME Sec III Div 1 thin-shell with 4 MPa He, R=60-100 cm, S=138 MPa SA-508 at 350°C, +3 mm corrosion (was 1, below ASME pressure-driven minimum)
    'Vessel Lower Plenum Height': 30,  # cm — GA MHTGR / HTR-PM-class flow distributor (was 2.848, unit-conv bug)
    'Vessel Upper Plenum Height': 47.152,       # cm — outlet plenum for hot-leg gas exit
    'Vessel Upper Gas Gap': 0,
    'Vessel Bottom Depth': 32.129,
    'Vessel Material': 'stainless_steel',
    # Guard vessel intentionally removed: He is inert, no chemical-leak hazard requiring secondary containment
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
    'Number of Operators': 2,
    'Levelization Period': 40,  # years
    'Refueling Period': 30,
    'Emergency Shutdowns Per Year': 0.2,
    'Startup Duration after Refueling': 2,
    'Startup Duration after Emergency Shutdown': 14,
    'Reactors Monitored Per Operator': 9,
    'Security Staff Per Shift': 1
})

# Based on https://digital.library.unt.edu/ark:/67531/metadc893980/m2/1/high_res_d/919556.pdf (tables 17 and 18):
# Estimated helium mass per MWt is 3.3 kg/MWt.
params['Onsite Coolant Inventory'] = 3.3 * params['Power MWt']  # kg
# According to https://www.nationalacademies.org/read/12844/chapter/6#69, the helium loss rate is 10% per year,
# so 1/10 of the initial inventory is replenished annually.
# Without purification, helium needs to be replaced more frequently.
params['Replacement Coolant Inventory'] = params['Onsite Coolant Inventory'] / 10
params['Annual Coolant Supply Frequency'] = 1 if params['Primary Loop Purification'] else 4

total_refueling_period = params['Fuel Lifetime'] + params['Refueling Period'] + params['Startup Duration after Refueling'] # days
total_refueling_period_yr = total_refueling_period/365
params['A75: Vessel Replacement Period (cycles)']        = np.floor(10/total_refueling_period_yr)
params['A75: Core Barrel Replacement Period (cycles)']   = np.floor(10/total_refueling_period_yr)
params['A75: Reflector Replacement Period (cycles)']     = np.floor(10/total_refueling_period_yr)
params['A75: Drum Replacement Period (cycles)']          = np.floor(10/total_refueling_period_yr)
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

# --- PTC (Production Tax Credit) ---
# The PTC is a per-MWh credit earned for every MWh of electricity produced and sold
# during the credit period. Under the IRA (Section 45Y), advanced nuclear facilities
# placed in service after Dec 31, 2024 may qualify for the Clean Electricity PTC.
# Note: ITC and PTC are mutually exclusive — only one can be selected per project.

# Base credit rate ($/MWh):
#   - $3/MWh  if prevailing wage requirements are NOT met
#   - $15/MWh if prevailing wage + apprenticeship requirements ARE met (5x multiplier)
# Assumed here: $15/MWh (prevailing wage requirements met)
# Units: $/MWh
# params['PTC credit value'] = 15.0  # $/MWh

# Duration of the PTC credit period.
# Under the IRA Section 45Y, the credit is available for 10 years after the facility
# is placed in service.
# Units: years
# Typical value: 10 years
# params['PTC credit period'] = 10  # years

# --- PTC Bonus Multipliers (optional, stackable) ---
# Under the IRA, additional bonus credits can be stacked on top of the base PTC
# if the project meets certain criteria. Each bonus is expressed as a fraction
# added to the base multiplier of 1.0.
# - domestic_content_bonus: +10% if the facility uses US-made iron, steel, and
#   manufactured products (Section 45Y domestic content adder)
#   Typical value: 0.10 (10%)
# - energy_community_bonus: +10% if the facility is sited in an "energy community"
#   (areas affected by coal plant closures or fossil fuel employment decline)
#   Typical value: 0.10 (10%)
# To disable bonuses, set both to 0.0 or remove them entirely.
# params['domestic_content_bonus'] = 0.10   # fraction — assumes domestic content standard is met
# params['energy_community_bonus'] = 0.10   # fraction — assumes facility is in an energy community

# --- Corporate Tax Rate ---
# The US federal corporate tax rate used to gross up the PTC tax credit to its
# before-tax revenue equivalent in the LCOE calculation. Since MOUSE uses a
# before-tax LCOE, the PTC must be converted to a before-tax equivalent.
# The current US federal corporate tax rate is 21% (as of 2024).
# Municipal utilities and non-profit cooperatives may use 0.0 (tax-exempt).
# Units: fraction (e.g. 0.21 for 21%)
# Typical values: 0.21 (federal only), 0.27 (federal + average state)
# params['Tax Rate'] = 0.21  # fraction

# --- IRA Sunset: Number of Units Claiming ITC/PTC ---
# Caps how many units in the deployment sequence may avail the credit. A unit
# is eligible only if its position is <= this cutoff. FOAK = unit 1; the NOAK
# column = unit 'NOAK Unit Number'. When a unit is past the cutoff, the
# ITC/PTC-adjusted outputs fall back to the un-subsidized values, producing
# a step in the LCOE-vs-deployment-scale curve at the sunset point.
# Only applies when ITC or PTC is enabled above. Defaults to effectively
# infinite when omitted (every unit claims the credit).
# params['Number of Units Claiming ITC/PTC'] = 10

# **************************************************************************************************************************
#                                           Sec. 11: Fleet Mode
# **************************************************************************************************************************
# Fleet Mode estimates shared manufacturing and servicing campuses for a fleet of
# GCMRs. Campus sizes and resource requirements are driven primarily by the annual
# reactor production rate.
params['Fleet Mode'] = True

## User Inputs

params['Production Rate'] = 100 #reactors per year
params['Deployment Period'] = 10 # years
params['Fleet'] = params['Deployment Period'] * params['Production Rate']
params['Average Distance From Serv to GenSite'] = 1000 #miles (statutory miles)
params['Servicing Total Time'] = 10 / 365  # years; 5.98 servicing days plus 78 testing hours, rounded up to 10 whole days.
params['Extra RPV Fraction'] = 0.015
params['Max Reactors Per Servicing Facility'] = 1000
params['Servicing Facility OCC Learning Rate'] = 0.30
params['Servicing Facility Learning Cap'] = 5
params['Cask Learning Rate'] = 0.15
params['Cask Learning Cap'] = 100

# ***
## END of User Inputs
# ***


## CONSTANTS

params['m3_to_kg_He_RT_atmospheric'] = 0.1661 #kg/m^3
params['Shift To Headcount'] = 5


## GENERATING SITE

params['CoolantInventoryRPV_Mass'] = 40.5125 #kg
params['Cycle Length'] = (params['Fuel Lifetime'] + params['Refueling Period']) / 365  # years
params['Fuel Mass In Core'] = params['Uranium Mass'] #[kgU] This is a MOUSE output.
params['Water Supply Frequency'] = 4 #[year^-1]
params['Maintenance Visit Frequency'] = 1 / (params['Cycle Length'] / 2)
params['GenSite Downtime'] = 1/12 # years (one month)

#To Discuss: Not Used in Cost Database at this time
# params['GenSite Operators Per Shift'] = 1
# params['GenSite Security Staff Per Shift'] = 1
# params['GenSite Shifts Per Day'] = 2
# params['GenSite Staff Rotation Working Fraction'] = 0.45


## FLEET



## MANUFACTURING CAMPUS

params['MFG Construction Duration'] = 36  # months
params['Manufacturing Campus Area'] = 56.1651 + 1.2066 * params['Production Rate'] ** 0.3919
params['MFG Campus Land Area'] = params['Manufacturing Campus Area']
params['MFG Emergency Generator Power'] = -10851.2764 + 10851.2773 * params['Production Rate'] ** 0.0001
params['MFG Warehouse Building Area'] = 5598.3987 + 70.3374 * params['Production Rate'] ** 0.5539
params['MFG Administration Building Area'] = 3947.3106 + 77.6167 * params['Production Rate'] ** 0.4768
params['MFG Warehouse Staff'] = np.ceil(0.8576 + 0.2539 * params['Production Rate'] ** 0.6533)
params['MFG Security Staff Per Shift'] = np.ceil(2.9981 + 0.3340 * params['Production Rate'] ** 0.4768)
params['MFG Maintenance Staff Headcount'] = np.ceil(2.9981 + 0.3340 * params['Production Rate'] ** 0.4768)
params['MFG Local Transport Vehicle Count'] = np.rint(-4340.5105 + 4340.5109 * params['Production Rate'] ** 0.0001)
params['MFG Utility Vehicle Count'] = np.rint(-4339.5105 + 4340.5109 * params['Production Rate'] ** 0.0001)
params['MFG Fire Station Count'] = 1  # REVIEW NEEDED: provisional assumption; absent from the source parameter file.
if params['Production Rate'] <= 500:
    params['MFG Guard Station Count'] = 1
elif params['Production Rate'] < 2000:
    params['MFG Guard Station Count'] = 2
else:
    params['MFG Guard Station Count'] = np.rint(params['Production Rate'] / 1000)

params['MFG Campus Power'] = 1.295614 + 0.00813495 * params['Production Rate'] ** 0.620945  # MWe
params['MFG Switchyard Rating'] = np.ceil(params['MFG Campus Power'] * 2)
params['MFG Testing Line Annual Rate'] = int(np.floor(8000 / 92))
params['MFG Testing Line Count'] = np.ceil(params['Production Rate'] / params['MFG Testing Line Annual Rate'])
params['MFG Testing Engineering Headcount'] = 6 + (params['MFG Testing Line Count'] - 1) * 2
params['MFG Testing Number of Operators'] = 5 * params['MFG Testing Line Count']  # REVIEW NEEDED: five total operators per testing line.
params['MFG Testing Coolant Allotment'] = 0.15 * (1 * 24.417 + 9 * 11.114) * 8.2402
params['MFG Local Control Building Area'] = 40 * params['MFG Testing Line Count']  # m^2; REVIEW NEEDED.
params['MFG Security Building Area'] = 8775 / (3.2808 ** 2)  # m^2; REVIEW NEEDED: uses servicing security-building area.

# ToDo: Road Length and Site Perimeter should not use the same relationship.
params['MFG Road Length'] = -1539388.0000975644 + 1540607.3153829477 * params['Production Rate'] ** 8.588964624690402e-05
params['MFG Site Perimeter'] = -1539388.0000975644 + 1540607.3153829477 * params['Production Rate'] ** 8.588964624690402e-05
params['MFG Protected Perimeter'] = 1421.7487567019996 + 60.41360643368738 * params['Production Rate'] ** 0.22340396387855505

params['MFG Security Camera Count'] = 200 + -157.0909090909111 * (params['Production Rate'] ** -0.037788560889399164)
params['MFG Motion Detector Count'] = 93.59999999999997 + 2.8444444444444628 * (params['Production Rate'] ** 0.35218251811136164)


## SERVICING CAMPUS

params['SER Construction Duration'] = 36  # months
params['Servicing Rate'] = params['Fleet'] / params['Cycle Length']  # reactors/year
(
    params['Servicing Facility Count'],
    params['Servicing Facility Reactor Counts'],
    params['Servicing Facility Design Capacity'],
    params['Servicing Rate Per Facility'],
) = servicing_facility_allocation(
    params['Fleet'],
    params['Servicing Rate'],
    params['Max Reactors Per Servicing Facility'],
)
params['Servicing Facility OCC Learning Multipliers'] = (
    servicing_facility_occ_learning_multipliers(
        params['Servicing Facility Count'],
        params['Servicing Facility OCC Learning Rate'],
        params['Servicing Facility Learning Cap'],
    )
)

params['SER Campus Area'] = (760 - 160.5) * (params['Servicing Rate Per Facility'] - 30) / (300 - 30) + 160.5
params['SER Campus Land Area'] = params['SER Campus Area']

scale_var_SER = np.rint(params['Servicing Rate Per Facility'] / 3) #Value of 3 should remain hardcoded. Couldn't let production rate be the variable for the Servicing Campus directly, for flexibility, but the relationships were built based on production rate with our basic assumption of a 3:1 ratio between servicing rate and production rate. Simple fix was to define this variable (scale_var_SER).

params['SER Switchyard Rating'] = 0 + 12 * scale_var_SER ** 0.698970
params['SER Switchyard Average Power'] = params['SER Switchyard Rating']/2

params['Servicing Hot Cell Annual Rate'] = int(np.floor(365* (11/12) / 3 ))
params['Servicing Hot Cell Count'] = np.ceil( params['Servicing Rate Per Facility'] / params['Servicing Hot Cell Annual Rate'] )
params['Radioactive Waste Processing Hot Cell Count'] = 1
params['He Gas Replenishment Per Hot Cell'] = ((3*3*5) * 2 * params['Servicing Hot Cell Annual Rate'] * 2 + 0.1 * (10*30*7)*12) * params['m3_to_kg_He_RT_atmospheric']
params['He Gas Replenishment'] = (params['Servicing Hot Cell Count'] * params['He Gas Replenishment Per Hot Cell'] + params['Radioactive Waste Processing Hot Cell Count'] * params['He Gas Replenishment Per Hot Cell'] + params['CoolantInventoryRPV_Mass'] * params['Servicing Rate Per Facility'])

params['SER Number of Operators Per Shift'] = np.ceil(0.0 + 33.33333333333329 * (scale_var_SER ** 0.4771212547196621))
params['SER Engineering Headcount'] = np.ceil(0.0 + 124.99999999999939 * (scale_var_SER ** 0.3010299956639814))
params['SER Maintenance Staff Per Shift'] = np.ceil(0.0 + 266.6666666666656 * (scale_var_SER ** 0.47712125471966266))
params['SER Security Staff Per Shift'] = np.ceil(0.0 + 3.1304347826086882 * (scale_var_SER ** 0.5835765856339492))

params['Roundtrip Time'] = 2*params['Average Distance From Serv to GenSite'] / 40 / (8760/2) #[years] based on average speed of 40 miles/hour, 12 hours of driving per day
params['Roundtrip Time Reactor Transport'] = 2*params['Average Distance From Serv to GenSite'] / 15 / (8760/2) #[years] based on average speed of 15 miles/hour, 12 hours of driving per day
params['Generating Sites Count'] = np.floor(
    params['Fleet']
    / (
        (
            (params['Roundtrip Time Reactor Transport'] + params['Servicing Total Time'] + params['Cycle Length'])
            / params['Cycle Length']
        )
        * (1 + params['Extra RPV Fraction'])
    )
)  # TimeOnSite = Cycle Length; TimeOffSite = reactor roundtrip transport time + total servicing time.
params['Dwell Time GenSite'] = 1/365 #[years]
params['Dwell Time Serv'] = 1/365 #[years]
params['Dwell Time Reactor Transport GenSite'] = 2/365 #[years]
params['Dwell Time Reactor Transport Serv'] = 2/365 #[years]

#ToDo: Road Length and Site Perimeter shouldn't be the same.
params['SER Road Length'] = -1539388.0000975644 + 1540607.3153829477 * (scale_var_SER ** 8.588964624690402e-05)
params['SER Site Perimeter'] = -1539388.0000975644 + 1540607.3153829477 * (scale_var_SER ** 8.588964624690402e-05)
params['SER Controlled Perimeter'] = 0
params['SER Protected Perimeter'] = 1421.7487567019996 + 60.41360643368738 * (scale_var_SER ** 0.22340396387855505)

params['Used Fuel Storage Lifetime Capacity'] = params['Servicing Facility Design Capacity'] * 50 / (params['Cycle Length'] + params['GenSite Downtime']) * params['Fuel Mass In Core']

params['SER Security Building Area'] = 8775 / (3.2808 ** 2)
params['SER Administration Building Area'] = 258000 / (3.2808 ** 2)

params['SER Helium Flowrate'] = params['He Gas Replenishment'] / (0.9 * 8766 * 3600)

params['SER Local Transport Vehicle Count'] = np.ceil(scale_var_SER * 5 / (0.9 * 365.25))  # Five movements per reactor serviced, adjusted for 90% vehicle capacity factor.
params['SER Utility Vehicle Count'] = np.ceil(0.0 + 3.1304347826086882 * (scale_var_SER ** 0.5835765856339492))
params['Reactor Transport Vehicle Count'] = np.ceil(params['Servicing Rate Per Facility'] * (params['Roundtrip Time Reactor Transport']+params['Dwell Time Reactor Transport GenSite']+params['Dwell Time Reactor Transport Serv']) + 1)
params['Helium Transport Truck Count'] = np.ceil(params['Servicing Facility Design Capacity'] * params['Annual Coolant Supply Frequency'] * ((params['Roundtrip Time']+params['Dwell Time GenSite']+params['Dwell Time Serv']) * 1.05))
params['Water Tanker Truck Count'] = np.ceil(params['Servicing Facility Design Capacity'] * params['Water Supply Frequency'] * ((params['Roundtrip Time']+params['Dwell Time GenSite']+params['Dwell Time Serv']) * 1.05))
params['Maintenance Truck Count'] = np.ceil(params['Servicing Facility Design Capacity'] * params['Maintenance Visit Frequency'] * ((params['Roundtrip Time']+params['Dwell Time GenSite']+params['Dwell Time Serv']) * 1.05))

params['SER Radiation Monitor Count'] = np.ceil(params['SER Site Perimeter'] / 1000 * 1.2)
params['SER Security Camera Count'] = np.ceil(params['SER Campus Area'] / params['Manufacturing Campus Area'] * params['MFG Security Camera Count'])
params['SER Motion Detector Count'] = np.ceil(params['SER Campus Area'] / params['Manufacturing Campus Area'] * params['MFG Motion Detector Count'])

params['Reactor Transport Cask Count'] = round(params['Servicing Facility Design Capacity'] * ((4 / 12) / params['Cycle Length']) * 1.5, -1)
params['Annual Used Fuel Cask Consumption'] = 1 * params['Servicing Rate Per Facility']
params['Annual Reactor Cask Replacement'] = np.ceil(0.05 * params['Reactor Transport Cask Count'])
params['Annual Radwaste Cask Consumption'] = 0.5 * params['Servicing Rate Per Facility']

params['Servicing Hot Cell Building Area'] = 0.0 + 411.428571 * (scale_var_SER ** 0.942)
params['Helium Purification and Storage Building Area'] = 0.0 + 72.0 * (scale_var_SER ** 0.6198)
params['SER Local Control Building Area'] = 40 * scale_var_SER
params['SER Remote Control Building Area'] = 40 * scale_var_SER
params['SER Rad Waste Management Building Area'] = 0.0 + 40.5 * (scale_var_SER ** 0.7447)
params['Radwaste Storage Warehouses Area'] = (0.0 + 15422.94 * (scale_var_SER ** 0.994))  # m^2
params['SER Emergency Generator Power'] = 0.0 + 3.266667 * (scale_var_SER ** 0.632)
params['Parts Service Center and Warehouse Building Area'] = 0.0 + 675.0 * (scale_var_SER ** 0.6021)
params['Service Air Water Building Count'] = np.rint( 0.0 + 1.333333 * (scale_var_SER ** 0.1761) )
params['SER Fire Station Count'] = np.rint( 0.0 + 1.333333 * (scale_var_SER ** 0.1761) )
params['SER Guard Station Count'] = np.rint( 0.0 + 1.0 * (scale_var_SER ** 0.301) )
params['Helicopter Count'] = np.rint( 0.0 + 0.5 * (scale_var_SER ** 0.301) )

# **************************************************************************************************************************
#                                           Sec. 12: Post Processing
# **************************************************************************************************************************
params['Number of Samples'] = 100  # number of samples for cost uncertainty analysis
# Estimate costs using the cost database file and save the output to an Excel file
estimate = detailed_bottom_up_cost_estimate('cost/Cost_Database.xlsx')
elapsed_time = (time.time() - time_start) / 60  # calculate execution time
print('Execution time:', np.round(elapsed_time, 1), 'minutes')
