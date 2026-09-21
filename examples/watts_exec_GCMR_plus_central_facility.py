# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

"""
This script performs a bottom-up cost estimate for a Gas Cooled Microreactor (GCMR)
WITH Central Facility cost estimation enabled.

This example extends the GCMR Design A example by adding central facility cost estimation.
A central facility is shared infrastructure (e.g., fuel handling, waste processing, control center)
that supports multiple reactor units deployed at the same site or region.

Key additions compared to watts_exec_GCMR_Design_A.py:
  - 'Estimate Central Facility': True  — enables central facility cost calculation
  - 'Maximum Number of Operating Reactors': 10  — number of reactors the facility supports
  - 'Central Facility Construction Duration': 24 months  — construction time for the facility

The output Excel file will contain an additional sheet "central facility cost estimate"
with the breakdown of central facility costs.

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
    'Enrichment': 0.1975,  # The enrichment is a fraction. It has to be between 0 and 1
    'UO2 atom fraction': 0.7,  # Mixing UO2 and UC by atom fraction
    'Radial Reflector': 'Graphite',
    'Axial Reflector': 'Graphite',
    'Matrix Material': 'Graphite',  # matrix material is the background material within the compact fuel element between TRISO particles
    'Moderator': 'Graphite',  # the moderator is outside the compact fuel region
    'Moderator Booster Materials': ['ZrH'],
    'Coolant': 'Helium',
    'Common Temperature': 850,  # Kelvins
    'Control Drum Absorber': 'B4C_enriched',  # The absorber material in the control drums
    'Control Drum Reflector': 'Graphite',  # The reflector material in the control drums
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
    'Packing Fraction': 0.3,
    
    # Coolant channel and booster dimensions
    'Coolant Channel Radius': 0.35,  # cm
    'Moderator Booster Radii': [0.55],  # cm
    'Lattice Pitch': 2.25,
    'Assembly Rings': 6,
    'Core Rings': 5,
})
params['Assembly FTF'] = params['Lattice Pitch']*(params['Assembly Rings']-1)*np.sqrt(3)
params['Radial Reflector Thickness'] = 27.393 # cm # radial reflector
params['Axial Reflector Thickness'] = params['Radial Reflector Thickness'] # cm
params['Core Radius'] = params['Assembly FTF']*params['Core Rings'] +  params['Radial Reflector Thickness']
params['Active Height'] = 250

# **************************************************************************************************************************
#                                           Sec. 3: Control Drums
# ************************************************************************************************************************** 
update_params({
    'Drum Radius': 9, # cm   
    'Drum Absorber Thickness': 1, # cm
    'Drum Height': params['Active Height'] + 2*params['Axial Reflector Thickness'],
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
    'Heat Flux Criteria': 0.9,  # MW/m^2 (needs review)
    'Burnup Steps': [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0,
                     30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0]  # MWd_per_Kg
    })

params['Power MWe'] = params['Power MWt'] * params['Thermal Efficiency'] 
params['Heat Flux'] = calculate_heat_flux_TRISO(params) # MW/m^2

# **************************************************************************************************************************
#                                           Sec. 5: Running OpenMC
# **************************************************************************************************************************
params['Particles'] = 2000
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
run_openmc(build_openmc_model_GCMR, heat_flux_monitor, params)
fuel_calculations(params)  # calculate the fuel mass and SWU

# **************************************************************************************************************************
#                                         Sec. 6: Primary Loop + Balance of Plant
# ************************************************************************************************************************** 
params.update({
    'Primary Loop Purification': True,
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

# Based on https://digital.library.unt.edu/ark:/67531/metadc893980/m2/1/high_res_d/919556.pdf (tables 17 and 18):
# Estimated helium mass per MWt is 3.3 kg/MWt.
params['Onsite Coolant Inventory'] = 3.3 * params['Power MWt']  # kg
# According to https://www.nationalacademies.org/read/12844/chapter/6#69, the helium loss rate is 10% per year,
# so 1/10 of the initial inventory is replenished annually.
# Without purification, helium needs to be replaced more frequently.
params['Replacement Coolant Inventory'] = params['Onsite Coolant Inventory'] / 10
params['Annual Coolant Supply Frequency'] = 1 if params['Primary Loop Purification'] else 6

total_refueling_period = params['Fuel Lifetime'] + params['Refueling Period'] + params['Startup Duration after Refueling'] # days
total_refueling_period_yr = total_refueling_period/365
params['A75: Outer Vessel Structure Replacement Period (years)'] = 20
params['A75: Inner Vessel Structure Replacement Period (cycles)'] = np.floor(10/total_refueling_period_yr)
params['A75: Reflector Replacement Period (cycles)']     = np.floor(10/total_refueling_period_yr)
params['A75: Reactor Control Devices Replacement Period (cycles)'] = np.floor(10/total_refueling_period_yr)
params['A75: Moderator Booster Replacement Period (cycles)'] = 1
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

# **************************************************************************************************************************
#                                           Sec. 11: Central Facility Costing
# **************************************************************************************************************************
# A central facility is shared infrastructure that supports multiple reactor units
# deployed at the same site or region. This includes:
#   - Servicing Facility: reactor refueling, defueling, and maintenance hot cells
#   - Manufacturing/Factory Facility: reactor component fabrication
#   - New Reactor Facility: fresh fuel storage, reactor fueling, and testing
#   - Radioactive Waste Management Facility: waste processing and storage
#   - Transportation infrastructure: vehicles and casks for reactor/fuel transport
#
# When 'Estimate Central Facility' is True, the cost estimation reads from the
# "Central Facility Database" sheet in Cost_Database.xlsx and produces a separate
# cost breakdown sheet in the output Excel file.
#
# All capacity/rate parameters for facilities and operations are assumed to be YEARLY
# unless otherwise specified.

# --- Overall Central Facility Parameters ---
update_params({
    'Estimate Central Facility': True,  # Enable central facility cost estimation

    # Maximum number of reactor units the central facility is designed to support.
    # Used to calculate fleet-wide metrics and normalize costs per kW.
    'Maximum Number of Operating Reactors': 100,  # total number of reactors served by all facilities

    # Construction duration for the central facility (may differ from reactor construction).
    # Used for calculating financing costs (interest during construction).
    'Central Facility Construction Duration': 120,  # months (this includes manufacturing facility and the service facility)

    # Total electrical capacity of the central facility itself (for its own operations).
    'Central Facility Power MWe': 50,  # MWe

    # Perimeter of the entire central facility site (for security fencing, etc.).
    'Site Perimeter': 20000,  # meters

    # Number of maintenance staff per shift at the central facility.
    'Maintenance Staff Per Shift': 40,  # FTEs per shift
})

# Derived parameter: total thermal power of the operating reactor fleet
params['Power Mwt of Operating Fleet'] = params['Power MWt'] * params['Maximum Number of Operating Reactors']

# --- Servicing Facility Parameters ---
# The servicing facility handles reactor refueling, defueling, inspection, and maintenance.
# It includes hot cells for handling irradiated components.
update_params({
    # Building volumes for servicing facility structures (concrete volumes in m^3)
    'Servicing Building Roof Volume': 500,  # m^3
    'Servicing Building Basement Volume': 500,  # m^3
    'Servicing Building Walls Volume': 300,  # m^3
    'Servicing Building Volume': 5000,  # m^3 (total enclosed volume)

    'Helium Purification and Storage Building Roof Volume': 100,  # m^3
    'Helium Purification and Storage Building Basement Volume': 100,  # m^3
    'Helium Purification and Storage Building Walls Volume': 80,  # m^3
    'Helium Purification and Storage Building Volume': 1000,  # m^3

    'Servicing Facility Integrated Control Room Roof Volume': 50,  # m^3 — mainly the hot cell
    'Servicing Facility Integrated Control Room Basement Volume': 50,  # m^3
    'Servicing Facility Integrated Control Room Walls Volume': 40,  # m^3
    'Servicing Facility Integrated Control Room Volume': 500,  # m^3

    'Servicing Facility Admin Building Roof Volume': 80,  # m^3
    'Servicing Facility Admin Building Basement Volume': 80,  # m^3
    'Servicing Facility Admin Building Walls Volume': 60,  # m^3
    'Servicing Facility Admin Building Volume': 800,  # m^3

    'Servicing Facility Security Building Roof Volume': 30,  # m^3
    'Servicing Facility Security Building Basement Volume': 30,  # m^3
    'Servicing Facility Security Building Walls Volume': 25,  # m^3
    'Servicing Facility Security Building Volume': 300,  # m^3

    # Perimeter of the servicing facility (for security fencing)
    'Servicing Facility Perimeter': 2000,  # meters

    # Number of reactors serviced per year
    'Total Servicing Rate': 50,  # reactors/year

    # Number of hot cells for reactor servicing operations
    'Servicing Hot Cell Count': 10,  # number of hot cells ('Total Servicing Rate' divided by hot cell capacity)

    # Volume of each servicing hot cell
    'Servicing Hot Cell Volume': 500,  # m^3 per hot cell (empty interior volume)

    # Number of defueling/refueling lines (typically equals hot cell count)
    'Defueling/Refueling Line Count': 10,  # number of lines
})

# Total volume of all servicing hot cells combined
params['Total Servicing Hot Cell Volume'] = params['Servicing Hot Cell Count'] * params['Servicing Hot Cell Volume']
# Thermal power processed by servicing facility (assumes 5% power for low-power testing per hot cell)
params['Power Mwt Processed by Servicing'] = 0.05 * params['Power MWt'] * params['Servicing Hot Cell Count']

# --- Manufacturing/Factory Facility Parameters ---
# The manufacturing facility fabricates reactor components and assembles new reactors.
update_params({
    # Building volumes for manufacturing facility structures
    'Fabrication Building Roof Volume': 400,  # m^3
    'Fabrication Building Basement Volume': 400,  # m^3
    'Fabrication Building Walls Volume': 300,  # m^3
    'Fabrication Building Volume': 4000,  # m^3

    'Feed and Product Warehouse Building Roof Volume': 200,  # m^3
    'Feed and Product Warehouse Building Basement Volume': 200,  # m^3
    'Feed and Product Warehouse Building Walls Volume': 150,  # m^3
    'Feed and Product Warehouse Building Volume': 2000,  # m^3

    'Manufacturing Facility Integrated Control Building Roof Volume': 50,  # m^3
    'Manufacturing Facility Integrated Control Building Basement Volume': 50,  # m^3
    'Manufacturing Facility Integrated Control Building Walls Volume': 40,  # m^3
    'Manufacturing Facility Integrated Control Building Volume': 500,  # m^3

    'Manufacturing Facility Admin Building Roof Volume': 80,  # m^3
    'Manufacturing Facility Admin Building Basement Volume': 80,  # m^3
    'Manufacturing Facility Admin Building Walls Volume': 60,  # m^3
    'Manufacturing Facility Admin Building Volume': 800,  # m^3

    'Manufacturing Facility Security Building Roof Volume': 30,  # m^3
    'Manufacturing Facility Security Building Basement Volume': 30,  # m^3
    'Manufacturing Facility Security Building Walls Volume': 25,  # m^3
    'Manufacturing Facility Security Building Volume': 300,  # m^3

    # Perimeter of the manufacturing/factory facility
    'Factory Perimeter': 3000,  # meters

    # Number of new reactors produced per year
    'New Reactor Production Rate': 20,  # reactors/year
})

# --- New Reactor Facility Parameters ---
# The new reactor facility handles fresh fuel storage, initial fueling, and reactor testing
# before deployment to field sites. # these are fresh reactors while the other servicing facility is for the used ones
update_params({
    # Building volumes for new reactor facility structures
    'Fresh Fuel Storage and Inspection Building Roof Volume': 150,  # m^3
    'Fresh Fuel Storage and Inspection Building Basement Volume': 150,  # m^3
    'Fresh Fuel Storage and Inspection Building Walls Volume': 120,  # m^3
    'Fresh Fuel Storage and Inspection Building Volume': 1500,  # m^3

    'Reactor Fueling Building Roof Volume': 200,  # m^3
    'Reactor Fueling Building Basement Volume': 200,  # m^3
    'Reactor Fueling Building Walls Volume': 150,  # m^3
    'Reactor Fueling Building Volume': 2000,  # m^3

    'Reactor Testing Building Roof Volume': 300,  # m^3
    'Reactor Testing Building Basement Volume': 300,  # m^3
    'Reactor Testing Building Walls Volume': 250,  # m^3
    'Reactor Testing Building Volume': 3000,  # m^3

    'New Reactor Fuel and Testing Facility Admin Building Roof Volume': 80,  # m^3
    'New Reactor Fuel and Testing Facility Admin Building Basement Volume': 80,  # m^3
    'New Reactor Fuel and Testing Facility Admin Building Walls Volume': 60,  # m^3
    'New Reactor Fuel and Testing Facility Admin Building Volume': 800,  # m^3

    'New Reactor Fuel and Testing Facility Security Building Roof Volume': 30,  # m^3
    'New Reactor Fuel and Testing Facility Security Building Basement Volume': 30,  # m^3
    'New Reactor Fuel and Testing Facility Security Building Walls Volume': 25,  # m^3
    'New Reactor Fuel and Testing Facility Security Building Volume': 300,  # m^3

    # Perimeter of the new reactor facility
    'New Reactor Facility Perimeter': 2500,  # meters

    # Number of fueling lines for new reactors
    'New Reactor Fueling Line Count': 5,  # number of lines

    # Number of testing lines/hot cells for new reactors
    'New Reactor Testing Line Count': 10,  # number of lines

    # Hot cell specifications for reactor testing
    'New Reactor Testing Hot Cell Count': 10,  # number of hot cells
    'New Reactor Testing Hot Cell Volume': 500,  # m^3 per hot cell
})

# Total volume of all new reactor testing hot cells
params['New Reactor Testing Hot Cell Volume'] = params['New Reactor Testing Hot Cell Count'] * params['New Reactor Testing Hot Cell Volume']

# Total electrical capacity processed by new reactor facility (production rate × reactor power)
params['Power Mwe Processed by New Reactor Facility'] = params['Power MWe'] * params['New Reactor Production Rate']

# --- Radioactive Waste Management Facility Parameters ---
# The radioactive waste management facility handles processing and storage of
# radioactive waste from reactor operations and servicing.
update_params({
    # Building volumes for waste management facility structures
    'Radioactive Waste Processing Building Roof Volume': 200,  # m^3
    'Radioactive Waste Processing Building Basement Volume': 200,  # m^3
    'Radioactive Waste Processing Building Walls Volume': 150,  # m^3
    'Radioactive Waste Processing Building Volume': 2000,  # m^3

    'Radioactive Waste Storage Building Roof Volume': 300,  # m^3
    'Radioactive Waste Storage Building Basement Volume': 300,  # m^3
    'Radioactive Waste Storage Building Walls Volume': 250,  # m^3
    'Radioactive Waste Storage Building Volume': 3000,  # m^3

    'Radioactive Waste Management Facility Integrated Control Room Roof Volume': 50,  # m^3
    'Radioactive Waste Management Facility Integrated Control Room Basement Volume': 50,  # m^3
    'Radioactive Waste Management Facility Integrated Control Room Walls Volume': 40,  # m^3
    'Radioactive Waste Management Facility Integrated Control Room Volume': 500,  # m^3

    'Radioactive Waste Management Facility Admin Building Roof Volume': 80,  # m^3
    'Radioactive Waste Management Facility Admin Building Basement Volume': 80,  # m^3
    'Radioactive Waste Management Facility Admin Building Walls Volume': 60,  # m^3
    'Radioactive Waste Management Facility Admin Building Volume': 800,  # m^3

    'Radioactive Waste Management Facility Security Building Roof Volume': 30,  # m^3
    'Radioactive Waste Management Facility Security Building Basement Volume': 30,  # m^3
    'Radioactive Waste Management Facility Security Building Walls Volume': 25,  # m^3
    'Radioactive Waste Management Facility Security Building Volume': 300,  # m^3

    # Perimeter of the radioactive waste management facility
    'Radioactive Waste Management Facility Perimeter': 2000,  # meters
})

# --- Transportation Parameters --- 
# Vehicles and equipment for transporting reactors, fuel, and waste between facilities.
update_params({
    # Local transport vehicles (within the central facility complex)
    'Local Transport Vehicle Count': 80,  # number of vehicles inside the facility

    # Vehicles for transporting complete reactor units to/from field sites
    'Reactor Transport Vehicle Count': 20,  # number of vehicles

    # Vehicles for transporting spare parts and radioactive waste
    'Spares/Waste Transport Vehicle Count': 50,  # number of vehicles

    # General purpose transport vehicles
    'General Transport Vehicle Count': 100,  # number of vehicles
})

# --- Cask Parameters ---
# Specialized containers for transporting spent fuel, reactors, and radioactive waste inside facility
# These are consumable items that need periodic replacement.
update_params({
    # Annual replacement rate for spent fuel transport casks
    'Annual Spent Fuel Cask Replacement': 20,  # casks/year — 20 casks are used and disposed of annually for spent fuel transport

    # Annual replacement rate for reactor transport casks
    'Annual Reactor Cask Replacement': 10,  # casks/year

    # Annual replacement rate for radioactive waste transport casks
    'Annual Rad Waste Cask Replacement': 50,  # casks/year
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
params['PTC credit value'] = 15.0  # $/MWh

# Duration of the PTC credit period.
# Under the IRA Section 45Y, the credit is available for 10 years after the facility
# is placed in service.
# Units: years
# Typical value: 10 years
params['PTC credit period'] = 10  # years

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
params['domestic_content_bonus'] = 0.10   # fraction — assumes domestic content standard is met
params['energy_community_bonus'] = 0.10   # fraction — assumes facility is in an energy community

# --- Corporate Tax Rate ---
# The US federal corporate tax rate used to gross up the PTC tax credit to its
# before-tax revenue equivalent in the LCOE calculation. Since MOUSE uses a
# before-tax LCOE, the PTC must be converted to a before-tax equivalent.
# The current US federal corporate tax rate is 21% (as of 2024).
# Municipal utilities and non-profit cooperatives may use 0.0 (tax-exempt).
# Units: fraction (e.g. 0.21 for 21%)
# Typical values: 0.21 (federal only), 0.27 (federal + average state)
params['Tax Rate'] = 0.21  # fraction

# **************************************************************************************************************************
#                                           Sec. 12: Post Processing
# **************************************************************************************************************************
params['Number of Samples'] = 100  # number of samples for cost uncertainty analysis
# Estimate costs using the cost database file and save the output to an Excel file
# Note: Output will include both reactor costs and central facility costs (separate sheets).
estimate = detailed_bottom_up_cost_estimate('cost/Cost_Database.xlsx')
elapsed_time = (time.time() - time_start) / 60  # calculate execution time
print('Execution time:', np.round(elapsed_time, 1), 'minutes')
