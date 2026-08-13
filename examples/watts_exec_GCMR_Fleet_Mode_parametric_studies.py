# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

"""
This script performs a parametric cost study for a Gas Cooled Microreactor (GCMR)
with Fleet Mode enabled. Each run records the standard reactor cost metrics and,
automatically, the manufacturing- and servicing-campus OCC, TCI, annual cost,
LCOE, per-reactor metrics, high-level accounts 10, 20, 30, 40, 60, and 70,
and their standard deviations.

Fleet Mode represents a reactor fleet supported by two shared campuses:
  - Manufacturing Campus: manufactures and tests new reactor units.
  - Servicing Campus: services operating reactors and manages fleet logistics.

The fleet parameters are derived primarily from the annual reactor production rate.
The cost-engine integration and separate campus outputs will be implemented
independently from this example configuration.

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
from cost.cost_estimation import parametric_studies

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
# run_openmc(build_openmc_model_GCMR, heat_flux_monitor, params)
params['Fuel Lifetime'] = 2003 # days
params['Mass U235'] = 80972 # grams
params['Mass U238'] = 327919 # grams
params['Uranium Mass'] = 409 # Kg

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
    'Number of Operators': 2,
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

# **************************************************************************************************************************
#                                           Sec. 11: Fleet Mode
# **************************************************************************************************************************
# Fleet Mode estimates shared manufacturing and servicing campuses for a fleet of
# GCMRs. Campus sizes and resource requirements are driven primarily by the annual
# reactor production rate.
def update_fleet_mode_params(production_rate):
    """Recalculate every Fleet Mode parameter for one production-rate case."""
    params['Fleet Mode'] = True

    ## User Inputs

    params['Production Rate'] = production_rate #reactors per year
    params['Deployment Period'] = 10 # years
    params['Fleet'] = params['Deployment Period'] * params['Production Rate']
    params['Generating Sites Count'] = params['Fleet']  # i think that there was something here to account for capacity factor or something but ..
    #my concern here is that we are double counting since the capacity factor is already included in the electricity generation estimation

    params['Average Distance From Serv to GenSite'] = 1000 #miles (statutory miles)

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

    params['Manufacturing Campus Area'] = 56.1651 +	1.2066 * params['Production Rate'] ** 0.3919

    # params['MFG Emergency Generator Power'] = -10851.2764 + 10851.2773 * (params['Production Rate'] ** 0.0001)
    # params['MFG Warehouse Building Area'] = 5598.3987 + 70.3374 * (params['Production Rate'] ** 0.5539)
    # params['MFG Administration Building Area'] = 3947.3106 + 77.6167 * (params['Production Rate'] ** 0.4768)
    # params['MFG Warehouse Staff'] = np.ceil( 0.8576 + 0.2539 * (params['Production Rate'] ** 0.6533) )
    # params['MFG Security Staff Per Shift'] = np.ceil( 2.9981 + 0.3340 * (params['Production Rate'] ** 0.4768) )
    # params['MFG Maintenance Staff Headcount'] = np.ceil( 2.9981 + 0.3340 * (params['Production Rate'] ** 0.4768) )
    # params['MFG Local Transport Vehicle Count'] = np.rint( -4340.5105 + 4340.5109 * (params['Production Rate'] ** 0.0001) )
    # params['MFG Utility Vehicle Count'] = np.rint( -4339.5105 + 4340.5109 * (params['Production Rate'] ** 0.0001) )
    # # params['MFG Guard Station Count'] = 1.0000 + 1e-15 * (params['Production Rate'] ** 5.0000)
    # if params['Production Rate'] <= 500:
    #     params['MFG Guard Station Count'] = 1

    # elif params['Production Rate'] < 2000:
    #     params['MFG Guard Station Count'] = 2

    # else:
    #     params['MFG Guard Station Count'] = np.rint(params['Production Rate']/1000)


    # params['MFG Campus Power'] = 1.295614 +	0.00813495 * params['Production Rate'] **  0.620945 #MWe
    # params['MFG Switchyard Rating'] = np.ceil(params['MFG Campus Power'] * 2)

    # params['MFG Testing Line Annual Rate'] = int(np.floor(8000/92))
    # params['MFG Testing Line Count'] = np.ceil(params['Production Rate'] / params['MFG Testing Line Annual Rate'])
    # params['MFG Testing Engineering Headcount'] = 6 + ((params['MFG Testing Line Count'] - 1)*2)
    # params['MFG Testing Coolant Allotment'] = 0.15 * (1 * 24.417 + 9 * 11.114) * 8.2402


    #ToDo: Road Length and Site Perimeter shouldn't be the same.
    # params['MFG Road Length'] =  -1539388.0000975644 + 1540607.3153829477 * (params['Production Rate'] ** 8.588964624690402e-05)
    # params['MFG Site Perimeter'] = -1539388.0000975644 + 1540607.3153829477 * (params['Production Rate'] ** 8.588964624690402e-05)
    # params['MFG Protected Perimeter'] = 1421.7487567019996 + 60.41360643368738 * (params['Production Rate'] ** 0.22340396387855505)

    params['MFG Security Camera Count'] = 200 + -157.0909090909111 * (params['Production Rate'] ** -0.037788560889399164)
    params['MFG Motion Detector Count'] = 93.59999999999997 + 2.8444444444444628 * (params['Production Rate'] ** 0.35218251811136164)

    # Manufacturing-campus inputs required by the cost database. These must be
    # recalculated inside the loop for every production-rate case.
    params['MFG Construction Duration'] = 120  # months; REVIEW NEEDED.
    params['MFG Campus Land Area'] = params['Manufacturing Campus Area']
    params['MFG Emergency Generator Power'] = -10851.2764 + 10851.2773 * params['Production Rate'] ** 0.0001
    params['MFG Warehouse Building Area'] = 5598.3987 + 70.3374 * params['Production Rate'] ** 0.5539
    params['MFG Administration Building Area'] = 3947.3106 + 77.6167 * params['Production Rate'] ** 0.4768
    params['MFG Warehouse Staff'] = np.ceil(0.8576 + 0.2539 * params['Production Rate'] ** 0.6533)
    params['MFG Security Staff Per Shift'] = np.ceil(2.9981 + 0.3340 * params['Production Rate'] ** 0.4768)
    params['MFG Maintenance Staff Headcount'] = np.ceil(2.9981 + 0.3340 * params['Production Rate'] ** 0.4768)
    params['MFG Local Transport Vehicle Count'] = np.rint(-4340.5105 + 4340.5109 * params['Production Rate'] ** 0.0001)
    params['MFG Utility Vehicle Count'] = np.rint(-4339.5105 + 4340.5109 * params['Production Rate'] ** 0.0001)
    params['MFG Fire Station Count'] = 1  # REVIEW NEEDED.
    if params['Production Rate'] <= 500:
        params['MFG Guard Station Count'] = 1
    elif params['Production Rate'] < 2000:
        params['MFG Guard Station Count'] = 2
    else:
        params['MFG Guard Station Count'] = np.rint(params['Production Rate'] / 1000)

    params['MFG Campus Power'] = 1.295614 + 0.00813495 * params['Production Rate'] ** 0.620945
    params['MFG Switchyard Rating'] = np.ceil(params['MFG Campus Power'] * 2)
    params['MFG Testing Line Annual Rate'] = int(np.floor(8000 / 92))
    params['MFG Testing Line Count'] = np.ceil(params['Production Rate'] / params['MFG Testing Line Annual Rate'])
    params['MFG Testing Engineering Headcount'] = 6 + (params['MFG Testing Line Count'] - 1) * 2
    params['MFG Testing Number of Operators'] = 5 * params['MFG Testing Line Count']  # REVIEW NEEDED.
    params['MFG Testing Coolant Allotment'] = 0.15 * (1 * 24.417 + 9 * 11.114) * 8.2402
    params['MFG Local Control Building Area'] = 40 * params['MFG Testing Line Count']  # m^2; REVIEW NEEDED.
    params['MFG Security Building Area'] = 8775 / (3.2808 ** 2)  # m^2; REVIEW NEEDED.
    params['MFG Road Length'] = -1539388.0000975644 + 1540607.3153829477 * params['Production Rate'] ** 8.588964624690402e-05
    params['MFG Site Perimeter'] = -1539388.0000975644 + 1540607.3153829477 * params['Production Rate'] ** 8.588964624690402e-05
    params['MFG Protected Perimeter'] = 1421.7487567019996 + 60.41360643368738 * params['Production Rate'] ** 0.22340396387855505


    ## SERVICING CAMPUS

    params['SER Construction Duration'] = 120 #months; carried over from the former central-facility assumption.
    params['Servicing Rate'] = params['Fleet'] / params['Cycle Length']  # reactors/year

    params['SER Campus Area'] = (760 - 160.5) * (params['Servicing Rate'] - 30) / (300 - 30) + 160.5
    params['SER Campus Land Area'] = params['SER Campus Area']

    scale_var_SER = np.rint(params['Servicing Rate'] / 3) #Value of 3 should remain hardcoded. Couldn't let production rate be the variable for the Servicing Campus directly, for flexibility, but the relationships were built based on production rate with our basic assumption of a 3:1 ratio between servicing rate and production rate. Simple fix was to define this variable (scale_var_SER).

    params['SER Switchyard Rating'] = 0 + 12 * scale_var_SER ** 0.698970
    params['SER Switchyard Average Power'] = params['SER Switchyard Rating']/2

    params['Servicing Hot Cell Annual Rate'] = int(np.floor(365* (11/12) / 3 ))
    params['Servicing Hot Cell Count'] = np.ceil( params['Servicing Rate'] / params['Servicing Hot Cell Annual Rate'] )
    params['Radioactive Waste Processing Hot Cell Count'] = 1
    params['He Gas Replenishment Per Hot Cell'] = ((3*3*5) * 2 * params['Servicing Hot Cell Annual Rate'] * 2 + 0.1 * (10*30*7)*12) * params['m3_to_kg_He_RT_atmospheric']
    params['He Gas Replenishment'] = (params['Servicing Hot Cell Count'] * params['He Gas Replenishment Per Hot Cell'] + params['Radioactive Waste Processing Hot Cell Count'] * params['He Gas Replenishment Per Hot Cell'] + params['CoolantInventoryRPV_Mass'] * params['Servicing Rate'])

    params['SER Number of Operators Per Shift'] = np.ceil( 0.0 + 5.625 * (scale_var_SER ** 0.426) )
    params['SER Engineering Headcount'] = np.ceil( 0.0 + 20.0 * (scale_var_SER ** 0.301) )
    params['SER Maintenance Staff Per Shift'] = 40  # REVIEW NEEDED: provisional value carried over from the old central-facility example.
    params['SER Security Staff Per Shift'] = 1  # REVIEW NEEDED: provisional value carried over from the old reactor-site staffing input.

    params['Roundtrip Time'] = 2*params['Average Distance From Serv to GenSite'] / 40 / (8760/2) #[years] based on average speed of 40 miles/hour, 12 hours of driving per day
    params['Roundtrip Time Reactor Transport'] = 2*params['Average Distance From Serv to GenSite'] / 15 / (8760/2) #[years] based on average speed of 15 miles/hour, 12 hours of driving per day
    params['Dwell Time GenSite'] = 1/365 #[years]
    params['Dwell Time Serv'] = 1/365 #[years]
    params['Dwell Time Reactor Transport GenSite'] = 2/365 #[years]
    params['Dwell Time Reactor Transport Serv'] = 2/365 #[years]

    #ToDo: Road Length and Site Perimeter shouldn't be the same.
    params['SER Road Length'] = -1539388.0000975644 + 1540607.3153829477 * (scale_var_SER ** 8.588964624690402e-05)
    params['SER Site Perimeter'] = -1539388.0000975644 + 1540607.3153829477 * (scale_var_SER ** 8.588964624690402e-05)
    params['SER Controlled Perimeter'] = 0
    params['SER Protected Perimeter'] = 1421.7487567019996 + 60.41360643368738 * (scale_var_SER ** 0.22340396387855505)

    params['Used Fuel Storage Lifetime Capacity'] = params['Generating Sites Count'] * 50 / (params['Cycle Length'] + params['GenSite Downtime']) * params['Fuel Mass In Core']

    params['SER Security Building Area'] = 8775 / (3.2808 ** 2)
    params['SER Administration Building Area'] = 258000 / (3.2808 ** 2)

    params['SER Helium Flowrate'] = params['He Gas Replenishment'] / (0.9 * 8766 * 3600)

    params['SER Local Transport Vehicle Count'] = 80  # REVIEW NEEDED: provisional value carried over from the old central-facility example.
    params['SER Utility Vehicle Count'] = 100  # REVIEW NEEDED: provisional value based on the old general transport vehicle count.
    params['Reactor Transport Vehicle Count'] = np.ceil((params['Fleet'] / params['Cycle Length']) * (params['Roundtrip Time Reactor Transport']+params['Dwell Time Reactor Transport GenSite']+params['Dwell Time Reactor Transport Serv']) + 1)
    params['Helium Transport Truck Count'] = np.ceil(params['Fleet'] * params['Annual Coolant Supply Frequency'] * ((params['Roundtrip Time']+params['Dwell Time GenSite']+params['Dwell Time Serv']) * 1.05))
    params['Water Tanker Truck Count'] = np.ceil(params['Fleet'] * params['Water Supply Frequency'] * ((params['Roundtrip Time']+params['Dwell Time GenSite']+params['Dwell Time Serv']) * 1.05))
    params['Maintenance Truck Count'] = np.ceil(params['Fleet'] * params['Maintenance Visit Frequency'] * ((params['Roundtrip Time']+params['Dwell Time GenSite']+params['Dwell Time Serv']) * 1.05))

    params['SER Radiation Monitor Count'] = np.ceil(params['SER Site Perimeter'] / 1000 * 1.2)
    params['SER Security Camera Count'] = np.ceil(params['SER Campus Area'] / params['Manufacturing Campus Area'] * params['MFG Security Camera Count'])
    params['SER Motion Detector Count'] = np.ceil(params['SER Campus Area'] / params['Manufacturing Campus Area'] * params['MFG Motion Detector Count'])

    params['Reactor Transport Cask Count'] = round(params['Fleet'] * ((4 / 12) / params['Cycle Length']) * 1.5, -1)
    params['Annual Used Fuel Cask Consumption'] = 1 * params['Servicing Rate']
    params['Annual Reactor Cask Replacement'] = np.ceil(0.05 * params['Reactor Transport Cask Count'])
    params['Annual Radwaste Cask Consumption'] = 0.5 * params['Servicing Rate']

    params['Servicing Hot Cell Building Area'] = 0.0 + 411.428571 * (scale_var_SER ** 0.942)
    params['Helium Purification and Storage Building Area'] = 0.0 + 72.0 * (scale_var_SER ** 0.6198)
    params['SER Local Control Building Area'] = 40 * scale_var_SER
    params['SER Remote Control Building Area'] = 40 * scale_var_SER
    params['SER Rad Waste Management Building Area'] = 0.0 + 40.5 * (scale_var_SER ** 0.7447)
    params['Radwaste Storage Warehouses Area'] = 0.0 + 15422.94 * (scale_var_SER ** 0.994)
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
# Inputs/design outputs to include alongside the automatically tracked reactor
# and manufacturing/servicing-campus cost results. Add or remove parameter names as needed.
tracked_params_list = [
    'Production Rate',
    'Deployment Period',
    'Fleet',
    'Fuel Lifetime',
    'Cycle Length',
    'Servicing Rate',
    'Generating Sites Count',
]

# Each case recalculates all Fleet Mode dependencies and appends one row to this
# example's output CSV.
for production_rate in [1, 10, 100, 500, 1000]:
    update_fleet_mode_params(production_rate)
    parametric_studies('cost/Cost_Database.xlsx', tracked_params_list)
elapsed_time = (time.time() - time_start) / 60  # calculate execution time
print('Execution time:', np.round(elapsed_time, 1), 'minutes')
