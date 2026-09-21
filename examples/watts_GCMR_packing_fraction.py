# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

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
    'plotting': "N",  # "Y" or "N": Yes or No
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
for params['Packing Fraction'] in np.linspace(0.25, 0.35, 2):
    update_params({
        'Fuel Pin Materials': ['UCO', 'buffer_graphite', 'PyC', 'SiC', 'PyC'],
        'Fuel Pin Radii': [0.0250, 0.0350, 0.0390, 0.0425, 0.0465],  # cm # https://art.inl.gov/NRC%20Training%202019/04_TRISO_Fuel.pdf
        'Compact Fuel Radius': 0.6225,  # cm
        'Coolant Channel Radius': 0.35,  # cm
        'Moderator Booster Radii': [0.55],  # cm
        'Lattice Pitch': 2.25,
        'Assembly Rings': 6,
        'Core Rings': 5,
    })
    params['Assembly FTF'] = params['Lattice Pitch']*(params['Assembly Rings']-1)*np.sqrt(3)
    params['Radial Reflector Thickness'] = 27.393 # cm # radial reflector
    params['Axial Reflector Thickness'] = params['Radial Reflector Thickness'] # cm
    params['Core Radius'] = params['Assembly FTF']*params['Core Rings'] + params['Radial Reflector Thickness']
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
        'Heat Flux Criteria': 0.9,  # MW/m^2
        'Burnup Steps': [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0,
                         30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0]  # MWd_per_Kg
        })

    params['Power MWe'] = params['Power MWt'] * params['Thermal Efficiency'] 
    params['Heat Flux'] = calculate_heat_flux_TRISO(params) # MW/m^2

    # **************************************************************************************************************************
    #                                           Sec. 5: Running OpenMC
    # **************************************************************************************************************************

    # --- Shutdown Margin (SDM) ---
    # Set to True to run an additional ARI simulation and calculate the shutdown margin (pcm).
    # Currently disabled to save computation time. Enable for final design verification.
    # See watts_exec_LTMR.py for an example with Shutdown Margin Calc = True.
    params['Shutdown Margin Calc'] = False  # True or False

    # --- Isothermal Temperature Coefficient ---
    # Set to True to calculate the temperature reactivity coefficient (pcm/K).
    # Currently disabled to save computation time. Enable for safety analysis.
    # See watts_exec_LTMR.py for an example with Isothermal Temperature Coefficients = True.
    params['Isothermal Temperature Coefficients'] = False  # True or False

    # --- Temperature Perturbation ---
    # Required ONLY when 'Isothermal Temperature Coefficients' is True.
    # Units: Kelvin. Typical range: 50-300 K.
    # Uncomment and set this parameter if enabling the temperature coefficient calculation above.
    # params['Temperature Perturbation'] = 100  # K

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
        'Primary Loop Count': 2,
        'Primary Loop per loop load fraction': 0.5,
        'Primary Loop Inlet Temperature': 300 + 273.15, # K
        'Primary Loop Outlet Temperature': 550 + 273.15, # K
        'Secondary Loop Inlet Temperature': 270 + 273.15, # K — cold-end PCHE pinch 30°C with 300°C primary inlet (was 290 -> 10°C pinch, below realistic PCHE design)
        'Secondary Loop Outlet Temperature': 500 + 273.15, # K,
        'Primary Loop Pressure Drop': 50e3, # Pa
    })
    params['Primary HX Mass'] = calculate_heat_exchanger_mass(params)  # Kg
    mass_flow_rate(params)
    compressor_power(params)

    params.update({
        'BoP Count': 2,
        'BoP per loop load fraction': 0.5,
        })
    params['BoP Power kWe'] = 1000 * params['Power MWe'] * params['BoP per loop load fraction']

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

    total_refueling_period = params['Fuel Lifetime'] + params['Refueling Period'] + params['Startup Duration after Refueling']
    total_refueling_period_yr = total_refueling_period/365
    params['A75: Outer Vessel Structure Replacement Period (years)'] = 20
    params['A75: Inner Vessel Structure Replacement Period (cycles)'] = np.floor(10/total_refueling_period_yr)
    params['A75: Reflector Replacement Period (cycles)']   = np.floor(10/total_refueling_period_yr)
    params['A75: Reactor Control Devices Replacement Period (cycles)'] = np.floor(10/total_refueling_period_yr)
    params['A75: Moderator Booster Replacement Period (cycles)'] = 1
params['A75: Moderator Replacement Period (cycles)'] = 1
    params['Maintenance to Direct Cost Ratio']              = 0.015
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

    # --- Tax Credits (ITC / PTC) ---
    # No tax credits are applied in this example.
    # To apply ITC (Investment Tax Credit), uncomment and set:
    #   params['ITC credit level'] = 0.30  # fraction (e.g. 0.06, 0.30, 0.40, 0.50)
    #   See watts_exec_LTMR.py for full documentation on ITC parameters.
    #
    # To apply PTC (Production Tax Credit), uncomment and set:
    #   params['PTC credit value'] = 15.0   # $/MWh
    #   params['PTC credit period'] = 10    # years
    #   params['Tax Rate'] = 0.21           # fraction
    #   params['domestic_content_bonus'] = 0.10  # fraction (optional)
    #   params['energy_community_bonus'] = 0.10  # fraction (optional)
    #   See watts_exec_GCMR_Design_A.py for full documentation on PTC parameters.
    #
    # Note: ITC and PTC are mutually exclusive — only one can be selected per project.

    # **************************************************************************************************************************
    #                                           Sec. 11: Post Processing
    # **************************************************************************************************************************
    params['Number of Samples'] = 100
    tracked_params_list = ["Packing Fraction", "Lattice Pitch", "Number Of TRISO Particles Per Compact Fuel",
                           "Total Number of TRISO Particles", "Core Radius", "Heat Flux", "Fuel Lifetime",
                           "Mass U235", "Mass U238", "Uranium Mass", "Max Peaking Factor"]

    parametric_studies('cost/Cost_Database.xlsx', tracked_params_list)

    elapsed_time = (time.time() - time_start) / 60
    print('Execution time:', np.round(elapsed_time, 2), 'minutes')
