# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

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
for params['Fuel'] in ['TRIGA_fuel', 'UO2']:
    update_params({
        'reactor type': "LTMR",
        'TRISO Fueled': "No",
        'Enrichment': 0.1975,  # Fraction between 0 and 1
        "H_Zr_ratio": 1.6,
        'U_met_wo': 0.3,
        'er_wo': 0,         # Erbium (burnable poison)
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
    #                                           Sec. 2: Geometry
    # **************************************************************************************************************************  

    update_params({
        'Fuel Pin Materials': ['Zr', None, params['Fuel'], None, 'SS304'],
        'Fuel Pin Radii': [0.28575, 0.3175, 1.5113, 1.5367, 1.5875],  # cm
        'Moderator Pin Materials': ['ZrH', 'SS304'],  
        'Moderator Pin Inner Radius': 1.5367,  # cm
        'Moderator Pin Radii': [1.5367, 1.5875],  # cm
        "Pin Gap Distance": 0.1,  # cm
        'Pins Arrangement': LTMR_pins_arrangement,
        'Number of Rings per Assembly': 12,
        'Radial Reflector Thickness': 14,  # cm
    })

    params['Lattice Apothem'] = calculate_hex_apothem(params)
    params['Assembly FTF'] = 2 * params['Lattice Radius']
    params['Active Height'] = 78.4
    params['Axial Reflector Thickness'] = params['Radial Reflector Thickness']
    params['Fuel Pin Count'] = calculate_pins_in_assembly(params, "FUEL")
    params['Moderator Pin Count'] = calculate_pins_in_assembly(params, "MODERATOR")
    params['Moderator Mass'] = calculate_moderator_mass(params)
    params['Core Radius'] = calculate_core_radius_from_hex(params)

    # **************************************************************************************************************************
    #                                           Sec. 3: Control Drums
    # ************************************************************************************************************************** 

    update_params({
        'Number of Drums': 12,
        # When the user does not specify the drum radius, the code automatically sets it to the largest allowable value that avoids drum overlap
        #'Drum Radius': 9.016,  # cm
        'Drum Absorber Thickness': 1,  # cm
        'Drum Absorber Arc Degrees': 120,
    })

    update_ltmr_reflector_geometry_from_drums(params)
    calculate_drums_volumes_and_masses(params)
    calculate_reflector_mass_LTMR(params)

    # **************************************************************************************************************************
    #                                           Sec. 4: Overall System
    # ************************************************************************************************************************** 

    update_params({
        'Power MWt': 20,  # MWt
        'Thermal Efficiency': 0.31,
        'Heat Flux Criteria': 0.9,  # MW/m^2
        'Burnup Steps': [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 140.0]
    })
    params['Power MWe'] = params['Power MWt'] * params['Thermal Efficiency']
    params['Heat Flux'] = calculate_heat_flux(params)

    # **************************************************************************************************************************
    #                                           Sec. 5: Running OpenMC
    # **************************************************************************************************************************

    # --- Shutdown Margin  ---
    # When True, additional OpenMC simulations are run to evaluate shutdown margin.
    # The model is evaluated in both operating (ARO) and shutdown (ARI) drum configurations.
    # Shutdown margin is then computed from the two configurations.
    # Recommended: True for final design verification; can be set to False to save
    # computation time during early design exploration.
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
    run_openmc(build_openmc_model_LTMR, heat_flux_monitor, params)
    fuel_calculations(params)

    # **************************************************************************************************************************
    #                                           Sec. 6: Primary Loop + Balance of Plant
    # ************************************************************************************************************************** 

    update_params({
        'Secondary HX Mass': 0,
        'Primary Pump': 'Yes',
        'Secondary Pump': 'No',
        'Pump Isentropic Efficiency': 0.8,
        'Primary Loop Inlet Temperature': 430 + 273.15,
        'Primary Loop Outlet Temperature': 520 + 273.15,
        'Secondary Loop Inlet Temperature': 395 + 273.15,
        'Secondary Loop Outlet Temperature': 495 + 273.15,
    })

    params['Primary HX Mass'] = calculate_heat_exchanger_mass(params)
    params.update({
        'BoP Count': 2,
        'BoP per loop load fraction': 0.5,
        })
    params['BoP Power kWe'] = 1000 * params['Power MWe'] * params['BoP per loop load fraction']
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
        'Out Of Vessel Shield Effective Density Factor': 0.5
    })

    params['In Vessel Shield Outer Radius'] = params['Core Radius'] + params['In Vessel Shield Thickness']

    # **************************************************************************************************************************
    #                                           Sec. 8: Vessels Calculations
    # ************************************************************************************************************************** 

    update_params({
        'Vessel Radius': params['Core Radius'] + params['In Vessel Shield Thickness'],
        'Vessel Thickness': 2,  # cm — ASME BPVC Sec III Div 5 minimum at 520°C; refs: Oklo Aurora (INL-22/68167), FFTF, EBR-II
        'Vessel Lower Plenum Height': 50,  # cm — IAEA TECDOC-1908; inlet manifold + flow distributor + grid plate (was 2.848, unit-conv bug)
        'Vessel Upper Plenum Height': 47.152,  # cm — includes cover-gas headspace (no separate Upper Gas Gap tracked)
        'Vessel Upper Gas Gap': 0,
        'Vessel Bottom Depth': 32.129,
        'Vessel Material': 'stainless_steel',
        'Gap Between Vessel And Guard Vessel': 5,  # cm — OECD/NEA SFR Vessel Design Guidelines 2017; ASME Sec III Div 5 NH-3000 (was 2)
        'Guard Vessel Thickness': 1,  # cm — ASME Sec III Class 3 minimum + IAEA TECDOC-1531 (was 0.5)
        'Guard Vessel Material': 'stainless_steel',
        'Gap Between Guard Vessel And Cooling Vessel': 5,  # cm
        'Cooling Vessel Thickness': 0.5,  # cm
        'Cooling Vessel Material': 'stainless_steel',
        'Gap Between Cooling Vessel And Intake Vessel': 5,  # cm — Hejzlar & Buongiorno 2007 NED RVACS minimum (was 3)
        'Intake Vessel Thickness': 0.5,  # cm
        'Intake Vessel Material': 'stainless_steel'
    })

    vessels_specs(params)
    calculate_shielding_masses(params)

    # **************************************************************************************************************************
    #                                           Sec. 9: Operation
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
    ## Based on https://www.edf.fr/sites/default/files/mediatheque/dp_creys_2017.pdf :
    ## 5,500 tonnes of sodium from the reactor vessel and secondary circuits at the Creys-Malville plant (France), which is 3,000 MWt.
    # This gives a rough estimate of 1833 kg/MWt.
    params['Onsite Coolant Inventory'] = 1833 * params['Power MWt']
    params['Replacement Coolant Inventory'] = 0  # NaK is assumed not to require replacement

    total_refueling_period = params['Fuel Lifetime'] + params['Refueling Period'] + params['Startup Duration after Refueling']
    total_refueling_period_yr = total_refueling_period/365
    params['A75: Vessel Replacement Period (cycles)']      = np.floor(10/total_refueling_period_yr)
    params['A75: Core Barrel Replacement Period (cycles)'] = np.floor(10/total_refueling_period_yr)
    params['A75: Reflector Replacement Period (cycles)']   = np.floor(10/total_refueling_period_yr)
    params['A75: Drum Replacement Period (cycles)']        = np.floor(10/total_refueling_period_yr)
    params['Maintenance to Direct Cost Ratio']              = 0.015
    params['A78: CAPEX to Decommissioning Cost Ratio'] = 0.15

    # **************************************************************************************************************************
    #                                           Sec. 12: Buildings & Economic Parameters
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
        'NOAK Unit Number': 100
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
    tracked_params_list = ["Fuel", "Fuel Lifetime", "Mass U235", "Mass U238", "Uranium Mass", "Max Peaking Factor"]
    parametric_studies('cost/Cost_Database.xlsx', tracked_params_list)

    elapsed_time = (time.time() - time_start) / 60
    print('Execution time:', np.round(elapsed_time, 2), 'minutes')
