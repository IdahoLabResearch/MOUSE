"""
This script performs a bottom-up cost estimate for a Liquid Metal Thermal Microreactor (LTMR).
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
#    'cross_sections_xml_location': '/home/hannbn/openmc_data/endfb-viii.0-hdf5/cross_sections.xml',
#    'simplified_chain_thermal_xml': '/home/hannbn/openmc_data/simplified_thermal_chain11.xml'
})

# **************************************************************************************************************************
#                                                Sec. 1: Materials
# **************************************************************************************************************************

update_params({
    'reactor type': "HPMR",
    'TRISO Fueled': "No",
    'Fuel': 'homog_TRISO',    # The fuel use is homogenized triso particles
    'Enrichment': 0.19985,
    'Reflector' : 'Be',
    'Coolant': 'homog_hp', 
    'Moderator': 'graphite2',
    'Gap': 'Helium',
    'Control Drum Absorber': 'B4C_natural',  # The absorber material in the control drums
    'Control Drum Reflector': 'Be',
    'Common Temperature': 1000, 
})

# **************************************************************************************************************************
#                                           Sec. 2: Geometry: Fuel Pins, Moderator Pins, Coolant, Hexagonal Lattice
# **************************************************************************************************************************

update_params({
    'Fuel Pin Materials': ['homog_TRISO', 'Helium'],
    'Fuel Pin Radii': [1.00, 1.05],
    'Heat Pipe Materials': ['homog_hp', 'Helium'],
    'Heat Pipe Radii': [1.10, 1.15],
    "Pin Gap Distance": 0.46,
    'Number of Rings per Assembly': 6,
    'Number of Rings per Core': 3,
    'Lattice Pitch': 3.4,
    'Fuel Pin Count per Assembly': 72,
    'Height': 200, # Total height with top and bottom reflector used in the correction factor 
})

params['Assembly FTF'] = (params['Lattice Pitch']*(params['Number of Rings per Assembly']-1)*np.sqrt(3)) + (2*(params['Fuel Pin Radii'])[-1]) + params['Pin Gap Distance']  #32 cm
params['Reflector Thickness'] = params['Assembly FTF'] / 2
params['Core Radius'] = params['Assembly FTF']* params['Number of Rings per Core'] +  params['Reflector Thickness']     #112 cm
params['Active Height'] = (10/7) * params['Core Radius']  #160 cm    # From the ratio between core radius and active height of HPMR
params['hexagonal Core Edge Length'] = (params['Assembly FTF'] * (params['Number of Rings per Core']-1)) + (params['Assembly FTF']/2) + 6.6  # The edge lenght is 86.6 as in the originial input so 6.6 is added based on this value
params['Fuel Assemblies Count'] =  (3 * params['Number of Rings per Core']**2) - (3 * params['Number of Rings per Core'])
params['Fuel Pin Count'] = params['Fuel Assemblies Count'] * params['Fuel Pin Count per Assembly']
params['Energy Groups'] = np.array([0.0, 0.625, 1e4, 1e6])  # Three energy groups defined according to the energy spectrum in HPMR report for correction factor calculations
# **************************************************************************************************************************
#                                           Sec. 3: Control Drums
# ************************************************************************************************************************** 

update_params({
    'Drum Radius': params['Core Radius'] / 7.442, 
    'Drum Absorber Thickness': 1,  # cm
    'Drum Height': params['Active Height']
})

calculate_drums_volumes_and_masses(params)
calculate_reflector_mass_HPMR(params)

# **************************************************************************************************************************
#                                           Sec. 4: Overall System
# ************************************************************************************************************************** 

update_params({
    'Power MWt': 5, 
    'Thermal Efficiency': 0.36,
    'Heat Flux Criteria': 0.9,  # MW/m^2 
    #'Burnup Steps': [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 140.0]  # MWd_per_Kg
    'Time Steps': [0.01*tf,   0.99*tf,   3.00*tf,   6.00*tf,  20.00*tf,  70.00*tf, 100.00*tf, 165.00*tf, 365.00*tf, 365.00*tf, 365.00*tf],
    'Power': [5.0e+06, 5.0e+06, 5.0e+06, 5.0e+06, 5.0e+06, 5.0e+06, 5.0e+06, 5.0e+06, 5.0e+06, 5.0e+06, 5.0e+06],
})
params['Power MWe'] = params['Power MWt'] * params['Thermal Efficiency']
params['Heat Flux'] =  calculate_heat_flux(params)

heat_flux_monitor = monitor_heat_flux(params)
run_openmc(build_openmc_model_HPMR, heat_flux_monitor, params)
fuel_calculations(params)  # calculate the fuel mass and SWU