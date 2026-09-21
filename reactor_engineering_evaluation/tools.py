# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
import numpy as np 

from core_design.openmc_materials_database import collect_materials_data

def ellipsoid_shell(a, b, c):
    return 4*np.pi*np.power(((a*b)**1.6 + (a*c)**1.6 + (b*c)**1.6)/3, 1/1.6)

def circle_area(r):
    return (np.pi) * r **2


def materials_densities(material, params=None):
    material_densities = {
    "stainless_steel": 8.0,  # Approximate density of stainless steel
    "SS316": 8.0,            # Approximate density of SS316
    "SS304": 7.93,           # Approximate density of SS304
    "low_alloy_steel": 7.85, # Approximate density of SA508 Gr3 Cls 1
    "SA508": 7.85,           # Approximate density of SA508 Gr3 Cls 1
    "B4C_enriched": 2.52,    # Approximate density of boron carbide
    "B4C_natural": 2.52,     # Approximate density of boron carbide
    "WEP": 1.1,              # WEP density (water extended polymer)
    }
    if material in material_densities:
        return material_densities[material]  # in gram/cm^3
    if params is None:
        raise KeyError(
            f"No engineering density is available for material {material!r}."
        )
    materials_database = collect_materials_data(params)
    if material not in materials_database:
        raise KeyError(
            f"Material {material!r} is not in the MOUSE materials database."
        )
    density = materials_database[material].density
    if density is None:
        raise ValueError(f"Material {material!r} has no density.")
    return float(density)

def material_specific_heat(material):
    material_cp = {
        "Helium": 5193,  # J/(kg·K)
        "NaK": 982.      # J/(kg·K)
    }
    return material_cp[material]  # J/(kg·K)

def cylinder_annulus_mass(
    outer_radius,
    inner_radius,
    height,
    material,
    params=None,
):

    volume = 3.14* (outer_radius**2 - inner_radius**2) * height
    mass = volume * materials_densities(material, params=params) / 1000  # kg
    return mass # in kg

def calculate_shielding_masses(params):
    params['In Vessel Shield Mass'] = cylinder_annulus_mass(params['In Vessel Shield Outer Radius'],\
    params['In Vessel Shield Inner Radius'], params['Vessel Height'], params['In Vessel Shield Material'], params=params)

    # Represent the out-of-vessel shield as an open-top cylindrical enclosure:
    # an annular side shield surrounding the vessel system plus a full circular
    # bottom shield.  A solid top shield is intentionally excluded so that the
    # reactor-service penetrations and RVACS intake/exhaust path remain open.
    shield_thickness = float(params['Out Of Vessel Shield Thickness'])
    inner_radius = float(params['Vessels Total Radius'])
    outer_radius = inner_radius + shield_thickness
    shield_height = float(params['Vessels Total Height'])

    params['Outer Shield Inner Radius'] = inner_radius
    params['Outer Shield Outer Radius'] = outer_radius

    side_volume = np.pi * (outer_radius**2 - inner_radius**2) * shield_height
    bottom_volume = np.pi * outer_radius**2 * shield_thickness
    total_shield_volume = side_volume + bottom_volume

    material_density = materials_densities(
        params['Out Of Vessel Shield Material'],
        params=params,
    )
    effective_density_factor = float(
        params['Out Of Vessel Shield Effective Density Factor']
    )
    params['Out Of Vessel Shield Mass'] = (
        total_shield_volume
        * material_density
        * effective_density_factor
        / 1000
    )

def mass_flow_rate(params):
    loop_factor = 1
    thermal_power_MW = params['Power MWt']
    if 'Primary Loop per loop load fraction' in params.keys():
        loop_factor = params['Primary Loop per loop load fraction']
        thermal_power_MW = params['Power MWt'] * loop_factor
        
    deltaT =  params['Primary Loop Outlet Temperature'] - params['Primary Loop Inlet Temperature']
    if params['reactor type'] == "HPMR":
        coolant = params['Secondary Coolant']
    else:    
        coolant = params['Coolant']
    coolant_specific_heat = material_specific_heat(coolant)
    m_dot = 1e6 * thermal_power_MW/ (deltaT * coolant_specific_heat)
    params['Coolant Mass Flow Rate']  = m_dot / loop_factor # For Reactor Mass Flow Rate
    params['Primary Loop Mass Flow Rate'] = m_dot # For individual Primary Loop Mass Flow Rate
    
def compressor_power(params):
    # Estimates the required compressor power based on a simplified
    # model using pressure drop and compressor isentropic efficiency

    rho_he = 3.3297  # kg/m3. TODO: Consider importing CoolProp to estimate density based on cold-leg temperature and pressure
    power = params['Primary Loop Pressure Drop']*params['Primary Loop Mass Flow Rate']/params['Compressor Isentropic Efficiency']/rho_he
    params['Primary Loop Compressor Power'] = power # W
    return

def compressor_wheel_diameter(params):
    # Estimates the approximate compressor size based on its specific
    # diameter, matched to the MIGHTR horizontal HTGR design.
    # Ref for specific diameter:
    #  https://www.dropbox.com/scl/fi/fnqdg2hyi6y4ozu9p7nyu/final-report-str-mech-ARDP-redacted-V3.pdf?rlkey=h97dii28tvf0bxtffo8q62tn5&st=zsls1bs2&dl=0
    ref_specific_diameter = 3.6  # dimensionless
    rho_He = 3.330  # kg/m3 for He at 4 MPa, 300 °C. TODO: use a He density correlation or CoolProp to estimate density based on cold-leg temperature and pressure
    Vdot_gcmr = params['Primary Loop Mass Flow Rate'] / rho_He  # m3/s — volumetric flow rate
    dP = params['Primary Loop Pressure Drop']
    diameter = ref_specific_diameter/1.054 / (dP/rho_He)**0.25 * np.sqrt(Vdot_gcmr) # m
    return diameter

def GCMR_integrated_heat_transfer_vessel(params):
    # Calculates the required parameters for the
    # GCMR Integrated Heat Transfer Vessel that houses:
    #   circulator, PCHE, piping, valves, insulation

    contingency = 0.3  # accounts for the volume/mass of valves, fittings, and connections
    PCHE_volume = (params['Primary HX Mass'] / (materials_densities(params['HX Material'])*1e3) / 0.4)  # accounts for assumed 60% coolant channel void fraction
    compressor_volume = (compressor_wheel_diameter(params))**3  # approximated as a cube with side = wheel diameter. TODO: improve compressor sizing

    vessel_inner_volume = (1+contingency)*(PCHE_volume + compressor_volume)  # assumes a cube-like structure
    vessel_outer_volume = (vessel_inner_volume**(1/3)+ 1e-2*params['Integrated Heat Transfer Vessel Thickness'])**3  # m3
    vessel_volume = vessel_outer_volume - vessel_inner_volume
    vessel_density = materials_densities(params['Integrated Heat Transfer Vessel Material'])*1e3
    params['Integrated Heat Transfer Vessel Outer Volume'] = vessel_outer_volume
    params['Integrated Heat Transfer Vessel Mass'] = vessel_volume * vessel_density
    
    if params['Integrated Heat Transfer Vessel Thickness'] == 0:
        params['Integrated Heat Transfer Vessel Outer Volume'] = 0
        params['Integrated Heat Transfer Vessel Mass'] = 0

    # Rough estimate of the mass supported by the support structure:
    # Primary HX + Integrated Heat Transfer Vessel + Compressor + Valves/Fittings/Bolts/etc.
    # params['Integrated Heat Transfer System Mass'] = params['Primary HX Mass'] + (vessel_volume * vessel_density) + compressor_volume*8000
