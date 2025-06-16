import numpy as np 

def ellipsoid_shell(a, b, c):
    return 4*np.pi*np.power(((a*b)**1.6 + (a*c)**1.6 + (b*c)**1.6)/3, 1/1.6)

def circle_area(r):
    return (np.pi) * r **2


def materials_densities(material):
    material_densities = {
    "stainless_steel": 8.0,       # Approximate density of stainless steel
    "B4C_natural": 2.52,        # Approximate density of boron carbide
    "WEP": 1.1 # WEP density (water extended polymer)
    }
    return material_densities[material] # in gram/cm^3

# def material_specific_heat(material):
#     material_cp= {
#     "Helium": 5193      # J/(Kg.K)

#     }
#     return material_cp[material] # in gram/cm^3    

def cylinder_annulus_mass(outer_radius , inner_radius,height, material ):
    volume = 3.14* (outer_radius**2 - inner_radius**2) * height
    mass = volume* materials_densities(material)/1000 # Kilograms
    return mass # in kg

def calculate_shielding_masses(params):
    params['In Vessel Shield Mass'] = cylinder_annulus_mass(params['In Vessel Shield Outer Radius Radius'],\
    params['In Vessel Shield Inner Radius'], params['Vessel Height'], params['In Vessel Shielding Material'] )

    params['Outer Shield Outer Radius'] = params['Out Of Vessel Shield Thickness']+ params['Vessels Total Radius']

    outer_shield_mass = cylinder_annulus_mass(params['Outer Shield Outer Radius'], params['Out Of Vessel Shield Thickness'],\
    params['Vessels Total Height'], params['Out Of Vessel Shield Material']) 
    params['Out Of Vessel Shield Mass'] = params['Out Of Vessel Shield Effective Density Factor'] * outer_shield_mass

# def mass_flow_rate(thermal_power_MW, deltaT, coolant):
#     coolant_specific_heat = material_specific_heat(coolant)
#     m_dot = 1000000 * thermal_power_MW/ (deltaT * coolant_specific_heat)
#     return m_dot
    
