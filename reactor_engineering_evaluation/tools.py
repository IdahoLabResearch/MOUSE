import numpy as np 

def ellipsoid_shell(a, b, c):
    return 4*np.pi*np.power(((a*b)**1.6 + (a*c)**1.6 + (b*c)**1.6)/3, 1/1.6)

def circle_area(r):
    return (np.pi) * r **2


def materials_densities(material):
    material_densities = {
    "stainless_steel": 8.0,       # Approximate density of stainless steel
    "boron_carbide": 2.52,        # Approximate density of boron carbide
    "water_extended_polymer": 1.1 # WEP density
    }
    return material_densities[material] # in gram/cm^3

def cylinder_annulus_mass(outer_radius , inner_radius,height, material ):
    volume = 3.14* (outer_radius**2 - inner_radius**2) * height
    mass = volume* materials_densities(material)/1000 # Kilograms
    return mass # in kg


    
