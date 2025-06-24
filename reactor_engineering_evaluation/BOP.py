# This is oversimplified and need to be modified
def calculate_heat_exchanger_mass(params):
    HX_mass_primary =1000 *  0.8 * params['Power MWt'] # Kg 
    HX_mass_seconary = 1000 * 1.365 * params['Power MWt'] # Kg
    return  HX_mass_primary, HX_mass_seconary,

# This is oversimplified and need to be modified

def calculate_pump_mechanical_power(params):
    primary_pump_mech_power   =  2.1514 * params['Power MWt'] # kw
    sec_pump_mech_power  = 1.984   * params['Power MWt']
    return primary_pump_mech_power, sec_pump_mech_power
