# This is oversimplified and need to be modified
def calculate_heat_exchanger_mass(params):
    HX_mass_intermediate =1000 *  0.8 * params['power_MW_th'] # Kg
    HX_mass_primary = 1000 * 1.365 * params['power_MW_th'] # Kg
    return  HX_mass_primary, HX_mass_intermediate,

# This is oversimplified and need to be modified

def calculate_pump_mechanical_power(params):
    primary_pump_mech_power   =  2.1514 * params['power_MW_th'] # kw
    sec_pump_mech_power  = 1.984   * params['power_MW_th']