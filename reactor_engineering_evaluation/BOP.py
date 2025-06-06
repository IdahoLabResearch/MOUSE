from tools import *
from math import log


# This is oversimplified and need to be modified
def calculate_heat_exchanger_mass(hx_thermal_load,
                                  primary_coolant_heat_capacity,
                                  th_in,
                                  th_out,
                                  tc_in,
                                  tc_out):
    """
      Assuming a printed circuit heat exchanger (PCHE).
      Input for this function are as follows
        - hx_thermal_load : thermal load of PCHE in [MW]
        - th_in           : PCHE hot side inlet temprture [K]
        - th_out          : PCHE hot side outlet temperature [K]
        - tc_in           : PCHE cold side inlet temperature [K]
        - tc_out          : PCHE cold side outlet temperature [K]
    """
    # Assumption on overall heat transfer coefficient
    U = 500  # [w/m2/K] This is an average value obtained by scanning of literature

    # PCHE channel dimensions assumptions
    hx_channel_diameter = 0.0015   # [m]
    hx_channel_pitch = 0.00225     # [m]
    hx_channel_length = 1          # [m]
    hx_void_fraction = 0.6
    hx_channel_thick = 0.003       # [m]

    rho_ss = 7.8   #  density of stainless steel :: ton/3
    ss_unit_cost = 100000 # unit cost of material :: $/ton
    unit_cost_per_energy = 0.1  # [$/W]

    hx_channel_perimeter = pi* hx_channel_diameter/2 + hx_channel_diameter
    hx_channel_ht_area = hx_channel_perimeter* hx_channel_length   # assuming a semi circular channel

    delta_t1 = th_in - tc_out
    delta_t2 = th_out - tc_in

    LMTD = (delta_t1 - delta_t2)/log(delta_t1/delta_t2)
    ht_area = hx_thermal_load*1e6/(U* LMTD)
    nchannels = ht_area/hx_channel_ht_area
    hx_alloy_volume = nchannels* hx_channel_pitch* hx_channel_thick - nchannels* pi/8* hx_channel_diameter**2
    
    hx_mass = hx_alloy_volume* rho_ss
    hx_cost = hx_mass* ss_unit_cost



    return hx_mass, hx_cost

# This is oversimplified and need to be modified

def calculate_primary_pump_mechanical_power(core_mass_flow_rate, core_active_height):
    """
      Pump electric power [kW] = mdot*g*h / 1000 
    """
    g    = 9.81                               # [m/s^2]
    h    = core_active_height                 # [m] 
    mdot = core_mass_flow_rate                # [kg/s]

    return  mdot* g* h / 1000     # [kWe]

def calculate_secondary_pump_mechanical_power(secondary_mass_flow_rate):
    """
      Pump electric power [kW] = mdot*g*h / 1000
    """
    g    = 9.81                               # [m/s^2]
    h    = 58.56* 0.3048                      # [m]  I found that this was an assumption based on literature of reactors with similar power
    mdot = secondary_mass_flow_rate           # [kg/s]

    return mdot* g* h / 1000     # [kWe]




def calculate_building_structure_cost(building):
    building_name       = building[0]
    inner_width         = building[1]
    inner_length        = building[2]
    inner_height        = building[3]
    wall_thickness      = building[4]
    slab_roof_thickness = building[5]
    basemat_thickness   = building[6]

    outer_width  = inner_width + 2*wall_thickness
    outer_length = inner_length + 2*wall_thickness
    outer_height = inner_height + slab_roof_thickness + basemat_thickness

    # --- Unit costs : From RSMeans
    slab_roof_unit_cost = 52.0466087   # [$/cf]
    basemat_unit_cost   = 40.851       # [$/cf]
    walls_unit_cost     = 31.25025     # [$/cf]

    cf_to_m3            = 35.3147      # [cf/m3]



    # --- Calc
    slab_roof_volume = outer_width* outer_length* slab_roof_thickness
    basemat_volume   = outer_width* outer_length* basemat_thickness
    walls_volume     = 2* inner_width* inner_height* wall_thickness +\
                       2* inner_length* inner_height* wall_thickness

    slab_roof_cost = slab_roof_volume* slab_roof_unit_cost* cf_to_m3
    basemat_cost   = basemat_volume* basemat_unit_cost* cf_to_m3
    walls_cost     = walls_volume* walls_unit_cost* cf_to_m3
    building_cost  = slab_roof_cost+ basemat_cost+ walls_cost
    return building_cost


def calculate_reactor_building_structure_volume():
    """
      * Reactor building considers that the internal walls have the dimensions
        of an ISO container.
      * The wall thickness is assumed to be 2 [m]
    """
    building_name        = 'Reactor Building'
    inner_width          = 2.6  # [m]
    inner_length         = 2.6  # [m]
    inner_height         = 11   # [m]
    wall_thickness       = 2.0  # [m]
    slab_roof_thickness  = 2.0  # [m]
    basemat_thickness    = 2.0  # [m]

    reactor_building_dimensions = [building_name, inner_width, inner_length,
                                   inner_height, wall_thickness, slab_roof_thickness,
                                   basemat_thickness]

    reactor_building_structure_cost = calculate_building_structure_cost(reactor_building_dimensions)
    return reactor_building_structure_cost


def calculate_energy_conversion_building_structure_volume():
    """
      * Energy conversion building considers that the internal walls have dimensions 
        of an ISO container placed horizontally.
    """
    building_name        = 'Energy Conversion Building'
    inner_width          = 2.6  # [m]
    inner_length         = 6.0  # [m]
    inner_height         = 5.6  # [m]
    wall_thickness       = 2.0  # [m]
    slab_roof_thickness  = 2.0  # [m]
    basemat_thickness    = 2.0  # [m]
 
    energy_conversion_building_structure_dimensions = [building_name, inner_width, inner_length,
                                                       inner_height, wall_thickness, slab_roof_thickness,
                                                       basemat_thickness]
    energy_conversion_building_structure_cost = calculate_building_structure_cost(energy_conversion_building_structure_dimensions)
    return energy_conversion_building_structure_cost
