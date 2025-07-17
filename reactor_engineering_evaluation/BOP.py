from reactor_engineering_evaluation.tools import *
from math import log

def calculate_heat_exchanger_mass(params):


    """
      Assuming a printed circuit heat exchanger (PCHE).
      Input for this function are as follows
        - hx_thermal_load : thermal load of PCHE in [MW]
        - th_in           : PCHE hot side inlet temprture [K]
        - th_out          : PCHE hot side outlet temperature [K]
        - tc_in           : PCHE cold side inlet temperature [K]
        - tc_out          : PCHE cold side outlet temperature [K]
    """
    hx_thermal_load = params['Power MWt']

    th_in  = params['Primary Loop Inlet Temperature']
    th_out = params['Primary Loop Outlet Temperature']
    tc_in  = params['Secondary Loop Inlet Temperature']
    tc_out = params['Secondary Loop Outlet Temperature']

    # Assumption on overall heat transfer coefficient
    U = 500  # [w/m2/K] This is an average value obtained by scanning of literature

    # PCHE channel dimensions assumptions
    hx_channel_diameter = 0.0015   # [m]
    hx_channel_pitch = 0.00225     # [m]
    hx_channel_length = 1          # [m]
    hx_void_fraction = 0.6
    hx_channel_thick = 0.003       # [m]

    rho_ss = 7850   #  density of stainless steel :: Kg/m^3

    hx_channel_perimeter = 3.14* hx_channel_diameter/2 + hx_channel_diameter
    hx_channel_ht_area = hx_channel_perimeter* hx_channel_length   # assuming a semi circular channel

    delta_t1 = abs(th_in - tc_out)
    delta_t2 = abs(th_out - tc_in)

    LMTD = (delta_t1 - delta_t2)/log(delta_t1/delta_t2)
    ht_area = hx_thermal_load*1e6/(U* LMTD)
    nchannels = ht_area/hx_channel_ht_area
    hx_alloy_volume = nchannels* hx_channel_pitch* hx_channel_thick - nchannels* 3.14/8* hx_channel_diameter**2 # m^3
    
    hx_mass = hx_alloy_volume* rho_ss
    return hx_mass


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

  

def calculate_building_structure_volumes(building):
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

    return slab_roof_volume, basemat_volume, walls_volume


  
def calculate_reactor_building_structure_volume(building_char):
    """
      * Reactor building considers that the internal walls have the dimensions
        of an ISO container.
      * The wall thickness is assumed to be 2 [m]
    """
    reactor_building_dimensions = building_char

    rb_slab_roof_vol, rb_basemat_vol, rb_walls_vol = calculate_building_structure_volumes(reactor_building_dimensions)
    return rb_slab_roof_vol, rb_basemat_vol, rb_walls_vol


  
def calculate_energy_conversion_building_structure_volume(building_char):
    """
      * Energy conversion building considers that the internal walls have dimensions 
        of an ISO container placed horizontally.
    """
    energy_conversion_building_structure_dimensions = building_char

    eb_slab_roof_vol, eb_basemat_vol, eb_walls_vol = calculate_building_structure_volumes(energy_conversion_building_structure_dimensions)
    return eb_slab_roof_vol, eb_basemat_vol, eb_walls_vol


def calculate_control_building_structure_volume(building_char):
    """
      * The control building dimensions are not completely based on assumptions.
      * The idea is that the control building is assumed to be occupied by two operators
        .. if needed.
    """
    control_building_dimensions = building_char

    cb_slab_roof_vol, cb_basemat_vol, cb_walls_vol = calculate_building_structure_volumes(control_building_dimensions)
    return cb_slab_roof_vol, cb_basemat_vol, cb_walls_vol


  
def calculate_refueling_building_strucutre_volume(building_char):
    """
      * Refueling building dimensions are entirely based on assumptions.
      * Some key notes prior to assumptions:
        - refueling building is larger than control building
        - Radioactive material needs to be handled so area needs to be large in
          order to account for existing shielding and equipment.
    """
    refueling_building_dimensions = building_char

    rb_slab_roof_vol, rb_basemat_vol, rb_walls_vol = calculate_building_structure_volumes(refueling_building_dimensions)
    return rb_slab_roof_vol, rb_basemat_vol, rb_walls_vol

  
def calculate_spent_fuel_building_structure_volume(building_char):
    """
      * It is expected that spent fuel building will have less equipment compared to refueling area
      * As a result, it has a smaller area
    """
    spent_fuel_building_dimensions = building_char

    sfb_slab_roof_vol, sfb_basemat_vol, sfb_walls_vol = calculate_building_structure_volumes(spent_fuel_building_dimensions)
    return sfb_slab_roof_vol, sfb_basemat_vol, sfb_walls_vol


def calculate_emergency_building_structure_volume(building_char):
    """
      * Dimensions are solely based on assumptions.
      * No supporting details for assumptions
    """
    emergency_building_dimensions = building_char

    eb_slab_roof_vol, eb_basemat_vol, eb_walls_vol = calculate_building_structure_volumes(emergency_building_dimensions)
    return eb_slab_roof_vol, eb_basemat_vol, eb_walls_vol

  
def calculate_storage_building_structure_volume(building_char):
    """
      * Dimensions are solely based on assumptinos
      * No supporting details for assumptions
    """
    storage_building_dimensions = building_char

    sb_slab_roof_vol, sb_basemat_vol, sb_walls_vol = calculate_building_structure_volumes(storage_building_dimensions)
    return sb_slab_roof_vol, sb_basemat_vol, sb_walls_vol

  
def calculate_radwaste_building_structure_volume(building_char):
    """
      * Dimensions are solely based on dimensions.
      * No elaboration on details of the assumptions
    """
    radwaste_storage_building_dimensions = building_char

    radb_slab_roof_vol, radb_basemat_vol, radb_walls_vol = calculate_building_structure_volumes(radwaste_storage_building_dimensions)
    return radb_slab_roof_vol, radb_basemat_vol, radb_walls_vol