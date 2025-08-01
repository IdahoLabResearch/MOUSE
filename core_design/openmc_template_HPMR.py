# Importing libraries
import numpy as np
import openmc
from core_design.openmc_materials_database import collect_materials_data
from core_design.utils import create_universe_plot, circle_area, create_cells


# **************************************************************************************************************************
#                                                Sec. 0 : Helper Functions
# **************************************************************************************************************************

def create_pin_regions(params,pin_type):
    if pin_type == 'fuel':
        pin_radii = {'fuel_meat': params['Fuel Pin Radii'][0],
                     'gap': params['Fuel Pin Radii'][1]}
    
        region_keys = ['fuel_meat', 'gap', 'moderator']

    elif pin_type == 'heat pipe':
        pin_radii = {'htpipe': params['Heat Pipe Radii'][0],
                     'gap': params['Heat Pipe Radii'][1]}
        
        region_keys = ['htpipe', 'gap', 'moderator']

    else:
        raise ValueError("Invalid pin type. Must be 'fuel' or 'heat pipe'.")  

    # Creating surfaces
    # Create cylindrical surfaces for each of the specified radii
    shells = [openmc.ZCylinder(r=r) for r in pin_radii.values()]

    # Define the regions within the pin using the created cylindrical surfaces
    regions = {}
    for i, key in enumerate(region_keys[:-1]):
        if i == 0:
            regions[key] = -shells[i]
        else:
            regions[key] = +shells[i-1] & -shells[i]
    regions[region_keys[-1]] = +shells[-1]

    return regions  


def create_fuel_pin(params, materials_database):
    #materials_database = collect_materials_data(params)
    fuel_pin_regions = create_pin_regions(params, 'fuel')

    # Creating fuel pin materials
    fuel_materials = [materials_database[mat] for mat in params['Fuel Pin Materials']]

    # Append the moderator to the list of fuel pin materials
    fuel_materials.append(materials_database[params['Moderator']])
    # creating the fuel pin universe
    fuel_cells = create_cells(fuel_pin_regions, fuel_materials)
    fuel_pin_universe = openmc.Universe(cells=fuel_cells.values())

    return fuel_pin_universe, fuel_materials


def create_htpipe_pin(params, materials_database):
    htpipe_pin_regions = create_pin_regions(params, 'heat pipe')

    # Creating heat pipe pin materials
    htpipe_materials = [materials_database[mat] for mat in params['Heat Pipe Materials']]

    # Append the moderator to the list of heat pipe materials
    htpipe_materials.append(materials_database[params['Moderator']])
    # creating the heat pipe pin universe
    htpipe_cells = create_cells(htpipe_pin_regions, htpipe_materials)
    htpipe_universe = openmc.Universe(cells=htpipe_cells.values())

    return htpipe_universe, htpipe_materials


def create_assembly(params, fuel_pin_universe, htpipe_universe, materials_database):

    # Define the fuel assembly planes
    rods_pitch    = params['Lattice Pitch']
    ass_rings     = params['Number of Rings per Assembly'] 
    l2            = params['Assembly FTF'] / np.sqrt(3.0)
    l2            = l2 - 0.4/np.sqrt(3.0)   
    c2            = np.sqrt(3.)/3.
    x2            = 0.0
    y2            = 0.0
    r_right       = openmc.YPlane(surface_id=5,y0= y2 + np.sqrt(3.)/2.*l2,)
    r_left        = openmc.YPlane(surface_id=6,y0= y2 - np.sqrt(3.)/2.*l2,)
    r_upper_right = openmc.Plane(surface_id=7,  b=c2,  a=1., d= l2+y2*c2+x2, name='r_upper_right')  # y = -x/sqrt(3) + a
    r_upper_left  = openmc.Plane(surface_id=8, b=-c2,  a=1., d= l2-y2*c2+x2, name='r_upper_left ')  # y = x/sqrt(3) + a
    r_lower_right = openmc.Plane(surface_id=9, b=-c2,  a=1., d=-l2-y2*c2+x2, name='r_lower_right')  # y = x/sqrt(3) - a
    r_lower_left  = openmc.Plane(surface_id=10,  b=c2,  a=1., d=-l2+y2*c2+x2, name='r_lower_left ')  # y = -x/sqrt(3) - a
    l2            = l2 + 0.4/np.sqrt(3.0)
    c2            = np.sqrt(3.)/3.
    x2            = 0.0
    y2            = 0.0
    g_right       = openmc.YPlane(surface_id=11,y0= y2 + np.sqrt(3.)/2.*l2,)
    g_left        = openmc.YPlane(surface_id=12,y0= y2 - np.sqrt(3.)/2.*l2,)
    g_upper_right = openmc.Plane(surface_id=13,  b=c2,  a=1., d= l2+y2*c2+x2, name='g_upper_right')  # y = -x/sqrt(3) + a
    g_upper_left  = openmc.Plane(surface_id=14, b=-c2,  a=1., d= l2-y2*c2+x2, name='g_upper_left ')  # y = x/sqrt(3) + a
    g_lower_right = openmc.Plane(surface_id=15, b=-c2,  a=1., d=-l2-y2*c2+x2, name='g_lower_right')  # y = x/sqrt(3) - a
    g_lower_left  = openmc.Plane(surface_id=16,  b=c2,  a=1., d=-l2+y2*c2+x2, name='g_lower_left ')  # y = -x/sqrt(3) - a

    # Define the fuel assembly cells 
    assembly_reg_1 = openmc.Cell(cell_id=332, name='assembly_reg_1')
    assembly_gap_12 = openmc.Cell(cell_id=334, name='assembly_gap_12')

    # Define graphite cell
    graphite_cell = openmc.Cell(cell_id=348, name='out_univ')

    # Define the fuel assembly regions
    assembly_reg_1.region   = +r_left & -r_right & -r_upper_right & -r_upper_left & +r_lower_right & +r_lower_left 
    assembly_gap_12.region  = (-r_left | +r_right | +r_upper_right | +r_upper_left | -r_lower_right | -r_lower_left) & +g_left & -g_right & -g_upper_right & -g_upper_left & +g_lower_right & +g_lower_left 

    # Define the assembly gap filling and graphite cell filling
    assembly_gap_12.fill    = materials_database[params['Moderator']]
    graphite_cell.fill  = materials_database[params['Moderator']]
    # Define the graphite universe
    graphite_universe = openmc.Universe(cells=[graphite_cell])

    # Create a hexagonal lattice for the fuel assembly
    lattice_hex_1             = openmc.HexLattice(lattice_id=55)
    lattice_hex_1.center      = (0., 0.,)
    lattice_hex_1.pitch       = (rods_pitch,)
    lattice_hex_1.n_rings     = (ass_rings)
    lattice_hex_1.orientation = 'x'
    lattice_hex_1.outer = graphite_universe

    all_rings = []

    for ring_idx in range(ass_rings):
        if ring_idx == 0:
            num_positions = 1
        else:
            num_positions = ring_idx * 6

        ring = []
        if ring_idx % 2 == 0:
            # Even ring (alternate heat pipe and fuel rod)
            for i in range(num_positions):
                if i % 2 == 0:
                    ring.append(htpipe_universe)
                else:
                    ring.append(fuel_pin_universe)
        else:
            # Odd ring (all fuel rods)
            for i in range(num_positions):
                ring.append(fuel_pin_universe) 

        all_rings.append(ring)

    lattice_hex_1.universes = all_rings[::-1]
    
    # Filling of fuel assembly with the lattice
    assembly_reg_1.fill = lattice_hex_1

    # Create fuel assembly universe
    fuel_assembly = openmc.Universe(cells=[assembly_reg_1, assembly_gap_12])

    # Define graphite assembly cell
    grp_cc_cnt = openmc.Cell(cell_id=342, name='grp_cc_cnt')

    # Define graphite assembly region
    grp_cc_cnt.region = +g_left & -g_right & -g_upper_right & -g_upper_left & +g_lower_right & +g_lower_left 

    # Define the graphite assembly filling
    grp_cc_cnt.fill = materials_database[params['Moderator']]
    # Define the graphite assembly universe
    graphite_assembly = openmc.Universe(cells=[grp_cc_cnt])

    return fuel_assembly, graphite_assembly, graphite_universe


def create_hex_core_geometry(params, fuel_assembly, graphite_assembly, graphite_universe, materials_database):

    # Define the hex core planes
    assembly_pitch = params['Assembly FTF']
    no_of_core_rings = params['Number of Rings per Core'] + 1  #The outermost ring in the original input is filled with graphite universes which is extra ring to the no of rings in the core
    lc             = params['hexagonal Core Edge Length']
    cr             = np.sqrt(3.)/3.
    x2             = 0.0
    y2             = 0.0
    c_right        = openmc.XPlane(surface_id=31,x0= x2 + np.sqrt(3.)/2.*lc,)
    c_left         = openmc.XPlane(surface_id=32,x0= x2 - np.sqrt(3.)/2.*lc,)
    c_upper_right  = openmc.Plane(surface_id=33,  a=cr,  b=1., d= lc+x2*cr+y2, name='c_upper_right') 
    c_upper_left   = openmc.Plane(surface_id=34, a=-cr,  b=1., d= lc-x2*cr+y2, name='c_upper_left')
    c_lower_right  = openmc.Plane(surface_id=35, a=-cr,  b=1., d=-lc-x2*cr+y2, name='c_lower_right')
    c_lower_left   = openmc.Plane(surface_id=36,  a=cr,  b=1., d=-lc+x2*cr+y2, name='c_lower_left')

    lc             = params['hexagonal Core Edge Length'] + 0.06  #From the original input
    cr             = np.sqrt(3.)/3.
    x2             = 0.0
    y2             = 0.0
    c_right_1      = openmc.XPlane(surface_id=301,x0= x2 + np.sqrt(3.)/2.*lc,)
    c_left_1       = openmc.XPlane(surface_id=302,x0= x2 - np.sqrt(3.)/2.*lc,)
    c_upper_right_1 = openmc.Plane(surface_id=303,  a=cr,  b=1., d= lc+x2*cr+y2, name='c_upper_right_1') 
    c_upper_left_1  = openmc.Plane(surface_id=304, a=-cr,  b=1., d= lc-x2*cr+y2, name='c_upper_left_1')
    c_lower_right_1 = openmc.Plane(surface_id=305, a=-cr,  b=1., d=-lc-x2*cr+y2, name='c_lower_right_1')
    c_lower_left_1  = openmc.Plane(surface_id=306,  a=cr,  b=1., d=-lc+x2*cr+y2, name='c_lower_left_1')

    # Define the hex core cells 
    core_reg = openmc.Cell(cell_id=345, name='core_reg')
    core_reg_out = openmc.Cell(cell_id=346, name='core_reg_out')

    # Define the hex core regions
    core_reg.region = +c_left & -c_right & -c_upper_right & -c_upper_left & +c_lower_right & +c_lower_left 
    core_reg_out.region = (-c_left | +c_right | +c_upper_right | +c_upper_left | -c_lower_right | -c_lower_left) & +c_left_1 & -c_right_1 & -c_upper_right_1 & -c_upper_left_1 & +c_lower_right_1 & +c_lower_left_1      

    # Create a hexagonal lattice for the hex core
    core_hex             = openmc.HexLattice(lattice_id=65)
    core_hex.center      = (0., 0.,)
    core_hex.pitch       = (assembly_pitch,)
    core_hex.n_rings     = (no_of_core_rings)
    core_hex.orientation = 'y'
    
    core_rings =[]
    
    for ring_idx in range(no_of_core_rings):
        if ring_idx == 0:
            num_positions = 1
        else:
            num_positions = ring_idx * 6

        if ring_idx == 0:
            ring = [graphite_assembly]

        elif ring_idx == no_of_core_rings - 1:
            ring = [graphite_universe] * num_positions

        else:
            ring = [fuel_assembly] * num_positions

        core_rings.append(ring)

    core_hex.universes = core_rings[::-1]

    # Filling of hex core with the lattice
    core_reg.fill = core_hex

    # Filling of hex core gap
    core_reg_out.fill = materials_database[params['Gap']]             


    return core_reg, core_reg_out   

def create_control_drums(params, materials_database):

    # Define the control drums planes
    theta        = np.pi/180.0
    cr_gap_radius = params['Drum Radius']
    cr_out_radius = cr_gap_radius - 0.05 
    cr_in_radius = cr_out_radius - params['Drum Absorber Thickness'] 
    core_radius = params['Core Radius']

    lc             = params['hexagonal Core Edge Length'] + 0.06  #From the original input
    cr             = np.sqrt(3.)/3.
    x2             = 0.0
    y2             = 0.0
    c_right_1      = openmc.XPlane(surface_id=307,x0= x2 + np.sqrt(3.)/2.*lc,)
    c_left_1       = openmc.XPlane(surface_id=308,x0= x2 - np.sqrt(3.)/2.*lc,)
    c_upper_right_1 = openmc.Plane(surface_id=309,  a=cr,  b=1., d= lc+x2*cr+y2, name='c_upper_right_1') 
    c_upper_left_1  = openmc.Plane(surface_id=310, a=-cr,  b=1., d= lc-x2*cr+y2, name='c_upper_left_1')
    c_lower_right_1 = openmc.Plane(surface_id=311, a=-cr,  b=1., d=-lc-x2*cr+y2, name='c_lower_right_1')
    c_lower_left_1  = openmc.Plane(surface_id=312,  a=cr,  b=1., d=-lc+x2*cr+y2, name='c_lower_left_1')

    cr_top = openmc.Plane(surface_id=41, a=-np.tan(60 * theta),b=1.0,  name='cr_top')
    cr_bot = openmc.Plane(surface_id=42, a=-np.tan(-60 * theta), b=1.0,  name='cr_bot')
    cr_in  = openmc.ZCylinder(surface_id=43, x0=0.0, y0=0.0, r=cr_in_radius,  name='cr_in')
    cr_out = openmc.ZCylinder(surface_id=44, x0=0.0, y0=0.0, r=cr_out_radius,  name='cr_out')
    cr_gap = openmc.ZCylinder(surface_id=45, x0=0.0, y0=0.0, r=cr_gap_radius,  name='cr_gap')

    CR_000_180 = openmc.YPlane(surface_id=70, y0= 0.0, name='CR_000')
    CR_030_210 = openmc.Plane(surface_id=71,   a=-np.tan(30 * theta),  b=1.0, d= 0, name='CR_030')
    CR_060_240 = openmc.Plane(surface_id=72,   a=-np.tan(60 * theta),  b=1.0, d= 0, name='CR_060')
    CR_090_270 = openmc.XPlane(surface_id=73, x0= 0.0, name='CR_090')
    CR_120_300 = openmc.Plane(surface_id=74,  a=-np.tan(120 * theta),  b=1.0, d= 0, name='CR_120')
    CR_150_330 = openmc.Plane(surface_id=75,  a=-np.tan(150 * theta),  b=1.0, d= 0, name='CR_150')

    core_out   = openmc.ZCylinder(surface_id=82, x0=0.0, y0=0.0, r=core_radius,  name='core_out_2')
    core_out.boundary_type= 'vacuum'

    # Define the control drums cells 
    # The cells for the one control drum
    cr_drum = openmc.Cell(cell_id=21, name='cr_drum')
    cr_refl = openmc.Cell(cell_id=22, name='cr_refl')
    cr_gpp  = openmc.Cell(cell_id=23, name='cr_gap')
    cr_ass  = openmc.Cell(cell_id=24, name='cr_ass')

    # The cells for control drums before translation
    cr_000 = openmc.Cell(cell_id=40, name='cr_000')
    cr_030 = openmc.Cell(cell_id=41, name='cr_030')
    cr_060 = openmc.Cell(cell_id=42, name='cr_060')
    cr_090 = openmc.Cell(cell_id=43, name='cr_090')
    cr_120 = openmc.Cell(cell_id=44, name='cr_120')
    cr_150 = openmc.Cell(cell_id=45, name='cr_150')
    cr_180 = openmc.Cell(cell_id=46, name='cr_180')
    cr_210 = openmc.Cell(cell_id=47, name='cr_210')
    cr_240 = openmc.Cell(cell_id=48, name='cr_240')
    cr_270 = openmc.Cell(cell_id=49, name='cr_270')
    cr_300 = openmc.Cell(cell_id=50, name='cr_300')
    cr_330 = openmc.Cell(cell_id=51, name='cr_330')

    # The cells for control drums after translation
    cr_01 = openmc.Cell(cell_id=61, name='cr_01')
    cr_02 = openmc.Cell(cell_id=62, name='cr_02')
    cr_03 = openmc.Cell(cell_id=63, name='cr_03')
    cr_04 = openmc.Cell(cell_id=64, name='cr_04')
    cr_05 = openmc.Cell(cell_id=65, name='cr_05')
    cr_06 = openmc.Cell(cell_id=66, name='cr_06')
    cr_07 = openmc.Cell(cell_id=67, name='cr_07')
    cr_08 = openmc.Cell(cell_id=68, name='cr_08')
    cr_09 = openmc.Cell(cell_id=69, name='cr_09')
    cr_10 = openmc.Cell(cell_id=70, name='cr_10')
    cr_11 = openmc.Cell(cell_id=71, name='cr_11')
    cr_12 = openmc.Cell(cell_id=72, name='cr_12')

    # Define the one control drum regions
    cr_drum.region =  +cr_in & -cr_out & +cr_bot & -cr_top 
    cr_refl.region =  (-cr_in | -cr_bot | +cr_top) & -cr_out 
    cr_gpp.region  =  +cr_out & -cr_gap 
    cr_ass.region  =  +cr_gap 

    cr_drum.fill            = materials_database[params['Control Drum Absorber']]
    cr_refl.fill            = materials_database[params['Control Drum Reflector']]
    cr_gpp.fill             = materials_database[params['Gap']]
    cr_ass.fill             = materials_database[params['Reflector']]
    # Define a universe for the control drum 
    control_drum_uni = openmc.Universe(cells=[cr_drum,cr_refl,cr_gpp,cr_ass])

    # Fill all drums cells at different angles with the control drum universe
    cr_000.fill             = control_drum_uni
    cr_030.fill             = control_drum_uni
    cr_060.fill             = control_drum_uni
    cr_090.fill             = control_drum_uni
    cr_120.fill             = control_drum_uni
    cr_150.fill             = control_drum_uni
    cr_180.fill             = control_drum_uni
    cr_210.fill             = control_drum_uni
    cr_240.fill             = control_drum_uni
    cr_270.fill             = control_drum_uni
    cr_300.fill             = control_drum_uni
    cr_330.fill             = control_drum_uni

    # Adjust the control drums rotation
    rotation_angle          = 0
    cr_000.rotation         = [0,  0,    0 + rotation_angle]
    cr_330.rotation         = [0,  0,    0 + rotation_angle]
    cr_030.rotation         = [0,  0,   60 + rotation_angle]
    cr_060.rotation         = [0,  0,   60 + rotation_angle]
    cr_090.rotation         = [0,  0,  120 + rotation_angle]
    cr_120.rotation         = [0,  0,  120 + rotation_angle]
    cr_150.rotation         = [0,  0,  180 + rotation_angle]
    cr_180.rotation         = [0,  0,  180 + rotation_angle]
    cr_210.rotation         = [0,  0, -120 + rotation_angle]
    cr_240.rotation         = [0,  0, -120 + rotation_angle]
    cr_270.rotation         = [0,  0,  -60 + rotation_angle]
    cr_300.rotation         = [0,  0,  -60 + rotation_angle]

    # Translate the control drums
    params['Drum Tube Radius'] = params['Drum Radius'] + params['Drum Radius'] / 90 # cm
    cd_distance = ((params['Number of Rings per Core'] - 1) * params['Assembly FTF']) + (params['Assembly FTF'] / 2) + params['Drum Tube Radius']

    r_0                  = cd_distance            #calculated in the original input as (78 + 112)/2.0
    r_angle              = (0.0+30.0)/2
    cr_000.translation   = [r_0*np.cos((r_angle + 30*0)*theta),  r_0*np.sin((r_angle + 30*0)*theta),  0]
    cr_030.translation   = [r_0*np.cos((r_angle + 30*1)*theta),  r_0*np.sin((r_angle + 30*1)*theta),  0]
    cr_060.translation   = [r_0*np.cos((r_angle + 30*2)*theta),  r_0*np.sin((r_angle + 30*2)*theta),  0]
    cr_090.translation   = [r_0*np.cos((r_angle + 30*3)*theta),  r_0*np.sin((r_angle + 30*3)*theta),  0]
    cr_120.translation   = [r_0*np.cos((r_angle + 30*4)*theta),  r_0*np.sin((r_angle + 30*4)*theta),  0]
    cr_150.translation   = [r_0*np.cos((r_angle + 30*5)*theta),  r_0*np.sin((r_angle + 30*5)*theta),  0]
    cr_180.translation   = [r_0*np.cos((r_angle + 30*6)*theta),  r_0*np.sin((r_angle + 30*6)*theta),  0]
    cr_210.translation   = [r_0*np.cos((r_angle + 30*7)*theta),  r_0*np.sin((r_angle + 30*7)*theta),  0]
    cr_240.translation   = [r_0*np.cos((r_angle + 30*8)*theta),  r_0*np.sin((r_angle + 30*8)*theta),  0]
    cr_270.translation   = [r_0*np.cos((r_angle + 30*9)*theta),  r_0*np.sin((r_angle + 30*9)*theta),  0]
    cr_300.translation   = [r_0*np.cos((r_angle + 30*10)*theta),  r_0*np.sin((r_angle + 30*10)*theta),  0]
    cr_330.translation   = [r_0*np.cos((r_angle + 30*11)*theta),  r_0*np.sin((r_angle + 30*11)*theta),  0]

    # Create a universe for each control drum
    a1  = openmc.Universe(cells=[cr_000])
    a2  = openmc.Universe(cells=[cr_030])
    a3  = openmc.Universe(cells=[cr_060])
    a4  = openmc.Universe(cells=[cr_090])
    a5  = openmc.Universe(cells=[cr_120])
    a6  = openmc.Universe(cells=[cr_150])
    a7  = openmc.Universe(cells=[cr_180])
    a8  = openmc.Universe(cells=[cr_210])
    a9  = openmc.Universe(cells=[cr_240])
    a10 = openmc.Universe(cells=[cr_270])
    a11 = openmc.Universe(cells=[cr_300])
    a12 = openmc.Universe(cells=[cr_330])

    # Create the regions for all control drums after translation
    cr_01.region  = +CR_000_180 & -CR_030_210 & +c_right_1 & -core_out  
    cr_02.region  = +CR_030_210 & -CR_060_240 & +c_upper_right_1 & -core_out  
    cr_03.region  = +CR_060_240 & +CR_090_270 & +c_upper_right_1 & -core_out  
    cr_04.region  = -CR_090_270 & +CR_120_300 & +c_upper_left_1 & -core_out  
    cr_05.region  = -CR_120_300 & +CR_150_330 & +c_upper_left_1 & -core_out  
    cr_06.region  = -CR_150_330 & +CR_000_180 & -c_left_1 & -core_out  
    cr_07.region  = +CR_030_210 & -CR_000_180 & -c_left_1 & -core_out  
    cr_08.region  = +CR_060_240 & -CR_030_210 & -c_lower_left_1 & -core_out  
    cr_09.region  = -CR_090_270 & -CR_060_240 & -c_lower_left_1 & -core_out  
    cr_10.region  = +CR_090_270 & -CR_120_300 & -c_lower_right_1 & -core_out  
    cr_11.region  = +CR_120_300 & -CR_150_330 & -c_lower_right_1 & -core_out  
    cr_12.region  = +CR_150_330 & -CR_000_180 & +c_right_1 & -core_out  

    # Fill the control drums cells with the universes after translation
    cr_01.fill = a1
    cr_02.fill = a2
    cr_03.fill = a3
    cr_04.fill = a4
    cr_05.fill = a5
    cr_06.fill = a6
    cr_07.fill = a7
    cr_08.fill = a8
    cr_09.fill = a9
    cr_10.fill = a10
    cr_11.fill = a11
    cr_12.fill = a12

    return cr_01, cr_02, cr_03, cr_04, cr_05, cr_06, cr_07, cr_08, cr_09, cr_10, cr_11, cr_12 


def create_core_geometry(core_reg, core_reg_out, cr_01, cr_02, cr_03, cr_04, cr_05, cr_06, cr_07, cr_08, cr_09, cr_10, cr_11, cr_12):

   # Create a universe for the whole core
   core = openmc.Universe(cells=[core_reg, core_reg_out, cr_01, cr_02, cr_03, cr_04, cr_05, cr_06, cr_07, cr_08, cr_09, cr_10, cr_11, cr_12])
   core_geometry = openmc.Geometry(core)

   return core_geometry


# **************************************************************************************************************************
#                                                Sec. 1 : OpenMC Model
# **************************************************************************************************************************

"""
An OpenMC function that accepts an instance of "parameters" 
and generates the necessary XMl files
"""
def build_openmc_model_HPMR(params):
    # **************************************************************************************************************************
    #                                                Sec. 1.1 : MATERIALS
    # **************************************************************************************************************************

    # Create the materials 
    #materials = create_materials(params)
    materials_database = collect_materials_data(params)

    fuel = materials_database[params['Fuel']]
    coolant = materials_database[params['Coolant']]
    reflector = materials_database[params['Reflector']]
    moderator = materials_database[params['Moderator']]
    gap = materials_database[params['Gap']]
    control_drum_absorber = materials_database[params['Control Drum Absorber']]
    control_drum_reflector = materials_database[params['Control Drum Reflector']]
    print(f"fuel {fuel}")
    print(f"coolant {coolant}")
    print(f"reflector {reflector}")
    print(f"moderator {moderator}")
    print(f"control_drum_absorber {control_drum_absorber}")
    print(f"control_drum_reflector {control_drum_reflector}")
    # **************************************************************************************************************************
    #                                                Sec. 1.2 : GEOMETRY
    # **************************************************************************************************************************


    # Create the fuel pin
    fuel_pin_universe, fuel_materials = create_fuel_pin(params, materials_database)

    # Create the heat pipe
    htpipe_universe, htpipe_materials = create_htpipe_pin(params, materials_database)

    # Create the fuel assembly and graphite assembly
    fuel_assembly, graphite_assembly, graphite_universe = create_assembly(params, fuel_pin_universe, htpipe_universe, materials_database)

    # Create the hexagonal core geometry
    core_reg, core_reg_out = create_hex_core_geometry(params, fuel_assembly, graphite_assembly, graphite_universe, materials_database) 

    # Create the control drums
    cr_01, cr_02, cr_03, cr_04, cr_05, cr_06, cr_07, cr_08, cr_09, cr_10, cr_11, cr_12 = create_control_drums(params, materials_database)

    # Create the whole core geometry
    core_geometry = create_core_geometry(core_reg, core_reg_out, cr_01, cr_02, cr_03, cr_04, cr_05, cr_06, cr_07, cr_08, cr_09, cr_10, cr_11, cr_12)

    # Export the geometry to .xml file
    core_geometry.export_to_xml()

    # **************************************************************************************************************************
    #                                                Sec. 1.3 : PLOTTING
    # **************************************************************************************************************************
    # reading all the materials properties
    materials_database = collect_materials_data(params)
    
    if params['plotting'] == "Y":
        create_universe_plot(materials_database, core_geometry, 
                        plot_width = 2.01 * params['Core Radius'],
                        num_pixels = 2000, 
                        font_size = 16,
                        title = "Reactor Core", 
                        fig_size = 8, 
                        output_file_name = "core.png")

   # # **************************************************************************************************************************
   # #                                                Sec. 1.4 : VOLUME INFO for Depletion
   # # **************************************************************************************************************************
    #find where the fuel is in the fuel pin
    fuel_index = params['Fuel Pin Materials'].index(params['Fuel'])

    fissile_area = np.pi * 1 **2
    fuel.volume = fissile_area * round(params['Active Height'],0) * params['Fuel Pin Count']

   
    all_materials = fuel_materials +\
        htpipe_materials + [coolant, reflector, moderator, gap, control_drum_absorber, control_drum_reflector]
    
    # removing "None" materials
    all_materials_cleaned_list = [item for item in all_materials if item is not None]
    materials = openmc.Materials(list(set(all_materials_cleaned_list)))
   
    openmc.Materials.cross_sections = params['cross_sections_xml_location']
    materials.export_to_xml()
         #=================================================================================================
    #                                     tallies.xml File
    #=================================================================================================
    tallies_file = openmc.Tallies()

    group_edges = params[Energy Groups]  # Three energy groups
    groups = openmc.mgxs.EnergyGroups(group_edges)

    mgxs_lib = openmc.mgxs.Library(core_geometry)
    mgxs_lib.energy_groups = groups
    mgxs_lib.legendre_order     = 1
    mgxs_lib.mgxs_types = ['absorption', 'diffusion-coefficient', 'transport', 'scatter matrix', 'total', 'scatter']
    mgxs_lib.domain_type = 'universe'
    mgxs_lib.domains = [core]
    mgxs_lib.build_library()
    mgxs_lib.add_to_tallies_file(tallies_file, merge=False)
    tallies_file.export_to_xml()


    # # **************************************************************************************************************************
    # #                                                Sec. 1.5 : SIMULATION
    # # **************************************************************************************************************************
  
    settings = openmc.Settings()
    settings.batches = 100
    settings.inactive = 20
    settings.particles = 1000
    settings.temperature   = {'method': 'interpolation'}
    settings.export_to_xml()