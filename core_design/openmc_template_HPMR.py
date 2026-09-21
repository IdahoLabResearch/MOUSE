# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
# Importing libraries
import numpy as np
import openmc
from core_design.openmc_materials_database import collect_materials_data
from core_design.utils import create_universe_plot, circle_area, create_cells


def _load_depleted_fuel_override(params):
    """Load the operating depleted HPMR fuel for a lifecycle static case."""
    materials_xml = params.get('_Depleted Fuel Materials XML')
    if not materials_xml:
        return None

    if '_Operating Fuel Material ID' not in params:
        raise ValueError(
            "Static lifecycle calculation is missing the operating fuel "
            "material ID. Run the operating depletion before static cases."
        )

    operating_fuel_id = int(params['_Operating Fuel Material ID'])
    depleted_materials = openmc.Materials.from_xml(materials_xml)
    matches = [
        material
        for material in depleted_materials
        if material.id == operating_fuel_id
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one depleted fuel material with ID "
            f"{operating_fuel_id}, but found {len(matches)} in "
            f"{materials_xml}."
        )

    depleted_fuel = matches[0]
    depleted_fuel.temperature = params['Common Temperature']
    return depleted_fuel


# **************************************************************************************************************************
#                                                Sec. 0 : Helper Functions
# **************************************************************************************************************************

def create_pin_regions(params, pin_type):
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
            regions[key] = +shells[i - 1] & -shells[i]
    regions[region_keys[-1]] = +shells[-1]

    return regions


def create_fuel_pin(params, materials_database):
    fuel_pin_regions = create_pin_regions(params, 'fuel')

    # Creating fuel pin materials
    fuel_materials = [materials_database[mat] for mat in params['Fuel Pin Materials']]

    # Append the moderator to the list of fuel pin materials
    fuel_materials.append(materials_database[params['Moderator']])
    # Create the fuel pin universe
    fuel_cells = create_cells(fuel_pin_regions, fuel_materials)
    # Fuel region cell (used in distribcell tally)
    fuel_cell = fuel_cells['fuel_meat']
    fuel_pin_universe = openmc.Universe(cells=fuel_cells.values())

    return fuel_pin_universe, fuel_materials, fuel_cell


def create_htpipe_pin(params, materials_database):
    htpipe_pin_regions = create_pin_regions(params, 'heat pipe')

    # Creating heat pipe pin materials
    htpipe_materials = [materials_database[mat] for mat in params['Heat Pipe Materials']]

    # Append the moderator to the list of heat pipe materials
    htpipe_materials.append(materials_database[params['Moderator']])
    # Create the heat pipe pin universe
    htpipe_cells = create_cells(htpipe_pin_regions, htpipe_materials)
    htpipe_universe = openmc.Universe(cells=htpipe_cells.values())

    return htpipe_universe, htpipe_materials


def create_assembly(params, fuel_pin_universe, htpipe_universe, materials_database):

    # Define the fuel assembly planes
    rods_pitch = params['Lattice Pitch']
    ass_rings = params['Number of Rings per Assembly']
    l2 = params['Assembly FTF'] / np.sqrt(3.0)
    l2 = l2 - 0.4 / np.sqrt(3.0)
    c2 = np.sqrt(3.) / 3.
    x2 = 0.0
    y2 = 0.0
    r_right = openmc.YPlane(surface_id=5, y0=y2 + np.sqrt(3.) / 2. * l2)
    r_left = openmc.YPlane(surface_id=6, y0=y2 - np.sqrt(3.) / 2. * l2)
    r_upper_right = openmc.Plane(surface_id=7, b=c2, a=1., d=l2 + y2 * c2 + x2, name='r_upper_right')
    r_upper_left = openmc.Plane(surface_id=8, b=-c2, a=1., d=l2 - y2 * c2 + x2, name='r_upper_left ')
    r_lower_right = openmc.Plane(surface_id=9, b=-c2, a=1., d=-l2 - y2 * c2 + x2, name='r_lower_right')
    r_lower_left = openmc.Plane(surface_id=10, b=c2, a=1., d=-l2 + y2 * c2 + x2, name='r_lower_left ')
    l2 = l2 + 0.4 / np.sqrt(3.0)
    c2 = np.sqrt(3.) / 3.
    x2 = 0.0
    y2 = 0.0
    g_right = openmc.YPlane(surface_id=11, y0=y2 + np.sqrt(3.) / 2. * l2)
    g_left = openmc.YPlane(surface_id=12, y0=y2 - np.sqrt(3.) / 2. * l2)
    g_upper_right = openmc.Plane(surface_id=13, b=c2, a=1., d=l2 + y2 * c2 + x2, name='g_upper_right')
    g_upper_left = openmc.Plane(surface_id=14, b=-c2, a=1., d=l2 - y2 * c2 + x2, name='g_upper_left ')
    g_lower_right = openmc.Plane(surface_id=15, b=-c2, a=1., d=-l2 - y2 * c2 + x2, name='g_lower_right')
    g_lower_left = openmc.Plane(surface_id=16, b=c2, a=1., d=-l2 + y2 * c2 + x2, name='g_lower_left ')

    # Define the fuel assembly cells
    assembly_reg_1 = openmc.Cell(cell_id=332, name='assembly_reg_1')
    assembly_gap_12 = openmc.Cell(cell_id=334, name='assembly_gap_12')

    # Define graphite cell
    graphite_cell = openmc.Cell(cell_id=348, name='out_univ')

    # Define the fuel assembly regions
    assembly_reg_1.region = +r_left & -r_right & -r_upper_right & -r_upper_left & +r_lower_right & +r_lower_left
    assembly_gap_12.region = (-r_left | +r_right | +r_upper_right | +r_upper_left | -r_lower_right | -r_lower_left) & +g_left & -g_right & -g_upper_right & -g_upper_left & +g_lower_right & +g_lower_left

    # Define the assembly gap filling and graphite cell filling
    assembly_gap_12.fill = materials_database[params['Moderator']]
    graphite_cell.fill = materials_database[params['Moderator']]
    # Define the graphite universe
    graphite_universe = openmc.Universe(cells=[graphite_cell])

    # Create a hexagonal lattice for the fuel assembly
    lattice_hex_1 = openmc.HexLattice(lattice_id=55)
    lattice_hex_1.center = (0., 0.)
    lattice_hex_1.pitch = (rods_pitch,)
    lattice_hex_1.n_rings = ass_rings
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
    no_of_core_rings = params['Number of Rings per Core'] + 1
    lc = params['hexagonal Core Edge Length']
    cr = np.sqrt(3.) / 3.
    x2 = 0.0
    y2 = 0.0
    c_right = openmc.XPlane(surface_id=31, x0=x2 + np.sqrt(3.) / 2. * lc)
    c_left = openmc.XPlane(surface_id=32, x0=x2 - np.sqrt(3.) / 2. * lc)
    c_upper_right = openmc.Plane(surface_id=33, a=cr, b=1., d=lc + x2 * cr + y2, name='c_upper_right')
    c_upper_left = openmc.Plane(surface_id=34, a=-cr, b=1., d=lc - x2 * cr + y2, name='c_upper_left')
    c_lower_right = openmc.Plane(surface_id=35, a=-cr, b=1., d=-lc - x2 * cr + y2, name='c_lower_right')
    c_lower_left = openmc.Plane(surface_id=36, a=cr, b=1., d=-lc + x2 * cr + y2, name='c_lower_left')

    lc = params['hexagonal Core Edge Length'] + 0.06
    cr = np.sqrt(3.) / 3.
    x2 = 0.0
    y2 = 0.0
    c_right_1 = openmc.XPlane(surface_id=301, x0=x2 + np.sqrt(3.) / 2. * lc)
    c_left_1 = openmc.XPlane(surface_id=302, x0=x2 - np.sqrt(3.) / 2. * lc)
    c_upper_right_1 = openmc.Plane(surface_id=303, a=cr, b=1., d=lc + x2 * cr + y2, name='c_upper_right_1')
    c_upper_left_1 = openmc.Plane(surface_id=304, a=-cr, b=1., d=lc - x2 * cr + y2, name='c_upper_left_1')
    c_lower_right_1 = openmc.Plane(surface_id=305, a=-cr, b=1., d=-lc - x2 * cr + y2, name='c_lower_right_1')
    c_lower_left_1 = openmc.Plane(surface_id=306, a=cr, b=1., d=-lc + x2 * cr + y2, name='c_lower_left_1')

    # Define the hex core cells
    core_reg = openmc.Cell(cell_id=345, name='core_reg')
    core_reg_out = openmc.Cell(cell_id=346, name='core_reg_out')

    # Define the hex core regions
    core_reg.region = +c_left & -c_right & -c_upper_right & -c_upper_left & +c_lower_right & +c_lower_left
    core_reg_out.region = (-c_left | +c_right | +c_upper_right | +c_upper_left | -c_lower_right | -c_lower_left) & +c_left_1 & -c_right_1 & -c_upper_right_1 & -c_upper_left_1 & +c_lower_right_1 & +c_lower_left_1

    # Create a hexagonal lattice for the hex core
    core_hex = openmc.HexLattice(lattice_id=65)
    core_hex.center = (0., 0.)
    core_hex.pitch = (assembly_pitch,)
    core_hex.n_rings = no_of_core_rings
    core_hex.orientation = 'y'

    core_rings = []

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
    core_reg_out.fill = materials_database[params['Coolant']]

    return core_reg, core_reg_out


def _point_in_positive_halfspace(surface, x, y):
    coeffs = surface.coefficients
    if 'x0' in coeffs:   # XPlane: positive side is x >= x0
        return x >= coeffs['x0']
    if 'y0' in coeffs:   # YPlane: positive side is y >= y0
        return y >= coeffs['y0']
    a = coeffs.get('a', 0.0)
    b = coeffs.get('b', 0.0)
    d = coeffs.get('d', 0.0)
    return (a * x + b * y - d) >= 0.0


def _halfspace_toward_point(surface, x, y):
    return +surface if _point_in_positive_halfspace(surface, x, y) else -surface


def _make_radial_plane(angle_deg, surface_id, name):
    # Always use a general Plane with normal (-sin θ, cos θ) so that
    # _point_in_positive_halfspace works correctly for all angles,
    # including 0°, 90°, 180°, 270°, and 360°.
    theta = np.deg2rad(angle_deg)
    return openmc.Plane(surface_id=surface_id,
                        a=-np.sin(theta), b=np.cos(theta),
                        d=0.0, name=name)


def create_control_drums(params, materials_database):

    num_drums = int(params['Drum Count'])
    allowed_counts = [6, 12, 18, 24]
    if num_drums not in allowed_counts:
        raise ValueError(f"Unsupported Number of Drums = {num_drums}. Allowed values are {allowed_counts}.")

    # Define the control drums planes
    theta = np.pi / 180.0
    cr_gap_radius = params['Drum Radius']
    cr_out_radius = cr_gap_radius - 0.05
    cr_in_radius = cr_out_radius - params['Drum Absorber Thickness']
    core_radius = params['Core Radius']

    lc = params['hexagonal Core Edge Length'] + 0.06
    cr = np.sqrt(3.) / 3.
    x2 = 0.0
    y2 = 0.0
    c_right_1 = openmc.XPlane(surface_id=307, x0=x2 + np.sqrt(3.) / 2. * lc)
    c_left_1 = openmc.XPlane(surface_id=308, x0=x2 - np.sqrt(3.) / 2. * lc)
    c_upper_right_1 = openmc.Plane(surface_id=309, a=cr, b=1., d=lc + x2 * cr + y2, name='c_upper_right_1')
    c_upper_left_1 = openmc.Plane(surface_id=310, a=-cr, b=1., d=lc - x2 * cr + y2, name='c_upper_left_1')
    c_lower_right_1 = openmc.Plane(surface_id=311, a=-cr, b=1., d=-lc - x2 * cr + y2, name='c_lower_right_1')
    c_lower_left_1 = openmc.Plane(surface_id=312, a=cr, b=1., d=-lc + x2 * cr + y2, name='c_lower_left_1')

    cr_top = openmc.Plane(surface_id=41, a=-np.tan(60 * theta), b=1.0, name='cr_top')
    cr_bot = openmc.Plane(surface_id=42, a=-np.tan(-60 * theta), b=1.0, name='cr_bot')
    cr_in = openmc.ZCylinder(surface_id=43, x0=0.0, y0=0.0, r=cr_in_radius, name='cr_in')
    cr_out = openmc.ZCylinder(surface_id=44, x0=0.0, y0=0.0, r=cr_out_radius, name='cr_out')
    cr_gap = openmc.ZCylinder(surface_id=45, x0=0.0, y0=0.0, r=cr_gap_radius, name='cr_gap')

    core_out = openmc.ZCylinder(
        surface_id=82,
        x0=0.0,
        y0=0.0,
        r=core_radius,
        name='shielding_source_boundary',
    )
    core_out.boundary_type = 'vacuum'

    # Define the one control drum cells
    cr_drum = openmc.Cell(cell_id=21, name='cr_drum')
    cr_refl = openmc.Cell(cell_id=22, name='cr_refl')
    cr_gpp = openmc.Cell(cell_id=23, name='cr_gap')
    cr_ass = openmc.Cell(cell_id=24, name='cr_ass')

    # Define the one control drum regions
    cr_drum.region = +cr_in & -cr_out & +cr_bot & -cr_top
    cr_refl.region = (-cr_in | -cr_bot | +cr_top) & -cr_out
    cr_gpp.region = +cr_out & -cr_gap
    cr_ass.region = +cr_gap

    cr_drum.fill = materials_database[params['Control Drum Absorber']]
    cr_refl.fill = materials_database[params['Control Drum Reflector']]
    cr_gpp.fill = materials_database[params['Coolant']]
    cr_ass.fill = materials_database[params['Radial Reflector']]

    # Define a universe for the control drum
    control_drum_uni = openmc.Universe(cells=[cr_drum, cr_refl, cr_gpp, cr_ass])

    # Drum placement
    params['Drum Tube Radius'] = params['Drum Radius'] + params['Drum Radius'] / 90  # cm
    cd_distance = ((params['Number of Rings per Core'] - 1) * params['Assembly FTF']) + (params['Assembly FTF'] / 2) + params['Drum Tube Radius']

    drum_step_deg = 360.0 / num_drums
    half_step_deg = drum_step_deg / 2.0
    rotation_angle = 180 if params['Shutdown Margin Calc'] else 0

    # Full outside-hex condition (De Morgan of "inside hex").
    # Using all 6 faces handles every sector width (6, 12, 18, 24 drums)
    # without needing to map each sector to its nearest single hex face.
    outside_hex = (  +c_right_1 | -c_left_1
                   | +c_upper_right_1 | +c_upper_left_1
                   | -c_lower_right_1 | -c_lower_left_1)

    drum_sector_cells = []

    for i in range(num_drums):
        center_angle_deg = half_step_deg + i * drum_step_deg
        center_angle_rad = np.deg2rad(center_angle_deg)

        # Cell holding a translated/rotated control drum universe
        drum_instance = openmc.Cell(cell_id=1000 + i, name=f'cr_inst_{i + 1:02d}')
        drum_instance.fill = control_drum_uni

        # Snap absorber orientation to the nearest hex face normal (0°, 60°, …, 300°)
        base_rotation = (60.0 * np.floor((center_angle_deg + 15.0) / 60.0)) % 360.0
        drum_instance.rotation = [0, 0, base_rotation + rotation_angle]
        drum_instance.translation = [
            cd_distance * np.cos(center_angle_rad),
            cd_distance * np.sin(center_angle_rad),
            0.0
        ]

        drum_universe = openmc.Universe(cells=[drum_instance])

        # Sector angular boundaries — always general Planes so that
        # _halfspace_toward_point works correctly at every angle.
        start_surface = _make_radial_plane(i * drum_step_deg,       2000 + 2 * i,     f'CR_START_{i + 1:02d}')
        end_surface   = _make_radial_plane((i + 1) * drum_step_deg, 2000 + 2 * i + 1, f'CR_END_{i + 1:02d}')

        test_r = max(core_radius * 0.95, 0.1)
        test_x = test_r * np.cos(center_angle_rad)
        test_y = test_r * np.sin(center_angle_rad)

        start_halfspace = _halfspace_toward_point(start_surface, test_x, test_y)
        end_halfspace   = _halfspace_toward_point(end_surface,   test_x, test_y)

        sector_cell = openmc.Cell(cell_id=3000 + i, name=f'cr_{i + 1:02d}')
        sector_cell.region = start_halfspace & end_halfspace & outside_hex & -core_out
        sector_cell.fill = drum_universe

        drum_sector_cells.append(sector_cell)

    return drum_sector_cells


def create_core_geometry(core_reg, core_reg_out, drum_sector_cells):

    # Create a universe for the full core
    core = openmc.Universe(cells=[core_reg, core_reg_out] + drum_sector_cells)
    core_geometry = openmc.Geometry(core)

    return core_geometry, core


# **************************************************************************************************************************
#                                                Sec. 1 : OpenMC Model
# **************************************************************************************************************************

"""
An OpenMC function that accepts an instance of "parameters"
and generates the necessary XML files.
"""
def build_openmc_model_HPMR(params):
    params.setdefault('Cold Shutdown Temperature', 300)
    params.setdefault('Shutdown Margin Calc', False)
    params.setdefault('Isothermal Temperature Coefficients', False)

    if int(params['Drum Count']) not in [6, 12, 18, 24]:
        raise ValueError("Number of Drums must be one of [6, 12, 18, 24].")

    # **************************************************************************************************************************
    #                                                Sec. 1.1 : MATERIALS
    # **************************************************************************************************************************

    # Create the materials
    materials_database = collect_materials_data(params)

    fuel = materials_database[params['Fuel']]
    depleted_fuel_override = _load_depleted_fuel_override(params)
    if depleted_fuel_override is None:
        params['_Operating Fuel Material ID'] = int(fuel.id)
    else:
        fuel = depleted_fuel_override
        materials_database[params['Fuel']] = fuel
    coolant = materials_database[params['Cooling Device']]
    reflector = materials_database[params['Radial Reflector']]
    moderator = materials_database[params['Moderator']]
    gap = materials_database[params['Coolant']]
    control_drum_absorber = materials_database[params['Control Drum Absorber']]
    control_drum_reflector = materials_database[params['Control Drum Reflector']]

    # **************************************************************************************************************************
    #                                                Sec. 1.2 : GEOMETRY
    # **************************************************************************************************************************

    # Create the fuel pin
    fuel_pin_universe, fuel_materials, fuel_cell = create_fuel_pin(params, materials_database)

    # Create the heat pipe
    htpipe_universe, htpipe_materials = create_htpipe_pin(params, materials_database)

    # Create the fuel assembly and graphite assembly
    fuel_assembly, graphite_assembly, graphite_universe = create_assembly(
        params, fuel_pin_universe, htpipe_universe, materials_database
    )

    # Create the hexagonal core geometry
    core_reg, core_reg_out = create_hex_core_geometry(
        params, fuel_assembly, graphite_assembly, graphite_universe, materials_database
    )

    # Create the control drums
    drum_sector_cells = create_control_drums(params, materials_database)

    # Create the whole core geometry
    core_geometry, core = create_core_geometry(core_reg, core_reg_out, drum_sector_cells)

    # Export the geometry to .xml file
    core_geometry.export_to_xml()

    # **************************************************************************************************************************
    #                                                Sec. 1.3 : PLOTTING
    # **************************************************************************************************************************

    if params['plotting'] == "Y":
        create_universe_plot(materials_database, fuel_pin_universe,
                             plot_width=2.2 * params['Fuel Pin Radii'][-1],
                             num_pixels=500,
                             font_size=32,
                             title="Fuel Pin Universe",
                             fig_size=8,
                             output_file_name="fuel_pin_universe.png")

        create_universe_plot(materials_database, htpipe_universe,
                             plot_width=2.2 * params['Heat Pipe Radii'][-1],
                             num_pixels=500,
                             font_size=32,
                             title="Heat Pipe Universe",
                             fig_size=8,
                             output_file_name="heatpipe_universe.png")

        create_universe_plot(materials_database, fuel_assembly,
                             plot_width=1.3 * params['Assembly FTF'],
                             num_pixels=500,
                             font_size=32,
                             title="Fuel Assembly",
                             fig_size=8,
                             output_file_name="fuel_assembly.png")

        create_universe_plot(materials_database, core_geometry,
                             plot_width=2.01 * params['Core Radius'],
                             num_pixels=2000,
                             font_size=32,
                             title="Reactor Core",
                             fig_size=8,
                             output_file_name="core.png")

    # **************************************************************************************************************************
    #                                                Sec. 1.4 : VOLUME INFO for Depletion
    # **************************************************************************************************************************
    fissile_area = np.pi * params['Fuel Pin Radii'][0] ** 2
    fuel.volume = fissile_area * params['Active Height'] * params['Fuel Pin Count']

    all_materials = fuel_materials + htpipe_materials + [
        coolant, reflector, moderator, gap, control_drum_absorber, control_drum_reflector
    ]

    # Remove None materials
    all_materials_cleaned_list = [item for item in all_materials if item is not None]
    materials = openmc.Materials(list(set(all_materials_cleaned_list)))

    openmc.Materials.cross_sections = params['cross_sections_xml_location']
    materials.export_to_xml()

    # =================================================================================================
    #                                     tallies.xml file
    # =================================================================================================
    tallies_file = openmc.Tallies()

    group_edges = np.array([1e-5, 6.7e-2, 3.2e-1, 1, 4, 9.88, 4.81e1, 4.54e2, 4.9e4, 1.83e5, 8.21e5, 4e7])
    groups = openmc.mgxs.EnergyGroups(group_edges)

    mgxs_lib = openmc.mgxs.Library(core_geometry)
    mgxs_lib.energy_groups = groups
    mgxs_lib.legendre_order = 1
    mgxs_lib.mgxs_types = ['absorption', 'diffusion-coefficient', 'transport', 'scatter matrix', 'total', 'scatter']
    mgxs_lib.domain_type = 'universe'
    mgxs_lib.domains = [core]
    mgxs_lib.build_library()
    mgxs_lib.add_to_tallies_file(tallies_file, merge=False)

    shielding_source_boundary = next(
        surface for surface in core_geometry.get_all_surfaces().values()
        if surface.name == 'shielding_source_boundary'
    )
    shielding_leakage = openmc.Tally(name='boc_shielding_leakage_current')
    shielding_leakage.filters = [
        openmc.SurfaceFilter(shielding_source_boundary),
        openmc.EnergyFilter(group_edges),
    ]
    shielding_leakage.scores = ['current']
    tallies_file.append(shielding_leakage)

    shielding_kappa = openmc.Tally(name='boc_total_kappa_fission')
    shielding_kappa.scores = ['kappa-fission']
    tallies_file.append(shielding_kappa)

    # Peaking factor tally (pin power)
    pin_filter = openmc.DistribcellFilter(fuel_cell)
    pin_power = openmc.Tally(name='pin_power_kappa')
    pin_power.scores = ['kappa-fission']
    pin_power.filters = [pin_filter]
    tallies_file.append(pin_power)
    tallies_file.export_to_xml()

    # **************************************************************************************************************************
    #                                                Sec. 1.5 : SIMULATION
    # **************************************************************************************************************************

    settings = openmc.Settings()
    settings.batches = 100
    settings.inactive = 20

    if 'Particles' in params.keys():
        settings.particles = int(params['Particles'])
    else:
        settings.particles = 1000

    settings.temperature = {
        'default': params['Common Temperature'],
        'method': 'interpolation',
        'tolerance': 50.0
    }

    settings.export_to_xml()
