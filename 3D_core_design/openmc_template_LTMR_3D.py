# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
# Importing libraries
import openmc
import numpy as np
from core_design.openmc_materials_database import collect_materials_data
from core_design.utils import (
    create_universe_plot,
    circle_area,
    create_cells,
    calculate_hex_edge_length,
    calculate_hex_apothem,
)
import copy


# **************************************************************************************************************************
#                                                Sec. 0 : Helper Functions
# **************************************************************************************************************************

"""
Helper functions are smaller, reusable functions defined to perform specific tasks,
and then used later to simplify and organize the code.
"""


def create_pin_regions(params, pin_type):
    """
    Creating the pin regions
    @ In, params, dict, The parameters that are used to "fill in" input files with placeholders.
    @ In, pin_type, str, The type of pin ('moderator' or 'fuel').
    @ out, regions, dict, Regions of the specified pin.
    """

    if pin_type == 'moderator':
        pin_radii = {
            'moderator': params['Moderator Pin Radii'][0],
            'cladding':  params['Moderator Pin Radii'][1]
        }
        region_keys = ['moderator', 'cladding', 'coolant']

    elif pin_type == 'fuel':
        pin_radii = {
            'insert':    params['Fuel Pin Radii'][0],
            'gap1':      params['Fuel Pin Radii'][1],
            'fuel_meat': params['Fuel Pin Radii'][2],
            'gap2':      params['Fuel Pin Radii'][3],
            'cladding':  params['Fuel Pin Radii'][4]
        }
        region_keys = ['insert', 'gap1', 'fuel_meat', 'gap2', 'cladding', 'coolant']

    else:
        raise ValueError("Invalid pin type. Must be 'moderator' or 'fuel'.")

    shells = [openmc.ZCylinder(r=r) for r in pin_radii.values()]

    regions = {}
    for i, key in enumerate(region_keys[:-1]):
        if i == 0:
            regions[key] = -shells[i]
        else:
            regions[key] = +shells[i - 1] & -shells[i]
    regions[region_keys[-1]] = +shells[-1]

    return regions


# ──────────────────────────────────────────────────────────────────────────────
#  Drum layout helpers  (ported from openmc_template_LTMR.py)
# ──────────────────────────────────────────────────────────────────────────────

def _get_valid_drum_counts():
    return [6, 12, 18, 24, 30, 36]


def _get_drum_layout_quantities(params, drum_radius):
    number_of_drums = params['Number of Drums']
    valid_drum_counts = _get_valid_drum_counts()
    if number_of_drums not in valid_drum_counts:
        raise ValueError(f"Number of Drums must be one of {valid_drum_counts}, got {number_of_drums}")

    drums_per_side   = number_of_drums // 6
    hex_edge_length  = calculate_hex_edge_length(params)
    apothem          = calculate_hex_apothem(params)
    drum_tube_radius = drum_radius + drum_radius / 90.0
    side_length      = hex_edge_length

    return drums_per_side, hex_edge_length, apothem, drum_tube_radius, side_length


def _drum_positions_for_radius(params, drum_radius):
    drums_per_side, _, apothem, drum_tube_radius, side_length = \
        _get_drum_layout_quantities(params, drum_radius)

    face_angles = [k * np.pi / 3 for k in range(6)]
    positions = []

    for face_angle in face_angles:
        along_x = -np.sin(face_angle)
        along_y =  np.cos(face_angle)

        radial_distance = apothem + drum_tube_radius
        face_center_x   = radial_distance * np.cos(face_angle)
        face_center_y   = radial_distance * np.sin(face_angle)

        for i in range(drums_per_side):
            offset = side_length * (i - (drums_per_side - 1) / 2.0) / drums_per_side
            x = face_center_x + offset * along_x
            y = face_center_y + offset * along_y
            positions.append((x, y, np.degrees(face_angle)))

    return positions, drum_tube_radius, side_length


def _drum_radius_is_feasible(params, drum_radius):
    positions, drum_tube_radius, side_length = _drum_positions_for_radius(params, drum_radius)
    drums_per_side = params['Number of Drums'] // 6

    # Same-face spacing check
    same_face_spacing = side_length / drums_per_side
    if 2.0 * drum_tube_radius > same_face_spacing:
        return False

    # Full pairwise overlap check
    min_center_dist = 2.0 * drum_tube_radius
    for i in range(len(positions)):
        x1, y1, _ = positions[i]
        for j in range(i + 1, len(positions)):
            x2, y2, _ = positions[j]
            dist = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
            if dist < min_center_dist:
                return False

    return True


def calculate_max_drum_radius(params, tol=1e-6, max_iter=100):
    """
    Compute the maximum feasible drum radius (cm) that avoids overlap for the
    current lattice/drum configuration.
    """
    number_of_drums = params['Number of Drums']
    valid_drum_counts = _get_valid_drum_counts()
    if number_of_drums not in valid_drum_counts:
        raise ValueError(f"Number of Drums must be one of {valid_drum_counts}, got {number_of_drums}")

    drums_per_side  = number_of_drums // 6
    hex_edge_length = calculate_hex_edge_length(params)
    side_length     = hex_edge_length

    upper_bound = (side_length / (2.0 * drums_per_side)) * 90.0 / 91.0
    lower_bound = 0.0

    for _ in range(max_iter):
        mid = 0.5 * (lower_bound + upper_bound)
        if _drum_radius_is_feasible(params, mid):
            lower_bound = mid
        else:
            upper_bound = mid
        if upper_bound - lower_bound < tol:
            break

    return lower_bound


def resolve_drum_radius(params):
    """
    If Drum Radius is not provided, set it to the maximum feasible value (cm).
    The resolved numeric value overwrites params['Drum Radius'] so all
    downstream calculations use a number.
    """
    if 'Drum Radius' not in params:
        params['Drum Radius'] = calculate_max_drum_radius(params)

    drum_radius = params['Drum Radius']
    if not isinstance(drum_radius, (int, float, np.floating)):
        raise ValueError(
            f"Drum Radius must be numeric if provided, got {drum_radius!r}"
        )

    params['Drum Radius'] = float(drum_radius)
    return params['Drum Radius']


def create_drums_universe(params, control_drum_absorber_material,
                          control_drum_reflector_material, drum_positions):
    """
    Creating the universe of control drums.
    @ In, params, dict
    @ In, control_drum_absorber_material, openmc.Material
    @ In, control_drum_reflector_material, openmc.Material
    @ In, drum_positions, list of (x, y, face_angle_deg) tuples
    @ out, list of openmc.Universe, one per drum position
    """
    number_of_drums = params['Number of Drums']
    valid_drum_counts = _get_valid_drum_counts()
    if number_of_drums not in valid_drum_counts:
        raise ValueError(f"Number of Drums must be one of {valid_drum_counts}, got {number_of_drums}")

    absorber_thickness = params['Drum Absorber Thickness']
    drum_radius        = resolve_drum_radius(params)
    absorber_arc       = np.deg2rad(params['Drum Absorber Arc Degrees'])

    # Shutdown Margin Calc = True  → ARI orientation (absorber facing core)
    # Shutdown Margin Calc = False → ARO orientation (absorber facing away)
    drum_state     = 'shutdown' if params['Shutdown Margin Calc'] else 'operation'
    rotation_angle = 0 if drum_state == 'shutdown' else 180

    cd_inner_shell  = openmc.ZCylinder(r=drum_radius - absorber_thickness)
    cd_outer_shell  = openmc.ZCylinder(r=drum_radius)
    cutting_plane_1 = openmc.Plane(a=np.sin(absorber_arc / 2),  b= np.cos(absorber_arc / 2))
    cutting_plane_2 = openmc.Plane(a=np.sin(absorber_arc / 2),  b=-np.cos(absorber_arc / 2))

    drum_absorber_region  = +cd_inner_shell & -cd_outer_shell & -cutting_plane_1 & -cutting_plane_2
    drum_reflector_region = -cd_outer_shell & ~drum_absorber_region
    drum_outside_region   = +cd_outer_shell

    drum_absorber  = openmc.Cell(name='drum_absorber',  fill=control_drum_absorber_material,  region=drum_absorber_region)
    drum_reflector = openmc.Cell(name='drum_reflector', fill=control_drum_reflector_material, region=drum_reflector_region)
    drum_exterior  = openmc.Cell(name='drum_outside',   region=drum_outside_region)

    drum_reference = openmc.Universe(cells=(drum_reflector, drum_absorber, drum_exterior))

    drum_cells = []
    for i, (_, _, face_angle_deg) in enumerate(drum_positions):
        dc = openmc.Cell(name=f'drum_{i}', fill=drum_reference)
        dc.rotation = [0, 0, face_angle_deg + rotation_angle]
        drum_cells.append(dc)

    drums = [openmc.Universe(cells=(dc,)) for dc in drum_cells]
    return drums


def create_assembly_universe(params, fuel_pin_universe, moderator_pin_universe,
                             pin_pitch, reflector_material, outer_coolant_universe):
    """
    Creating the universe of the fuel assembly.
    """
    assembly        = openmc.HexLattice()
    assembly.center = (0., 0.)
    assembly.pitch  = (pin_pitch,)
    assembly.outer  = outer_coolant_universe

    rings = copy.deepcopy(params['Pins Arrangement'])
    rings = rings[-params['Number of Rings per Assembly']:]

    for i in range(len(rings)):
        for j in range(len(rings[i])):
            if rings[i][j] == 'FUEL':
                rings[i][j] = fuel_pin_universe
            elif rings[i][j] == 'MODERATOR':
                rings[i][j] = moderator_pin_universe

    assembly.universes = rings

    # Use the same hex edge length calculation as the 2D template
    hex_edge_length  = calculate_hex_edge_length(params)
    assembly_boundary = openmc.model.hexagonal_prism(
        edge_length=hex_edge_length,
        corner_radius=params['Fuel Pin Radii'][-1] + params["Pin Gap Distance"]
    )

    fuel_assembly_cell = openmc.Cell(fill=assembly,           region= assembly_boundary)
    reflector_cell     = openmc.Cell(fill=reflector_material, region=~assembly_boundary)

    assembly_universe = openmc.Universe(cells=[fuel_assembly_cell, reflector_cell])
    return assembly_universe


def create_control_drums_positions(params):
    """
    Place N/6 drums along each of the 6 flat faces of the hexagonal lattice,
    touching each face from outside, evenly spaced along the face length.
    Returns list of (x, y, face_angle_deg) tuples.
    """
    number_of_drums   = params['Number of Drums']
    valid_drum_counts = _get_valid_drum_counts()
    if number_of_drums not in valid_drum_counts:
        raise ValueError(
            f"Number of Drums must be one of {valid_drum_counts}, got {number_of_drums}"
        )

    drums_per_side  = number_of_drums // 6
    hex_edge_length = calculate_hex_edge_length(params)
    apothem         = calculate_hex_apothem(params)
    drum_radius     = resolve_drum_radius(params)
    drum_tube_radius = drum_radius + drum_radius / 90.0
    side_length     = hex_edge_length

    # Same-face neighbor spacing check
    same_face_spacing = side_length / drums_per_side
    if 2.0 * drum_tube_radius > same_face_spacing:
        max_drum_radius = (same_face_spacing / 2.0) * 90.0 / 91.0
        raise ValueError(
            f"Drums on the same hex face will overlap. "
            f"For Number of Drums = {number_of_drums}, the maximum Drum Radius is about "
            f"{max_drum_radius:.3f} cm, but got {drum_radius:.3f} cm."
        )

    face_angles = [k * np.pi / 3 for k in range(6)]
    positions   = []

    for face_angle in face_angles:
        along_x = -np.sin(face_angle)
        along_y =  np.cos(face_angle)

        radial_distance = apothem + drum_tube_radius
        face_center_x   = radial_distance * np.cos(face_angle)
        face_center_y   = radial_distance * np.sin(face_angle)

        for i in range(drums_per_side):
            offset = side_length * (i - (drums_per_side - 1) / 2.0) / drums_per_side
            x = face_center_x + offset * along_x
            y = face_center_y + offset * along_y
            positions.append((x, y, np.degrees(face_angle)))

    # Full pairwise overlap check
    min_center_dist = 2.0 * drum_tube_radius
    for i in range(len(positions)):
        x1, y1, _ = positions[i]
        for j in range(i + 1, len(positions)):
            x2, y2, _ = positions[j]
            dist = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
            if dist < min_center_dist:
                overlap = min_center_dist - dist
                raise ValueError(
                    f"Drum placement overlap detected between drums {i} and {j}. "
                    f"Overlap = {overlap:.3f} cm. "
                    f"Try reducing Drum Radius or reducing Number of Drums."
                )

    return positions


def update_ltmr_reflector_geometry_from_drums(params):
    """
    Derive Core Radius, Radial Reflector Thickness, Axial Reflector Thickness,
    and Drum Height from the actual drum layout implied by Number of Drums and
    Drum Radius. Mirrors the same function in openmc_template_LTMR.py.
    """
    drum_radius      = resolve_drum_radius(params)
    drum_tube_radius = drum_radius + drum_radius / 90.0

    drum_positions = create_control_drums_positions(params)

    max_outer_radius = max(
        np.sqrt(x ** 2 + y ** 2) + drum_tube_radius
        for x, y, _ in drum_positions
    )

    hex_apothem = calculate_hex_apothem(params)

    params['Core Radius']               = max_outer_radius
    params['Radial Reflector Thickness'] = params['Core Radius'] - hex_apothem
    params['Axial Reflector Thickness']  = params['Radial Reflector Thickness']
    params['Drum Height']               = params['Active Height'] + 2 * params['Axial Reflector Thickness']

    return drum_positions


def create_core_geometry(params, drums, drums_positions, assembly_universes, reflector_material):

    N_AXIAL     = len(assembly_universes)
    half_h      = params['Active Height'] / 2.0
    axial_thick = params['Axial Reflector Thickness']
    zone_h      = params['Active Height'] / N_AXIAL

    z_planes      = [openmc.ZPlane(z0=-half_h + i * zone_h) for i in range(N_AXIAL + 1)]
    z_bottom_refl = openmc.ZPlane(z0=-(half_h + axial_thick), boundary_type='vacuum')
    z_top_refl    = openmc.ZPlane(z0= (half_h + axial_thick), boundary_type='vacuum')

    drum_radius      = resolve_drum_radius(params)
    drum_tube_radius = drum_radius + drum_radius / 90.0
    params['Drum Tube Radius'] = drum_tube_radius

    axial_bounds = +z_planes[0] & -z_planes[-1]

    drum_shells = []
    drum_cells  = []
    for (x, y, _), du in zip(drums_positions, drums):
        drum_shell = openmc.ZCylinder(x0=x, y0=y, r=drum_tube_radius)
        drum_shells.append(drum_shell)
        drum_cell = openmc.Cell(fill=du, region=-drum_shell & axial_bounds)
        drum_cell.translation = (x, y, 0)
        drum_cells.append(drum_cell)

    drums_outside = +drum_shells[0]
    for d in drum_shells[1:]:
        drums_outside = drums_outside & +d

    outer_surface = openmc.ZCylinder(r=params['Core Radius'], boundary_type='vacuum')

    core_cells = []
    for i, au in enumerate(assembly_universes):
        zone_region = -outer_surface & drums_outside & +z_planes[i] & -z_planes[i + 1]
        core_cells.append(openmc.Cell(fill=au, region=zone_region))

    top_refl_cell = openmc.Cell(name='axial_reflector_top',
                                fill=reflector_material,
                                region=-outer_surface & +z_planes[-1] & -z_top_refl)
    bot_refl_cell = openmc.Cell(name='axial_reflector_bot',
                                fill=reflector_material,
                                region=-outer_surface & +z_bottom_refl & -z_planes[0])

    core          = openmc.Universe(cells=core_cells + [top_refl_cell, bot_refl_cell] + drum_cells)
    core_geometry = openmc.Geometry(core)
    return core_geometry, core


# **************************************************************************************************************************
#                                                Sec. 1 : OpenMC Model
# **************************************************************************************************************************

def build_openmc_model_LTMR_3D(params):
    """
    OpenMC Model — 3D with 1/4 rotational symmetry and 15 axial fuel zones.
    @ In, params, watts.parameters.Parameters
    """

    params.setdefault('Shutdown Margin Calc', False)
    params.setdefault('Isothermal Temperature Coefficients', False)

    # Ensure Drum Radius is always numeric before any downstream drum geometry use
    resolve_drum_radius(params)

    # **************************************************************************************************************************
    #                                                Sec. 1.1 : MATERIALS
    # **************************************************************************************************************************

    materials_database    = collect_materials_data(params)
    fuel                  = materials_database[params['Fuel']]
    coolant               = materials_database[params['Coolant']]
    reflector             = materials_database[params['Radial Reflector']]
    control_drum_absorber = materials_database[params['Control Drum Absorber']]
    control_drum_reflector= materials_database[params['Control Drum Reflector']]

    # **************************************************************************************************************************
    #                                                Sec. 1.2 : GEOMETRY: Fuel Pins, Moderator Pins, Coolant
    # **************************************************************************************************************************

    fuel_pin_regions = create_pin_regions(params, 'fuel')

    N_AXIAL = 15

    # One fuel material clone per axial zone
    fuel_materials_axial = [fuel.clone() for _ in range(N_AXIAL)]
    for i, fm in enumerate(fuel_materials_axial):
        fm.name = f'{params["Fuel"]}_axial_{i}'

    # One fuel pin universe per axial zone
    fuel_pin_universes_axial = []
    fuel_cells_axial         = []
    for i in range(N_AXIAL):
        zone_mats = []
        for mat in params['Fuel Pin Materials']:
            if mat is None:
                zone_mats.append(None)
            elif mat == params['Fuel']:
                zone_mats.append(fuel_materials_axial[i])
            else:
                zone_mats.append(materials_database[mat])
        zone_mats.append(coolant)
        zone_cells = create_cells(fuel_pin_regions, zone_mats)
        fuel_cells_axial.append(zone_cells['fuel_meat'])
        fuel_pin_universes_axial.append(openmc.Universe(cells=zone_cells.values()))

    fuel_cell = fuel_cells_axial[0]  # representative cell for distribcell tally

    if params['plotting'] == "Y":
        create_universe_plot(
            materials_database, fuel_pin_universes_axial[0],
            plot_width=2.2 * params['Fuel Pin Radii'][-1],
            num_pixels=500, font_size=32,
            title="Fuel Pin Universe", fig_size=8,
            output_file_name="fuel_pin_universe.png"
        )

    # Moderator pins
    moderator_pin_regions = create_pin_regions(params, 'moderator')
    moderator_materials   = []
    for mat in params['Moderator Pin Materials']:
        if mat is None:
            moderator_materials.append(None)
        else:
            moderator_materials.append(materials_database[mat])
    moderator_materials.append(coolant)

    if len(moderator_pin_regions) != len(moderator_materials):
        raise ValueError(
            f"The number of moderator pin regions ({len(moderator_pin_regions)}) must match "
            f"the number of introduced materials ({len(moderator_materials)})."
        )

    moderator_cells       = create_cells(moderator_pin_regions, moderator_materials)
    moderator_pin_universe = openmc.Universe(cells=moderator_cells.values())

    if params['plotting'] == "Y":
        create_universe_plot(
            materials_database, moderator_pin_universe,
            plot_width=2.2 * params['Moderator Pin Radii'][-1],
            num_pixels=500, font_size=32,
            title="Moderator Pin Universe", fig_size=8,
            output_file_name="moderator_pin_universe.png"
        )

    # Coolant universe
    coolant_cell    = openmc.Cell(fill=coolant)
    coolant_universe = openmc.Universe(cells=(coolant_cell,))

    # **************************************************************************************************************************
    #                                                Sec. 1.3 : Fuel Assembly Universes (one per axial zone)
    # **************************************************************************************************************************

    pin_pitch = 2 * params['Fuel Pin Radii'][-1] + params["Pin Gap Distance"]

    assembly_universes = []
    for i in range(N_AXIAL):
        assembly_universes.append(
            create_assembly_universe(
                params,
                fuel_pin_universes_axial[i],
                moderator_pin_universe,
                pin_pitch,
                reflector,
                coolant_universe
            )
        )

    # **************************************************************************************************************************
    #                                                Sec. 1.4 : Control Drums
    # **************************************************************************************************************************

    # Derive Core Radius and reflector thicknesses from the actual drum layout,
    # then build the drum positions and universes.
    control_drum_positions = update_ltmr_reflector_geometry_from_drums(params)
    drums = create_drums_universe(
        params, control_drum_absorber, control_drum_reflector, control_drum_positions
    )

    # **************************************************************************************************************************
    #                                                Sec. 1.5 : Volume Info for Depletion
    # **************************************************************************************************************************

    fuel_index   = params['Fuel Pin Materials'].index(params['Fuel'])
    fissile_area = circle_area(params['Fuel Pin Radii'][fuel_index]) \
                 - circle_area(params['Fuel Pin Radii'][fuel_index - 1])
    zone_height  = params['Active Height'] / N_AXIAL
    for fm in fuel_materials_axial:
        fm.volume = fissile_area * zone_height * params['Fuel Pin Count']

    # **************************************************************************************************************************
    #                                                Sec. 1.6 : Core Geometry
    # **************************************************************************************************************************

    core_geometry, core = create_core_geometry(
        params,
        drums,
        drums_positions=control_drum_positions,
        assembly_universes=assembly_universes,
        reflector_material=reflector
    )

    core_geometry.export_to_xml()

    # Export materials AFTER geometry is built so traversal catches every clone
    materials = openmc.Materials(list(core_geometry.get_all_materials().values()))
    openmc.Materials.cross_sections = params['cross_sections_xml_location']
    materials.export_to_xml()

    if params['plotting'] == "Y":
        create_universe_plot(
            materials_database, core_geometry,
            plot_width=2.01 * params['Core Radius'],
            num_pixels=2000, font_size=32,
            title="Reactor Core", fig_size=8,
            output_file_name="core_sym.png"
        )

    # **************************************************************************************************************************
    #                                                Sec. 1.7 : TALLIES
    # **************************************************************************************************************************

    tallies_file = openmc.Tallies()

    group_edges = np.array([1e-5, 6.7e-2, 3.2e-1, 1, 4, 9.88, 4.81e1, 4.54e2, 4.9e4, 1.83e5, 8.21e5, 4e7])
    groups      = openmc.mgxs.EnergyGroups(group_edges)

    mgxs_lib                = openmc.mgxs.Library(core_geometry)
    mgxs_lib.energy_groups  = groups
    mgxs_lib.legendre_order = 1
    mgxs_lib.mgxs_types     = ['absorption', 'diffusion-coefficient', 'transport',
                                'scatter matrix', 'total', 'scatter']
    mgxs_lib.domain_type    = 'universe'
    mgxs_lib.domains        = [core]
    mgxs_lib.build_library()
    mgxs_lib.add_to_tallies_file(tallies_file, merge=False)

    pin_filter = openmc.DistribcellFilter(fuel_cell)
    pin_power  = openmc.Tally(name='pin_power_kappa')
    pin_power.scores  = ['kappa-fission']
    pin_power.filters = [pin_filter]
    tallies_file.append(pin_power)
    tallies_file.export_to_xml()

    # **************************************************************************************************************************
    #                                                Sec. 1.8 : SIMULATION SETTINGS
    # **************************************************************************************************************************

    settings = openmc.Settings()
    source   = openmc.Source()
    source.space = openmc.stats.Box(
        (-params['Core Radius'], -params['Core Radius'],
         -(params['Active Height']/2 + params['Axial Reflector Thickness'])),
        ( params['Core Radius'],  params['Core Radius'],
           params['Active Height']/2 + params['Axial Reflector Thickness'])
    )
    settings.source = source
    settings.batches  = 100
    settings.inactive = 50

    if 'Particles' in params.keys():
        settings.particles = int(params['Particles'])
    else:
        settings.particles = 1000

    if params['Isothermal Temperature Coefficients']:
        settings.temperature = {
            'default':   params['Common Temperature'],
            'method':    'interpolation',
            'tolerance': 50.0
        }
    else:
        settings.temperature = {'method': 'interpolation'}

    settings.export_to_xml()