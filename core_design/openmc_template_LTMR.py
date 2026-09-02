# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
# Importing libraries
import openmc
import numpy as np
from core_design.openmc_materials_database import collect_materials_data
from core_design.pins_arrangement import (
    SUPPORTED_LTMR_SHUTDOWN_ROD_COUNTS,
)
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
Helper functions are smaller, reusable functions that are defined to perform specific tasks,
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
        # Extract the radii values for different regions of the moderator pin from the input parameters
        pin_radii = {
            'moderator': params['Moderator Pin Radii'][0],
            'cladding': params['Moderator Pin Radii'][1]
        }
        region_keys = ['moderator', 'cladding', 'coolant']

    elif pin_type == 'fuel':
        # Extract the radii values for different regions of the fuel pin from the input parameters
        pin_radii = {
            'insert': params['Fuel Pin Radii'][0],
            'gap1': params['Fuel Pin Radii'][1],
            'fuel_meat': params['Fuel Pin Radii'][2],
            'gap2': params['Fuel Pin Radii'][3],
            'cladding': params['Fuel Pin Radii'][4]
        }
        region_keys = ['insert', 'gap1', 'fuel_meat', 'gap2', 'cladding', 'coolant']

    else:
        raise ValueError("Invalid pin type. Must be 'moderator' or 'fuel'.")

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

def create_shutdown_rod_universe(params, materials_database):
    """
    Simplified shutdown channel:

    ARO:
        B4C absorber and its SS304 cladding are withdrawn together;
        the entire channel contains NaK.

    ARI:
        Enriched B4C absorber and its SS304 cladding are inserted.
    """

    absorber_radius = params['Shutdown Rod Absorber Radius']
    clad_radius = params['Shutdown Rod Clad Radius']

    if absorber_radius <= 0.0:
        raise ValueError(
            "Shutdown Rod Absorber Radius must be greater than zero."
        )

    if absorber_radius >= clad_radius:
        raise ValueError(
            "Shutdown Rod Absorber Radius must be smaller than "
            "Shutdown Rod Clad Radius."
        )

    if clad_radius >= params['Fuel Pin Radii'][-1]:
        raise ValueError(
            "Shutdown Rod Clad Radius must be smaller than "
            "the existing pin outer radius."
        )

    absorber_surface = openmc.ZCylinder(
        r=absorber_radius,
        name='shutdown_absorber_surface'
    )

    clad_surface = openmc.ZCylinder(
        r=clad_radius,
        name='shutdown_clad_surface'
    )

    absorber_material = materials_database[
        params['Shutdown Rod Absorber']
    ]

    clad_material = materials_database[
        params['Shutdown Rod Cladding']
    ]

    coolant_material = materials_database[
        params['Coolant']
    ]

    if params['Shutdown Margin Calc']:
        print(">>> SHUTDOWN RODS: INSERTED (B4C AND CLADDING)")
        absorber_cell = openmc.Cell(
            name='shutdown_absorber_inserted',
            fill=absorber_material,
            region=-absorber_surface
        )

        clad_cell = openmc.Cell(
            name='shutdown_rod_cladding_inserted',
            fill=clad_material,
            region=+absorber_surface & -clad_surface
        )

        outside_coolant_cell = openmc.Cell(
            name='shutdown_rod_outer_coolant',
            fill=coolant_material,
            region=+clad_surface
        )

        cells = [
            absorber_cell,
            clad_cell,
            outside_coolant_cell,
        ]
    else:
        print(">>> SHUTDOWN RODS: WITHDRAWN (B4C AND CLADDING)")
        withdrawn_channel_cell = openmc.Cell(
            name='shutdown_channel_withdrawn',
            fill=coolant_material
        )
        cells = [withdrawn_channel_cell]

    shutdown_rod_universe = openmc.Universe(
        name='shutdown_rod_universe',
        cells=cells
    )

    return shutdown_rod_universe

def _get_valid_drum_counts():
    return [6, 12, 18, 24, 30, 36]


def _get_drum_layout_quantities(params, drum_radius):
    number_of_drums = params['Number of Drums']
    valid_drum_counts = _get_valid_drum_counts()
    if number_of_drums not in valid_drum_counts:
        raise ValueError(f"Number of Drums must be one of {valid_drum_counts}, got {number_of_drums}")

    drums_per_side = number_of_drums // 6
    hex_edge_length = calculate_hex_edge_length(params)
    apothem = calculate_hex_apothem(params)

    drum_tube_radius = drum_radius + drum_radius / 90.0
    side_length = hex_edge_length

    return drums_per_side, hex_edge_length, apothem, drum_tube_radius, side_length


def _drum_positions_for_radius(params, drum_radius):
    drums_per_side, _, apothem, drum_tube_radius, side_length = _get_drum_layout_quantities(params, drum_radius)

    face_angles = [k * np.pi / 3 for k in range(6)]
    positions = []

    for face_angle in face_angles:
        along_x = -np.sin(face_angle)
        along_y = np.cos(face_angle)

        radial_distance = apothem + drum_tube_radius
        face_center_x = radial_distance * np.cos(face_angle)
        face_center_y = radial_distance * np.sin(face_angle)

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

    drums_per_side = number_of_drums // 6
    hex_edge_length = calculate_hex_edge_length(params)
    side_length = hex_edge_length

    # Same-face spacing gives a safe upper bound for the physical drum radius
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


def create_drums_universe(params, control_drum_absorber_material, control_drum_reflector_material, drum_positions):
    number_of_drums = params['Number of Drums']
    valid_drum_counts = _get_valid_drum_counts()
    if number_of_drums not in valid_drum_counts:
        raise ValueError(f"Number of Drums must be one of {valid_drum_counts}, got {number_of_drums}")

    absorber_thickness = params['Drum Absorber Thickness']
    drum_radius = resolve_drum_radius(params)
    absorber_arc = np.deg2rad(params['Drum Absorber Arc Degrees'])

    # Shutdown Margin Calc = True means build drums in shutdown (ARI) orientation.
    # Shutdown Margin Calc = False means build drums in operating (ARO) orientation.
    drum_state = 'shutdown' if params['Shutdown Margin Calc'] else 'operation'
    rotation_angle = 0 if drum_state == 'shutdown' else 180

    cd_inner_shell = openmc.ZCylinder(r=drum_radius - absorber_thickness)
    cd_outer_shell = openmc.ZCylinder(r=drum_radius)
    cutting_plane_1 = openmc.Plane(a=np.sin(absorber_arc / 2), b=np.cos(absorber_arc / 2))
    cutting_plane_2 = openmc.Plane(a=np.sin(absorber_arc / 2), b=-np.cos(absorber_arc / 2))

    drum_absorber_region = +cd_inner_shell & -cd_outer_shell & -cutting_plane_1 & -cutting_plane_2
    drum_reflector_region = -cd_outer_shell & ~drum_absorber_region
    drum_outside_region = +cd_outer_shell

    drum_absorber = openmc.Cell(
        name='drum_absorber',
        fill=control_drum_absorber_material,
        region=drum_absorber_region
    )
    drum_reflector = openmc.Cell(
        name='drum_reflector',
        fill=control_drum_reflector_material,
        region=drum_reflector_region
    )
    drum_exterior = openmc.Cell(name='drum_outside', region=drum_outside_region)
    drum_reference = openmc.Universe(cells=(drum_reflector, drum_absorber, drum_exterior))

    drum_cells = []
    for i, (_, _, face_angle_deg) in enumerate(drum_positions):
        dc = openmc.Cell(name=f'drum_{i}', fill=drum_reference)
        dc.rotation = [0, 0, face_angle_deg + rotation_angle]
        drum_cells.append(dc)

    drums = [openmc.Universe(cells=(dc,)) for dc in drum_cells]
    return drums


# def create_assembly_universe(params, fuel_pin_universe, moderator_pin_universe, pin_pitch, reflector_material, outer_coolant_universe):
#     """
#     Creating the universe of the fuel assembly
#     @ In, params, dict, The parameters that are used to "fill in" input files with placeholders.
#     @ In, fuel_pin_universe, openmc.universe.Universe
#     @ In, moderator_pin_universe, openmc.universe.Universe
#     @ In, pin_pitch, float, the center-to-center distance between adjacent fuel/moderator pins
#     @ In, reflector_material, openmc.material.Material, the material of the outer radial reflector
#     @ In, outer_coolant_universe, openmc.universe.Universe, the OpenMC universe of the coolant in the assembly
#     @ out, assembly_universe, openmc.universe.Universe, the fuel assembly universe
#     """

#     assembly = openmc.HexLattice()
#     assembly.center = (0., 0.)
#     assembly.pitch = (pin_pitch,)
#     assembly.outer = outer_coolant_universe

#     rings = copy.deepcopy(params['Pins Arrangement'])
#     rings = rings[-params['Number of Rings per Assembly']:]

#     for i in range(len(rings)):
#         for j in range(len(rings[i])):
#             if rings[i][j] == 'FUEL':
#                 rings[i][j] = fuel_pin_universe
#             elif rings[i][j] == 'MODERATOR':
#                 rings[i][j] = moderator_pin_universe

#     assembly.universes = rings

#     hex_edge_length = calculate_hex_edge_length(params)
#     assembly_boundary = openmc.model.hexagonal_prism(
#         edge_length=hex_edge_length,
#         corner_radius=params['Fuel Pin Radii'][-1] + params["Pin Gap Distance"]
#     )

#     fuel_assembly_cell = openmc.Cell(fill=assembly, region=assembly_boundary)
#     reflector_cell = openmc.Cell(fill=reflector_material, region=~assembly_boundary)

#     assembly_universe = openmc.Universe(cells=[fuel_assembly_cell, reflector_cell])

#     return assembly_universe


def _validate_shutdown_rod_layout(params, rings):
    """Validate the selected active LTMR shutdown-channel configuration."""
    actual_count = sum(row.count('SHUTDOWN') for row in rings)
    requested_count = params.get('Number of Shutdown Rods', actual_count)

    if requested_count not in SUPPORTED_LTMR_SHUTDOWN_ROD_COUNTS:
        raise ValueError(
            "Number of Shutdown Rods must be one of "
            f"{list(SUPPORTED_LTMR_SHUTDOWN_ROD_COUNTS)}, got "
            f"{requested_count!r}."
        )

    if actual_count != requested_count:
        raise ValueError(
            "The selected Pins Arrangement contains "
            f"{actual_count} shutdown channels, but Number of Shutdown "
            f"Rods is {requested_count}. Select the matching predefined "
            "LTMR pin arrangement."
        )

    for ring_index, ring in enumerate(rings):
        if len(ring) == 1:
            continue
        if len(ring) % 6 != 0:
            raise ValueError(
                f"LTMR ring {ring_index} has {len(ring)} positions; "
                "a sixfold-symmetric hexagonal ring must be divisible by 6."
            )

        positions_per_sector = len(ring) // 6
        reference_sector = [
            value == 'SHUTDOWN'
            for value in ring[:positions_per_sector]
        ]
        for sector_index in range(1, 6):
            start = sector_index * positions_per_sector
            sector = [
                value == 'SHUTDOWN'
                for value in ring[start:start + positions_per_sector]
            ]
            if sector != reference_sector:
                raise ValueError(
                    "Shutdown channels in LTMR ring "
                    f"{ring_index} are not sixfold symmetric."
                )

    return actual_count

def create_assembly_universe(params, fuel_pin_universe, moderator_pin_universe, shutdown_rod_universe, pin_pitch, reflector_material, outer_coolant_universe):
    """
    Creating the universe of the fuel assembly.
    """

    assembly = openmc.HexLattice()
    assembly.center = (0.0, 0.0)
    assembly.pitch = (pin_pitch,)
    assembly.outer = outer_coolant_universe

    rings = copy.deepcopy(params['Pins Arrangement'])
    rings = rings[-params['Number of Rings per Assembly']:]
    _validate_shutdown_rod_layout(params, rings)

    for i in range(len(rings)):
        for j in range(len(rings[i])):
            if rings[i][j] == 'FUEL':
                rings[i][j] = fuel_pin_universe

            elif rings[i][j] == 'MODERATOR':
                rings[i][j] = moderator_pin_universe

            elif rings[i][j] == 'SHUTDOWN':
                rings[i][j] = shutdown_rod_universe

            else:
                raise ValueError(
                    f"Unknown pin label at ring {i}, position {j}: "
                    f"{rings[i][j]!r}"
                )

    assembly.universes = rings

    hex_edge_length = calculate_hex_edge_length(params)

    assembly_boundary = openmc.model.hexagonal_prism(
        edge_length=hex_edge_length,
        corner_radius=(
            params['Fuel Pin Radii'][-1]
            + params['Pin Gap Distance']
        )
    )

    fuel_assembly_cell = openmc.Cell(
        fill=assembly,
        region=assembly_boundary
    )

    reflector_cell = openmc.Cell(
        fill=reflector_material,
        region=~assembly_boundary
    )

    assembly_universe = openmc.Universe(
        cells=[
            fuel_assembly_cell,
            reflector_cell
        ]
    )

    return assembly_universe


def create_control_drums_positions(params):
    """
    Place N/6 drums along each of the 6 flat faces of the hexagonal lattice,
    touching each face from outside, evenly spaced along the face length.
    """
    number_of_drums = params['Number of Drums']
    valid_drum_counts = _get_valid_drum_counts()
    if number_of_drums not in valid_drum_counts:
        raise ValueError(
            f"Number of Drums must be one of {valid_drum_counts}, got {number_of_drums}"
        )

    drums_per_side = number_of_drums // 6

    hex_edge_length = calculate_hex_edge_length(params)
    apothem = calculate_hex_apothem(params)

    drum_radius = resolve_drum_radius(params)
    drum_tube_radius = drum_radius + drum_radius / 90.0
    side_length = hex_edge_length

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

    positions = []
    for face_angle in face_angles:
        along_x = -np.sin(face_angle)
        along_y = np.cos(face_angle)

        # Drum centers are placed just outside each hex face
        radial_distance = apothem + drum_tube_radius
        face_center_x = radial_distance * np.cos(face_angle)
        face_center_y = radial_distance * np.sin(face_angle)

        for i in range(drums_per_side):
            offset = side_length * (i - (drums_per_side - 1) / 2.0) / drums_per_side
            x = face_center_x + offset * along_x
            y = face_center_y + offset * along_y
            positions.append((x, y, np.degrees(face_angle)))

    # Full pairwise overlap check at the placement stage
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
    For LTMR, derive:
    - Core Radius
    - Radial Reflector Thickness
    - Axial Reflector Thickness
    - Drum Height

    from the actual drum layout implied by Number of Drums and Drum Radius.
    """
    drum_radius = resolve_drum_radius(params)
    drum_tube_radius = drum_radius + drum_radius / 90.0

    drum_positions = create_control_drums_positions(params)

    max_outer_radius = max(
        np.sqrt(x ** 2 + y ** 2) + drum_tube_radius
        for x, y, _ in drum_positions
    )

    hex_apothem = calculate_hex_apothem(params)

    params['Core Radius'] = max_outer_radius
    params['Radial Reflector Thickness'] = params['Core Radius'] - hex_apothem
    params['Axial Reflector Thickness'] = params['Radial Reflector Thickness']
    params['Drum Height'] = params['Active Height'] + 2 * params['Axial Reflector Thickness']

    return drum_positions


def create_core_geometry(params, drums, drums_positions, assembly_universe):
    """
    Build the full 2D radial core geometry with reflector-embedded control drums.
    Includes overlap checks and places each drum at its prescribed position.
    """
    drum_radius = resolve_drum_radius(params)
    params['Drum Tube Radius'] = drum_radius + drum_radius / 90.0
    drum_tube_radius = params['Drum Tube Radius']

    # Outer vacuum boundary
    outer_surface = openmc.ZCylinder(r=params['Core Radius'], boundary_type='vacuum')

    # Pairwise overlap check
    for i in range(len(drums_positions)):
        x1, y1, _ = drums_positions[i]
        for j in range(i + 1, len(drums_positions)):
            x2, y2, _ = drums_positions[j]
            dist = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
            if dist < 2.0 * drum_tube_radius:
                overlap = 2.0 * drum_tube_radius - dist
                raise ValueError(
                    f"Drums {i} and {j} overlap by {overlap:.3f} cm! "
                    f"Reduce Drum Radius or reduce Number of Drums."
                )

    drum_shells = []
    drum_cells = []

    for (x, y, _), du in zip(drums_positions, drums):
        drum_shell = openmc.ZCylinder(x0=x, y0=y, r=drum_tube_radius)
        drum_shells.append(drum_shell)

        # Limit each drum cell to the outer vacuum boundary too
        drum_cell = openmc.Cell(fill=du, region=-drum_shell & -outer_surface)
        drum_cell.translation = (x, y, 0)
        drum_cells.append(drum_cell)

    drums_outside = +drum_shells[0]
    for d in drum_shells[1:]:
        drums_outside = drums_outside & +d

    core_cell = openmc.Cell(fill=assembly_universe, region=-outer_surface & drums_outside)

    core = openmc.Universe(cells=[core_cell] + drum_cells)
    core_geometry = openmc.Geometry(core)

    return core_geometry, core


# **************************************************************************************************************************
#                                                Sec. 1 : OpenMC Model
# **************************************************************************************************************************

"""
An OpenMC function that accepts an instance of "parameters"
and generates the necessary XML files.
"""


def _load_depleted_fuel_override(params):
    """Load the operating depletion fuel composition for a static case."""
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


def build_openmc_model_LTMR(params):
    """
    OpenMC Model
    @ In, params, watts.parameters.Parameters, The parameters that are used to "fill in"
    input files with placeholders. params mostly behaves like a Python dictionary with a few extra capabilities.
    """

    params.setdefault('Shutdown Margin Calc', False)
    params.setdefault('Isothermal Temperature Coefficients', False)

    # Ensure Drum Radius is always numeric before any downstream drum geometry use
    resolve_drum_radius(params)

    # **************************************************************************************************************************
    #                                                Sec. 1.1 : MATERIALS
    # **************************************************************************************************************************

    # Read all material properties
    materials_database = collect_materials_data(params)

    # Read material properties for the fuel, coolant, reflector, and control drums
    # (each drum contains two materials: absorber and reflector)
    fuel = materials_database[params['Fuel']]
    depleted_fuel_override = _load_depleted_fuel_override(params)
    if depleted_fuel_override is None:
        params['_Operating Fuel Material ID'] = int(fuel.id)
    else:
        fuel = depleted_fuel_override
        materials_database[params['Fuel']] = fuel
    coolant = materials_database[params['Coolant']]
    reflector = materials_database[params['Radial Reflector']]
    control_drum_absorber = materials_database[params['Control Drum Absorber']]
    control_drum_reflector = materials_database[params['Control Drum Reflector']]
    shutdown_rod_absorber = materials_database[params['Shutdown Rod Absorber']]
    shutdown_rod_cladding = materials_database[params['Shutdown Rod Cladding']]

    # **************************************************************************************************************************
    #                                                Sec. 1.2 : Pin Cell Universes and Coolant
    # **************************************************************************************************************************

    # Create fuel pin regions
    fuel_pin_regions = create_pin_regions(params, 'fuel')

    # Create fuel pin materials
    fuel_materials = []
    for mat in params['Fuel Pin Materials']:
        if mat is None:
            fuel_materials.append(None)
        else:
            material_1 = materials_database[mat]
            fuel_materials.append(material_1)
    fuel_materials.append(coolant)

    # Give the user an error message if the number of materials is not the same as the number of regions
    if len(fuel_pin_regions) != len(fuel_materials):
        raise ValueError(
            f"The number of fuel pin regions ({len(fuel_pin_regions)}) must match "
            f"the number of introduced materials ({len(fuel_materials)})."
        )

    # Create the fuel pin universe
    fuel_cells = create_cells(fuel_pin_regions, fuel_materials)

    # Fuel region cell (used in distribcell tally)
    fuel_cell = fuel_cells['fuel_meat']
    fuel_pin_universe = openmc.Universe(cells=fuel_cells.values())

    if params['plotting'] == "Y":
        create_universe_plot(
            materials_database,
            fuel_pin_universe,
            plot_width=2.2 * params['Fuel Pin Radii'][-1],
            num_pixels=500,
            font_size=32,
            title="Fuel Pin Universe",
            fig_size=8,
            output_file_name="fuel_pin_universe.png"
        )

    # Create moderator pin regions
    moderator_pin_regions = create_pin_regions(params, 'moderator')

    # Create moderator pin materials
    moderator_materials = []
    for mat in params['Moderator Pin Materials']:
        if mat is None:
            moderator_materials.append(None)
        else:
            material_1 = materials_database[mat]
            moderator_materials.append(material_1)
    moderator_materials.append(coolant)

    # Give the user an error message if the number of materials is not the same as the number of regions
    if len(moderator_pin_regions) != len(moderator_materials):
        raise ValueError(
            f"The number of moderator pin regions ({len(moderator_pin_regions)}) must match "
            f"the number of introduced materials ({len(moderator_materials)})."
        )

    # Create the moderator pin universe
    moderator_cells = create_cells(moderator_pin_regions, moderator_materials)
    moderator_pin_universe = openmc.Universe(cells=moderator_cells.values())

    if params['plotting'] == "Y":
        create_universe_plot(
            materials_database,
            moderator_pin_universe,
            plot_width=2.2 * params['Moderator Pin Radii'][-1],
            num_pixels=500,
            font_size=32,
            title="Moderator Pin Universe",
            fig_size=8,
            output_file_name="moderator_pin_universe.png"
        )

    # Create shutdown-rod universe
    shutdown_rod_universe = create_shutdown_rod_universe(params, materials_database)

    if params['plotting'] == "Y":
        create_universe_plot(
            materials_database,
            shutdown_rod_universe,
            plot_width=2.2 * params['Fuel Pin Radii'][-1],
            num_pixels=500,
            font_size=32,
            title="Shutdown Rod Universe",
            fig_size=8,
            output_file_name="shutdown_rod_universe.png"
        )

    # Coolant universe
    coolant_cell = openmc.Cell(fill=coolant)
    coolant_universe = openmc.Universe(cells=(coolant_cell,))

    # **************************************************************************************************************************
    #                                                Sec. 1.3 : Fuel Assembly Universe
    # **************************************************************************************************************************

    # Center-to-center distance between adjacent fuel/moderator pins
    pin_pitch = 2 * params['Fuel Pin Radii'][-1] + params["Pin Gap Distance"]

    assembly_universe = create_assembly_universe(
        params,
        fuel_pin_universe,
        moderator_pin_universe,
        shutdown_rod_universe,
        pin_pitch,
        reflector,
        coolant_universe
    )

    # **************************************************************************************************************************
    #                                                Sec. 1.4 : Material Volumes and Materials XML
    # **************************************************************************************************************************
    # Find the fuel region index within the fuel pin
    fuel_index = params['Fuel Pin Materials'].index(params['Fuel'])

    fissile_area = (
        circle_area(params['Fuel Pin Radii'][fuel_index])
        - circle_area(params['Fuel Pin Radii'][fuel_index - 1])
    )
    fuel.volume = fissile_area * params['Active Height'] * params['Fuel Pin Count']

    all_materials = (
        fuel_materials
        + moderator_materials
        + [coolant, reflector, control_drum_absorber, control_drum_reflector, shutdown_rod_absorber, shutdown_rod_cladding]
    )

    # Remove None materials while preserving their deterministic order
    all_materials_cleaned_list = [item for item in all_materials if item is not None]
    materials = openmc.Materials(list(dict.fromkeys(all_materials_cleaned_list)))

    openmc.Materials.cross_sections = params['cross_sections_xml_location']
    materials.export_to_xml()

    # **************************************************************************************************************************
    #                                                Sec. 1.5 : Control Drum Placement and Core Geometry
    # **************************************************************************************************************************

    control_drum_positions = update_ltmr_reflector_geometry_from_drums(params)
    drums = create_drums_universe(
        params,
        control_drum_absorber,
        control_drum_reflector,
        control_drum_positions
    )

    core_geometry, core = create_core_geometry(
        params,
        drums,
        drums_positions=control_drum_positions,
        assembly_universe=assembly_universe
    )

    core_geometry.export_to_xml()

    # if params['plotting'] == "Y":
    #     drum_state_label = "shutdown" if params['Shutdown Margin Calc'] else "operation"
    #     create_universe_plot(
    #         materials_database,
    #         core_geometry,
    #         plot_width=2.01 * params['Core Radius'],
    #         num_pixels=2000,
    #         font_size=32,
    #         title="Reactor Core",
    #         fig_size=8,
    #         output_file_name=f"core_{drum_state_label}.png"
    #     )

    if params['plotting'] == "Y":

        print(
            "CORE PLOT: Shutdown Margin Calc =",
            params['Shutdown Margin Calc']
        )

        drum_state_label = (
            "shutdown"
            if params['Shutdown Margin Calc']
            else "operation"
        )

        print(
            "CORE PLOT FILE =",
            f"core_{drum_state_label}.png"
        )

        create_universe_plot(
            materials_database,
            core_geometry,
            plot_width=2.01 * params['Core Radius'],
            num_pixels=2000,
            font_size=32,
            title=f"Reactor Core - {drum_state_label}",
            fig_size=8,
            output_file_name=f"core_{drum_state_label}.png"
        )

    # **************************************************************************************************************************
    #                                                Sec. 1.6 : TALLIES
    # **************************************************************************************************************************

    tallies_file = openmc.Tallies()

    # 11 energy groups from HPMR report table 5 in eV
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

    # Peaking factor tally (pin power)
    pin_filter = openmc.DistribcellFilter(fuel_cell)
    pin_power = openmc.Tally(name='pin_power_kappa')
    pin_power.scores = ['kappa-fission']
    pin_power.filters = [pin_filter]
    tallies_file.append(pin_power)
    tallies_file.export_to_xml()

    # **************************************************************************************************************************
    #                                                Sec. 1.7 : OpenMC Settings
    # **************************************************************************************************************************

    point = openmc.stats.Point((0, 0, 0))
    source = openmc.Source(space=point)
    settings = openmc.Settings()
    settings.source = source
    settings.batches = 100
    settings.inactive = 50

    if 'Particles' in params.keys():
        settings.particles = int(params['Particles'])
    else:
        settings.particles = 2000

    if '_OpenMC Seed' in params:
        settings.seed = int(params['_OpenMC Seed'])

    settings.temperature = {
        'default': params['Common Temperature'],
        'method': 'interpolation',
        'tolerance': 50.0
    }

    settings.export_to_xml()
