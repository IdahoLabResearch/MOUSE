# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
import numpy as np

from core_design.openmc_materials_database import collect_materials_data
from core_design.utils import (
    circle_area,
    cylinder_volume,
    calculate_number_of_rings,
    calculate_hex_edge_length,
    calculate_hex_apothem,
)


def _get_valid_ltmr_drum_counts():
    return [6, 12, 18, 24, 30, 36]


def _get_ltmr_drum_layout_quantities(params, drum_radius):
    number_of_drums = params['Number of Drums']
    valid_drum_counts = _get_valid_ltmr_drum_counts()
    if number_of_drums not in valid_drum_counts:
        raise ValueError(f"Number of Drums must be one of {valid_drum_counts}, got {number_of_drums}")

    drums_per_side = number_of_drums // 6
    hex_edge_length = calculate_hex_edge_length(params)
    apothem = calculate_hex_apothem(params)

    drum_tube_radius = drum_radius + drum_radius / 90.0
    side_length = hex_edge_length

    return drums_per_side, hex_edge_length, apothem, drum_tube_radius, side_length


def _ltmr_drum_positions_for_radius(params, drum_radius):
    drums_per_side, _, apothem, drum_tube_radius, side_length = _get_ltmr_drum_layout_quantities(params, drum_radius)

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


def _ltmr_drum_radius_is_feasible(params, drum_radius):
    positions, drum_tube_radius, side_length = _ltmr_drum_positions_for_radius(params, drum_radius)
    drums_per_side = params['Number of Drums'] // 6

    same_face_spacing = side_length / drums_per_side
    if 2.0 * drum_tube_radius > same_face_spacing:
        return False

    min_center_dist = 2.0 * drum_tube_radius
    for i in range(len(positions)):
        x1, y1, _ = positions[i]
        for j in range(i + 1, len(positions)):
            x2, y2, _ = positions[j]
            dist = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
            if dist < min_center_dist:
                return False

    return True


def _calculate_max_ltmr_drum_radius(params, tol=1e-6, max_iter=100):
    number_of_drums = params['Number of Drums']
    valid_drum_counts = _get_valid_ltmr_drum_counts()
    if number_of_drums not in valid_drum_counts:
        raise ValueError(f"Number of Drums must be one of {valid_drum_counts}, got {number_of_drums}")

    drums_per_side = number_of_drums // 6
    side_length = calculate_hex_edge_length(params)

    upper_bound = (side_length / (2.0 * drums_per_side)) * 90.0 / 91.0
    lower_bound = 0.0

    for _ in range(max_iter):
        mid = 0.5 * (lower_bound + upper_bound)
        if _ltmr_drum_radius_is_feasible(params, mid):
            lower_bound = mid
        else:
            upper_bound = mid

        if upper_bound - lower_bound < tol:
            break

    return lower_bound


def _calculate_max_gcmr_drum_radius(params):
    """
    Return the maximum allowable GCMR drum radius (cm).

    Each drum sits in a hexagonal cell of the outermost core ring with
    flat-to-flat size Assembly_FTF.  The drum tube (drum_radius * 46/45)
    must fit within the cell's inscribed-circle radius (apothem = FTF/2),
    giving:  drum_radius_max = Assembly_FTF / 2 * (45/46)
    """
    return params['Assembly FTF'] / 2 * (45 / 46)


def _calculate_max_hpmr_drum_radius(params):
    """
    Return the maximum HPMR drum radius (cm) such that no two adjacent
    drums overlap on the placement ring.

    Drums are placed on a ring of radius
        cd_distance = r0 + drum_tube_radius
    where r0 = (N_rings - 1)*Assembly_FTF + Assembly_FTF/2
    and   drum_tube_radius = drum_radius * (91/90).

    No-overlap condition (chord >= 2 * drum_tube_radius):
        2 * cd_distance * sin(pi/n) >= 2 * drum_tube_radius
    Solving the self-referential equation gives:
        max_drum_tube = r0 * sin(pi/n) / (1 - sin(pi/n))
    """
    n_drums = int(params.get('Drum Count', 12))
    r0 = ((params['Number of Rings per Core'] - 1) * params['Assembly FTF']
          + params['Assembly FTF'] / 2)
    sin_half = np.sin(np.pi / n_drums)
    max_drum_tube = r0 * sin_half / (1.0 - sin_half)
    return max_drum_tube * (90.0 / 91.0)


def _resolve_drum_radius(params):
    """
    If Drum Radius is not provided, set it to the maximum feasible value (cm).
    Supported for LTMR, GCMR, and HPMR reactor types.
    """
    if params.get('reactor type') == "LTMR":
        if 'Drum Radius' not in params:
            params['Drum Radius'] = _calculate_max_ltmr_drum_radius(params)

    elif params.get('reactor type') == "GCMR":
        if 'Drum Radius' not in params:
            params['Drum Radius'] = _calculate_max_gcmr_drum_radius(params)

    elif params.get('reactor type') == "HPMR":
        if 'Drum Radius' not in params:
            params['Drum Radius'] = _calculate_max_hpmr_drum_radius(params)

    if 'Drum Radius' not in params:
        raise KeyError("Drum Radius is required for this reactor type.")

    drum_radius = params['Drum Radius']
    if not isinstance(drum_radius, (int, float, np.floating)):
        raise ValueError(f"Drum Radius must be numeric if provided, got {drum_radius!r}")

    params['Drum Radius'] = float(drum_radius)
    return params['Drum Radius']


def calculate_drums_volumes_and_masses(params):
    drum_radius = _resolve_drum_radius(params)
    absorber_thickness = params['Drum Absorber Thickness']

    # --- GCMR: auto-resolve dependent geometry parameters ---
    # Drum ring cells sit at Core_Rings × Assembly_FTF from the core center.
    # The drum tube (radius = drum_radius × 46/45) extends that far beyond the ring,
    # so the reflector must be at least that thick to fully enclose the drums.
    if params.get('reactor type') == 'GCMR':
        drum_tube_radius = drum_radius * (46 / 45)

        if 'Radial Reflector Thickness' not in params:
            params['Radial Reflector Thickness'] = drum_tube_radius

        if 'Axial Reflector Thickness' not in params:
            params['Axial Reflector Thickness'] = params['Radial Reflector Thickness']

        # Core Radius is always kept consistent with Radial Reflector Thickness
        params['Core Radius'] = (params['Assembly FTF'] * params['Core Rings']
                                 + params['Radial Reflector Thickness'])

        if 'Drum Height' not in params:
            params['Drum Height'] = params['Active Height'] + 2 * params['Axial Reflector Thickness']

    # --- HPMR: auto-resolve dependent geometry parameters ---
    if params.get('reactor type') == 'HPMR':
        drum_tube_radius = drum_radius + drum_radius / 90.0
        n_drums = int(params.get('Drum Count', 12))
        r0 = ((params['Number of Rings per Core'] - 1) * params['Assembly FTF']
               + params['Assembly FTF'] / 2)
        cd_distance = r0 + drum_tube_radius
        hex_apothem = 0.5 * np.sqrt(3) * params['hexagonal Core Edge Length']

        if 'Radial Reflector Thickness' not in params:
            # Minimum reflector thickness that fully encloses all drums
            params['Radial Reflector Thickness'] = cd_distance + drum_radius - hex_apothem

        # Core Radius is always kept consistent with Radial Reflector Thickness
        params['Core Radius'] = hex_apothem + params['Radial Reflector Thickness']
        
        if 'Axial Reflector Thickness' not in params:
            params['Axial Reflector Thickness'] = params['Radial Reflector Thickness']

        if 'Drum Height' not in params:
            params['Drum Height'] = params['Active Height'] + 2 * params['Axial Reflector Thickness']

        # Validate that drums fit inside the reflector
        drum_outer = cd_distance + drum_radius
        if drum_outer > params['Core Radius']:
            raise ValueError(
                f"\n\n--- HPMR DRUM RADIUS ERROR ---\n"
                f"Drum outer extent ({drum_outer:.4f} cm) exceeds Core Radius "
                f"({params['Core Radius']:.4f} cm).\n"
                f"Reduce Drum Radius or increase Radial Reflector Thickness.\n"
                f"Maximum allowable Drum Radius: "
                f"{_calculate_max_hpmr_drum_radius(params):.4f} cm\n"
            )

        if drum_radius <= absorber_thickness:
            raise ValueError(
                f"\n\n--- HPMR DRUM RADIUS ERROR ---\n"
                f"Drum Radius ({drum_radius:.4f} cm) must be greater than "
                f"Drum Absorber Thickness ({absorber_thickness:.4f} cm).\n"
            )

        # Validate no-overlap on the placement ring
        chord = 2.0 * cd_distance * np.sin(np.pi / n_drums)
        if 2.0 * drum_tube_radius > chord:
            raise ValueError(
                f"\n\n--- HPMR DRUM OVERLAP ERROR ---\n"
                f"Drums overlap: tube diameter ({2*drum_tube_radius:.4f} cm) "
                f"exceeds the chord between adjacent drum centres ({chord:.4f} cm).\n"
                f"Maximum allowable Drum Radius: "
                f"{_calculate_max_hpmr_drum_radius(params):.4f} cm\n"
            )

    drum_height = params['Drum Height']

    # --- GCMR drum radius validation ---
    # Each drum sits inside a hexagonal cell of the outermost core ring.
    # The cell's inscribed-circle radius (apothem) is Assembly_FTF / 2.
    # The drum tube (drum + small clearance gap) must fit within this apothem.
    # The tube radius is drum_radius * (46/45) due to the 1/45 gap factor.
    if params.get('reactor type') == 'GCMR' and 'Assembly FTF' in params:
        apothem = params['Assembly FTF'] / 2
        drum_tube_radius = drum_radius * (46 / 45)
        max_drum_radius = apothem * (45 / 46)

        if drum_tube_radius > apothem:
            raise ValueError(
                f"\n\n--- DRUM RADIUS ERROR ---\n"
                f"Drum radius {drum_radius:.4f} cm is too large.\n"
                f"The drum tube radius ({drum_tube_radius:.4f} cm) exceeds the hex cell "
                f"apothem (Assembly FTF / 2 = {apothem:.4f} cm), causing the drum to "
                f"overlap into adjacent cells.\n"
                f"Maximum allowable drum radius: {max_drum_radius:.4f} cm\n"
            )

        if drum_radius <= absorber_thickness:
            raise ValueError(
                f"\n\n--- DRUM RADIUS ERROR ---\n"
                f"Drum radius ({drum_radius:.4f} cm) must be greater than the absorber "
                f"thickness ({absorber_thickness:.4f} cm).\n"
                f"The inner drum shell radius would be zero or negative, which is not physical.\n"
                f"Minimum allowable drum radius: > {absorber_thickness:.4f} cm\n"
            )

    drum_volume = np.pi * drum_radius * drum_radius * drum_height
    if 'coating_angle' in params:
        drum_absorp_vol = (
            np.pi * (drum_radius * drum_radius)
            - np.pi / 180 * params['coating_angle'] * (drum_radius - absorber_thickness) * (drum_radius - absorber_thickness)
        ) * drum_height / 3
    else:
        drum_absorp_vol = (
            np.pi * (
                drum_radius * drum_radius
                - (drum_radius - absorber_thickness) * (drum_radius - absorber_thickness)
            ) * drum_height
        ) / 3

    drum_refl_vol = drum_volume - drum_absorp_vol

    if params['reactor type'] == "LTMR":
        number_of_drums = params['Number of Drums']
        params['Drum Count'] = number_of_drums

    elif params['reactor type'] == "GCMR":
        if 'Drum Count' in params:
            number_of_drums = params['Drum Count']
        else:
            number_of_drums = 6 * (params['Core Rings'] - 1)
            params['Drum Count'] = number_of_drums

    elif params['reactor type'] == "HPMR":
        number_of_drums = int(params.get('Drum Count', 12))
        valid_drum_counts = [6, 12, 18, 24]
        if number_of_drums not in valid_drum_counts:
            raise ValueError(f"Number of Drums must be one of {valid_drum_counts}, got {number_of_drums}")
        params['Drum Count'] = number_of_drums

    else:
        raise ValueError(f"Unsupported reactor type: {params['reactor type']}")

    all_drums_volume = drum_volume * number_of_drums
    drum_absorp_vol_all = drum_absorp_vol * number_of_drums
    drum_refl_vol_all = drum_refl_vol * number_of_drums

    materials_database = collect_materials_data(params)
    drums_absorber_density = materials_database[params['Control Drum Absorber']].density
    drums_reflector_density = materials_database[params['Control Drum Reflector']].density

    drum_absorp_all_mass = drum_absorp_vol_all * drums_absorber_density / 1000  # kg
    drum_refl_all_mass = drum_refl_vol_all * drums_reflector_density / 1000  # kg

    control_drums_mass = drum_absorp_all_mass + drum_refl_all_mass
    params['All Drums Volume'] = all_drums_volume
    params['Control Drum Absorber Mass'] = drum_absorp_all_mass
    params['Control Drum Reflector Mass'] = drum_refl_all_mass
    params['Control Drums Mass'] = control_drums_mass
    params['All Drums Area'] = params['All Drums Volume'] / params['Drum Height']


def hexagonal_area_from_ftf(ftf_distance):
    """Calculate hexagonal area directly from flat-to-flat distance."""
    return (np.sqrt(3) / 2) * ftf_distance ** 2


def calculate_reflector_mass_LTMR(params):
    _resolve_drum_radius(params)

    hex_area = hexagonal_area_from_ftf(params['Assembly FTF'])
    core_radius = params['Core Radius']
    area_of_all_drums = params['All Drums Area']
    drum_height = params['Drum Height']

    # Assumes all drums lie fully inside the reflector region.
    area_reflector = np.pi * core_radius * core_radius - hex_area - area_of_all_drums  # cm^2
    if area_reflector < 0:
        raise ValueError(
            "LTMR radial reflector area is negative. "
            f"Core Radius ({core_radius:.4f} cm) is too small for the active "
            f"hex area ({hex_area:.4f} cm^2) and drum area "
            f"({area_of_all_drums:.4f} cm^2). "
            "Update the LTMR reflector geometry from the drum layout before "
            "calculating reflector mass, or increase Radial Reflector Thickness."
        )

    vol_reflector = area_reflector * drum_height  # cm^3

    materials_database = collect_materials_data(params)
    rad_reflector_density = materials_database[params['Radial Reflector']].density
    ax_reflector_density = materials_database[params['Axial Reflector']].density

    mass_reflector_rad = vol_reflector * rad_reflector_density / 1000  # kg
    params['Radial Reflector Mass'] = mass_reflector_rad
    params['Axial Reflector Mass'] = (1 / 1000) * ax_reflector_density * cylinder_volume(
        core_radius,
        params['Axial Reflector Thickness']
    )


def calculate_reflector_mass_GCMR(params):
    materials_database = collect_materials_data(params)
    tot_number_assemblies = calculate_number_of_rings(params['Core Rings'])
    reflector_height = params['Active Height']
    reflector_volume = reflector_height * (
        circle_area(params['Core Radius'])
        - tot_number_assemblies * hexagonal_area_from_ftf(params['Assembly FTF'])
        - params['All Drums Area']
    )

    rad_reflector_density = materials_database[params['Radial Reflector']].density
    rad_reflector_mass = rad_reflector_density * reflector_volume / 1000  # kg
    params['Radial Reflector Mass'] = rad_reflector_mass
    params['Axial Reflector Mass'] = 2 * (1 / 1000) * materials_database[params['Axial Reflector']].density * cylinder_volume(
        params['Core Radius'],
        params['Axial Reflector Thickness']
    )


def calculate_moderator_mass_GCMR(params):
    materials_database = collect_materials_data(params)
    AR = params['Assembly Rings']
    CR = params['Core Rings']
    tot_number_assemblies = calculate_number_of_rings(CR)

    # Area of one hexagonal assembly cell in the core lattice
    hex_area = hexagonal_area_from_ftf(params['Assembly FTF'])

    # Fuel compact area per assembly (packing fraction only — graphite matrix inside the compact counts as moderator)
    num_fuel_regions_per_hex = calculate_number_of_rings(AR - 1)
    area_fuel_per_hex = params['Packing Fraction'] * circle_area(params['Compact Fuel Radius']) * num_fuel_regions_per_hex

    # Every hex cell in the assembly (fuel, booster, or coolant lattice hex) has 6 coolant channels
    # at its vertices. Vertex sharing gives 2 effective channels per cell across all calculate_number_of_rings(AR)
    # cells (inner fuel rings + outer booster/coolant ring).
    area_coolant_per_hex = 2 * calculate_number_of_rings(AR) * circle_area(params['Coolant Channel Radius'])

    # Booster pin count: each booster_lattice_hex has its pin at the CENTER (not a vertex), so it
    # is never shared between adjacent cells — 1 full pin per booster cell.
    # Standard assemblies have 6*(AR-1) booster cells in the outer ring.
    # Corner and edge assemblies have partial booster outer rings (some cells replaced by coolant lattice hexes).
    booster_materials = params['Moderator Booster Materials']
    booster_radii = params['Moderator Booster Radii']
    if len(booster_materials) != len(booster_radii):
        raise ValueError(
            f"'Moderator Booster Materials' (len={len(booster_materials)}) and "
            f"'Moderator Booster Radii' (len={len(booster_radii)}) must have the same length."
        )

    n_standard = calculate_number_of_rings(CR - 1)
    n_corners  = 6
    n_edges    = 6 * (CR - 2)
    total_booster_pins = (
        n_standard * 6 * (AR - 1)
        + n_corners * ((AR - 1) * 3 - 1)
        + n_edges   * ((AR - 1) * 4 - 1)
    )

    # Per-assembly average booster footprint (used to compute the moderator area displacement)
    area_moderator_booster_per_hex = (total_booster_pins / tot_number_assemblies) * circle_area(booster_radii[-1])

    # Per-region booster mass (annular regions from innermost to outermost)
    tot_booster_mass = 0.0
    for i, (mat_name, r_outer) in enumerate(zip(booster_materials, booster_radii)):
        r_inner = booster_radii[i - 1] if i > 0 else 0.0
        annular_area = circle_area(r_outer) - circle_area(r_inner)
        density = materials_database[mat_name].density
        region_mass = (
            total_booster_pins
            * annular_area
            * params['Active Height']
            * density
            / 1000
        )
        params[f'Moderator Booster Mass {mat_name}'] = region_mass
        tot_booster_mass += region_mass

    moderator_area = hex_area - area_fuel_per_hex - area_coolant_per_hex - area_moderator_booster_per_hex
    tot_moderator_mass = (
        tot_number_assemblies
        * moderator_area
        * params['Active Height']
        * materials_database[params['Moderator']].density
        / 1000
    )
    params['Moderator Mass'] = tot_moderator_mass
    params['Moderator Booster Mass'] = tot_booster_mass


def calculate_reflector_and_moderator_mass_HPMR(params):
    materials_database = collect_materials_data(params)
    assembly_long_diag = 1.1547 * params['Assembly FTF']
    assembly_side_length = params['Assembly FTF'] / np.sqrt(3)
    big_hex_ftf = (
        params['Number of Rings per Core'] * assembly_long_diag
        + (params['Number of Rings per Core'] - 1) * assembly_side_length
    )
    big_hex_area = hexagonal_area_from_ftf(big_hex_ftf)
    reflector_volume = (circle_area(params['Core Radius']) - big_hex_area) * params['Active Height']

    rad_reflector_density = materials_database[params['Radial Reflector']].density
    rad_reflector_mass = rad_reflector_density * reflector_volume / 1000  # kg
    params['Radial Reflector Mass'] = rad_reflector_mass
    params['Axial Reflector Mass'] = 2 * (1 / 1000) * materials_database[params['Axial Reflector']].density * cylinder_volume(
        params['Core Radius'],
        params['Axial Reflector Thickness']
    )

    # Moderator = big hex minus the fuel and heatpipes
    fuel_area = params['Fuel Pin Count'] * circle_area(params['Fuel Pin Radii'][-1])
    heatpipe_area = params['Number of Heatpipes'] * circle_area(params['Heat Pipe Radii'][-1])
    params['Moderator Total Area'] = big_hex_area - fuel_area - heatpipe_area
    params['Moderator Mass'] = (
        params['Moderator Total Area']
        * params['Active Height']
        * materials_database[params['Moderator']].density
        / 1000
    )


def calculate_reflector_and_moderator_mass_HPMR_vtb(params):
    materials_database = collect_materials_data(params)
    assembly_long_diag = 1.1547 * params['Assembly FTF']
    assembly_side_length = params['Assembly FTF'] / np.sqrt(3)
    big_hex_ftf = (
        params['Number of Rings per Core'] * assembly_long_diag
        + (params['Number of Rings per Core'] - 1) * assembly_side_length
    )
    big_hex_area = hexagonal_area_from_ftf(big_hex_ftf)
    rad_reflector_volume = (circle_area(params['Core Radius']) - big_hex_area) * params['Active Height']
    rad_reflector_density = materials_database[params['Radial Reflector']].density
    rad_reflector_mass = rad_reflector_density * rad_reflector_volume / 1000
    params['Radial Reflector Mass'] = rad_reflector_mass
    params['Axial Reflector Mass'] = 2 * (1 / 1000) * materials_database[params['Axial Reflector']].density * cylinder_volume(
        params['Core Radius'],
        params['Axial Reflector Thickness']
    )

    fuel_area = params['Fuel Pin Count'] * circle_area(params['Fuel Pin Radii'][-1])
    heatpipe_area = params['Number of Heatpipes'] * circle_area(params['Heat Pipe Radii'][-1])
    moderator_booster_area = params['Number of Moderator Booster'] * circle_area(params['Moderator Booster Raddi'])
    params['Moderator Booster Mass'] = (
        moderator_booster_area
        * params['Active Height']
        * materials_database[params['Moderator Booster']].density
        / 1000
    )

    drum_radius = params['Drum Radius']
    drum_area = np.pi * drum_radius * drum_radius
    number_of_drums = 12
    params['Moderator Total Area'] = big_hex_area - fuel_area - heatpipe_area - moderator_booster_area - drum_area * number_of_drums
    params['Moderator Mass'] = (
        params['Moderator Total Area']
        * params['Active Height']
        * materials_database[params['Moderator']].density
        / 1000
    )


def calculate_moderator_mass(params):
    # For the moderator pins
    materials_database = collect_materials_data(params)
    moderator_volume = params['Moderator Pin Count'] * circle_area(params['Moderator Pin Radii'][0]) * params['Active Height']
    moderator_mass = (1 / 1000) * moderator_volume * materials_database[params['Moderator Pin Materials'][0]].density
    return moderator_mass
