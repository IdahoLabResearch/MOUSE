# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
"""Design-dependent shutdown shielding estimate for MOUSE.

This module implements a deliberately compact rigorous-two-step (R2S)
calculation suitable for cost-estimation studies:

1. Read the BOL, 11-group core leakage spectrum from the normal MOUSE
   eigenvalue statepoint.
2. Transport that spectrum through a one-dimensional cylindrical model of
   the in-vessel shield, steel vessels, and vacuum gaps.  The resulting source
   is therefore located immediately outside the intake vessel even though the
   vessels are not added to the detailed core model.
3. For each trial WEP thickness, transport the outside-intake source through
   an axisymmetric WEP/pit model, activate graded ordinary-concrete regions for
   60 years, decay them for 90 days, remove the reactor and WEP, and calculate
   photon dose 30 cm above the pit opening.

The geometry is intentionally coarse (three wall and three base activation
zones by default) and the
result is intended for cost estimation, not licensing or detailed shielding
design.  Material definitions live in openmc_materials_database.py so another
shield material can later be substituted by database name.
"""

from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import sys
from typing import Any, Mapping, MutableMapping, Sequence
import xml.etree.ElementTree as ET

import numpy as np
import openmc
import openmc.deplete

from core_design.correction_factor import corrected_keff_static
from core_design.openmc_materials_database import collect_materials_data
from reactor_engineering_evaluation.vessels_calcs import vessels_specs


GROUP_EDGES_EV = np.array([
    1.0e-5,
    6.7e-2,
    3.2e-1,
    1.0,
    4.0,
    9.88,
    4.81e1,
    4.54e2,
    4.9e4,
    1.83e5,
    8.21e5,
    4.0e7,
])
ELECTRON_VOLT_J = 1.602176634e-19
PSV_PER_SECOND_TO_MREM_PER_HOUR = 3.6e-4
SHIELDING_IMPLEMENTATION_VERSION = (
    "2026-09-21-fast-r2s-v8-conservative-zero-score-extrapolation"
)


def _shielding_status(*args, **kwargs) -> None:
    """Write shielding progress where WATTS immediately mirrors it."""
    kwargs.setdefault("file", sys.stderr)
    kwargs.setdefault("flush", True)
    print(*args, **kwargs)


def _create_unique_run_directory(
    working_directory: str | os.PathLike[str],
) -> Path:
    """Atomically create a unique shielding directory for this MOUSE run."""
    parent = Path(working_directory).resolve()
    parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")

    # mkdir(exist_ok=False) is the atomic claim operation.  If two independent
    # jobs reach shielding during the same second, one receives the base name
    # and the other advances to _001, _002, and so on without sharing output.
    for collision_index in range(10_000):
        suffix = "" if collision_index == 0 else f"_{collision_index:03d}"
        candidate = parent / f"run_{timestamp}{suffix}"
        try:
            candidate.mkdir(exist_ok=False)
        except FileExistsError:
            continue
        return candidate

    raise RuntimeError(
        "Could not allocate a unique shielding run directory after 10,000 "
        f"attempts under {parent}."
    )


def _require_supported_openmc() -> None:
    version_text = str(getattr(openmc, "__version__", "0.0"))
    match = re.match(r"(\d+)\.(\d+)", version_text)
    version = tuple(map(int, match.groups())) if match else (0, 0)
    if version < (0, 14):
        raise RuntimeError(
            "Dynamic shielding requires OpenMC 0.14 or newer because it uses "
            "openmc.deplete.get_microxs_and_flux. Installed version: "
            f"{version_text}"
        )


@contextmanager
def _working_directory(path: Path):
    """Temporarily run OpenMC in *path* without deleting existing outputs."""
    path.mkdir(parents=True, exist_ok=True)
    original = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(original)


def _source_class():
    """Return the non-deprecated source class when the installed version has it."""
    return getattr(openmc, "IndependentSource", openmc.Source)


def _latest_statepoint(directory: Path) -> Path:
    candidates = list(directory.glob("statepoint.*.h5"))
    if not candidates:
        raise FileNotFoundError(f"No OpenMC statepoint was written in {directory}")
    return max(candidates, key=lambda item: item.stat().st_mtime)


def _required(params: Mapping[str, Any], names: Sequence[str]) -> None:
    missing = [name for name in names if name not in params]
    if missing:
        raise KeyError(
            "Dynamic shielding needs these parameters before run_openmc(): "
            + ", ".join(missing)
        )


def prepare_shielding_geometry(params: MutableMapping[str, Any]) -> None:
    """Update derived shield/vessel dimensions after the core has been built."""
    _required(
        params,
        [
            "Core Radius",
            "In Vessel Shield Thickness",
            "Active Height",
            "Axial Reflector Thickness",
            "Vessel Thickness",
            "Vessel Lower Plenum Height",
            "Vessel Upper Plenum Height",
            "Vessel Upper Gas Gap",
            "Vessel Bottom Depth",
            "Vessel Material",
            "Gap Between Vessel And Guard Vessel",
            "Guard Vessel Thickness",
            "Guard Vessel Material",
            "Gap Between Guard Vessel And Cooling Vessel",
            "Cooling Vessel Thickness",
            "Cooling Vessel Material",
            "Gap Between Cooling Vessel And Intake Vessel",
            "Intake Vessel Thickness",
            "Intake Vessel Material",
        ],
    )
    params["In Vessel Shield Inner Radius"] = float(params["Core Radius"])
    params["In Vessel Shield Outer Radius"] = (
        float(params["Core Radius"])
        + float(params["In Vessel Shield Thickness"])
    )
    params["Vessel Radius"] = params["In Vessel Shield Outer Radius"]
    vessels_specs(params)


def _group_midpoints(edges: np.ndarray = GROUP_EDGES_EV) -> np.ndarray:
    return np.sqrt(edges[:-1] * edges[1:])


def _energy_distribution(probabilities: Sequence[float]):
    values = np.asarray(probabilities, dtype=float)
    values = np.clip(values, 0.0, None)
    total = float(values.sum())
    if total <= 0.0:
        raise ValueError("The leakage energy spectrum contains no positive weight.")
    return openmc.stats.Discrete(_group_midpoints(), values / total)


def extract_bol_core_leakage(
    params: MutableMapping[str, Any],
    statepoint_file: str | os.PathLike[str],
) -> dict[str, Any]:
    """Convert the BOL current tally into a physical core-leakage source."""
    with openmc.StatePoint(statepoint_file) as statepoint:
        leakage_tally = statepoint.get_tally(
            name="boc_shielding_leakage_current"
        )
        current = np.asarray(leakage_tally.mean, dtype=float).reshape(
            -1, len(GROUP_EDGES_EV) - 1
        ).sum(axis=0)
        current_std = np.sqrt(
            np.square(
                np.asarray(leakage_tally.std_dev, dtype=float).reshape(
                    -1, len(GROUP_EDGES_EV) - 1
                )
            ).sum(axis=0)
        )
        kappa_fission_ev_per_source = float(
            np.asarray(
                statepoint.get_tally(name="boc_total_kappa_fission").mean,
                dtype=float,
            ).sum()
        )

    current = np.clip(current, 0.0, None)
    radial_leakage_fraction = float(current.sum())
    if radial_leakage_fraction <= 0.0:
        raise RuntimeError(
            "The BOL shielding leakage tally is zero. Check the named outer "
            "vacuum boundary and boc_shielding_leakage_current tally."
        )
    if kappa_fission_ev_per_source <= 0.0:
        raise RuntimeError("The BOL kappa-fission tally is zero.")

    total_height = (
        float(params["Active Height"])
        + 2.0 * float(params["Axial Reflector Thickness"])
    )
    total_leakage_fraction = radial_leakage_fraction
    axial_non_leakage_probability = 1.0
    diffusion_total_leakage_fraction = None
    leakage_method = "measured radial current"
    try:
        correction = corrected_keff_static(
            str(statepoint_file),
            total_height,
            core_radius=float(params["Core Radius"]),
        )
        axial_non_leakage_probability = float(
            correction["axial_non_leakage_probability"]
        )
        diffusion_total_non_leakage = float(
            correction["total_non_leakage_probability"]
        )
        if np.isfinite(diffusion_total_non_leakage):
            diffusion_total_leakage_fraction = float(
                np.clip(1.0 - diffusion_total_non_leakage, 0.0, 1.0)
            )
        if np.isfinite(axial_non_leakage_probability):
            axial_non_leakage_probability = float(
                np.clip(axial_non_leakage_probability, 0.0, 1.0)
            )
            radial_non_leakage_probability = 1.0 - radial_leakage_fraction
            total_leakage_fraction = 1.0 - (
                radial_non_leakage_probability
                * axial_non_leakage_probability
            )
            total_leakage_fraction = float(
                np.clip(total_leakage_fraction, 0.0, 1.0)
            )
            leakage_method = (
                "measured radial leakage combined with 1D diffusion axial "
                "non-leakage"
            )
    except Exception as error:  # The measured current remains a usable fallback.
        _shielding_status(
            "[MOUSE SHIELDING] Total-leakage correction unavailable; using "
            f"radial current ({error})."
        )

    fission_source_rate = (
        float(params["Power MWt"]) * 1.0e6
        / (kappa_fission_ev_per_source * ELECTRON_VOLT_J)
    )
    leakage_rate = fission_source_rate * total_leakage_fraction
    spectrum = current / radial_leakage_fraction

    params["BOL Radial Leakage Fraction"] = radial_leakage_fraction
    params["BOL Axial Non-Leakage Probability Used For Shielding"] = (
        axial_non_leakage_probability
    )
    params["BOL Total Leakage Fraction Used For Shielding"] = (
        total_leakage_fraction
    )
    params["BOL Core Leakage Source Rate (n/s)"] = leakage_rate

    return {
        "group_edges_eV": GROUP_EDGES_EV.tolist(),
        "probabilities": spectrum.tolist(),
        "probability_std_dev": (
            current_std / radial_leakage_fraction
        ).tolist(),
        "source_rate_n_per_s": leakage_rate,
        "radial_leakage_fraction_per_source": radial_leakage_fraction,
        "axial_non_leakage_probability": axial_non_leakage_probability,
        "total_leakage_fraction_per_source": total_leakage_fraction,
        "diffusion_total_leakage_fraction_diagnostic": (
            diffusion_total_leakage_fraction
        ),
        "kappa_fission_eV_per_source": kappa_fission_ev_per_source,
        "normalization_method": leakage_method,
        "location": "detailed-core outer boundary",
    }


def _unique_materials(materials: Sequence[openmc.Material]) -> openmc.Materials:
    by_id = {material.id: material for material in materials}
    result = openmc.Materials(list(by_id.values()))
    return result


def _fixed_source_settings(
    params: Mapping[str, Any],
    source,
    *,
    photon_transport: bool = False,
    batches: int | None = None,
    particles: int | None = None,
) -> openmc.Settings:
    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.batches = int(
        params.get("Shielding Batches", 10) if batches is None else batches
    )
    settings.particles = int(
        params.get("Shielding Particles", 5_000)
        if particles is None
        else particles
    )
    settings.source = source
    settings.photon_transport = photon_transport
    settings.temperature = {
        "method": "interpolation",
        "range": (250.0, max(1200.0, float(params.get("Common Temperature", 300.0)))),
    }
    return settings


def transport_source_to_intake_vessel(
    params: MutableMapping[str, Any],
    core_source: Mapping[str, Any],
    output_dir: str | os.PathLike[str],
) -> dict[str, Any]:
    """Run a 1D cylindrical transfer calculation through all vessel layers."""
    materials_database = collect_materials_data(params)
    layers = [
        (
            float(params["In Vessel Shield Thickness"]),
            params["In Vessel Shield Material"],
            "in-vessel shield",
        ),
        (float(params["Vessel Thickness"]), params["Vessel Material"], "vessel"),
        (
            float(params["Gap Between Vessel And Guard Vessel"]),
            None,
            "vessel/guard vacuum gap",
        ),
        (
            float(params["Guard Vessel Thickness"]),
            params["Guard Vessel Material"],
            "guard vessel",
        ),
        (
            float(params["Gap Between Guard Vessel And Cooling Vessel"]),
            None,
            "guard/cooling vacuum gap",
        ),
        (
            float(params["Cooling Vessel Thickness"]),
            params["Cooling Vessel Material"],
            "cooling vessel",
        ),
        (
            float(params["Gap Between Cooling Vessel And Intake Vessel"]),
            None,
            "cooling/intake vacuum gap",
        ),
        (
            float(params["Intake Vessel Thickness"]),
            params["Intake Vessel Material"],
            "intake vessel",
        ),
    ]

    inner_radius = float(params["Core Radius"])
    inner_surface = openmc.ZCylinder(r=inner_radius)
    cells = [openmc.Cell(name="core leakage source region", region=-inner_surface)]
    used_materials: list[openmc.Material] = []

    for thickness, material_name, layer_name in layers:
        if thickness < 0.0:
            raise ValueError(f"Negative thickness for {layer_name}: {thickness}")
        if thickness == 0.0:
            continue
        outer_radius = inner_radius + thickness
        outer_surface = openmc.ZCylinder(r=outer_radius)
        material = None
        if material_name is not None:
            if material_name not in materials_database:
                raise KeyError(
                    f"Shielding material {material_name!r} is not in "
                    "openmc_materials_database.py"
                )
            material = materials_database[material_name]
            used_materials.append(material)
        cells.append(
            openmc.Cell(
                name=layer_name,
                fill=material,
                region=+inner_surface & -outer_surface,
            )
        )
        inner_radius = outer_radius
        inner_surface = outer_surface

    inner_surface.boundary_type = "vacuum"
    geometry = openmc.Geometry(openmc.Universe(cells=cells))

    source_space = openmc.stats.CylindricalIndependent(
        r=openmc.stats.PowerLaw(0.0, float(params["Core Radius"]) * 0.999999, 1.0),
        phi=openmc.stats.Uniform(0.0, 2.0 * math.pi),
        z=openmc.stats.Discrete([0.0], [1.0]),
    )
    source = _source_class()(
        space=source_space,
        angle=openmc.stats.Isotropic(),
        energy=_energy_distribution(core_source["probabilities"]),
        particle="neutron",
    )
    settings = _fixed_source_settings(params, source)

    current_tally = openmc.Tally(name="outside_intake_leakage_current")
    current_tally.filters = [
        openmc.SurfaceFilter(inner_surface),
        openmc.EnergyFilter(GROUP_EDGES_EV),
    ]
    current_tally.scores = ["current"]
    model = openmc.Model(
        geometry=geometry,
        materials=_unique_materials(used_materials),
        settings=settings,
        tallies=openmc.Tallies([current_tally]),
    )
    model.materials.cross_sections = params["cross_sections_xml_location"]

    run_dir = Path(output_dir).resolve()
    with _working_directory(run_dir):
        model.run()
    statepoint_path = _latest_statepoint(run_dir)
    with openmc.StatePoint(statepoint_path) as statepoint:
        tally = statepoint.get_tally(name="outside_intake_leakage_current")
        current = np.asarray(tally.mean, dtype=float).reshape(
            -1, len(GROUP_EDGES_EV) - 1
        ).sum(axis=0)
        current_std = np.sqrt(
            np.square(
                np.asarray(tally.std_dev, dtype=float).reshape(
                    -1, len(GROUP_EDGES_EV) - 1
                )
            ).sum(axis=0)
        )

    current = np.clip(current, 0.0, None)
    transmission = float(current.sum())
    if transmission <= 0.0:
        raise RuntimeError(
            "No neutrons reached the outside of the intake vessel in the "
            "1D transfer calculation. Increase Shielding Particles."
        )
    source_rate = float(core_source["source_rate_n_per_s"]) * transmission
    outside_source = {
        "group_edges_eV": GROUP_EDGES_EV.tolist(),
        "probabilities": (current / transmission).tolist(),
        "probability_std_dev": (current_std / transmission).tolist(),
        "source_rate_n_per_s": source_rate,
        "vessel_transmission_per_core_leakage_neutron": transmission,
        "radius_cm": inner_radius,
        "location": "immediately outside intake vessel",
    }
    params["Leakage Source Radius Outside Intake Vessel (cm)"] = inner_radius
    params["Leakage Source Rate Outside Intake Vessel (n/s)"] = source_rate
    params["Vessel Neutron Transmission Fraction"] = transmission
    return outside_source


def _concrete_activation_depth_edges(
    params: Mapping[str, Any],
) -> list[float]:
    """Return concrete-zone depth edges measured from the pit surface."""
    concrete_thickness = float(
        params.get("Shielding Concrete Thickness", 200.0)
    )
    if concrete_thickness <= 0.0:
        raise ValueError("Shielding Concrete Thickness must be positive.")

    raw_boundaries = params.get(
        "Shielding Concrete Activation Zone Boundaries",
        [0.025 * concrete_thickness, 0.125 * concrete_thickness],
    )
    if not isinstance(raw_boundaries, Sequence) or isinstance(
        raw_boundaries, (str, bytes)
    ):
        raise TypeError(
            "Shielding Concrete Activation Zone Boundaries must be a "
            "sequence of depths in cm."
        )
    boundaries = [float(value) for value in raw_boundaries]
    if any(not np.isfinite(value) for value in boundaries):
        raise ValueError("Concrete activation-zone boundaries must be finite.")
    if any(value <= 0.0 or value >= concrete_thickness for value in boundaries):
        raise ValueError(
            "Concrete activation-zone boundaries must be greater than zero "
            "and less than Shielding Concrete Thickness."
        )
    if boundaries != sorted(set(boundaries)):
        raise ValueError(
            "Concrete activation-zone boundaries must be unique and strictly "
            "increasing."
        )
    return [0.0, *boundaries, concrete_thickness]


def _make_candidate_neutron_model(
    params: Mapping[str, Any],
    leakage_source: Mapping[str, Any],
    thickness_cm: float,
):
    """Build the axisymmetric WEP and graded concrete activation model."""
    materials_database = collect_materials_data(params)
    material_name = str(params.get("Out Of Vessel Shield Material", "WEP"))
    if material_name not in materials_database:
        raise KeyError(
            f"Out-of-vessel shield material {material_name!r} is not in "
            "openmc_materials_database.py"
        )
    shield_material = materials_database[material_name].clone()
    shield_material.name = f"{material_name} effective shielding material"
    density_factor = float(
        params.get("Out Of Vessel Shield Effective Density Factor", 1.0)
    )
    if density_factor <= 0.0:
        raise ValueError(
            "Out Of Vessel Shield Effective Density Factor must be positive."
        )
    shield_material.set_density(
        "g/cm3",
        float(shield_material.density) * density_factor,
    )
    concrete = materials_database["ordinary_concrete"]

    intake_radius = float(leakage_source["radius_cm"])
    source_height = float(params["Vessels Total Height"])
    source_half_height = 0.5 * source_height
    shield_outer_radius = intake_radius + thickness_cm
    shield_outer_half_height = source_half_height + thickness_cm
    clearance = float(params.get("Shielding Pit Clearance", 10.0))
    concrete_thickness = float(params.get("Shielding Concrete Thickness", 200.0))
    concrete_depth_edges = _concrete_activation_depth_edges(params)

    center_z = -clearance - shield_outer_half_height
    inner_bottom = center_z - source_half_height
    inner_top = center_z + source_half_height
    shield_bottom = center_z - shield_outer_half_height
    shield_top = center_z + shield_outer_half_height
    pit_bottom_value = shield_bottom - clearance
    pit_top_value = 0.0
    world_bottom_value = pit_bottom_value - concrete_thickness
    detector_distance = float(params.get("Shielding Detector Distance", 30.0))
    detector_radius = float(params.get("Shielding Detector Radius", 10.0))
    world_top_value = detector_distance + detector_radius + 10.0

    inner_cylinder = openmc.ZCylinder(r=intake_radius)
    shield_cylinder = openmc.ZCylinder(r=shield_outer_radius)
    pit_inner_cylinder = openmc.ZCylinder(r=shield_outer_radius + clearance)
    world_cylinder = openmc.ZCylinder(
        r=shield_outer_radius + clearance + concrete_thickness,
        boundary_type="vacuum",
    )
    inner_bottom_plane = openmc.ZPlane(z0=inner_bottom)
    inner_top_plane = openmc.ZPlane(z0=inner_top)
    shield_bottom_plane = openmc.ZPlane(z0=shield_bottom)
    shield_top_plane = openmc.ZPlane(z0=shield_top)
    pit_bottom_plane = openmc.ZPlane(z0=pit_bottom_value)
    pit_top_plane = openmc.ZPlane(z0=pit_top_value)
    world_bottom_plane = openmc.ZPlane(
        z0=world_bottom_value,
        boundary_type="vacuum",
    )
    world_top_plane = openmc.ZPlane(
        z0=world_top_value,
        boundary_type="vacuum",
    )

    inner_region = (
        -inner_cylinder & +inner_bottom_plane & -inner_top_plane
    )
    shield_outer_region = (
        -shield_cylinder & +shield_bottom_plane & -shield_top_plane
    )
    shield_region = shield_outer_region & ~inner_region
    pit_region = (
        -pit_inner_cylinder & +pit_bottom_plane & -pit_top_plane
    )
    above_region = (
        -world_cylinder & +pit_top_plane & -world_top_plane
    )

    shield_cell = None
    if thickness_cm > 0.0:
        shield_cell = openmc.Cell(
            name=f"{material_name} shield",
            fill=shield_material,
            region=shield_region,
        )
    pit_void_cell = openmc.Cell(
        name="pit void",
        region=(pit_region & ~shield_region) if shield_cell else pit_region,
    )
    above_cell = openmc.Cell(name="air/void above pit", region=above_region)

    side_cylinders = [pit_inner_cylinder]
    side_cylinders.extend(
        openmc.ZCylinder(r=pit_inner_cylinder.r + depth)
        for depth in concrete_depth_edges[1:-1]
    )
    side_cylinders.append(world_cylinder)

    base_planes = [pit_bottom_plane]
    base_planes.extend(
        openmc.ZPlane(z0=pit_bottom_value - depth)
        for depth in concrete_depth_edges[1:-1]
    )
    base_planes.append(world_bottom_plane)

    activation_cells: list[openmc.Cell] = []
    activation_materials: list[openmc.Material] = []
    concrete_zones: list[dict[str, Any]] = []
    concrete_height = pit_top_value - world_bottom_value

    for index, (inner_cylinder_zone, outer_cylinder_zone) in enumerate(
        zip(side_cylinders[:-1], side_cylinders[1:]),
        start=1,
    ):
        inner_depth = concrete_depth_edges[index - 1]
        outer_depth = concrete_depth_edges[index]
        material = concrete.clone()
        material.name = (
            "activated ordinary concrete - side "
            f"{inner_depth:g}-{outer_depth:g} cm"
        )
        material.depletable = True
        region = (
            +inner_cylinder_zone
            & -outer_cylinder_zone
            & +world_bottom_plane
            & -pit_top_plane
        )
        volume = (
            math.pi
            * (outer_cylinder_zone.r ** 2 - inner_cylinder_zone.r ** 2)
            * concrete_height
        )
        cell = openmc.Cell(
            name=f"ordinary concrete side activation zone {index}",
            fill=material,
            region=region,
        )
        cell.volume = volume
        material.volume = volume
        activation_cells.append(cell)
        activation_materials.append(material)
        concrete_zones.append({
            "kind": "side",
            "index": index,
            "inner_depth_cm": inner_depth,
            "outer_depth_cm": outer_depth,
            "r_inner_cm": float(inner_cylinder_zone.r),
            "r_outer_cm": float(outer_cylinder_zone.r),
            "z_lower_cm": world_bottom_value,
            "z_upper_cm": pit_top_value,
        })

    for index, (upper_plane, lower_plane) in enumerate(
        zip(base_planes[:-1], base_planes[1:]),
        start=1,
    ):
        inner_depth = concrete_depth_edges[index - 1]
        outer_depth = concrete_depth_edges[index]
        material = concrete.clone()
        material.name = (
            "activated ordinary concrete - base "
            f"{inner_depth:g}-{outer_depth:g} cm"
        )
        material.depletable = True
        region = -pit_inner_cylinder & +lower_plane & -upper_plane
        z_lower = float(lower_plane.z0)
        z_upper = float(upper_plane.z0)
        volume = math.pi * pit_inner_cylinder.r ** 2 * (z_upper - z_lower)
        cell = openmc.Cell(
            name=f"ordinary concrete base activation zone {index}",
            fill=material,
            region=region,
        )
        cell.volume = volume
        material.volume = volume
        activation_cells.append(cell)
        activation_materials.append(material)
        concrete_zones.append({
            "kind": "base",
            "index": index,
            "inner_depth_cm": inner_depth,
            "outer_depth_cm": outer_depth,
            "r_inner_cm": 0.0,
            "r_outer_cm": float(pit_inner_cylinder.r),
            "z_lower_cm": z_lower,
            "z_upper_cm": z_upper,
        })

    cells = [pit_void_cell, *activation_cells, above_cell]
    if shield_cell is not None:
        cells.insert(0, shield_cell)
    geometry = openmc.Geometry(openmc.Universe(cells=cells))
    source_radius = max(1.0e-6, intake_radius * (1.0 - 1.0e-7))
    source_space = openmc.stats.CylindricalIndependent(
        r=openmc.stats.Discrete([source_radius], [1.0]),
        phi=openmc.stats.Uniform(0.0, 2.0 * math.pi),
        z=openmc.stats.Uniform(inner_bottom, inner_top),
    )
    source = _source_class()(
        space=source_space,
        angle=openmc.stats.Isotropic(),
        energy=_energy_distribution(leakage_source["probabilities"]),
        particle="neutron",
    )
    settings = _fixed_source_settings(params, source)
    model = openmc.Model(
        geometry=geometry,
        materials=_unique_materials(
            [shield_material, *activation_materials]
        ),
        settings=settings,
    )
    model.materials.cross_sections = params["cross_sections_xml_location"]

    geometry_data = {
        "pit_inner_radius": pit_inner_cylinder.r,
        "world_radius": world_cylinder.r,
        "pit_bottom": pit_bottom_value,
        "pit_top": pit_top_value,
        "world_bottom": world_bottom_value,
        "world_top": world_top_value,
        "detector_distance": detector_distance,
        "detector_radius": detector_radius,
        "concrete_activation_depth_edges_cm": concrete_depth_edges,
        "concrete_zones": concrete_zones,
    }
    return model, activation_cells, activation_materials, geometry_data


def _activation_schedule(params: Mapping[str, Any]) -> tuple[list[float], list[float]]:
    irradiation_days = float(params.get("Shielding Irradiation Years", 60.0)) * 365.25
    decay_days = float(params.get("Shielding Decay Days", 90.0))
    if irradiation_days <= 0.0 or decay_days < 0.0:
        raise ValueError("The shielding irradiation/decay schedule is invalid.")

    # The independent activation operator uses fixed fluxes and microscopic
    # cross sections, so its depletion matrix is constant during irradiation.
    # Subdividing 60 years into hundreds of identical 90-day steps adds runtime
    # without improving this deliberately approximate shielding calculation.
    timesteps = [irradiation_days]
    source_flags = [1.0]
    if decay_days > 0.0:
        timesteps.append(decay_days)
        source_flags.append(0.0)
    return timesteps, source_flags


def _activate_concrete(
    params: Mapping[str, Any],
    model: openmc.Model,
    activation_cells: Sequence[openmc.Cell],
    activation_materials: Sequence[openmc.Material],
    leakage_rate: float,
    output_dir: Path,
) -> list[openmc.Material]:
    chain_file = Path(str(params["simplified_chain_thermal_xml"]))
    if not chain_file.is_file():
        raise FileNotFoundError(
            f"Shielding depletion chain not found: {chain_file}"
        )
    chain_root = ET.parse(chain_file).getroot()
    chain_nuclides = set()
    photon_source_nuclides = set()
    for node in chain_root.iter("nuclide"):
        name = node.attrib.get("name")
        if name is None:
            continue
        chain_nuclides.add(name)
        if any(
            source.attrib.get("particle") == "photon"
            for source in node.iter("source")
        ):
            photon_source_nuclides.add(name)

    if not photon_source_nuclides:
        raise RuntimeError(
            f"Depletion chain {chain_file} contains no decay-photon spectra. "
            "Use a chain containing <source particle=\"photon\" ...> entries."
        )

    diagnostic_products = {"Na24", "Mn54", "Fe59", "Co60"}
    retained_products = sorted(chain_nuclides & diagnostic_products)
    photon_products = sorted(photon_source_nuclides & diagnostic_products)
    if not retained_products or not photon_products:
        _shielding_status(
            "[MOUSE SHIELDING] WARNING: the selected shielding chain does not "
            "retain decay-photon data for any of Na24, Mn54, Fe59, or Co60. "
            "Concrete activation dose may be strongly underestimated."
        )
    _shielding_status(
        f"[MOUSE SHIELDING] Activation chain: {chain_file} "
        f"({len(photon_source_nuclides)} nuclides with photon spectra; "
        f"diagnostic products: {', '.join(photon_products) or 'none'})"
    )
    openmc.config["chain_file"] = str(chain_file)
    openmc.config["cross_sections"] = str(params["cross_sections_xml_location"])

    neutron_dir = output_dir / "neutron_activation_transport"
    _shielding_status(
        "[MOUSE SHIELDING] Calculating concrete activation fluxes and "
        "microscopic cross sections...",
        flush=True,
    )
    with _working_directory(neutron_dir):
        fluxes, micro_xs = openmc.deplete.get_microxs_and_flux(
            model,
            list(activation_cells),
            energies=GROUP_EDGES_EV,
            chain_file=str(chain_file),
        )
    _shielding_status(
        "[MOUSE SHIELDING] Concrete activation transport complete.",
        flush=True,
    )

    activation_dir = output_dir / "activation"
    activation_dir.mkdir(parents=True, exist_ok=True)
    operator = openmc.deplete.IndependentOperator(
        openmc.Materials(list(activation_materials)),
        fluxes,
        micro_xs,
        chain_file=str(chain_file),
        normalization_mode="source-rate",
    )
    operator.output_dir = activation_dir
    timesteps, source_flags = _activation_schedule(params)
    source_rates = [leakage_rate * flag for flag in source_flags]
    integrator = openmc.deplete.PredictorIntegrator(
        operator,
        timesteps,
        source_rates=source_rates,
        timestep_units="d",
    )
    _shielding_status(
        "[MOUSE SHIELDING] Activating concrete for "
        f"{float(params.get('Shielding Irradiation Years', 60.0)):.3g} years "
        "using one constant-flux depletion step, followed by "
        f"{float(params.get('Shielding Decay Days', 90.0)):.3g} days of decay.",
        flush=True,
    )
    with _working_directory(activation_dir):
        integrator.integrate()
    _shielding_status(
        "[MOUSE SHIELDING] Concrete activation and decay complete.",
        flush=True,
    )

    results_path = activation_dir / "depletion_results.h5"
    results = openmc.deplete.Results(results_path)
    final_step = results[-1]
    return [
        final_step.get_material(str(material.id))
        for material in activation_materials
    ]


def _photon_dose_rate(
    params: Mapping[str, Any],
    activated_materials: Sequence[openmc.Material],
    geometry_data: Mapping[str, Any],
    output_dir: Path,
) -> tuple[float, float, float]:
    pit_inner_radius = float(geometry_data["pit_inner_radius"])
    world_radius = float(geometry_data["world_radius"])
    pit_bottom = float(geometry_data["pit_bottom"])
    pit_top = float(geometry_data["pit_top"])
    world_bottom = float(geometry_data["world_bottom"])
    world_top = float(geometry_data["world_top"])
    detector_distance = float(geometry_data["detector_distance"])
    detector_radius = float(geometry_data["detector_radius"])

    # The activated depletion materials can contain thousands of vanishingly
    # small trace densities when a complete chain is used.  Exporting those
    # traces to materials.xml can make the OpenMC C++ reader's std::stod call
    # underflow.  Activation changes the bulk concrete composition negligibly,
    # so use clean ordinary concrete for photon attenuation while retaining the
    # depleted materials below solely to construct the decay-photon sources.
    concrete = collect_materials_data(params)["ordinary_concrete"]
    concrete_zones = list(geometry_data["concrete_zones"])
    if len(activated_materials) != len(concrete_zones):
        raise RuntimeError(
            "Activated concrete materials do not match the concrete-zone "
            "geometry."
        )

    pit_inner_cylinder = openmc.ZCylinder(r=pit_inner_radius)
    world_cylinder = openmc.ZCylinder(r=world_radius, boundary_type="vacuum")
    pit_bottom_plane = openmc.ZPlane(z0=pit_bottom)
    pit_top_plane = openmc.ZPlane(z0=pit_top)
    world_bottom_plane = openmc.ZPlane(z0=world_bottom, boundary_type="vacuum")
    world_top_plane = openmc.ZPlane(z0=world_top, boundary_type="vacuum")
    detector_sphere = openmc.Sphere(
        x0=0.0,
        y0=0.0,
        z0=pit_top + detector_distance,
        r=detector_radius,
    )

    central_void_region = (
        -pit_inner_cylinder
        & +pit_bottom_plane
        & -world_top_plane
        & +detector_sphere
    )
    upper_annulus_region = (
        +pit_inner_cylinder
        & -world_cylinder
        & +pit_top_plane
        & -world_top_plane
    )

    detector_cell = openmc.Cell(name="shutdown dose detector", region=-detector_sphere)
    detector_volume = 4.0 / 3.0 * math.pi * detector_radius ** 3
    detector_cell.volume = detector_volume
    central_void_cell = openmc.Cell(
        name="empty pit after reactor removal",
        region=central_void_region,
    )
    upper_annulus_cell = openmc.Cell(
        name="void above concrete side wall",
        region=upper_annulus_region,
    )
    radial_surfaces: dict[float, openmc.ZCylinder] = {
        pit_inner_radius: pit_inner_cylinder,
        world_radius: world_cylinder,
    }
    axial_surfaces: dict[float, openmc.ZPlane] = {
        pit_bottom: pit_bottom_plane,
        world_bottom: world_bottom_plane,
    }
    for zone in concrete_zones:
        if zone["kind"] == "side":
            for radius in (zone["r_inner_cm"], zone["r_outer_cm"]):
                radius = float(radius)
                radial_surfaces.setdefault(radius, openmc.ZCylinder(r=radius))
        elif zone["kind"] == "base":
            for height in (zone["z_lower_cm"], zone["z_upper_cm"]):
                height = float(height)
                axial_surfaces.setdefault(height, openmc.ZPlane(z0=height))
        else:
            raise ValueError(f"Unknown concrete-zone kind: {zone['kind']!r}")

    transport_materials: list[openmc.Material] = []
    concrete_cells: list[openmc.Cell] = []
    spaces = []
    for zone in concrete_zones:
        material = concrete.clone()
        material.name = (
            "ordinary concrete photon transport - "
            f"{zone['kind']} {zone['index']}"
        )
        transport_materials.append(material)
        if zone["kind"] == "side":
            inner_radius = float(zone["r_inner_cm"])
            outer_radius = float(zone["r_outer_cm"])
            region = (
                +radial_surfaces[inner_radius]
                & -radial_surfaces[outer_radius]
                & +world_bottom_plane
                & -pit_top_plane
            )
            space = openmc.stats.CylindricalIndependent(
                r=openmc.stats.PowerLaw(inner_radius, outer_radius, 1.0),
                phi=openmc.stats.Uniform(0.0, 2.0 * math.pi),
                z=openmc.stats.Uniform(world_bottom, pit_top),
            )
        else:
            z_lower = float(zone["z_lower_cm"])
            z_upper = float(zone["z_upper_cm"])
            region = (
                -pit_inner_cylinder
                & +axial_surfaces[z_lower]
                & -axial_surfaces[z_upper]
            )
            space = openmc.stats.CylindricalIndependent(
                r=openmc.stats.PowerLaw(0.0, pit_inner_radius, 1.0),
                phi=openmc.stats.Uniform(0.0, 2.0 * math.pi),
                z=openmc.stats.Uniform(z_lower, z_upper),
            )
        concrete_cells.append(
            openmc.Cell(
                name=(
                    "activated concrete photon region - "
                    f"{zone['kind']} {zone['index']}"
                ),
                fill=material,
                region=region,
            )
        )
        spaces.append(space)

    geometry = openmc.Geometry(
        openmc.Universe(
            cells=[
                *concrete_cells,
                detector_cell,
                central_void_cell,
                upper_annulus_cell,
            ]
        )
    )

    source_cls = _source_class()
    sources = []
    strengths = []
    for material, space in zip(activated_materials, spaces):
        photon_energy = material.get_decay_photon_energy()
        if photon_energy is None:
            continue
        strength = float(photon_energy.integral())
        if strength <= 0.0:
            continue
        sources.append(
            source_cls(
                space=space,
                angle=openmc.stats.Isotropic(),
                energy=photon_energy,
                particle="photon",
                strength=strength,
            )
        )
        strengths.append(strength)
    if not sources:
        raise RuntimeError(
            "The simplified depletion chain produced no decay-photon source. "
            "It must retain decay photon spectra for the activated concrete "
            "nuclides."
        )

    photon_batches = int(
        params.get(
            "Shielding Photon Batches",
            params.get("Shielding Batches", 10),
        )
    )
    photon_particles = int(params.get("Shielding Photon Particles", 50_000))
    retry_multiplier = int(params.get("Shielding Photon Retry Multiplier", 4))
    maximum_retries = int(params.get("Shielding Photon Maximum Retries", 1))
    if photon_batches <= 0 or photon_particles <= 0:
        raise ValueError("Shielding photon batches and particles must be positive.")
    if retry_multiplier < 2:
        raise ValueError("Shielding Photon Retry Multiplier must be at least 2.")
    if maximum_retries < 0:
        raise ValueError("Shielding Photon Maximum Retries cannot be negative.")

    settings = _fixed_source_settings(
        params,
        sources,
        photon_transport=True,
        batches=photon_batches,
        particles=photon_particles,
    )
    dose_energy, dose_coefficients = openmc.data.dose_coefficients(
        "photon", "AP"
    )
    dose_tally = openmc.Tally(name="shutdown_photon_effective_dose")
    dose_tally.filters = [
        openmc.CellFilter(detector_cell),
        openmc.ParticleFilter(["photon"]),
        openmc.EnergyFunctionFilter(
            dose_energy,
            dose_coefficients,
            interpolation="cubic",
        ),
    ]
    dose_tally.scores = ["flux"]
    model = openmc.Model(
        geometry=geometry,
        materials=_unique_materials(transport_materials),
        settings=settings,
        tallies=openmc.Tallies([dose_tally]),
    )
    model.materials.cross_sections = params["cross_sections_xml_location"]

    total_photon_rate = float(sum(strengths))
    # The decay-energy distributions carry their absolute intensities and the
    # fixed-source tally is therefore already weighted by the decay-photon
    # source rate.  Convert the volume-integrated pSv/s response to mrem/h;
    # multiplying by total_photon_rate here would count the source rate twice.
    scale = PSV_PER_SECOND_TO_MREM_PER_HOUR / detector_volume
    for attempt in range(maximum_retries + 1):
        attempt_particles = photon_particles * retry_multiplier ** attempt
        model.settings.particles = attempt_particles
        photon_dir_name = (
            "photon_transport"
            if attempt == 0
            else f"photon_transport_retry_{attempt}"
        )
        photon_dir = output_dir / photon_dir_name
        _shielding_status(
            "[MOUSE SHIELDING] Running shutdown-photon dose transport with "
            f"{total_photon_rate:.5e} photon/s, {attempt_particles} "
            f"particles/batch, and {photon_batches} batches"
            + ("..." if attempt == 0 else f" (retry {attempt})..."),
            flush=True,
        )
        with _working_directory(photon_dir):
            model.run()
        _shielding_status("[MOUSE SHIELDING] Photon dose transport complete.")
        with openmc.StatePoint(_latest_statepoint(photon_dir)) as statepoint:
            tally = statepoint.get_tally(name="shutdown_photon_effective_dose")
            dose_per_source = float(np.asarray(tally.mean, dtype=float).sum())
            dose_std_per_source = float(
                np.sqrt(np.square(np.asarray(tally.std_dev, dtype=float)).sum())
            )

        dose = dose_per_source * scale
        dose_std = dose_std_per_source * scale
        if dose > 0.0:
            return dose, dose_std, total_photon_rate
        if attempt < maximum_retries:
            _shielding_status(
                "[MOUSE SHIELDING] Photon tally was zero; retrying only the "
                f"photon calculation with {retry_multiplier}x more particles.",
                flush=True,
            )

    _shielding_status(
        "[MOUSE SHIELDING] WARNING: The photon calculation still scored "
        "zero dose after "
        f"{maximum_retries + 1} attempt(s); the final attempt used "
        f"{model.settings.particles} particles/batch and {photon_batches} "
        "batches. The zero-score point will be marked unresolved and excluded "
        "from the conservative exponential tail extrapolation.",
        flush=True,
    )
    return float("nan"), float("nan"), total_photon_rate


def evaluate_shielding_candidate(
    params: Mapping[str, Any],
    leakage_source: Mapping[str, Any],
    thickness_cm: float,
    output_dir: str | os.PathLike[str],
) -> dict[str, float]:
    """Evaluate one WEP thickness with neutron activation and photon transport."""
    thickness_cm = float(thickness_cm)
    if thickness_cm < 0.0:
        raise ValueError("Shield thickness cannot be negative.")
    candidate_dir = Path(output_dir).resolve()
    candidate_dir.mkdir(parents=True, exist_ok=True)
    model, cells, materials, geometry_data = _make_candidate_neutron_model(
        params,
        leakage_source,
        thickness_cm,
    )
    activated_materials = _activate_concrete(
        params,
        model,
        cells,
        materials,
        float(leakage_source["source_rate_n_per_s"]),
        candidate_dir,
    )
    dose, dose_std, photon_rate = _photon_dose_rate(
        params,
        activated_materials,
        geometry_data,
        candidate_dir,
    )
    result = {
        "thickness_cm": thickness_cm,
        "dose_mrem_per_h": dose,
        "dose_std_dev_mrem_per_h": dose_std,
        "decay_photon_rate_per_s": photon_rate,
        "dose_resolved": bool(np.isfinite(dose) and dose > 0.0),
    }
    (candidate_dir / "result.json").write_text(
        json.dumps(result, indent=2) + "\n",
        encoding="utf-8",
    )
    return result


def _serializable_shielding_inputs(params: Mapping[str, Any]) -> dict[str, Any]:
    keys = [
        "reactor type",
        "Common Temperature",
        "cross_sections_xml_location",
        "simplified_chain_thermal_xml",
        "Out Of Vessel Shield Material",
        "Out Of Vessel Shield Effective Density Factor",
        "WEP Nominal Density",
        "Vessels Total Height",
        "Shielding Batches",
        "Shielding Particles",
        "Shielding Photon Batches",
        "Shielding Photon Particles",
        "Shielding Photon Retry Multiplier",
        "Shielding Photon Maximum Retries",
        "Shielding Dose Limit",
        "Shielding Minimum Thickness",
        "Shielding Initial Upper Thickness",
        "Shielding Maximum Thickness",
        "Shielding Thickness Growth Factor",
        "Shielding Thickness Tolerance",
        "Shielding Confirmation Increment",
        "Shielding Maximum Search Iterations",
        "Shielding Extrapolation Confidence Multiplier",
        "Shielding Working Directory",
        "Shielding Temperature",
        "Shielding Concrete Thickness",
        "Shielding Concrete Activation Zone Boundaries",
        "Shielding Pit Clearance",
        "Shielding Detector Distance",
        "Shielding Detector Radius",
        "Shielding Irradiation Years",
        "Shielding Activation Step Days",
        "Shielding Decay Days",
        "Enrichment",
        "UO2 atom fraction",
        "H_Zr_ratio",
        "U_met_wo",
        "er_wo",
        "B10_at_frac_B",
    ]
    result = {}
    for key in keys:
        if key not in params:
            continue
        value = params[key]
        if isinstance(value, (str, int, float, bool)):
            result[key] = value
        elif isinstance(value, (list, tuple)) and all(
            isinstance(item, (str, int, float, bool)) for item in value
        ):
            result[key] = list(value)
    return result


def _exponential_crossing_thickness(
    lower_result: Mapping[str, float],
    upper_result: Mapping[str, float],
    target: float,
) -> float:
    """Interpolate a dose-limit crossing linearly in log-dose space."""
    lower_thickness = float(lower_result["thickness_cm"])
    upper_thickness = float(upper_result["thickness_cm"])
    lower_dose = float(lower_result["dose_mrem_per_h"])
    upper_dose = float(upper_result["dose_mrem_per_h"])
    if not (
        lower_thickness < upper_thickness
        and lower_dose > target
        and 0.0 < upper_dose <= target
    ):
        raise ValueError(
            "Exponential interpolation requires a thinner failing point and "
            "a thicker passing point with positive doses."
        )
    log_slope = (
        math.log(upper_dose) - math.log(lower_dose)
    ) / (upper_thickness - lower_thickness)
    if not np.isfinite(log_slope) or log_slope >= 0.0:
        raise ValueError("The local dose/thickness slope is not decreasing.")
    predicted = lower_thickness + (
        math.log(target) - math.log(lower_dose)
    ) / log_slope
    return float(np.clip(predicted, lower_thickness, upper_thickness))


def _conservative_exponential_tail_extrapolation(
    results: Sequence[Mapping[str, Any]],
    target: float,
    confidence_multiplier: float = 1.645,
) -> dict[str, Any]:
    """Extrapolate a dose crossing from the last two resolved tail points.

    Zero-score and otherwise unresolved points are excluded.  Each resolved
    dose is replaced by ``mean + confidence_multiplier * standard deviation``
    before fitting in log-dose space.  This makes the scoping estimate more
    conservative than extrapolating the raw Monte Carlo means.
    """
    if target <= 0.0:
        raise ValueError("The dose target must be positive.")
    if confidence_multiplier < 0.0:
        raise ValueError(
            "Shielding Extrapolation Confidence Multiplier cannot be negative."
        )

    valid_points = []
    for result in results:
        thickness = float(result["thickness_cm"])
        dose = float(result["dose_mrem_per_h"])
        dose_std = float(result["dose_std_dev_mrem_per_h"])
        if not (
            np.isfinite(thickness)
            and np.isfinite(dose)
            and dose > 0.0
            and np.isfinite(dose_std)
            and dose_std >= 0.0
        ):
            continue
        conservative_dose = dose + confidence_multiplier * dose_std
        if not np.isfinite(conservative_dose) or conservative_dose <= 0.0:
            continue
        valid_points.append(
            {
                "thickness_cm": thickness,
                "dose_mrem_per_h": dose,
                "dose_std_dev_mrem_per_h": dose_std,
                "conservative_dose_mrem_per_h": conservative_dose,
            }
        )

    valid_points.sort(key=lambda point: point["thickness_cm"])
    if len(valid_points) < 2:
        raise ValueError(
            "At least two resolved positive-dose points are required for "
            "exponential tail extrapolation."
        )

    thinner, thicker = valid_points[-2:]
    thinner_thickness = thinner["thickness_cm"]
    thicker_thickness = thicker["thickness_cm"]
    thinner_dose = thinner["conservative_dose_mrem_per_h"]
    thicker_dose = thicker["conservative_dose_mrem_per_h"]
    if thicker_thickness <= thinner_thickness:
        raise ValueError("Tail-fit thicknesses must be strictly increasing.")
    if thicker_dose <= target:
        raise ValueError(
            "The final resolved conservative dose already meets the target; "
            "a fail/pass interpolation should be used instead."
        )

    log_slope = (
        math.log(thicker_dose) - math.log(thinner_dose)
    ) / (thicker_thickness - thinner_thickness)
    if not np.isfinite(log_slope) or log_slope >= 0.0:
        raise ValueError(
            "The last two resolved conservative dose points do not have a "
            "decreasing exponential tail."
        )

    predicted = thicker_thickness + (
        math.log(target) - math.log(thicker_dose)
    ) / log_slope
    if not np.isfinite(predicted) or predicted <= thicker_thickness:
        raise ValueError(
            "The exponential tail fit did not predict a thicker dose crossing."
        )

    return {
        "predicted_thickness_cm": float(predicted),
        "log_slope_per_cm": float(log_slope),
        "confidence_multiplier": float(confidence_multiplier),
        "points_used": [thinner, thicker],
    }


def _round_up_to_increment(value: float, increment: float) -> float:
    """Round a nonnegative thickness upward to a practical increment."""
    if increment <= 0.0:
        raise ValueError("Shielding Confirmation Increment must be positive.")
    return math.ceil((float(value) - 1.0e-12) / increment) * increment


def estimate_dynamic_shielding(
    params: MutableMapping[str, Any],
    statepoint_file: str | os.PathLike[str],
) -> dict[str, Any]:
    """Estimate and store the minimum shield thickness meeting the dose target."""
    _require_supported_openmc()
    prepare_shielding_geometry(params)
    concrete_depth_edges = _concrete_activation_depth_edges(params)
    params["Shielding Concrete Activation Zone Boundaries"] = (
        concrete_depth_edges[1:-1]
    )
    if "Shielding Confirmation Increment" not in params:
        params["Shielding Confirmation Increment"] = 0.5 * float(
            params.get("Shielding Thickness Tolerance", 10.0)
        )
    root = _create_unique_run_directory(
        str(params.get("Shielding Working Directory", "shielding_runs"))
    )

    _shielding_status("\n" + "=" * 78)
    _shielding_status(
        "[MOUSE SHIELDING] Starting early design-dependent shielding estimate"
    )
    _shielding_status(
        "[MOUSE SHIELDING] Implementation version: "
        f"{SHIELDING_IMPLEMENTATION_VERSION}"
    )
    _shielding_status(
        "[MOUSE SHIELDING] Concrete activation depth zones: "
        + ", ".join(
            f"{inner:g}-{outer:g} cm"
            for inner, outer in zip(
                concrete_depth_edges[:-1],
                concrete_depth_edges[1:],
            )
        )
        + " for both the side wall and base"
    )
    _shielding_status("[MOUSE SHIELDING] Fuel depletion has NOT started yet")
    _shielding_status("=" * 78)

    core_source = extract_bol_core_leakage(params, statepoint_file)
    _shielding_status(
        "[MOUSE SHIELDING] BOL core leakage source: "
        f"{core_source['source_rate_n_per_s']:.5e} n/s"
    )
    _shielding_status(
        "[MOUSE SHIELDING] Leakage normalization: radial "
        f"{100.0 * core_source['radial_leakage_fraction_per_source']:.3f}%, "
        "axial non-leakage "
        f"{core_source['axial_non_leakage_probability']:.5f}, combined total "
        f"{100.0 * core_source['total_leakage_fraction_per_source']:.3f}%"
    )
    leakage_source = transport_source_to_intake_vessel(
        params,
        core_source,
        root / "vessel_transfer",
    )
    _shielding_status(
        "[MOUSE SHIELDING] Source immediately outside intake vessel: "
        f"{leakage_source['source_rate_n_per_s']:.5e} n/s"
    )

    source_record = {
        "core_source": core_source,
        "outside_intake_source": leakage_source,
    }
    (root / "shielding_leakage_source.json").write_text(
        json.dumps(source_record, indent=2) + "\n",
        encoding="utf-8",
    )
    (root / "shielding_inputs.json").write_text(
        json.dumps(_serializable_shielding_inputs(params), indent=2) + "\n",
        encoding="utf-8",
    )

    target = float(params.get("Shielding Dose Limit", 0.5))
    tolerance = float(params.get("Shielding Thickness Tolerance", 10.0))
    minimum = max(0.0, float(params.get("Shielding Minimum Thickness", 0.0)))
    initial_upper = max(
        minimum + tolerance,
        float(params.get("Shielding Initial Upper Thickness", 50.0)),
    )
    maximum = float(params.get("Shielding Maximum Thickness", 300.0))
    growth_factor = float(params.get("Shielding Thickness Growth Factor", 1.5))
    configured_refinements = int(
        params.get("Shielding Maximum Search Iterations", 3)
    )
    pre_fit_midpoint_limit = min(configured_refinements, 1)
    confirmation_increment = float(
        params.get("Shielding Confirmation Increment", 0.5 * tolerance)
    )
    if target <= 0.0:
        raise ValueError("Shielding Dose Limit must be positive.")
    if tolerance <= 0.0:
        raise ValueError("Shielding Thickness Tolerance must be positive.")
    if maximum < initial_upper:
        raise ValueError(
            "Shielding Maximum Thickness must be at least the initial upper "
            "thickness."
        )
    if growth_factor <= 1.0:
        raise ValueError("Shielding Thickness Growth Factor must exceed 1.0.")
    if configured_refinements < 0:
        raise ValueError("Shielding Maximum Search Iterations cannot be negative.")
    if confirmation_increment <= 0.0:
        raise ValueError("Shielding Confirmation Increment must be positive.")
    extrapolation_confidence_multiplier = float(
        params.get("Shielding Extrapolation Confidence Multiplier", 1.645)
    )
    if extrapolation_confidence_multiplier < 0.0:
        raise ValueError(
            "Shielding Extrapolation Confidence Multiplier cannot be negative."
        )
    evaluated: dict[float, dict[str, Any]] = {}

    def evaluate(value: float) -> dict[str, Any]:
        rounded = round(float(value), 5)
        if rounded not in evaluated:
            _shielding_status(
                f"[MOUSE SHIELDING] Evaluating {rounded:.2f} cm "
                f"{params.get('Out Of Vessel Shield Material', 'WEP')}..."
            )
            result = evaluate_shielding_candidate(
                params,
                leakage_source,
                rounded,
                root / f"candidate_{rounded:.2f}_cm",
            )
            result_dose = float(result["dose_mrem_per_h"])
            result_std = float(result["dose_std_dev_mrem_per_h"])
            thinner_values = [
                value
                for value in evaluated
                if value < rounded
                and np.isfinite(evaluated[value]["dose_mrem_per_h"])
            ]
            if thinner_values:
                nearest_thinner = max(thinner_values)
                thinner_dose = evaluated[nearest_thinner]["dose_mrem_per_h"]
                if np.isfinite(result_dose) and result_dose > thinner_dose:
                    _shielding_status(
                        "[MOUSE SHIELDING] WARNING: Mean dose increased from "
                        f"{thinner_dose:.5g} mrem/h at {nearest_thinner:.2f} cm "
                        f"to {result['dose_mrem_per_h']:.5g} mrem/h at "
                        f"{rounded:.2f} cm. This is likely Monte Carlo noise; "
                        "the scoping calculation will continue without an "
                        "uncertainty-driven rerun."
                    )
            thicker_values = [
                value
                for value in evaluated
                if value > rounded
                and np.isfinite(evaluated[value]["dose_mrem_per_h"])
            ]
            if thicker_values:
                nearest_thicker = min(thicker_values)
                thicker_dose = evaluated[nearest_thicker]["dose_mrem_per_h"]
                if np.isfinite(result_dose) and thicker_dose > result_dose:
                    _shielding_status(
                        "[MOUSE SHIELDING] WARNING: Mean dose increases from "
                        f"{result['dose_mrem_per_h']:.5g} mrem/h at "
                        f"{rounded:.2f} cm to {thicker_dose:.5g} mrem/h at "
                        f"{nearest_thicker:.2f} cm. This is likely Monte Carlo "
                        "noise; the raw values are retained, while the fitted "
                        "search uses only its local fail/pass bracket."
                    )
            evaluated[rounded] = result
            if np.isfinite(result_dose):
                _shielding_status(
                    f"[MOUSE SHIELDING] {rounded:.2f} cm -> "
                    f"{result_dose:.5g} +/- {result_std:.2g} mrem/h"
                )
            else:
                _shielding_status(
                    f"[MOUSE SHIELDING] {rounded:.2f} cm -> unresolved "
                    "zero-score photon tally; excluded from fitting"
                )
        return evaluated[rounded]

    fit_prediction = None
    extrapolation_details = None
    extrapolated_prediction = None
    extrapolation_error = None
    thickness_extrapolated = False
    estimated_criterion_met = False
    unresolved_thickness = None
    confirmation_thickness = None
    confirmation_performed = False
    pre_fit_midpoint_evaluations = 0

    # Start at the expected neighborhood.  The zero/minimum-shield case is
    # evaluated only when the initial point already passes and a lower bracket
    # is therefore needed.  This avoids one full activation calculation in the
    # common case where the initial guess is too thin.
    trial_thickness = initial_upper
    trial_result = evaluate(trial_thickness)
    trial_dose = float(trial_result["dose_mrem_per_h"])
    if not np.isfinite(trial_dose):
        criterion_met = False
        selected = trial_result
        low_result = None
        high_result = None
        unresolved_thickness = trial_thickness
    elif trial_dose <= target:
        if trial_thickness == minimum:
            criterion_met = True
            selected = trial_result
            low_result = None
            high_result = trial_result
        else:
            minimum_result = evaluate(minimum)
            minimum_dose = float(minimum_result["dose_mrem_per_h"])
            if np.isfinite(minimum_dose) and minimum_dose <= target:
                criterion_met = True
                selected = minimum_result
                low_result = None
                high_result = minimum_result
            elif np.isfinite(minimum_dose):
                criterion_met = True
                low_result = minimum_result
                high_result = trial_result
                selected = high_result
            else:
                criterion_met = False
                selected = minimum_result
                low_result = None
                high_result = None
                unresolved_thickness = minimum
    else:
        low_result = trial_result
        high_result = None
        while trial_thickness < maximum:
            next_thickness = min(
                maximum,
                max(
                    trial_thickness + tolerance,
                    trial_thickness * growth_factor,
                ),
            )
            trial_thickness = next_thickness
            trial_result = evaluate(trial_thickness)
            trial_dose = float(trial_result["dose_mrem_per_h"])
            if not np.isfinite(trial_dose):
                unresolved_thickness = trial_thickness
                break
            if trial_dose <= target:
                high_result = trial_result
                break
            low_result = trial_result

        criterion_met = high_result is not None
        selected = high_result if criterion_met else low_result

    if criterion_met and low_result is not None:
        lower = float(low_result["thickness_cm"])
        upper = float(high_result["thickness_cm"])

        # One midpoint is allowed to stabilize a wide bracket before fitting.
        # This is deliberately capped at one regardless of the legacy input so
        # the new search cannot revert to a long bisection sequence.
        if pre_fit_midpoint_limit and upper - lower > tolerance:
            midpoint = 0.5 * (lower + upper)
            midpoint_result = evaluate(midpoint)
            pre_fit_midpoint_evaluations = 1
            midpoint_dose = float(midpoint_result["dose_mrem_per_h"])
            if not np.isfinite(midpoint_dose):
                _shielding_status(
                    "[MOUSE SHIELDING] Midpoint photon tally was unresolved; "
                    "retaining the existing resolved fail/pass bracket."
                )
            elif midpoint_dose <= target:
                high_result = midpoint_result
            else:
                low_result = midpoint_result
            lower = float(low_result["thickness_cm"])
            upper = float(high_result["thickness_cm"])

        try:
            fit_prediction = _exponential_crossing_thickness(
                low_result,
                high_result,
                target,
            )
            proposed = _round_up_to_increment(
                fit_prediction,
                confirmation_increment,
            )
            proposed = float(np.clip(proposed, lower, upper))
            _shielding_status(
                "[MOUSE SHIELDING] Local exponential interpolation predicts "
                f"{fit_prediction:.2f} cm; bounded confirmation candidate is "
                f"{proposed:.2f} cm."
            )
        except ValueError as error:
            proposed = 0.5 * (lower + upper)
            _shielding_status(
                "[MOUSE SHIELDING] Exponential interpolation unavailable; "
                f"using one midpoint confirmation ({error})."
            )

        proposed_key = round(proposed, 5)
        if proposed_key not in evaluated and lower < proposed < upper:
            confirmation_thickness = proposed
            confirmation_performed = True
            confirmation_result = evaluate(proposed)
            if confirmation_result["dose_mrem_per_h"] <= target:
                selected = confirmation_result
            else:
                # Stop after one confirmation.  The existing upper point is
                # already known to pass and is retained without another run.
                selected = high_result
                _shielding_status(
                    "[MOUSE SHIELDING] The single confirmation point failed; "
                    "selecting the existing passing upper bracket without "
                    "another calculation."
                )
        else:
            selected = high_result
            _shielding_status(
                "[MOUSE SHIELDING] The fitted/rounded value coincides with "
                "an existing bracket point; no confirmation run is needed."
            )

    resolved_failing_results = [
        result
        for result in evaluated.values()
        if np.isfinite(result["dose_mrem_per_h"])
        and result["dose_mrem_per_h"] > target
    ]
    last_direct_failure_thickness = (
        max(result["thickness_cm"] for result in resolved_failing_results)
        if resolved_failing_results
        else None
    )

    if not criterion_met:
        try:
            extrapolation_details = _conservative_exponential_tail_extrapolation(
                list(evaluated.values()),
                target,
                confidence_multiplier=extrapolation_confidence_multiplier,
            )
            extrapolated_prediction = float(
                extrapolation_details["predicted_thickness_cm"]
            )
            applied_thickness = min(
                maximum,
                _round_up_to_increment(
                    extrapolated_prediction,
                    confirmation_increment,
                ),
            )
            tail_point = extrapolation_details["points_used"][-1]
            fitted_dose = tail_point[
                "conservative_dose_mrem_per_h"
            ] * math.exp(
                extrapolation_details["log_slope_per_cm"]
                * (applied_thickness - tail_point["thickness_cm"])
            )
            thickness_extrapolated = True
            estimated_criterion_met = extrapolated_prediction <= maximum
            selected = {
                "thickness_cm": float(applied_thickness),
                "dose_mrem_per_h": float(fitted_dose),
                "dose_std_dev_mrem_per_h": float("nan"),
                "decay_photon_rate_per_s": float("nan"),
                "dose_resolved": False,
                "dose_estimated": True,
            }
            if estimated_criterion_met:
                _shielding_status(
                    "[MOUSE SHIELDING] No directly scored passing point was "
                    "available. Conservative exponential extrapolation of the "
                    "last two resolved positive-dose points predicts "
                    f"{extrapolated_prediction:.2f} cm; the applied cost-model "
                    f"thickness is {applied_thickness:.2f} cm."
                )
            else:
                _shielding_status(
                    "[MOUSE SHIELDING] Conservative exponential extrapolation "
                    f"predicts {extrapolated_prediction:.2f} cm, above the "
                    f"configured {maximum:.2f} cm limit. The configured maximum "
                    "is retained and the estimated criterion remains unmet."
                )
        except ValueError as error:
            extrapolation_error = str(error)
            selected = {
                "thickness_cm": float(maximum),
                "dose_mrem_per_h": float("nan"),
                "dose_std_dev_mrem_per_h": float("nan"),
                "decay_photon_rate_per_s": float("nan"),
                "dose_resolved": False,
                "dose_estimated": False,
            }
            _shielding_status(
                "[MOUSE SHIELDING] WARNING: Conservative exponential "
                f"extrapolation was unavailable ({error}). The configured "
                f"maximum thickness of {maximum:.2f} cm is retained and the "
                "dose criterion is not claimed."
            )

    params["Out Of Vessel Shield Thickness"] = selected["thickness_cm"]
    params["Estimated Shutdown Dose Rate (mrem/h)"] = selected[
        "dose_mrem_per_h"
    ]
    params["Estimated Shutdown Dose Rate Std Dev (mrem/h)"] = selected[
        "dose_std_dev_mrem_per_h"
    ]
    params["Shielding Dose Criterion Met"] = criterion_met
    params["Shielding Dose Criterion Estimated Met"] = bool(
        criterion_met or estimated_criterion_met
    )
    params["Shielding Thickness Extrapolated"] = thickness_extrapolated
    params["Shielding Extrapolated Required Thickness (cm)"] = (
        extrapolated_prediction
        if extrapolated_prediction is not None
        else float("nan")
    )
    params["Shielding Calculation Directory"] = str(root)

    summary = {
        "selected": selected,
        "dose_limit_mrem_per_h": target,
        "criterion_met": criterion_met,
        "criterion_estimated_met": bool(
            criterion_met or estimated_criterion_met
        ),
        "thickness_extrapolated": thickness_extrapolated,
        "required_thickness_found": criterion_met,
        "required_thickness_estimated": estimated_criterion_met,
        "required_thickness_lower_bound_cm": (
            None if criterion_met else last_direct_failure_thickness
        ),
        "search_settings": {
            "minimum_thickness_cm": float(
                params.get("Shielding Minimum Thickness", 0.0)
            ),
            "initial_upper_thickness_cm": float(
                params.get("Shielding Initial Upper Thickness", 50.0)
            ),
            "maximum_thickness_cm": maximum,
            "growth_factor": growth_factor,
            "thickness_tolerance_cm": tolerance,
            "confirmation_increment_cm": confirmation_increment,
            "configured_maximum_search_iterations": configured_refinements,
            "effective_pre_fit_midpoint_limit": pre_fit_midpoint_limit,
            "maximum_fitted_confirmation_calculations": 1,
            "search_method": (
                "geometric bracket, at most one midpoint, local exponential "
                "interpolation, at most one confirmation; unresolved zero-score "
                "points are excluded and trigger conservative exponential tail "
                "extrapolation from the last two resolved points"
            ),
            "selection_uses_mean_dose": True,
            "uncertainty_driven_retries": False,
            "tail_extrapolation_confidence_multiplier": (
                extrapolation_confidence_multiplier
            ),
        },
        "search_diagnostics": {
            "pre_fit_midpoint_evaluations": pre_fit_midpoint_evaluations,
            "exponential_prediction_cm": fit_prediction,
            "confirmation_performed": confirmation_performed,
            "confirmation_thickness_cm": confirmation_thickness,
            "unresolved_zero_score_thickness_cm": unresolved_thickness,
            "extrapolated_prediction_cm": extrapolated_prediction,
            "extrapolation_details": extrapolation_details,
            "extrapolation_error": extrapolation_error,
        },
        "evaluations": [evaluated[key] for key in sorted(evaluated)],
        "calculation_directory": str(root),
    }
    (root / "shielding_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    if criterion_met:
        status = "MEETS"
    elif estimated_criterion_met:
        status = "IS ESTIMATED TO MEET"
    else:
        status = "DOES NOT MEET"
    _shielding_status("-" * 78)
    if criterion_met:
        _shielding_status(
            "[MOUSE SHIELDING] Estimated required thickness: "
            f"{selected['thickness_cm']:.2f} cm "
            f"({selected['thickness_cm'] / 2.54:.2f} in)"
        )
    elif thickness_extrapolated and estimated_criterion_met:
        _shielding_status(
            "[MOUSE SHIELDING] Conservative extrapolated required thickness: "
            f"{selected['thickness_cm']:.2f} cm "
            f"({selected['thickness_cm'] / 2.54:.2f} in). This thickness was "
            "not directly verified by a nonzero photon tally."
        )
    elif thickness_extrapolated:
        _shielding_status(
            "[MOUSE SHIELDING] Conservative extrapolated required thickness "
            f"({extrapolated_prediction:.2f} cm) exceeds the configured "
            f"{maximum:.2f} cm limit. The configured maximum is retained for "
            "costing and the dose criterion remains unmet."
        )
    else:
        _shielding_status(
            "[MOUSE SHIELDING] No directly passing or defensibly extrapolated "
            f"thickness was found through {maximum:.2f} cm. The configured "
            "maximum is retained for costing and the dose criterion is not "
            "claimed."
        )
    _shielding_status(
        "[MOUSE SHIELDING] Estimated shutdown dose: "
        f"{selected['dose_mrem_per_h']:.5g} mrem/h; {status} "
        f"the {target:g} mrem/h criterion"
    )
    _shielding_status(f"[MOUSE SHIELDING] Detailed outputs: {root}")
    _shielding_status("[MOUSE SHIELDING] Main fuel depletion will start next")
    _shielding_status("=" * 78 + "\n")
    return summary
