# Copyright 2026 Battelle Energy Alliance, LLC
# Released under the MIT License.
"""Screening calculations for transport of an irradiated microreactor module.

The runtime model deliberately avoids OpenMC and depletion calculations. It uses:

* the MOUSE full-power fuel lifetime;
* a finite-irradiation Way-Wigner decay-heat approximation;
* a gamma-energy fraction and representative photon energy documented by ANL;
* NIST photon attenuation and air energy-absorption coefficients; and
* a closed cylindrical external transport shield around the MOUSE reactor module.

This is an early-design screening model. It does not credit attenuation from the
fuel, reflector, vessels, coolant, or existing reactor shielding and does not yet
include activation gamma rays, shutdown neutrons, photon buildup, penetrations,
impact limiters, or a certified transport-package structure.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
import math
from typing import Dict, Tuple

import pandas as pd

_REPO_ROOT = Path(__file__).resolve().parent.parent
_IRR_DIR = _REPO_ROOT / "assets" / "irradiated_transport"

_GAMMA_FRACTION = 0.50
_REPRESENTATIVE_PHOTON_ENERGY_MEV = 0.70
_AIR_MASS_ENERGY_ABSORPTION_CM2_G = 0.02915
_WAY_WIGNER_COEFFICIENT = 0.066
_DAYS_PER_MONTH = 30.0


@lru_cache(maxsize=1)
def load_shield_material_properties() -> pd.DataFrame:
    """Load and validate the candidate shield-material screening inputs."""
    path = _IRR_DIR / "shield_material_properties.csv"
    df = pd.read_csv(path)
    required = {
        "material",
        "density_kg_m3",
        "raw_material_cost_2025_usd_per_kg",
        "mass_attenuation_cm2_g_at_0p7_mev",
        "screening_scope",
        "notes",
    }
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(
            f"{path.name} is missing columns: {', '.join(sorted(missing))}"
        )
    numeric_columns = (
        "density_kg_m3",
        "raw_material_cost_2025_usd_per_kg",
        "mass_attenuation_cm2_g_at_0p7_mev",
    )
    for column in numeric_columns:
        df[column] = pd.to_numeric(df[column], errors="raise")
        if (df[column] <= 0.0).any():
            raise ValueError(f"{path.name} contains a nonpositive {column}.")
    return df


@lru_cache(maxsize=1)
def load_irradiated_cost_inputs() -> pd.DataFrame:
    return pd.read_csv(
        _IRR_DIR / "irradiated_transport_cost_inputs_2025.csv"
    ).set_index("id", drop=False)


def available_shield_materials() -> Tuple[str, ...]:
    return tuple(load_shield_material_properties()["material"].astype(str))


def calculate_decay_heat_w(
    thermal_power_mwt: float,
    fuel_lifetime_days: float,
    cooldown_months: float,
) -> float:
    """Finite-irradiation Way-Wigner screening estimate of decay heat."""
    power_w = max(float(thermal_power_mwt), 0.0) * 1.0e6
    irradiation_s = max(float(fuel_lifetime_days), 0.0) * 86400.0
    cooldown_s = max(float(cooldown_months), 1.0e-6) * _DAYS_PER_MONTH * 86400.0
    if power_w <= 0 or irradiation_s <= 0:
        return 0.0
    return float(
        power_w
        * _WAY_WIGNER_COEFFICIENT
        * (cooldown_s ** -0.2 - (cooldown_s + irradiation_s) ** -0.2)
    )


def _dose_mrem_h(
    thickness_cm: float,
    distance_m: float,
    gamma_power_w: float,
    density_g_cm3: float,
    mass_attenuation_cm2_g: float,
) -> float:
    """Point-source air-kerma approximation for the representative photon."""
    attenuation = math.exp(
        -float(mass_attenuation_cm2_g)
        * float(density_g_cm3)
        * max(float(thickness_cm), 0.0)
    )
    distance = max(float(distance_m), 1.0e-9)
    energy_fluence_w_m2 = (
        max(float(gamma_power_w), 0.0) * attenuation
        / (4.0 * math.pi * distance**2)
    )
    air_mu_en_m2_kg = _AIR_MASS_ENERGY_ABSORPTION_CM2_G * 0.1
    # For photons, 1 Gy is approximated as 1 Sv for this screening conversion.
    return float(energy_fluence_w_m2 * air_mu_en_m2_kg * 3600.0 * 100000.0)


def _solve_thickness_cm(
    target_dose_mrem_h: float,
    distance_m: float,
    gamma_power_w: float,
    density_g_cm3: float,
    mass_attenuation_cm2_g: float,
) -> Tuple[float, float]:
    target = max(float(target_dose_mrem_h), 1.0e-12)
    unshielded = _dose_mrem_h(
        0.0,
        distance_m,
        gamma_power_w,
        density_g_cm3,
        mass_attenuation_cm2_g,
    )
    if unshielded <= target:
        return 0.0, unshielded
    low, high = 0.0, 5.0
    while (
        _dose_mrem_h(
            high,
            distance_m,
            gamma_power_w,
            density_g_cm3,
            mass_attenuation_cm2_g,
        )
        > target
    ):
        high *= 2.0
        if high > 500.0:
            raise ValueError("Required shielding exceeds the 500 cm screening limit.")
    for _ in range(80):
        mid = 0.5 * (low + high)
        if (
            _dose_mrem_h(
                mid,
                distance_m,
                gamma_power_w,
                density_g_cm3,
                mass_attenuation_cm2_g,
            )
            <= target
        ):
            high = mid
        else:
            low = mid
    return float(high), unshielded


def estimate_irradiated_transport_shield(
    reactor_type: str,
    params: Dict[str, object],
    cooldown_months: float,
    target_dose_mrem_h: float,
    dose_distance_m: float,
    shield_material: str,
    module_mass_kg: float,
    module_length_m: float,
    module_width_m: float,
    module_height_m: float,
) -> Dict[str, float | str]:
    """Estimate added shield thickness, mass, dimensions, and raw material cost."""
    reactor_type = str(reactor_type).upper()
    if reactor_type not in {"LTMR", "GCMR", "HPMR"}:
        raise ValueError(f"Unsupported reactor type: {reactor_type}")

    properties = load_shield_material_properties()
    selected = properties[properties["material"].astype(str) == str(shield_material)]
    if selected.empty:
        raise ValueError(f"Unsupported shielding material: {shield_material}")
    material = selected.iloc[0]

    fuel_lifetime_days = float(params.get("Fuel Lifetime", 0.0))
    thermal_power_mwt = float(params.get("Power MWt", 0.0))
    decay_heat_w = calculate_decay_heat_w(
        thermal_power_mwt, fuel_lifetime_days, cooldown_months
    )
    gamma_power_w = _GAMMA_FRACTION * decay_heat_w

    density_kg_m3 = float(material["density_kg_m3"])
    density_g_cm3 = density_kg_m3 / 1000.0
    mass_attenuation_cm2_g = float(
        material["mass_attenuation_cm2_g_at_0p7_mev"]
    )
    thickness_cm, unshielded_dose = _solve_thickness_cm(
        target_dose_mrem_h,
        dose_distance_m,
        gamma_power_w,
        density_g_cm3,
        mass_attenuation_cm2_g,
    )
    thickness_m = thickness_cm / 100.0

    length_m = max(float(module_length_m), 0.001)
    diameter_m = max(float(module_width_m), float(module_height_m), 0.001)
    radius_m = 0.5 * diameter_m
    outer_volume_m3 = math.pi * (radius_m + thickness_m) ** 2 * (
        length_m + 2.0 * thickness_m
    )
    inner_volume_m3 = math.pi * radius_m**2 * length_m
    shield_mass_kg = max(
        0.0, density_kg_m3 * (outer_volume_m3 - inner_volume_m3)
    )
    shield_cost_2025 = shield_mass_kg * float(
        material["raw_material_cost_2025_usd_per_kg"]
    )

    return {
        "shield_material": str(shield_material),
        "material_screening_scope": str(material["screening_scope"]),
        "material_notes": str(material["notes"]),
        "cooldown_months": float(cooldown_months),
        "target_dose_mrem_h": float(target_dose_mrem_h),
        "dose_distance_m": float(dose_distance_m),
        "decay_heat_w": decay_heat_w,
        "gamma_power_w": gamma_power_w,
        "gamma_fraction": _GAMMA_FRACTION,
        "representative_photon_energy_mev": _REPRESENTATIVE_PHOTON_ENERGY_MEV,
        "mass_attenuation_cm2_g": mass_attenuation_cm2_g,
        "unshielded_dose_mrem_h": unshielded_dose,
        "shield_thickness_cm": thickness_cm,
        "shield_mass_kg": shield_mass_kg,
        "shield_raw_material_cost_2025_usd": shield_cost_2025,
        "shielded_module_mass_kg": float(module_mass_kg) + shield_mass_kg,
        "shielded_module_length_m": length_m + 2.0 * thickness_m,
        "shielded_module_width_m": diameter_m + 2.0 * thickness_m,
        "shielded_module_height_m": diameter_m + 2.0 * thickness_m,
    }
