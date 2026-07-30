# Copyright 2026 Battelle Energy Alliance, LLC
# Released under the MIT License.
"""Screening calculations for transport of an irradiated microreactor module.

The runtime model deliberately avoids OpenMC and depletion calculations. It uses:

* the MOUSE full-power fuel lifetime;
* a finite-irradiation Way-Wigner decay-heat approximation;
* a fixed 50% decay-gamma energy fraction;
* a cooldown-dependent normalized 48-group photon spectrum from the Manit
  Shah workbook;
* a point-source dose/lead attenuation treatment evaluated outward from the
  shield surface; and
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

import numpy as np
import pandas as pd

_REPO_ROOT = Path(__file__).resolve().parent.parent
_IRR_DIR = _REPO_ROOT / "assets" / "irradiated_transport"
_J_PER_MEV = 1.602176634e-13
_GAMMA_FRACTION = 0.50
_WAY_WIGNER_COEFFICIENT = 0.066
_DAYS_PER_MONTH = 30.0


@lru_cache(maxsize=1)
def _load_gamma_reference() -> pd.DataFrame:
    path = _IRR_DIR / "decay_gamma_reference.csv"
    df = pd.read_csv(path)
    required = {
        "source_family", "reference_power_mwt", "cooldown_months",
        "energy_group", "energy_mid_mev", "gamma_energy_fraction",
        "lead_mass_attenuation_cm2_g",
    }
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"{path.name} is missing columns: {', '.join(sorted(missing))}")
    return df


@lru_cache(maxsize=1)
def load_shield_material_properties() -> pd.DataFrame:
    return pd.read_csv(_IRR_DIR / "shield_material_properties.csv")


@lru_cache(maxsize=1)
def load_irradiated_cost_inputs() -> pd.DataFrame:
    return pd.read_csv(
        _IRR_DIR / "irradiated_transport_cost_inputs_2025.csv"
    ).set_index("id", drop=False)


def available_shield_materials() -> Tuple[str, ...]:
    return tuple(load_shield_material_properties()["material"].astype(str))


def _reference_spectrum(
    reactor_type: str,
    power_mwt: float,
    cooldown_months: float,
) -> Tuple[pd.DataFrame, str, float, float]:
    """Return a cooldown-dependent normalized spectrum and lead coefficients.

    LTMR and GCMR use the Shah HTPM family as a spectrum-shape analogue;
    HPMR uses the Shah HPMR family. LTMR retains the 10-MWt HTPM shape used in
    the preliminary screening calculations, while GCMR and HPMR use the closest
    available Shah reference power. Spectrum shape is linearly
    interpolated between integer months through month 36 and held at the
    36-month shape for longer cooldown periods.
    """
    family = "HPMR" if reactor_type == "HPMR" else "HTPM"
    df = _load_gamma_reference()
    family_df = df[df["source_family"] == family]
    powers = np.sort(family_df["reference_power_mwt"].unique().astype(float))
    if reactor_type == "LTMR":
        # The preliminary LTMR screening used the 10-MWt HTPM photon-spectrum
        # shape as a fixed thermal-spectrum analogue across LTMR powers.
        ref_power = 10.0
    else:
        ref_power = float(powers[np.argmin(np.abs(powers - float(power_mwt)))])

    m = max(0.0, min(float(cooldown_months), 36.0))
    m0 = int(math.floor(m))
    m1 = int(math.ceil(m))
    alpha = m - m0

    def at_month(month: int) -> pd.DataFrame:
        part = family_df[
            (family_df["reference_power_mwt"] == ref_power)
            & (family_df["cooldown_months"] == month)
        ].sort_values("energy_group")
        if len(part) != 48:
            raise ValueError(
                f"Incomplete {family} reference spectrum at {ref_power} MWt, month {month}."
            )
        return part.reset_index(drop=True)

    s0 = at_month(m0)
    if m1 == m0:
        spectrum = s0.copy()
    else:
        s1 = at_month(m1)
        spectrum = s0.copy()
        spectrum["gamma_energy_fraction"] = (
            (1.0 - alpha) * s0["gamma_energy_fraction"].to_numpy(dtype=float)
            + alpha * s1["gamma_energy_fraction"].to_numpy(dtype=float)
        )
    fractions = spectrum["gamma_energy_fraction"].to_numpy(dtype=float)
    fractions = np.clip(fractions, 0.0, None)
    if fractions.sum() <= 0:
        raise ValueError("Reference gamma spectrum has zero total energy fraction.")
    spectrum["gamma_energy_fraction"] = fractions / fractions.sum()
    return spectrum, family, ref_power, m


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
    spectrum: pd.DataFrame,
    density_g_cm3: float,
) -> float:
    energy = spectrum["energy_mid_mev"].to_numpy(dtype=float)
    energy_fraction = spectrum["gamma_energy_fraction"].to_numpy(dtype=float)
    mu_rho = spectrum["lead_mass_attenuation_cm2_g"].to_numpy(dtype=float)
    photon_rate = gamma_power_w * energy_fraction / (energy * _J_PER_MEV)
    attenuation = np.exp(-mu_rho * density_g_cm3 * max(float(thickness_cm), 0.0))
    # Shah workbook coefficient: 0.14 gives uSv/h; multiply by 0.1 to mrem/h.
    return float(np.sum(
        0.014 * (photon_rate / 1.0e6) * energy * attenuation
        / max(float(distance_m), 1.0e-9) ** 2
    ))


def _solve_thickness_cm(
    target_dose_mrem_h: float,
    distance_from_shield_surface_m: float,
    module_radius_m: float,
    gamma_power_w: float,
    spectrum: pd.DataFrame,
    density_g_cm3: float,
) -> Tuple[float, float, float]:
    """Solve thickness with the detector outside the evolving shield surface."""
    target = max(float(target_dose_mrem_h), 1.0e-12)
    surface_offset_m = max(float(distance_from_shield_surface_m), 0.0)
    inner_radius_m = max(float(module_radius_m), 0.0)

    def dose(thickness_cm: float) -> float:
        source_to_detector_m = (
            inner_radius_m + max(float(thickness_cm), 0.0) / 100.0 + surface_offset_m
        )
        return _dose_mrem_h(
            thickness_cm,
            source_to_detector_m,
            gamma_power_w,
            spectrum,
            density_g_cm3,
        )

    unshielded = dose(0.0)
    if unshielded <= target:
        return 0.0, unshielded, inner_radius_m + surface_offset_m
    low, high = 0.0, 5.0
    while dose(high) > target:
        high *= 2.0
        if high > 500.0:
            raise ValueError("Required shielding exceeds the 500 cm screening limit.")
    for _ in range(80):
        mid = 0.5 * (low + high)
        if dose(mid) <= target:
            high = mid
        else:
            low = mid
    source_to_detector_m = inner_radius_m + high / 100.0 + surface_offset_m
    return float(high), unshielded, float(source_to_detector_m)


def _closed_shell_metrics(
    thickness_cm: float,
    density_kg_m3: float,
    raw_material_cost_usd_per_kg: float,
    module_mass_kg: float,
    module_length_m: float,
    module_diameter_m: float,
) -> Dict[str, float]:
    """Return closed cylindrical shell mass, cost, and outer dimensions."""
    thickness_m = max(float(thickness_cm), 0.0) / 100.0
    length_m = max(float(module_length_m), 0.001)
    diameter_m = max(float(module_diameter_m), 0.001)
    radius_m = 0.5 * diameter_m
    outer_volume_m3 = math.pi * (radius_m + thickness_m) ** 2 * (
        length_m + 2.0 * thickness_m
    )
    inner_volume_m3 = math.pi * radius_m ** 2 * length_m
    shield_mass_kg = max(
        0.0, float(density_kg_m3) * (outer_volume_m3 - inner_volume_m3)
    )
    return {
        "shield_thickness_cm": float(thickness_cm),
        "shield_mass_kg": shield_mass_kg,
        "shield_raw_material_cost_2025_usd": (
            shield_mass_kg * float(raw_material_cost_usd_per_kg)
        ),
        "shielded_module_mass_kg": float(module_mass_kg) + shield_mass_kg,
        "shielded_module_length_m": length_m + 2.0 * thickness_m,
        "shielded_module_width_m": diameter_m + 2.0 * thickness_m,
        "shielded_module_height_m": diameter_m + 2.0 * thickness_m,
    }


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
    """Estimate shielding with the Manit Shah 48-group photon spectrum.

    The requested dose distance is measured radially outward from the outer
    shield surface. The resulting thickness is the transport design basis used
    for shield mass, cost, package dimensions, and transportation screening.
    """
    reactor_type = str(reactor_type).upper()
    properties = load_shield_material_properties()
    selected = properties[properties["material"].astype(str) == str(shield_material)]
    if selected.empty:
        raise ValueError(f"Unsupported shielding material: {shield_material}")
    material = selected.iloc[0]
    if str(shield_material) != "Lead":
        raise ValueError("Only lead is validated in the current screening model.")

    fuel_lifetime_days = float(params.get("Fuel Lifetime", 0.0))
    thermal_power_mwt = float(params.get("Power MWt", 0.0))
    decay_heat_w = calculate_decay_heat_w(
        thermal_power_mwt, fuel_lifetime_days, cooldown_months
    )
    gamma_power_w = _GAMMA_FRACTION * decay_heat_w
    spectrum, source_family, ref_power, used_spectrum_month = _reference_spectrum(
        reactor_type, thermal_power_mwt, cooldown_months
    )
    density_kg_m3 = float(material["density_kg_m3"])
    density_g_cm3 = density_kg_m3 / 1000.0
    length_m = max(float(module_length_m), 0.001)
    diameter_m = max(float(module_width_m), float(module_height_m), 0.001)
    radius_m = 0.5 * diameter_m
    material_cost = float(material["raw_material_cost_2025_usd_per_kg"])

    thickness_cm, unshielded_dose, source_to_detector_m = _solve_thickness_cm(
        target_dose_mrem_h,
        dose_distance_m,
        radius_m,
        gamma_power_w,
        spectrum,
        density_g_cm3,
    )
    metrics = _closed_shell_metrics(
        thickness_cm,
        density_kg_m3,
        material_cost,
        module_mass_kg,
        length_m,
        diameter_m,
    )

    return {
        "shield_material": str(shield_material),
        "cooldown_months": float(cooldown_months),
        "target_dose_mrem_h": float(target_dose_mrem_h),
        "dose_distance_m": float(dose_distance_m),
        "dose_distance_from_shield_surface_m": float(dose_distance_m),
        "decay_heat_w": decay_heat_w,
        "gamma_power_w": gamma_power_w,
        "gamma_fraction": _GAMMA_FRACTION,
        "unshielded_dose_mrem_h": unshielded_dose,
        "source_to_detector_distance_m": source_to_detector_m,
        **metrics,
        "photon_model": "multigroup_shah",
        "photon_model_label": "Manit Shah 48-group",
        "photon_energy_group_count": 48,
        "source_spectrum_family": source_family,
        "source_spectrum_reference_power_mwt": ref_power,
        "source_spectrum_month_used": used_spectrum_month,
        "source_spectrum_capped_at_36_months": bool(float(cooldown_months) > 36.0),
    }
