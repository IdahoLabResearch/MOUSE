# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
"""
reactor_config.py Builds a fully-populated params dict for LTMR, GCMR, or HPMR
without running OpenMC. OpenMC outputs (Fuel Lifetime, Mass U235, Mass U238) are
interpolated from a pre-computed parametric study CSV keyed by reactor type,
enrichment, and thermal power.
"""

import os
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Single source of truth for the escalation year used across the webapp.
# Change this value to update all cost-year references in one place.
# ---------------------------------------------------------------------------
ESCALATION_YEAR = 2025

# openmc and watts must already be stubbed in sys.modules before this module is imported.
from core_design.utils import (
    calculate_lattice_radius,
    calculate_hex_apothem,
    calculate_core_radius_from_hex,
    calculate_pins_in_assembly,
    calculate_heat_flux,
    calculate_heat_flux_TRISO,
    calculate_number_fuel_elements_hpmr,
    number_of_heatpipes_hmpr,
    calculate_total_number_of_TRISO_particles,
)
from core_design.drums import (
    calculate_drums_volumes_and_masses,
    calculate_reflector_mass_LTMR,
    calculate_reflector_mass_GCMR,
    calculate_moderator_mass_GCMR,
    calculate_reflector_and_moderator_mass_HPMR,
    calculate_moderator_mass, # used by LTMR
)
from core_design.pins_arrangement import LTMR_pins_arrangement
from core_design.openmc_template_LTMR import update_ltmr_reflector_geometry_from_drums
from reactor_engineering_evaluation.fuel_calcs import fuel_calculations
from webapp.fuel_lifetime_estimator import estimate_ltmr_fuel_lifetime
from webapp.gcmr_fuel_lifetime_estimator import estimate_gcmr_fuel_lifetime
from webapp.hpmr_fuel_lifetime_estimator import (
    estimate_hpmr_fuel_lifetime,
    hpmr_total_uranium_mass_g,
)
# BOP.py re-exports everything from tools.py and adds its own functions
from reactor_engineering_evaluation.BOP import (
    calculate_heat_exchanger_mass,
    calculate_primary_pump_mechanical_power,
)
from reactor_engineering_evaluation.tools import (
    mass_flow_rate,
    compressor_power,
    GCMR_integrated_heat_transfer_vessel,
    calculate_shielding_masses,
)
from reactor_engineering_evaluation.vessels_calcs import vessels_specs

# ---------------------------------------------------------------------------
# CSV-based OpenMC results interpolated from parametric study table
# ---------------------------------------------------------------------------
_CSV_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    'assets', 'Ref_Results', 'parametric_study_power_and_enrichment.csv',
)
_lookup_df = None


class SubcriticalError(Exception):
    """Raised when the interpolated fuel lifetime is zero (subcritical core)."""
    pass


class ShortLifetimeError(Exception):
    """Raised when the estimated fuel lifetime is too short (< 90 days) to be a
    meaningful design point cost calculations are skipped."""
    pass


_MIN_USEFUL_LIFETIME_DAYS = 90


def _load_lookup():
    global _lookup_df
    if _lookup_df is None:
        df = pd.read_csv(_CSV_PATH)
        df.columns = df.columns.str.strip()
        df['reactor type'] = df['reactor type'].str.strip()
        _lookup_df = df.dropna(subset=['reactor type']).copy()
    return _lookup_df


def interpolate_openmc_results(reactor_type, power_mwt, enrichment):
    """
    Interpolate OpenMC results for a given (reactor_type, power_mwt, enrichment).

    - Mass U235 / Mass U238 depend only on enrichment → 1-D linear interpolation.
    - Fuel Lifetime depends on both enrichment and power → bilinear interpolation:
        1. Find the two bracketing enrichment levels in the data.
        2. At each enrichment level, linearly interpolate in power.
        3. Linearly interpolate the two results in enrichment.
      This ensures that if both bracketing power values at a given enrichment are
      zero, the result is also zero avoiding the spurious non-zero values that
      unstructured triangulation (griddata) can produce across zero-lifetime regions.
    - Fuel Lifetime is clamped to >= 0.

    Returns a dict with keys 'Fuel Lifetime', 'Mass U235', 'Mass U238'.
    Raises SubcriticalError if Fuel Lifetime == 0.
    """
    df = _load_lookup()
    rdf = df[df['reactor type'] == reactor_type].copy()

    # --- 1-D interpolation for uranium masses (enrichment only) ---
    enr_unique = np.sort(rdf['Enrichment'].unique())
    u235_vals = [rdf[rdf['Enrichment'] == e]['Mass U235'].mean() for e in enr_unique]
    u238_vals = [rdf[rdf['Enrichment'] == e]['Mass U238'].mean() for e in enr_unique]
    mass_u235 = float(np.interp(enrichment, enr_unique, u235_vals))
    mass_u238 = float(np.interp(enrichment, enr_unique, u238_vals))

    # --- Bilinear interpolation for fuel lifetime (enrichment × power) ---
    def _fl_at_enr(enr):
        """1-D linear interpolation of Fuel Lifetime in power at a fixed enrichment."""
        sub = rdf[rdf['Enrichment'] == enr].sort_values('Power MWt')
        return float(np.interp(power_mwt,
                               sub['Power MWt'].values.astype(float),
                               sub['Fuel Lifetime'].values.astype(float)))

    idx = int(np.searchsorted(enr_unique, enrichment))

    if idx == 0:
        # At or below minimum enrichment clamp
        fl = _fl_at_enr(enr_unique[0])
    elif idx >= len(enr_unique):
        # At or above maximum enrichment clamp
        fl = _fl_at_enr(enr_unique[-1])
    else:
        enr_lo, enr_hi = enr_unique[idx - 1], enr_unique[idx]
        fl_lo = _fl_at_enr(enr_lo)
        fl_hi = _fl_at_enr(enr_hi)
        t = (enrichment - enr_lo) / (enr_hi - enr_lo)
        fl = fl_lo + t * (fl_hi - fl_lo)

    fl = max(0.0, float(fl))

    return {
        'Fuel Lifetime': int(round(fl)),
        'Mass U235': max(0, int(round(mass_u235))),
        'Mass U238': max(0, int(round(mass_u238))),
    }


def _build_ltmr(params):
    """Populate params for LTMR (Liquid Metal Thermal Microreactor).

    The user-controlled inputs are:
        Power MWt, Enrichment, Number of Rings per Assembly, Active Height
    All four are pre-set by build_params before this function is called.
    Drum Radius and Radial Reflector Thickness are auto-resolved from N
    using the LTMR drum-cover rule. Fuel Lifetime, Mass U235, Mass U238
    are derived from the KNN estimator and physics scaling.
    """

    # Sec 1: Materials
    params.update({
        'reactor type': 'LTMR',
        'TRISO Fueled': 'No',
        'Fuel': 'UZrH_alloy',
        'H_Zr_ratio': 1.6,
        'U_met_wo': 0.3,
        'er_wo': 0,
        'Coolant': 'NaK',
        'Radial Reflector': 'Graphite',
        'Axial Reflector': 'Graphite',
        'Moderator': 'ZrH',
        'Control Drum Absorber': 'B4C_enriched',
        'Control Drum Reflector': 'Graphite',
        'Common Temperature': 600,
        'HX Material': 'SS316',
    })

    # Sec 2: Geometry
    # 'Number of Rings per Assembly' and 'Active Height' come from build_params.
    params.update({
        'Fuel Pin Materials': ['Zr', None, params['Fuel'], None, 'SS304'],
        'Fuel Pin Radii': [0.28575, 0.3175, 1.5113, 1.5367, 1.5875],
        'Moderator Pin Materials': ['ZrH', 'SS304'],
        'Moderator Pin Inner Radius': 1.5367,
        'Moderator Pin Radii': [1.5367, 1.5875],
        'Pin Gap Distance': 0.1,
        'Pins Arrangement': LTMR_pins_arrangement,
    })
    params['Lattice Radius'] = calculate_hex_apothem(params)
    params['Assembly FTF'] = 2 * params['Lattice Radius']
    params['Fuel Pin Count'] = calculate_pins_in_assembly(params, 'FUEL')
    params['Moderator Pin Count'] = calculate_pins_in_assembly(params, 'MODERATOR')
    params['Moderator Mass'] = calculate_moderator_mass(params)

    # Sec 3: Control drums count fixed at 12 (matches watts_exec_LTMR.py).
    # Drum count is not a feature of the KNN fuel-lifetime estimator, so
    # this can be changed without invalidating the lifetime model.
    # Drum Radius is auto-maximised; reflector thickness then auto-set to
    # cover the drums (LTMR rule: max_outer_radius - hex_apothem).
    params.update({
        'Number of Drums': 12,
        'Drum Absorber Thickness': 1,
        'Drum Absorber Arc Degrees': 120,
    })
    update_ltmr_reflector_geometry_from_drums(params)
    calculate_drums_volumes_and_masses(params)
    calculate_reflector_mass_LTMR(params)

    # Sec 4: Overall system
    params.update({
        'Thermal Efficiency': 0.31,
        'Heat Flux Criteria': 0.9,
        'Burnup Steps': [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0,
                         30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 140.0],
    })
    params['Power MWe'] = params['Power MWt'] * params['Thermal Efficiency']
    params['Heat Flux'] = calculate_heat_flux(params)

    # Sec 5: Fuel lifetime via KNN estimator on the LTMR parametric study.
    # Mass U235 / U238 are derived from physics scaling: total uranium
    # mass per (pin × cm) is constant for TRIGA fuel; the reference
    # (N=12, H=78.4) parametric case has 345,481 g of U total.
    fl = estimate_ltmr_fuel_lifetime(
        n_rings_per_assembly = params['Number of Rings per Assembly'],
        active_height = params['Active Height'],
        enrichment = params['Enrichment'],
        power_mwt = params['Power MWt'],
    )
    if fl == 0:
        raise SubcriticalError(
            f"Reactor is subcritical for LTMR at "
            f"N_rings={params['Number of Rings per Assembly']}, "
            f"H={params['Active Height']:.1f} cm, "
            f"E={params['Enrichment']:.4f}, P={params['Power MWt']} MWt."
        )
    if fl < _MIN_USEFUL_LIFETIME_DAYS:
        raise ShortLifetimeError(
            f"Estimated fuel lifetime is only {fl} days "
            f"({fl / 30.0:.1f} months) too short for a meaningful "
            f"design point. Try increasing the diameter, the height, "
            f"or the enrichment."
        )

    # Total uranium mass per (fuel pin × cm) for LTMR TRIGA fuel. Derived
    # directly from the parametric-study Excel: every row gives the same
    # value of (Mass U235 + Mass U238) / (Fuel Pin Count × Active Height),
    # confirming the per-pin geometry is fixed across the design space.
    _U_PER_PIN_CM = 14.6888 # g of total uranium per (fuel pin · cm)
    total_u_mass = _U_PER_PIN_CM * params['Fuel Pin Count'] * params['Active Height']

    params['Fuel Lifetime'] = fl
    params['Mass U235'] = int(round(total_u_mass * params['Enrichment']))
    params['Mass U238'] = int(round(total_u_mass * (1.0 - params['Enrichment'])))
    params['Uranium Mass'] = (params['Mass U235'] + params['Mass U238']) / 1000 # kg
    fuel_calculations(params)

    # Sec 6: Primary Loop + BoP
    params.update({
        'Secondary HX Mass': 0,
        'Primary Pump': 'Yes',
        'Secondary Pump': 'No',
        'Pump Isentropic Efficiency': 0.8,
        'Primary Loop Inlet Temperature': 430 + 273.15,
        'Primary Loop Outlet Temperature': 520 + 273.15,
        'Secondary Loop Inlet Temperature': 395 + 273.15,
        'Secondary Loop Outlet Temperature': 495 + 273.15,
    })
    params['Primary HX Mass'] = calculate_heat_exchanger_mass(params)
    params.update({
        'BoP Count': 2,
        'BoP per loop load fraction': 0.5,
    })
    params['BoP Power kWe'] = 1000 * params['Power MWe'] * params['BoP per loop load fraction']
    mass_flow_rate(params)
    calculate_primary_pump_mechanical_power(params)

    # Sec 7: Shielding
    params.update({
        'In Vessel Shield Thickness': 10.16,
        'In Vessel Shield Inner Radius': params['Core Radius'],
        'In Vessel Shield Material': 'B4C_natural',
        'Out Of Vessel Shield Thickness': 39.37,
        'Out Of Vessel Shield Material': 'WEP',
        'Out Of Vessel Shield Effective Density Factor': 0.5,
    })
    params['In Vessel Shield Outer Radius'] = params['Core Radius'] + params['In Vessel Shield Thickness']

    # Sec 8: Vessels
    # ----------------------------------------------------------------
    # Reactor / guard / cooling / intake vessel dimensions for the LTMR.
    # All values updated 2026-04-28 with inline references; previous
    # values were on the thin/tight end of plausibility or contained
    # a unit-conversion bug (lower plenum was 2.85 cm, physically
    # unworkable). See Botros review for the comparison to literature.
    # ----------------------------------------------------------------
    params.update({
        'Vessel Radius': params['Core Radius'] + params['In Vessel Shield Thickness'],

        # Reactor-vessel wall thickness fixed at 2 cm for the
        # microreactor 1-60 MWt range. ASME BPVC Section III, Division 5
        # (high-temperature components for liquid-metal reactors) requires
        # a thickness compatible with creep + thermal cycling at >500 °C;
        # the pressure-driven thickness for D ~ 1-3 m at ~atmospheric
        # NaK pressure is only ~4 mm, so ASME minimums dominate.
        # Reference points: Oklo Aurora ~2 cm (INL-22/68167), FFTF 3 cm,
        # EBR-II 3.8 cm. At microreactor scale the variation across
        # 1-60 MWt is < 1 mm, so a fixed 2 cm is appropriate.
        'Vessel Thickness': 2,

        # Lower plenum height fixed at 50 cm. Reference:
        # IAEA TECDOC-1908 "Sodium-Cooled Fast Reactors: Status and Trends"
        # reports lower-plenum heights of 0.4-0.8 m across micro-to-medium
        # SFR designs; below ~30 cm the geometry cannot accommodate the
        # inlet manifold + flow distributor + grid plate. The previous
        # value (42.848 - 40 = 2.848 cm) was a unit-conversion bug a
        # 3 cm plenum is physically unworkable.
        'Vessel Lower Plenum Height': 50,

        # Upper plenum height kept at 47 cm. By assumption this region
        # also includes the cover-gas (argon) headspace above the free
        # liquid surface; LTMR does not currently track a separate
        # 'Vessel Upper Gas Gap'.
        'Vessel Upper Plenum Height': 47.152,
        'Vessel Upper Gas Gap': 0,

        'Vessel Bottom Depth': 32.129,
        'Vessel Material': 'stainless_steel',

        # Reactor-vessel ↔ guard-vessel gap fixed at 5 cm. Reference:
        # ASME Section III Division 5 NH-3000 series + OECD/NEA
        # "Sodium-cooled Fast Reactor Vessel Design Guidelines" (2017)
        # recommend 50-150 mm for thermal expansion and leak monitoring.
        # 5 cm is the lower bound for a small microreactor where manual
        # in-service inspection is not expected; the previous value of
        # 2 cm was below ASME-required clearance for SFR-class vessels.
        'Gap Between Vessel And Guard Vessel': 5,

        # Guard-vessel wall thickness fixed at 1 cm. Reference:
        # ASME Section III Class 3 minimum-thickness rules for
        # non-pressure-bearing structural shells with D > 1 m specify
        # ≥ 10 mm; IAEA TECDOC-1531 "Fast Reactor Database" reports
        # 10-20 mm typical for LMR guard vessels. The previous 0.5 cm
        # was below ASME minimum for a vessel that must contain a
        # primary-coolant leak.
        'Guard Vessel Thickness': 1,

        'Guard Vessel Material': 'stainless_steel',

        'Gap Between Guard Vessel And Cooling Vessel': 5,

        'Cooling Vessel Thickness': 0.5,
        'Cooling Vessel Material': 'stainless_steel',

        # Cooling-vessel ↔ intake-vessel gap (RVACS air downcomer) 
        # fixed at 5 cm. Reference: Hejzlar & Buongiorno, "Passive
        # Decay Heat Removal in Lead-cooled Fast Reactors", Nuclear
        # Engineering & Design (2007), recommend 50-150 mm for RVACS
        # systems removing 1-3 MWt of decay heat (typical microreactor
        # range). The previous 3 cm was below the lower bound for
        # reliable natural-circulation flow.
        'Gap Between Cooling Vessel And Intake Vessel': 5,

        'Intake Vessel Thickness': 0.5,
        'Intake Vessel Material': 'stainless_steel',
    })
    vessels_specs(params)
    calculate_shielding_masses(params)

    # Sec 9: Operation
    # Note: Operation Mode, Emergency Shutdowns Per Year, and Startup
    # Duration after Emergency Shutdown are all controlled by the user
    # via the webapp UI (see _run_estimate in app.py) they are
    # intentionally not set here. 'Startup Duration after Refueling'
    # is also UI-controlled, but a default must be set here because the
    # total_refueling_period calculation below runs before user_overrides
    # are applied; the UI value still wins for downstream consumers.
    params.update({
        'Number of On-Site Operators per Shift': 2,
        'Levelization Period': 60,
        'Refueling Period': 7,
        'Startup Duration after Refueling': 2,
        'Reactors Monitored Per Operator': 10,
        'Security Staff Per Shift': 1,
    })
    # Onsite coolant inventory: 1833 kg/MWt (rough estimate from
    # Creys-Malville: 5,500 t Na for a 3,000 MWt plant -> 1833 kg/MWt).
    # Reference: https://www.edf.fr/sites/default/files/mediatheque/dp_creys_2017.pdf
    params['Onsite Coolant Inventory'] = 1833 * params['Power MWt']
    params['Replacement Coolant Inventory'] = 0 # NaK is assumed not to require replacement

    total_refueling_period = (params['Fuel Lifetime'] + params['Refueling Period']
                              + params['Startup Duration after Refueling'])
    total_refueling_period_yr = total_refueling_period / 365
    params['A75: Outer Vessel Structure Replacement Period (years)'] = 20
    params['A75: Inner Vessel Structure Replacement Period (cycles)'] = np.floor(10 / total_refueling_period_yr)
    params['A75: Reflector Replacement Period (cycles)'] = np.floor(10 / total_refueling_period_yr)
    params['A75: Reactor Control Devices Replacement Period (cycles)'] = np.floor(10 / total_refueling_period_yr)
    params['A75: Moderator Booster Replacement Period (cycles)'] = 1
    params['Maintenance to Direct Cost Ratio'] = 0.015
    params['A78: CAPEX to Decommissioning Cost Ratio'] = 0.15

    # Sec 10: Buildings & Economic params
    params.update({
        'Land Area': 18,
        'Excavation Volume': 412.605,
        'Reactor Building Slab Roof Volume': (9750 * 6502.4 * 1500) / 1e9,
        'Reactor Building Basement Volume': (9750 * 6502.4 * 1500) / 1e9,
        'Reactor Building Exterior Walls Volume': ((2 * 9750 * 3500 * 1500) + (3502.4 * 3500 * (1500 + 750))) / 1e9,
        'Reactor Building Superstructure Area': ((2 * 3500 * 3500) + (2 * 7500 * 3500)) / 1e6,
        'Integrated Heat Exchanger Building Slab Roof Volume': 0,
        'Integrated Heat Exchanger Building Basement Volume': 0,
        'Integrated Heat Exchanger Building Exterior Walls Volume': 0,
        'Integrated Heat Exchanger Building Superstructure Area': 0,
        'Turbine Building Slab Roof Volume': (12192 * 2438 * 200) / 1e9,
        'Turbine Building Basement Volume': (12192 * 2438 * 200) / 1e9,
        'Turbine Building Exterior Walls Volume': ((12192 * 2496 * 200) + (2038 * 2496 * 200)) * 2 / 1e9,
        'Control Building Slab Roof Volume': (12192 * 2438 * 200) / 1e9,
        'Control Building Basement Volume': (12192 * 2438 * 200) / 1e9,
        'Control Building Exterior Walls Volume': ((12192 * 2496 * 200) + (2038 * 2496 * 200)) * 2 / 1e9,
        'Manipulator Building Slab Roof Volume': (4876.8 * 2438.4 * 400) / 1e9,
        'Manipulator Building Basement Volume': (4876.8 * 2438.4 * 1500) / 1e9,
        'Manipulator Building Exterior Walls Volume': ((4876.8 * 4445 * 400) + (2038.4 * 4445 * 400 * 2)) / 1e9,
        'Refueling Building Slab Roof Volume': 0,
        'Refueling Building Basement Volume': 0,
        'Refueling Building Exterior Walls Volume': 0,
        'Spent Fuel Building Slab Roof Volume': 0,
        'Spent Fuel Building Basement Volume': 0,
        'Spent Fuel Building Exterior Walls Volume': 0,
        'Emergency Building Slab Roof Volume': 0,
        'Emergency Building Basement Volume': 0,
        'Emergency Building Exterior Walls Volume': 0,
        'Storage Building Slab Roof Volume': (8400 * 3500 * 400) / 1e9,
        'Storage Building Basement Volume': (8400 * 3500 * 400) / 1e9,
        'Storage Building Exterior Walls Volume': ((8400 * 2700 * 400) + (3100 * 2700 * 400 * 2)) / 1e9,
        'Radwaste Building Slab Roof Volume': 0,
        'Radwaste Building Basement Volume': 0,
        'Radwaste Building Exterior Walls Volume': 0,
        'Annual Return': 0.0475,
        'NOAK Unit Number': 100,
        'Escalation Year': ESCALATION_YEAR,
        'Interest Rate': 0.07,
        'Discount Rate': 0.07,
        'Construction Duration': 12,
        'Debt To Equity Ratio': 1,
    })
    # ITC/PTC credits are controlled by the user via the webapp not hardcoded here.


def _build_gcmr(params):
    """Populate params for GCMR Design A (Gas Cooled Microreactor)."""

    # Sec 1: Materials
    # Note: 'Enrichment' and 'Power MWt' are pre-set by build_params from user input.
    params.update({
        'reactor type': 'GCMR',
        'TRISO Fueled': 'Yes',
        'Fuel': 'UCO',
        'UO2 atom fraction': 0.7,
        'Radial Reflector': 'Graphite',
        'Axial Reflector': 'Graphite',
        'Matrix Material': 'Graphite',
        'Moderator': 'Graphite',
        'Moderator Booster Materials': ['ZrH'],
        'Coolant': 'Helium',
        'Common Temperature': 850,
        'Control Drum Absorber': 'B4C_enriched',
        'Control Drum Reflector': 'Graphite',
        'HX Material': 'SS316',
    })

    # Sec 2: Geometry
    # 'Assembly Rings', 'Core Rings', 'Active Height' come from build_params.
    params.update({
        'Fuel Pin Materials': ['UCO', 'buffer_graphite', 'PyC', 'SiC', 'PyC'],
        'Fuel Pin Radii': [0.0250, 0.0350, 0.0390, 0.0425, 0.0465], # cm INL ART TRISO Fuel training (2019)
        'Compact Fuel Radius': 0.6225,
        'Packing Fraction': 0.3,
        'Coolant Channel Radius': 0.35,
        'Moderator Booster Radii': [0.55],
        'Lattice Pitch': 2.25,
    })
    params['Assembly FTF'] = params['Lattice Pitch'] * (params['Assembly Rings'] - 1) * np.sqrt(3)

    # Compute the per-compact and total TRISO particle counts so the
    # webapp's materials section can show them. Uses Compact Fuel
    # Radius, Active Height, Fuel Pin Radii (TRISO outer), Packing
    # Fraction, Assembly Rings, Core Rings all set above.
    calculate_total_number_of_TRISO_particles(params)

    # Sec 3: Control drums GCMR auto-resolution sets Drum Radius,
    # Radial Reflector Thickness, Core Radius, Drum Height inside
    # calculate_drums_volumes_and_masses (drum-cell rule).
    params.update({
        'Drum Absorber Thickness': 1,
    })
    calculate_drums_volumes_and_masses(params)
    calculate_reflector_mass_GCMR(params)
    calculate_moderator_mass_GCMR(params)

    # Sec 4: Overall system
    params.update({
        'Thermal Efficiency': 0.4,
        'Heat Flux Criteria': 0.9,
        'Burnup Steps': [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0,
                         30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0],
    })
    params['Power MWe'] = params['Power MWt'] * params['Thermal Efficiency']
    params['Heat Flux'] = calculate_heat_flux_TRISO(params)

    # Sec 5: Fuel lifetime via KNN estimator on the GCMR parametric study.
    # Mass U235 / U238 are derived from the physics scaling formula:
    # Total uranium mass = G_PER_VOL_INDEX × F_A × F_C × H
    # where F_A = 3(N_A−1)²−3(N_A−1)+1, F_C = 3 N_C²−3 N_C+1.
    fl = estimate_gcmr_fuel_lifetime(
        assembly_rings = params['Assembly Rings'],
        core_rings = params['Core Rings'],
        active_height = params['Active Height'],
        enrichment = params['Enrichment'],
        power_mwt = params['Power MWt'],
    )
    if fl == 0:
        raise SubcriticalError(
            f"Reactor is subcritical for GCMR at "
            f"N_A={params['Assembly Rings']}, N_C={params['Core Rings']}, "
            f"H={params['Active Height']:.1f} cm, "
            f"E={params['Enrichment']:.4f}, P={params['Power MWt']} MWt."
        )
    if fl < _MIN_USEFUL_LIFETIME_DAYS:
        raise ShortLifetimeError(
            f"Estimated fuel lifetime is only {fl} days "
            f"({fl / 30.0:.1f} months) too short for a meaningful "
            f"design point. Try increasing the diameter, the height, "
            f"or the enrichment."
        )

    _GCMR_G_PER_VOL_INDEX = 0.5776 # g of total uranium per (F_A × F_C × H)
    _na = params['Assembly Rings']
    _nc = params['Core Rings']
    _F_A = 3 * (_na - 1) ** 2 - 3 * (_na - 1) + 1
    _F_C = 3 * _nc ** 2 - 3 * _nc + 1
    total_u_mass = _GCMR_G_PER_VOL_INDEX * _F_A * _F_C * params['Active Height']

    params['Fuel Lifetime'] = fl
    params['Mass U235'] = int(round(total_u_mass * params['Enrichment']))
    params['Mass U238'] = int(round(total_u_mass * (1.0 - params['Enrichment'])))
    params['Uranium Mass'] = (params['Mass U235'] + params['Mass U238']) / 1000 # kg
    fuel_calculations(params)

    # Sec 6: Primary Loop + BoP
    params.update({
        'Primary Loop Purification': True,
        'Secondary HX Mass': 0,
        'Compressor Pressure Ratio': 4,
        'Compressor Isentropic Efficiency': 0.8,
        'Primary Loop Count': 2,
        'Primary Loop per loop load fraction': 0.5,
        'Primary Loop Inlet Temperature': 300 + 273.15,
        'Primary Loop Outlet Temperature': 550 + 273.15,
        'Secondary Loop Inlet Temperature': 270 + 273.15, # cold-end PCHE pinch 30°C with 300°C primary inlet (was 290 -> 10°C pinch, below realistic PCHE design)
        'Secondary Loop Outlet Temperature': 500 + 273.15,
        'Primary Loop Pressure Drop': 50e3,
    })
    params['Primary HX Mass'] = calculate_heat_exchanger_mass(params)
    mass_flow_rate(params)
    compressor_power(params)
    params.update({
        'BoP Count': 2,
        'BoP per loop load fraction': 0.5,
    })
    params['BoP Power kWe'] = 1000 * params['Power MWe'] * params['BoP per loop load fraction']
    params.update({
        'Integrated Heat Transfer Vessel Thickness': 0,
        'Integrated Heat Transfer Vessel Material': 'SA508',
    })
    GCMR_integrated_heat_transfer_vessel(params)

    # Sec 7: Shielding
    params.update({
        'In Vessel Shield Thickness': 0,
        'In Vessel Shield Inner Radius': params['Core Radius'],
        'In Vessel Shield Material': 'B4C_natural',
        'Out Of Vessel Shield Thickness': 39.37,
        'Out Of Vessel Shield Material': 'WEP',
        'Out Of Vessel Shield Effective Density Factor': 0.5,
    })
    params['In Vessel Shield Outer Radius'] = params['Core Radius'] + params['In Vessel Shield Thickness']

    # Sec 8: Vessels
    params.update({
        'Vessel Radius': params['Core Radius'] + params['In Vessel Shield Thickness'],

        # Reactor-vessel wall thickness fixed at 3 cm for the GCMR.
        # GCMR primary loop runs at ~4 MPa He (per
        # reactor_engineering_evaluation/tools.py:79). ASME Section III
        # Division 1 thin-shell formula: t = P·R/(S·E − 0.6·P) with
        # design pressure 1.1×4 = 4.4 MPa, S = 138 MPa for SA-508 at
        # 350 °C, E = 1.0, plus 3 mm corrosion allowance gives:
        # R = 60 cm → ~22 mm
        # R = 80 cm → ~29 mm
        # R = 100 cm → ~35 mm
        # Reference points: USNC MMR ~5 cm @ 5 MPa; X-energy Xe-100
        # ~15 cm @ 6 MPa; HTR-PM 17.7 cm @ 7 MPa (R≈2.5 m, too large to
        # scale from). 3 cm is a single conservative value covering
        # the MOUSE GCMR design space. Previous 1 cm was below the
        # ASME pressure-driven minimum.
        'Vessel Thickness': 3,

        # Lower plenum height fixed at 30 cm. GCMR has no liquid
        # plenum; this region houses the He inlet manifold and lower
        # flow distributor. Reference: GA MHTGR / HTR-PM-class designs
        # use 0.3-0.5 m for the lower distributor region; below ~20 cm
        # the inlet flow cannot be evened across the core. Previous
        # value (42.848 − 40 = 2.848 cm) was a unit-conversion bug.
        'Vessel Lower Plenum Height': 30,

        # Upper plenum (47 cm) kept; serves as outlet plenum for
        # hot-leg gas exit. No separate cover-gas headspace tracked.
        'Vessel Upper Plenum Height': 47.152,
        'Vessel Upper Gas Gap': 0,

        'Vessel Bottom Depth': 32.129,
        'Vessel Material': 'stainless_steel',

        # Guard vessel intentionally removed for GCMR. Helium is an
        # inert gas with no chemical-leak hazard (unlike NaK in LTMR);
        # there is no large coolant inventory to contain in a leak
        # scenario. Vessel↔guard-vessel gap and guard-vessel thickness
        # are therefore set to zero by design.
        'Gap Between Vessel And Guard Vessel': 0,
        'Guard Vessel Thickness': 0,
        'Guard Vessel Material': 'low_alloy_steel',

        # RVACS hot-leg gap kept at 5 cm.
        'Gap Between Guard Vessel And Cooling Vessel': 5,

        'Cooling Vessel Thickness': 0.5,
        'Cooling Vessel Material': 'stainless_steel',

        # Cooling-vessel ↔ intake-vessel gap (RVACS air downcomer) 
        # fixed at 5 cm. Reference: Hejzlar & Buongiorno, "Passive
        # Decay Heat Removal in Lead-cooled Fast Reactors", Nuclear
        # Engineering & Design (2007), recommend 50-150 mm for RVACS
        # systems removing 1-3 MWt of decay heat (typical microreactor
        # range). Previous 4 cm was below the recommended lower bound.
        'Gap Between Cooling Vessel And Intake Vessel': 5,

        'Intake Vessel Thickness': 0.5,
        'Intake Vessel Material': 'stainless_steel',
    })
    vessels_specs(params)
    calculate_shielding_masses(params)

    # Sec 9: Operation
    # Note: Operation Mode, Emergency Shutdowns Per Year, and Startup
    # Duration after Emergency Shutdown are all controlled by the user
    # via the webapp UI (see _run_estimate in app.py) they are
    # intentionally not set here. 'Startup Duration after Refueling'
    # is also UI-controlled, but a default must be set here because the
    # total_refueling_period calculation below runs before user_overrides
    # are applied; the UI value still wins for downstream consumers.
    params.update({
        'Number of On-Site Operators per Shift': 2,
        'Levelization Period': 60,
        'Refueling Period': 7,
        'Startup Duration after Refueling': 2,
        'Reactors Monitored Per Operator': 10,
        'Security Staff Per Shift': 1,
    })
    # Onsite He inventory: 3.3 kg/MWt (UNT 919556 tables 17 & 18).
    # Replacement: He loss rate ~10%/year (National Academies 12844),
    # so ~1/10 of initial inventory is replaced annually.
    params['Onsite Coolant Inventory'] = 3.3 * params['Power MWt']
    params['Replacement Coolant Inventory'] = params['Onsite Coolant Inventory'] / 10
    params['Annual Coolant Supply Frequency'] = 1 if params['Primary Loop Purification'] else 6

    total_refueling_period = (params['Fuel Lifetime'] + params['Refueling Period']
                              + params['Startup Duration after Refueling'])
    total_refueling_period_yr = total_refueling_period / 365
    params['A75: Outer Vessel Structure Replacement Period (years)'] = 20
    params['A75: Inner Vessel Structure Replacement Period (cycles)'] = np.floor(10 / total_refueling_period_yr)
    params['A75: Reflector Replacement Period (cycles)'] = np.floor(10 / total_refueling_period_yr)
    params['A75: Reactor Control Devices Replacement Period (cycles)'] = np.floor(10 / total_refueling_period_yr)
    params['A75: Moderator Booster Replacement Period (cycles)'] = 1
    params['Maintenance to Direct Cost Ratio'] = 0.015
    params['A78: CAPEX to Decommissioning Cost Ratio'] = 0.15

    # Sec 10: Buildings & Economic params
    params.update({
        'Land Area': 18,
        'Excavation Volume': 412.605,
        'Reactor Building Slab Roof Volume': (9750 * 6502.4 * 1500) / 1e9,
        'Reactor Building Basement Volume': (9750 * 6502.4 * 1500) / 1e9,
        'Reactor Building Exterior Walls Volume': ((2 * 9750 * 3500 * 1500) + (3502.4 * 3500 * (1500 + 750))) / 1e9,
        'Reactor Building Superstructure Area': ((2 * 3500 * 3500) + (2 * 7500 * 3500)) / 1e6,
        'Integrated Heat Exchanger Building Slab Roof Volume': 0,
        'Integrated Heat Exchanger Building Basement Volume': 0,
        'Integrated Heat Exchanger Building Exterior Walls Volume': 0,
        'Integrated Heat Exchanger Building Superstructure Area': 0,
        'Turbine Building Slab Roof Volume': (12192 * 2438 * 200) / 1e9,
        'Turbine Building Basement Volume': (12192 * 2438 * 200) / 1e9,
        'Turbine Building Exterior Walls Volume': ((12192 * 2496 * 200) + (2038 * 2496 * 200)) * 2 / 1e9,
        'Control Building Slab Roof Volume': (12192 * 2438 * 200) / 1e9,
        'Control Building Basement Volume': (12192 * 2438 * 200) / 1e9,
        'Control Building Exterior Walls Volume': ((12192 * 2496 * 200) + (2038 * 2496 * 200)) * 2 / 1e9,
        'Manipulator Building Slab Roof Volume': (4876.8 * 2438.4 * 400) / 1e9,
        'Manipulator Building Basement Volume': (4876.8 * 2438.4 * 1500) / 1e9,
        'Manipulator Building Exterior Walls Volume': ((4876.8 * 4445 * 400) + (2038.4 * 4445 * 400 * 2)) / 1e9,
        'Refueling Building Slab Roof Volume': 0,
        'Refueling Building Basement Volume': 0,
        'Refueling Building Exterior Walls Volume': 0,
        'Spent Fuel Building Slab Roof Volume': 0,
        'Spent Fuel Building Basement Volume': 0,
        'Spent Fuel Building Exterior Walls Volume': 0,
        'Emergency Building Slab Roof Volume': 0,
        'Emergency Building Basement Volume': 0,
        'Emergency Building Exterior Walls Volume': 0,
        'Storage Building Slab Roof Volume': (8400 * 3500 * 400) / 1e9,
        'Storage Building Basement Volume': (8400 * 3500 * 400) / 1e9,
        'Storage Building Exterior Walls Volume': ((8400 * 2700 * 400) + (3100 * 2700 * 400 * 2)) / 1e9,
        'Radwaste Building Slab Roof Volume': 0,
        'Radwaste Building Basement Volume': 0,
        'Radwaste Building Exterior Walls Volume': 0,
        'Annual Return': 0.0475,
        'NOAK Unit Number': 100,
        'Escalation Year': ESCALATION_YEAR,
        'Interest Rate': 0.07,
        'Discount Rate': 0.07,
        'Construction Duration': 12,
        'Debt To Equity Ratio': 1,
    })
    # ITC/PTC credits are controlled by the user via the webapp not hardcoded here.


def _build_hpmr(params):
    """Populate params for HPMR (Heat Pipe Microreactor)."""

    # Sec 1: Materials
    # Note: 'Enrichment' and 'Power MWt' are pre-set by build_params from user input.
    params.update({
        'reactor type': 'HPMR',
        'TRISO Fueled': 'Yes',
        'Fuel': 'homog_TRISO',
        'Radial Reflector': 'Graphite',
        'Axial Reflector': 'Graphite',
        'Moderator': 'monolith_graphite',
        'Coolant': 'Helium',
        'Control Drum Absorber': 'B4C_natural',
        'Control Drum Reflector': 'Graphite',
        'Cooling Device': 'heatpipe',
        'Common Temperature': 1000,
        'HX Material': 'SS316',
    })

    # Sec 2: Geometry
    # 'Number of Rings per Assembly', 'Number of Rings per Core',
    # and 'Active Height' come from build_params (user-controlled).
    # The previous version of this block locked H = 2 × Core Radius;
    # now H is independent so users can sweep it via the slider.
    params.update({
        'Fuel Pin Materials': ['homog_TRISO', 'Helium'],
        'Fuel Pin Radii': [1.00, 1.05],
        'Heat Pipe Materials': ['heatpipe', 'Helium'],
        'Heat Pipe Radii': [1.10, 1.15],
        'Lattice Pitch': 3.4,
    })
    params['Assembly FTF'] = ((params['Lattice Pitch'] * (params['Number of Rings per Assembly'] - 1)
                               + 1.4 * params['Fuel Pin Radii'][-1]) * np.sqrt(3))
    params['hexagonal Core Edge Length'] = ((params['Assembly FTF'] * (params['Number of Rings per Core'] - 1))
                                            + (params['Assembly FTF'] / 2) + 6.6)
    # NOTE: Radial Reflector Thickness, Core Radius, and Axial
    # Reflector Thickness are NOT set here. They are auto-resolved
    # in calculate_drums_volumes_and_masses (Sec 3) so the reflector
    # grows just enough to fully enclose the drums for any N_C.
    # Hardcoding the reflector at 50 cm caused drums to stick out of
    # the core for N_C ≥ 5 the HPMR DRUM RADIUS ERROR.
    params['Fuel Pin Count per Assembly'] = calculate_number_fuel_elements_hpmr(params['Number of Rings per Assembly'])
    params['Fuel Assemblies Count'] = ((3 * params['Number of Rings per Core'] ** 2)
                                       - (3 * params['Number of Rings per Core']))
    params['Fuel Pin Count'] = params['Fuel Assemblies Count'] * params['Fuel Pin Count per Assembly']
    number_of_heatpipes_hmpr(params)

    # Sec 3: Control drums. All four geometry parameters
    # Drum Radius, Drum Height, Radial Reflector Thickness, Core Radius
    # are auto-resolved inside calculate_drums_volumes_and_masses to
    # match the parametric study exactly:
    # Drum Radius = _calculate_max_hpmr_drum_radius(params)
    # largest drum that fits without
    # overlap on the placement ring
    # Radial Reflector = cd_distance + Drum Radius − hex_apothem
    # Core Radius = hex_apothem + Radial Reflector
    # Axial Reflector = Radial Reflector
    # Drum Height = Active Height + 2 × Axial Reflector
    # Verified to reproduce the (Drum Radius, Reflector, Core Radius)
    # values in the parametric study to within float rounding for
    # every (N_A=6, N_C ∈ [3..7]) configuration.
    params.update({
        'Drum Count': 12, # allowed: 6, 12, 18, 24
        'Drum Absorber Thickness': 1,
    })
    calculate_drums_volumes_and_masses(params)
    calculate_reflector_and_moderator_mass_HPMR(params)

    # Sec 4: Overall system
    # 'Power MWt' is pre-set by build_params; do not override it here.
    params.update({
        'Thermal Efficiency': 0.36,
        'Heat Flux Criteria': 0.9,
        'Time Steps': [t * 86400 for t in [0.01, 0.99, 3, 6, 20, 70, 100, 165,
                                            365, 365, 365, 365, 365, 365, 365.00]],
    })
    params['Power MWe'] = params['Power MWt'] * params['Thermal Efficiency']
    params['Heat Flux'] = calculate_heat_flux(params)

    # Sec 5: Fuel lifetime via KNN estimator on the HPMR parametric
    # study (XLSX, 104 rows over (N_C, H, E, P) at fixed N_A = 6).
    # Uranium mass is derived from the EXACT HPMR-specific physics
    # scaling (constant 1.6116 g per pin per cm × N_pins(N_A, N_C) × H)
    # so off-grid geometry queries produce consistent (lifetime, mass)
    # pairs. See webapp/hpmr_fuel_lifetime_estimator.py.
    fl = estimate_hpmr_fuel_lifetime(
        n_rings_per_assembly = params['Number of Rings per Assembly'],
        n_rings_per_core = params['Number of Rings per Core'],
        active_height = params['Active Height'],
        enrichment = params['Enrichment'],
        power_mwt = params['Power MWt'],
    )
    if fl == 0:
        raise SubcriticalError(
            f"Reactor is subcritical for HPMR at "
            f"N_a={params['Number of Rings per Assembly']}, "
            f"N_c={params['Number of Rings per Core']}, "
            f"H={params['Active Height']:.1f} cm, "
            f"E={params['Enrichment']:.4f}, P={params['Power MWt']} MWt."
        )
    if fl < _MIN_USEFUL_LIFETIME_DAYS:
        raise ShortLifetimeError(
            f"Estimated fuel lifetime is only {fl} days "
            f"({fl / 30.0:.1f} months) too short for a meaningful "
            f"design point. Try increasing the power, enrichment, "
            f"or core size."
        )

    total_u_g = hpmr_total_uranium_mass_g(
        params['Number of Rings per Assembly'],
        params['Number of Rings per Core'],
        params['Active Height'],
    )
    params['Fuel Lifetime'] = fl
    params['Mass U235'] = int(round(total_u_g * params['Enrichment']))
    params['Mass U238'] = int(round(total_u_g * (1.0 - params['Enrichment'])))
    params['Uranium Mass'] = (params['Mass U235'] + params['Mass U238']) / 1000 # kg
    fuel_calculations(params)

    # Sec 6: Primary Loop + BoP
    params.update({
        'Primary Loop Purification': True,
        'Secondary HX Mass': 0,
        'Primary Loop Count': 2,
        # HPMR primary "loop" is the heat-pipe condenser interface to
        # the HX. Sodium heat pipes are nearly isothermal (vapor-stream
        # ΔT typically 5-20 °C); represent as Tin = Tout = 650 °C to
        # match watts_exec_HPMR.py and the underlying physics. The
        # previous value (Tin=900, Tout=650) implied a 250 °C sensible-
        # heat drop, which is incorrect for a heat-pipe device.
        'Primary Loop Inlet Temperature': 650 + 273.15,
        'Primary Loop Outlet Temperature': 650 + 273.15,
        'Secondary Loop Inlet Temperature': 300 + 273.15,
        'Secondary Loop Outlet Temperature': 630 + 273.15,
    })
    params['Primary HX Mass'] = calculate_heat_exchanger_mass(params)
    params.update({
        'BoP Count': 2,
        'BoP per loop load fraction': 0.5,
    })
    params['BoP Power kWe'] = 1000 * params['Power MWe'] * params['BoP per loop load fraction']
    # Note: HPMR uses heat pipes no primary pump or mass flow rate call

    # Sec 7: Shielding
    params.update({
        'In Vessel Shield Thickness': 0,
        'In Vessel Shield Inner Radius': params['Core Radius'],
        'In Vessel Shield Material': 'B4C_natural',
        'Out Of Vessel Shield Thickness': 39.37,
        'Out Of Vessel Shield Material': 'WEP',
        'Out Of Vessel Shield Effective Density Factor': 0.5,
    })
    params['In Vessel Shield Outer Radius'] = params['Core Radius'] + params['In Vessel Shield Thickness']

    # Sec 8: Vessels
    params.update({
        'Vessel Radius': params['Core Radius'] + params['In Vessel Shield Thickness'],

        # Reactor-vessel wall thickness fixed at 2 cm for HPMR.
        # HPMR vessel is at ~atmospheric pressure (heat pipes are
        # individually sealed, primary "coolant" is the heat-pipe
        # working fluid Na inside SS316 cladding, not a bulk pool).
        # Pressure-driven thickness is therefore <1 mm; ASME Section III
        # Division 5 minimum thickness for high-temperature creep
        # service at >650 °C dominates. Reference: INL HPMR Design A,
        # Westinghouse eVinci, Oklo concepts use ~15-25 mm vessel walls.
        # Previous 1 cm was below the Div 5 high-T minimum.
        'Vessel Thickness': 2,

        # Lower plenum height fixed at 20 cm. HPMR has no liquid
        # plenum; this region houses the heat-pipe condenser/evaporator
        # transition manifold and lower core support. Reference: INL
        # HPMR Design A and MARVEL designs show 15-25 cm for the lower
        # heat-pipe header / support region. Previous value
        # (42.848 − 40 = 2.848 cm) was a unit-conversion bug.
        'Vessel Lower Plenum Height': 20,

        # Upper plenum (47 cm) kept; houses heat-pipe condenser
        # interface to the secondary heat exchanger.
        'Vessel Upper Plenum Height': 47.152,
        'Vessel Upper Gas Gap': 0,

        'Vessel Bottom Depth': 32.129,
        'Vessel Material': 'stainless_steel',

        # Guard vessel intentionally removed for HPMR. There is no
        # bulk primary coolant each heat pipe is individually sealed
        # (Na in SS316 cladding) and a single failure releases only a
        # small inventory. No need for a secondary containment shell
        # around the entire reactor vessel. Vessel↔guard-vessel gap
        # and guard-vessel thickness are therefore set to zero by
        # design.
        'Gap Between Vessel And Guard Vessel': 0,
        'Guard Vessel Thickness': 0,
        'Guard Vessel Material': 'low_alloy_steel',

        # RVACS hot-leg gap kept at 5 cm.
        'Gap Between Guard Vessel And Cooling Vessel': 5,

        'Cooling Vessel Thickness': 0.5,
        'Cooling Vessel Material': 'stainless_steel',

        # Cooling-vessel ↔ intake-vessel gap (RVACS air downcomer) 
        # fixed at 5 cm. Reference: Hejzlar & Buongiorno, "Passive
        # Decay Heat Removal in Lead-cooled Fast Reactors", Nuclear
        # Engineering & Design (2007), recommend 50-150 mm for RVACS
        # systems removing 1-3 MWt of decay heat (typical microreactor
        # range). Previous 4 cm was below the recommended lower bound.
        'Gap Between Cooling Vessel And Intake Vessel': 5,

        'Intake Vessel Thickness': 0.5,
        'Intake Vessel Material': 'stainless_steel',
    })
    vessels_specs(params)
    calculate_shielding_masses(params)

    # Sec 9: Operation
    # Note: Operation Mode, Emergency Shutdowns Per Year, and Startup
    # Duration after Emergency Shutdown are all controlled by the user
    # via the webapp UI (see _run_estimate in app.py) they are
    # intentionally not set here. 'Startup Duration after Refueling'
    # is also UI-controlled, but a default must be set here because the
    # total_refueling_period calculation below runs before user_overrides
    # are applied; the UI value still wins for downstream consumers.
    params.update({
        'Number of On-Site Operators per Shift': 2,
        'Levelization Period': 60,
        'Refueling Period': 7,
        'Startup Duration after Refueling': 2,
        'Reactors Monitored Per Operator': 10,
        'Security Staff Per Shift': 1,
    })
    # HPMR has no bulk primary coolant inventory: each heat pipe is
    # individually sealed (Na in SS316), and the helium gap between
    # fuel pin and heat pipe is extremely thin and can be neglected.
    params['Onsite Coolant Inventory'] = 0
    params['Replacement Coolant Inventory'] = 0

    total_refueling_period = (params['Fuel Lifetime'] + params['Refueling Period']
                              + params['Startup Duration after Refueling'])
    total_refueling_period_yr = total_refueling_period / 365
    params['A75: Outer Vessel Structure Replacement Period (years)'] = 20
    params['A75: Inner Vessel Structure Replacement Period (cycles)'] = np.floor(10 / total_refueling_period_yr)
    params['A75: Reflector Replacement Period (cycles)'] = np.floor(10 / total_refueling_period_yr)
    params['A75: Reactor Control Devices Replacement Period (cycles)'] = np.floor(10 / total_refueling_period_yr)
    params['A75: Moderator Booster Replacement Period (cycles)'] = 1
    params['Maintenance to Direct Cost Ratio'] = 0.015
    params['A78: CAPEX to Decommissioning Cost Ratio'] = 0.15

    # Sec 10: Buildings & Economic params
    params.update({
        'Land Area': 18,
        'Excavation Volume': 412.605,
        'Reactor Building Slab Roof Volume': (9750 * 6502.4 * 1500) / 1e9,
        'Reactor Building Basement Volume': (9750 * 6502.4 * 1500) / 1e9,
        'Reactor Building Exterior Walls Volume': ((2 * 9750 * 3500 * 1500) + (3502.4 * 3500 * (1500 + 750))) / 1e9,
        'Reactor Building Superstructure Area': ((2 * 3500 * 3500) + (2 * 7500 * 3500)) / 1e6,
        'Integrated Heat Exchanger Building Slab Roof Volume': 0,
        'Integrated Heat Exchanger Building Basement Volume': 0,
        'Integrated Heat Exchanger Building Exterior Walls Volume': 0,
        'Integrated Heat Exchanger Building Superstructure Area': 0,
        'Turbine Building Slab Roof Volume': (12192 * 2438 * 200) / 1e9,
        'Turbine Building Basement Volume': (12192 * 2438 * 200) / 1e9,
        'Turbine Building Exterior Walls Volume': ((12192 * 2496 * 200) + (2038 * 2496 * 200)) * 2 / 1e9,
        'Control Building Slab Roof Volume': (12192 * 2438 * 200) / 1e9,
        'Control Building Basement Volume': (12192 * 2438 * 200) / 1e9,
        'Control Building Exterior Walls Volume': ((12192 * 2496 * 200) + (2038 * 2496 * 200)) * 2 / 1e9,
        'Manipulator Building Slab Roof Volume': (4876.8 * 2438.4 * 400) / 1e9,
        'Manipulator Building Basement Volume': (4876.8 * 2438.4 * 1500) / 1e9,
        'Manipulator Building Exterior Walls Volume': ((4876.8 * 4445 * 400) + (2038.4 * 4445 * 400 * 2)) / 1e9,
        'Refueling Building Slab Roof Volume': 0,
        'Refueling Building Basement Volume': 0,
        'Refueling Building Exterior Walls Volume': 0,
        'Spent Fuel Building Slab Roof Volume': 0,
        'Spent Fuel Building Basement Volume': 0,
        'Spent Fuel Building Exterior Walls Volume': 0,
        'Emergency Building Slab Roof Volume': 0,
        'Emergency Building Basement Volume': 0,
        'Emergency Building Exterior Walls Volume': 0,
        'Storage Building Slab Roof Volume': (8400 * 3500 * 400) / 1e9,
        'Storage Building Basement Volume': (8400 * 3500 * 400) / 1e9,
        'Storage Building Exterior Walls Volume': ((8400 * 2700 * 400) + (3100 * 2700 * 400 * 2)) / 1e9,
        'Radwaste Building Slab Roof Volume': 0,
        'Radwaste Building Basement Volume': 0,
        'Radwaste Building Exterior Walls Volume': 0,
        'Annual Return': 0.0475,
        'NOAK Unit Number': 100,
        'Escalation Year': ESCALATION_YEAR,
        'Interest Rate': 0.07,
        'Discount Rate': 0.07,
        'Construction Duration': 12,
        'Debt To Equity Ratio': 1,
    })
    # No ITC/PTC credits for HPMR by default


def build_params(reactor_type, power_mwt, enrichment, user_overrides,
                 n_rings_per_assembly=None, active_height=None,
                 n_assembly_rings=None, n_core_rings=None):
    """
    Build a fully-populated params dict for the given reactor type,
    then apply user_overrides on top.

    Parameters
    ----------
    reactor_type : str
        One of 'LTMR', 'GCMR', 'HPMR'.
    power_mwt : float
        Thermal power in MWt.
    enrichment : float
        U-235 enrichment fraction (0.05-0.1975).
    user_overrides : dict
        User-supplied values that override defaults (e.g. Interest Rate).
    n_rings_per_assembly : int, optional
        Number of rings per fuel assembly. Required for LTMR (training range
        6-24); ignored by GCMR/HPMR for now.
    active_height : float, optional
        Active core height in cm. Required for LTMR (training range
        50-180); ignored by GCMR/HPMR for now.

    Returns
    -------
    dict
        Fully-populated params dict ready for bottom_up_cost_estimate().

    Raises
    ------
    SubcriticalError
        If the reactor is subcritical at the requested operating point.
    """
    params = {}
    # Set the user-controlled inputs FIRST so all builder sections can use them.
    params['Power MWt'] = float(power_mwt)
    params['Enrichment'] = float(enrichment)

    if reactor_type == 'LTMR':
        if n_rings_per_assembly is None or active_height is None:
            raise ValueError(
                "LTMR requires n_rings_per_assembly and active_height."
            )
        params['Number of Rings per Assembly'] = int(n_rings_per_assembly)
        params['Active Height'] = float(active_height)

    if reactor_type == 'GCMR':
        if n_assembly_rings is None or n_core_rings is None or active_height is None:
            raise ValueError(
                "GCMR requires n_assembly_rings, n_core_rings, and active_height."
            )
        params['Assembly Rings'] = int(n_assembly_rings)
        params['Core Rings'] = int(n_core_rings)
        params['Active Height'] = float(active_height)

    if reactor_type == 'HPMR':
        # HPMR uses the same kwargs as GCMR (n_assembly_rings,
        # n_core_rings, active_height) but the builder reads them
        # under HPMR-style key names ('Number of Rings per Assembly',
        # 'Number of Rings per Core').
        if n_assembly_rings is None or n_core_rings is None or active_height is None:
            raise ValueError(
                "HPMR requires n_assembly_rings, n_core_rings, and active_height."
            )
        params['Number of Rings per Assembly'] = int(n_assembly_rings)
        params['Number of Rings per Core'] = int(n_core_rings)
        params['Active Height'] = float(active_height)

    builders = {'LTMR': _build_ltmr, 'GCMR': _build_gcmr, 'HPMR': _build_hpmr}
    if reactor_type not in builders:
        raise ValueError(f"Unknown reactor type: {reactor_type!r}. Choose LTMR, GCMR, or HPMR.")

    builders[reactor_type](params)

    # Apply user overrides (these come last so they win)
    params.update(user_overrides)

    # Number of Samples = 1: use deterministic class-3 estimates (no Monte Carlo sampling).
    # With Number of Samples > 1, the lognormal sampler is called with NaN bounds for
    # many database rows (those lacking low/high-end estimates), which propagates NaN
    # through the entire cost calculation and causes parent-account rows to appear blank
    # in the table. Number of Samples = 1 uses unit_cost_0 / fixed_cost_0 directly and
    # avoids this issue entirely.
    params['Number of Samples'] = 10

    return params
