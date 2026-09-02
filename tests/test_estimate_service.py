"""Smoke/regression guardrails for the Streamlit estimate service.

These tests intentionally use the same lightweight dependency strategy as
webapp/app.py: OpenMC/WATTS/matplotlib are stubbed, and material densities come
from webapp/materials_densities.json. They protect the app-facing calculation
path without requiring the full Streamlit app to launch.
"""

from __future__ import annotations

import json
import math
import os
import sys
import unittest
import warnings
from unittest.mock import MagicMock


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WEBAPP_DIR = os.path.join(REPO_ROOT, 'webapp')
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)
if WEBAPP_DIR not in sys.path:
    sys.path.insert(0, WEBAPP_DIR)


class _MaterialStub:
    def __init__(self, name=None, temperature=None):
        self.name = name
        self.temperature = temperature
        self.density = 0.0

    def set_density(self, units, value):
        self.density = value

    def add_nuclide(self, *args, **kwargs):
        pass

    def add_element(self, *args, **kwargs):
        pass

    def add_s_alpha_beta(self, *args, **kwargs):
        pass

    @staticmethod
    def mix_materials(materials, fractions, method, name=None):
        result = _MaterialStub(name=name)
        try:
            result.density = sum(m.density * f for m, f in zip(materials, fractions))
        except Exception:
            result.density = 0.0
        return result


class _MaterialsStub:
    def append(self, mat):
        pass

    def extend(self, mats):
        pass


class _ThinMaterial:
    __slots__ = ('name', 'density')

    def __init__(self, name, density):
        self.name = name
        self.density = density


def _install_runtime_stubs():
    openmc_stub = MagicMock()
    openmc_stub.Material = _MaterialStub
    openmc_stub.Materials = _MaterialsStub

    for mod in ['openmc', 'openmc.deplete', 'openmc.mgxs']:
        sys.modules[mod] = openmc_stub
    sys.modules['watts'] = MagicMock()

    mpl_stub = MagicMock()
    for mod in ['matplotlib', 'matplotlib.pyplot', 'matplotlib.patches',
                'matplotlib.colors']:
        sys.modules[mod] = mpl_stub


def _install_material_density_lookup():
    import core_design.openmc_materials_database as materials_db_mod

    densities_path = os.path.join(WEBAPP_DIR, 'materials_densities.json')
    with open(densities_path) as f:
        raw = json.load(f)

    thin_materials_by_reactor = {
        rtype: {
            name: _ThinMaterial(name, float(density))
            for name, density in mats.items()
        }
        for rtype, mats in raw.items()
    }

    def collect_materials_data(params):
        rtype = params.get('reactor type', 'LTMR')
        return thin_materials_by_reactor.get(
            rtype, thin_materials_by_reactor['LTMR']
        )

    materials_db_mod.collect_materials_data = collect_materials_data
    return raw


_install_runtime_stubs()
MATERIAL_DENSITIES_RAW = _install_material_density_lookup()

from webapp.estimate_service import (  # noqa: E402
    EstimateInputs,
    LcoeAtNoakInputs,
    run_estimate,
    run_lcoe_at_noak_unit,
)
from core_design.drums import (  # noqa: E402
    calculate_reflector_mass_LTMR,
    calculate_shutdown_rods_volumes_and_masses,
)


BASE_INPUTS = {
    'enrichment': 0.1975,
    'interest_rate': 0.05,
    'discount_rate': 0.05,
    'construction_duration': 12,
    'debt_to_equity': 0.5,
    'operation_mode': 'Remotely Monitored',
    'emergency_shutdowns': 2.0,
    'startup_duration': 21,
    'startup_duration_refueling': 30,
    'tax_credit_type': 'ITC',
    'tax_credit_value': 0.30,
    'plant_lifetime': 60,
    'tax_credit_units': 10,
}

REACTOR_CASES = {
    'LTMR': {
        'inputs': {
            'reactor_type': 'LTMR',
            'power_mwt': 20,
            'n_rings_per_assembly': 12,
            'active_height': 95,
        },
        'expected_lifetime_days': 2196,
        'expected_power_mwe': 6.2,
    },
    'GCMR': {
        'inputs': {
            'reactor_type': 'GCMR',
            'power_mwt': 15,
            'n_assembly_rings': 6,
            'n_core_rings': 5,
            'active_height': 192,
        },
        'expected_lifetime_days': 2186,
        'expected_power_mwe': 6.0,
    },
    'HPMR': {
        'inputs': {
            'reactor_type': 'HPMR',
            'power_mwt': 5,
            'n_assembly_rings': 6,
            'n_core_rings': 5,
            'active_height': 109,
        },
        'expected_lifetime_days': 5376,
        'expected_power_mwe': 1.8,
    },
}


class EstimateServiceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.filterwarnings('ignore')
        cls.results = {}
        for reactor_type, case in REACTOR_CASES.items():
            inputs = EstimateInputs(**BASE_INPUTS, **case['inputs'])
            cls.results[reactor_type] = run_estimate(inputs)

    def test_material_density_lookup_contains_all_app_reactors(self):
        self.assertEqual({'LTMR', 'GCMR', 'HPMR'},
                         set(MATERIAL_DENSITIES_RAW.keys()))
        for reactor_type, materials in MATERIAL_DENSITIES_RAW.items():
            self.assertGreater(
                len(materials), 10,
                f'{reactor_type} material-density lookup looks too small',
            )

    def test_all_reactors_return_non_empty_app_results(self):
        for reactor_type, result in self.results.items():
            with self.subTest(reactor_type=reactor_type):
                self.assertFalse(result.display_df.empty)
                self.assertFalse(result.enriched_df.empty)
                self.assertFalse(result.detailed_sorted_df.empty)
                self.assertEqual(reactor_type, result.params['reactor type'])
                self.assertIn('Account', result.display_df.columns)
                self.assertIn('Account Title', result.display_df.columns)
                self.assertIn('FOAK LCOE', result.enriched_df.columns)
                self.assertIn('NOAK LCOE', result.enriched_df.columns)

    def test_all_reactors_keep_expected_headline_physics(self):
        for reactor_type, case in REACTOR_CASES.items():
            with self.subTest(reactor_type=reactor_type):
                params = self.results[reactor_type].params
                self.assertEqual(case['expected_lifetime_days'],
                                 params['Fuel Lifetime'])
                self.assertAlmostEqual(case['expected_power_mwe'],
                                       params['Power MWe'],
                                       places=6)
                self.assertGreater(params['Capacity Factor'], 0)

    def test_cost_outputs_are_finite_and_positive(self):
        for reactor_type, result in self.results.items():
            with self.subTest(reactor_type=reactor_type):
                foak_lcoe = result.enriched_df['FOAK LCOE']
                noak_lcoe = result.enriched_df['NOAK LCOE']
                self.assertTrue((foak_lcoe.dropna() >= 0).all())
                self.assertTrue((noak_lcoe.dropna() >= 0).all())
                self.assertGreater(float(foak_lcoe.sum()), 0.0)
                self.assertGreater(float(noak_lcoe.sum()), 0.0)

    def test_reflector_masses_are_positive(self):
        for reactor_type, result in self.results.items():
            with self.subTest(reactor_type=reactor_type):
                params = result.params
                self.assertGreater(params['Radial Reflector Mass'], 0.0)
                self.assertGreater(params['Axial Reflector Mass'], 0.0)

    def test_ltmr_shutdown_rod_mass_uses_count_height_and_radii(self):
        params = {
            'reactor type': 'LTMR',
            'Number of Shutdown Rods': 6,
            'Shutdown Rod Height': 120.0,
            'Shutdown Rod Absorber Radius': 1.30,
            'Shutdown Rod Clad Radius': 1.50,
            'Shutdown Rod Absorber': 'B4C_enriched',
            'Shutdown Rod Cladding': 'SS304',
        }

        calculate_shutdown_rods_volumes_and_masses(params)

        absorber_volume = 6 * math.pi * 1.30 ** 2 * 120.0
        cladding_volume = 6 * math.pi * (1.50 ** 2 - 1.30 ** 2) * 120.0
        self.assertAlmostEqual(
            absorber_volume,
            params['Shutdown Rod Absorber Volume'],
        )
        self.assertAlmostEqual(
            cladding_volume,
            params['Shutdown Rod Cladding Volume'],
        )
        self.assertAlmostEqual(
            absorber_volume * 2.52 / 1000,
            params['Shutdown Rod Absorber Mass'],
        )
        self.assertAlmostEqual(
            cladding_volume * 7.98 / 1000,
            params['Shutdown Rod Cladding Mass'],
        )
        self.assertAlmostEqual(
            params['Shutdown Rod Absorber Mass']
            + params['Shutdown Rod Cladding Mass'],
            params['Shutdown Rods Mass'],
        )

    def test_ltmr_reflector_volumes_exclude_drum_tubes_without_overlap(self):
        core_radius = 83.11277148455888
        assembly_ftf = 77.14554296911777
        params = {
            'reactor type': 'LTMR',
            'Drum Radius': 9.016,
            'Number of Drums': 12,
            'Core Radius': core_radius,
            'Assembly FTF': assembly_ftf,
            'Active Height': 120.0,
            'Axial Reflector Thickness': 44.54,
            'Radial Reflector': 'Graphite',
            'Axial Reflector': 'Graphite',
        }

        calculate_reflector_mass_LTMR(params)

        tube_radius = 9.016 * 91 / 90
        tube_area = 12 * math.pi * tube_radius ** 2
        hex_area = math.sqrt(3) / 2 * assembly_ftf ** 2
        radial_volume = (
            math.pi * core_radius ** 2 - hex_area - tube_area
        ) * 120.0
        axial_volume = 2 * (
            math.pi * core_radius ** 2 - tube_area
        ) * 44.54

        self.assertAlmostEqual(tube_area, params['All Drum Tubes Area'])
        self.assertAlmostEqual(
            radial_volume,
            params['Radial Reflector Volume'],
        )
        self.assertAlmostEqual(
            axial_volume,
            params['Axial Reflector Volume'],
        )
        self.assertAlmostEqual(
            radial_volume * 1.60 / 1000,
            params['Radial Reflector Mass'],
        )
        self.assertAlmostEqual(
            axial_volume * 1.60 / 1000,
            params['Axial Reflector Mass'],
        )

    def test_noak_lcoe_anchor_returns_diagnostics(self):
        anchor = run_lcoe_at_noak_unit(LcoeAtNoakInputs(
            **BASE_INPUTS,
            **REACTOR_CASES['LTMR']['inputs'],
            noak_unit_number=10,
            include_diagnostics=True,
        ))

        self.assertTrue(math.isfinite(anchor.mean))
        self.assertGreater(anchor.mean, 0)
        self.assertTrue(math.isfinite(anchor.std))
        self.assertGreaterEqual(anchor.std, 0)
        self.assertIsNotNone(anchor.diag_df)
        self.assertFalse(anchor.diag_df.empty)
        self.assertEqual(10, anchor.diag_params['NOAK Unit Number'])

    def test_noak_lcoe_anchor_omits_diagnostics_by_default(self):
        anchor = run_lcoe_at_noak_unit(LcoeAtNoakInputs(
            **BASE_INPUTS,
            **REACTOR_CASES['LTMR']['inputs'],
            noak_unit_number=10,
        ))

        self.assertTrue(math.isfinite(anchor.mean))
        self.assertIsNone(anchor.diag_df)
        self.assertEqual(10, anchor.diag_params['NOAK Unit Number'])


if __name__ == '__main__':
    unittest.main()
