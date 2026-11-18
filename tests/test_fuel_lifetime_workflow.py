"""Exercise lifetime handling at real workflow boundaries without transport runs."""

from contextlib import ExitStack, redirect_stderr, redirect_stdout
from io import StringIO
import sys
import unittest
from unittest.mock import MagicMock, patch

import numpy as np

from core_design.fuel_lifetime import (
    EXTRAPOLATED,
    INTERPOLATED,
    LOWER_BOUND,
    SUBCRITICAL,
    InvalidFuelLifetimeError,
    estimate_fuel_lifetime,
    record_fuel_lifetime,
)


# Import the production module using the same dependency boundary as the app.
# The tests patch its OpenMC/WATTS objects individually, whether it was already
# imported by another regression module or first imported here.
_openmc = MagicMock()
_mpl = MagicMock()
with patch.dict(sys.modules, {
    'openmc': _openmc,
    'openmc.deplete': _openmc.deplete,
    'openmc.mgxs': _openmc.mgxs,
    'watts': MagicMock(),
    'matplotlib': _mpl,
    'matplotlib.pyplot': _mpl.pyplot,
    'matplotlib.patches': _mpl.patches,
}):
    from core_design import utils


class FuelLifetimeWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.stack = ExitStack()
        self.addCleanup(self.stack.close)
        self.output = StringIO()
        self.stack.enter_context(redirect_stdout(self.output))
        self.stack.enter_context(redirect_stderr(StringIO()))
        self.openmc = MagicMock()
        self.openmc.config = {}
        self.watts = MagicMock()
        self.stack.enter_context(patch.object(utils, 'openmc', self.openmc))
        self.stack.enter_context(patch.object(utils, 'watts', self.watts))

    def params(self, shutdown=False):
        return {
            'cross_sections_xml_location': 'cross_sections.xml',
            'simplified_chain_thermal_xml': 'chain.xml',
            'Power MWt': 15.0,
            'Time Steps': [86400.0, 86400.0],
            'Active Height': 100.0,
            'Axial Reflector Thickness': 20.0,
            'Core Radius': 50.0,
            'Common Temperature': 900.0,
            'Cold Shutdown Temperature': 300.0,
            'Shutdown Margin Calc': shutdown,
            'Isothermal Temperature Coefficients': False,
        }

    def record(self, params, times, keff):
        record_fuel_lifetime(params, estimate_fuel_lifetime(times, keff))

    def test_requested_shutdown_analysis_does_not_bypass_operating_guard(self):
        for values in ([0.99, 0.98], [1.04, 1.02], [1.02, 1.03]):
            with self.subTest(values=values):
                params = self.params()

                def depletion(target):
                    self.record(target, [0.0, 100.0], values)

                with patch.object(utils, 'run_depletion_analysis', side_effect=depletion), \
                     patch.object(utils, '_run_lifecycle_snapshot_calculations') as snapshots:
                    with self.assertRaises(InvalidFuelLifetimeError):
                        utils._run_operating_depletion_and_lifecycle(
                            MagicMock(), params, True, False,
                        )
                    snapshots.assert_not_called()
                self.assertTrue(params['Shutdown Margin Calc'])
                self.assertIn('Fuel Lifetime Status', params)

    def test_2d_invalid_lifetime_stops_before_peaking_or_mass_exports(self):
        for values, expected_status in (
            ([0.99, 0.98, 0.97], SUBCRITICAL),
            ([1.06, 1.04, 1.02], EXTRAPOLATED),
            ([1.02, 1.02, 1.03], LOWER_BOUND),
        ):
            with self.subTest(status=expected_status):
                params = self.params(shutdown=True)
                times = [0.0, 100.0, 200.0]

                def corrected(_results, _height, core_radius, lifetime_params):
                    self.assertIs(lifetime_params, params)
                    self.record(lifetime_params, times, values)
                    return (params['Fuel Lifetime'], times, values, values, .9, 10, .8, 20)

                results = self.openmc.deplete.Results.return_value
                results.export_to_materials.reset_mock()
                with patch.object(utils, 'corrected_keff_2d', side_effect=corrected), \
                     patch.object(utils, 'compute_pin_peaking_factors') as peaking:
                    with self.assertRaises(InvalidFuelLifetimeError):
                        utils.openmc_depletion(params, MagicMock(), MagicMock())
                    peaking.assert_not_called()
                    results.export_to_materials.assert_not_called()
                self.assertEqual(params['Fuel Lifetime Status'], expected_status)
                self.assertEqual(params['Fuel Lifetime Simulated Endpoint (days)'], 200.0)
                if expected_status == SUBCRITICAL:
                    self.assertEqual(params['Fuel Lifetime'], 0.0)

    def results_3d(self, values):
        results = MagicMock()
        results.get_keff.return_value = (
            np.array([0.0, 100.0, 200.0]) * 86400.0,
            np.array([[value, 0.001] for value in values]),
        )
        material = MagicMock()
        material.get_mass.side_effect = lambda nuclide: {'U235': 1000.0, 'U238': 9000.0}[nuclide]
        results.export_to_materials.return_value = [material]
        return results

    def test_3d_subcritical_shutdown_then_operating_crossing(self):
        params = self.params(shutdown=True)
        shutdown = self.results_3d([0.90, 0.89, 0.88])
        operating = self.results_3d([1.04, 1.02, 0.98])
        self.openmc.deplete.Results.side_effect = [shutdown, operating]

        lifetime, mass235, mass238 = utils.openmc_depletion_3d(params, MagicMock(), MagicMock())
        self.assertEqual(lifetime, 0.0)
        self.assertEqual((mass235, mass238), (1000.0, 9000.0))
        self.assertEqual(params['keff 3D'], [0.90, 0.89, 0.88])
        self.assertNotIn('Fuel Lifetime Status', params)
        self.assertNotIn('Fuel Lifetime', params)
        self.assertNotIn('ERROR:', self.output.getvalue())

        params['Shutdown Margin Calc'] = False
        lifetime, _, _ = utils.openmc_depletion_3d(params, MagicMock(), MagicMock())
        self.assertEqual(lifetime, 150.0)
        self.assertEqual(params['Fuel Lifetime'], 150.0)
        self.assertEqual(params['Fuel Lifetime Status'], INTERPOLATED)
        self.assertFalse(params['Fuel Lifetime Extrapolated'])
        self.assertEqual(params['Fuel Lifetime Simulated Endpoint (days)'], 200.0)

    def test_3d_invalid_operating_lifetime_blocks_before_material_export(self):
        for values, expected_status in (
            ([0.99, 0.98, 0.97], SUBCRITICAL),
            ([1.06, 1.04, 1.02], EXTRAPOLATED),
            ([1.02, 1.02, 1.03], LOWER_BOUND),
        ):
            with self.subTest(status=expected_status):
                params = self.params()
                results = self.results_3d(values)
                self.openmc.deplete.Results.return_value = results
                with self.assertRaises(InvalidFuelLifetimeError):
                    utils.openmc_depletion_3d(params, MagicMock(), MagicMock())
                results.export_to_materials.assert_not_called()
                self.assertEqual(params['Fuel Lifetime Status'], expected_status)

    def test_run_openmc_preserves_diagnostics_and_restores_flags_after_error(self):
        params = self.params(shutdown=True)
        params['Isothermal Temperature Coefficients'] = True
        params['Temperature Perturbation'] = 100.0

        def depletion(target):
            self.record(target, [0.0, 100.0, 200.0], [1.06, 1.04, 1.02])
            target['_Depleted Fuel Materials XML'] = 'temporary-materials.xml'
            target['_Operating Fuel Material ID'] = 42
            target['Common Temperature'] = 300.0

        def plugin(target, function):
            self.assertIs(target, params)
            return function()

        self.watts.PluginOpenMC.return_value.side_effect = plugin
        with patch.object(utils, 'run_depletion_analysis', side_effect=depletion), \
             patch.object(utils, '_run_lifecycle_snapshot_calculations') as snapshots:
            with self.assertRaisesRegex(InvalidFuelLifetimeError, 'Extend the depletion schedule'):
                utils.run_openmc(MagicMock(), MagicMock(), params)
            snapshots.assert_not_called()

        self.assertTrue(params['Shutdown Margin Calc'])
        self.assertTrue(params['Isothermal Temperature Coefficients'])
        self.assertEqual(params['Common Temperature'], 900.0)
        self.assertNotIn('_Depleted Fuel Materials XML', params)
        self.assertNotIn('_Operating Fuel Material ID', params)
        self.assertEqual(params['Fuel Lifetime Status'], EXTRAPOLATED)
        self.assertEqual(params['Fuel Lifetime Simulated Endpoint (days)'], 200.0)
        self.assertAlmostEqual(params['Fuel Lifetime'], 300.0)


if __name__ == '__main__':
    unittest.main()
