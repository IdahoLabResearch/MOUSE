import contextlib
import io
import math
import unittest
from collections import UserDict

from core_design.fuel_lifetime import (
    EXTRAPOLATED,
    INTERPOLATED,
    LOWER_BOUND,
    SUBCRITICAL,
    InvalidFuelLifetimeError,
    estimate_fuel_lifetime,
    record_fuel_lifetime,
    require_cost_eligible_lifetime,
)


class FuelLifetimeEstimationTests(unittest.TestCase):
    def test_crossing_is_interpolated_without_rounding_to_cycle_years(self):
        result = estimate_fuel_lifetime([0, 100, 200], [1.1, 1.02, 0.98])
        self.assertAlmostEqual(result.lifetime_days, 150)
        self.assertEqual(result.status, INTERPOLATED)
        self.assertFalse(result.extrapolated)
        self.assertEqual(result.simulated_endpoint_days, 200)

    def test_first_downward_crossing_is_used_even_if_keff_recovers(self):
        result = estimate_fuel_lifetime([0, 100, 200, 300], [1.1, 0.9, 1.1, 0.9])
        self.assertAlmostEqual(result.lifetime_days, 50)
        self.assertEqual(result.status, INTERPOLATED)
        self.assertEqual(result.simulated_endpoint_days, 300)

    def test_exact_critical_endpoint_is_a_simulated_crossing(self):
        result = estimate_fuel_lifetime([0, 100], [1.02, 1.0])
        self.assertEqual(result.lifetime_days, 100)
        self.assertEqual(result.status, INTERPOLATED)

    def test_extrapolation_records_prediction_separately_from_endpoint(self):
        result = estimate_fuel_lifetime([0, 100, 200], [1.1, 1.04, 1.02])
        self.assertAlmostEqual(result.lifetime_days, 300)
        self.assertEqual(result.status, EXTRAPOLATED)
        self.assertTrue(result.extrapolated)
        self.assertEqual(result.simulated_endpoint_days, 200)

    def test_flat_or_rising_final_trend_is_only_a_lower_bound(self):
        for values in ([1.1, 1.02, 1.02], [1.1, 1.02, 1.04]):
            with self.subTest(values=values):
                result = estimate_fuel_lifetime([0, 100, 200], values)
                self.assertEqual(result.lifetime_days, 200)
                self.assertEqual(result.status, LOWER_BOUND)
                self.assertFalse(result.extrapolated)

    def test_one_supercritical_point_cannot_establish_a_lifetime(self):
        result = estimate_fuel_lifetime([0], [1.1])
        self.assertEqual(result.lifetime_days, 0)
        self.assertEqual(result.status, LOWER_BOUND)

    def test_subcritical_start_is_zero_even_if_later_points_rise(self):
        for start in (0.98, 1.0):
            with self.subTest(start=start):
                result = estimate_fuel_lifetime([0, 100], [start, 1.1])
                self.assertEqual(result.lifetime_days, 0)
                self.assertEqual(result.status, SUBCRITICAL)
                self.assertFalse(result.extrapolated)

    def test_invalid_histories_fail_before_interpretation(self):
        cases = [
            ([], []), ([0], [1.1, 1.0]),
            ([0, 0], [1.1, 1.0]), ([1, 0], [1.1, 1.0]),
            ([-1, 0], [1.1, 1.0]), ([0, math.inf], [1.1, 1.0]),
            ([0, math.nan], [1.1, 1.0]), ([0, 1], [1.1, math.nan]),
            ([0, 1], [1.1, math.inf]), ([0, 1], [1.1, 0]),
            ([0, 1], [1.1, -0.1]), ([0, 'bad'], [1.1, 1.0]),
        ]
        for times, values in cases:
            with self.subTest(times=times, values=values):
                with self.assertRaises(ValueError):
                    estimate_fuel_lifetime(times, values)


class FuelLifetimeReportingAndEligibilityTests(unittest.TestCase):
    def record(self, times, keff):
        params = UserDict()
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            record_fuel_lifetime(params, estimate_fuel_lifetime(times, keff))
        return params, output.getvalue()

    def test_record_preserves_day_rounding_and_simulated_endpoint(self):
        params, message = self.record([0, 201], [1.01, 0.99])
        self.assertEqual(params['Fuel Lifetime'], 100)
        self.assertEqual(params['Fuel Lifetime Simulated Endpoint (days)'], 201)
        self.assertEqual(params['Fuel Lifetime Status'], INTERPOLATED)
        self.assertIs(params['Fuel Lifetime Extrapolated'], False)
        self.assertIn('Simulated endpoint: 201 days', message)
        require_cost_eligible_lifetime(params)

    def test_subcritical_report_is_a_clear_error_and_blocks_costs(self):
        params, message = self.record([0, 100], [0.98, 0.97])
        self.assertEqual(params['Fuel Lifetime'], 0)
        self.assertIn('ERROR:', message)
        self.assertIn('0 days', message)
        self.assertIn('Cost estimation will not run', message)
        with self.assertRaisesRegex(InvalidFuelLifetimeError, '0 days'):
            require_cost_eligible_lifetime(params)

    def test_explicit_3d_reporting_can_preserve_fractional_days(self):
        result = estimate_fuel_lifetime([0, 201], [1.01, 0.99])
        params = UserDict()
        with contextlib.redirect_stdout(io.StringIO()):
            record_fuel_lifetime(params, result, round_days=False)
        self.assertAlmostEqual(params['Fuel Lifetime'], 100.5)
        self.assertEqual(params['Fuel Lifetime'], result.lifetime_days)
        self.assertEqual(params['Fuel Lifetime Status'], INTERPOLATED)
        require_cost_eligible_lifetime(params)

    def test_provisional_predictions_and_bounds_require_extended_depletion(self):
        for values, status in (([1.04, 1.02], EXTRAPOLATED), ([1.02, 1.03], LOWER_BOUND)):
            with self.subTest(status=status):
                params, message = self.record([0, 100], values)
                self.assertEqual(params['Fuel Lifetime Status'], status)
                self.assertIn('WARNING:', message)
                self.assertIn('100 days', message)
                self.assertIn('Extend the depletion', message)
                with self.assertRaisesRegex(InvalidFuelLifetimeError, 'Extend the depletion schedule'):
                    require_cost_eligible_lifetime(params)

    def test_positive_legacy_app_lifetimes_remain_eligible(self):
        require_cost_eligible_lifetime({'Fuel Lifetime': 1786})
        require_cost_eligible_lifetime({})

    def test_invalid_lifetimes_cannot_bypass_guard_without_status(self):
        for lifetime in (0, -1, math.nan, math.inf, None, 'bad'):
            with self.subTest(lifetime=lifetime):
                with self.assertRaisesRegex(InvalidFuelLifetimeError, 'Cost estimation will not run'):
                    require_cost_eligible_lifetime({'Fuel Lifetime': lifetime})

    def test_explicit_invalid_or_inconsistent_metadata_is_rejected(self):
        cases = [
            {'Fuel Lifetime': 100, 'Fuel Lifetime Status': 'unknown'},
            {'Fuel Lifetime': 100, 'Fuel Lifetime Status': None},
            {'Fuel Lifetime': 100, 'Fuel Lifetime Extrapolated': True},
            {'Fuel Lifetime': 100, 'Fuel Lifetime Status': INTERPOLATED, 'Fuel Lifetime Extrapolated': True},
            {'Fuel Lifetime Status': INTERPOLATED},
            {'Fuel Lifetime Extrapolated': False},
            {'Fuel Lifetime Simulated Endpoint (days)': 100},
        ]
        for params in cases:
            with self.subTest(params=params):
                with self.assertRaises(InvalidFuelLifetimeError):
                    require_cost_eligible_lifetime(params)


if __name__ == '__main__':
    unittest.main()
