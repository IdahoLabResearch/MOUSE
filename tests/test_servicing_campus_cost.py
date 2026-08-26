import copy
import runpy
import unittest
from pathlib import Path

import numpy as np

from cost.cost_estimation import bottom_up_cost_estimate_servicing_campus
from cost.fleet_mode import (
    cumulative_unit_learning_multiplier,
    servicing_facility_allocation,
    servicing_facility_occ_learning_multipliers,
)
from cost.cost_scaling import non_standard_cost_scale


ROOT = Path(__file__).resolve().parents[1]
DATABASE = ROOT / "cost" / "Cost_Database.xlsx"
FLEET_PARAMS = ROOT / "assets" / "Cost_Database_Separate_Campuses define variables params.py"


def servicing_params():
    params = {
        "Power MWt": 15.0,
        "Thermal Efficiency": 0.4,
        "Annual Electricity Production": 15.0 * 0.4 * 0.9 * 365 * 24,
        "Annual Coolant Supply Frequency": 1,
        "Reactors Monitored Per Operator": 10,
        "Escalation Year": 2025,
        "Number of Samples": 1,
        "Interest Rate": 0.07,
        "Debt To Equity Ratio": 1.0,
        "Levelization Period": 60,
        "Discount Rate": 0.07,
        "Fleet Mode": True,
        "FTEs Per Onsite Operator Per Year": 1.0,
    }
    runpy.run_path(str(FLEET_PARAMS), init_globals={"params": params})
    return params


class ServicingCampusCostTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.params = servicing_params()
        cls.result = bottom_up_cost_estimate_servicing_campus(str(DATABASE), cls.params)

    def value(self, account):
        return self.result.loc[self.result["Account"] == account, "Value"].iloc[0]

    def test_servicing_rate_is_fleet_divided_by_cycle_length(self):
        self.assertAlmostEqual(
            self.params["Servicing Rate"],
            self.params["Fleet"] / self.params["Cycle Length"],
        )

    def test_servicing_facility_count_increases_above_1000_reactors(self):
        self.assertEqual(
            servicing_facility_allocation(1000, 100.0),
            (1, (1000,), 1000, 100.0),
        )
        self.assertEqual(
            servicing_facility_allocation(1001, 100.1),
            (2, (501, 500), 501, 50.05),
        )
        self.assertEqual(
            servicing_facility_allocation(10000, 1000.0),
            (10, (1000,) * 10, 1000, 100.0),
        )

    def test_servicing_facility_learning_propagates_to_occ_tci_and_lcoe(self):
        one_facility_params = copy.deepcopy(self.params)
        one_facility_params['Generating Sites Count'] = 1000
        one_facility_params['Servicing Facility Design Capacity'] = 1000
        one_facility = bottom_up_cost_estimate_servicing_campus(
            str(DATABASE), one_facility_params, sample_seeds=[12345]
        )

        two_facility_params = copy.deepcopy(one_facility_params)
        two_facility_params['Fleet'] = 2000
        two_facility_params['Generating Sites Count'] = 2000
        two_facility_params['Servicing Facility Count'] = 2
        two_facility_params['Servicing Facility Design Capacity'] = 1000
        two_facility_params['Servicing Rate'] *= 2
        two_facility_params['Servicing Rate Per Facility'] = (
            one_facility_params['Servicing Rate Per Facility']
        )
        two_facilities = bottom_up_cost_estimate_servicing_campus(
            str(DATABASE), two_facility_params, sample_seeds=[12345]
        )

        def result_value(result, account):
            return result.loc[result['Account'] == account, 'Value'].iloc[0]

        aggregate_learning_multiplier = sum(
            servicing_facility_occ_learning_multipliers(2, 0.30, 5)
        )
        for total_account in ['OCC', 'TCI']:
            self.assertAlmostEqual(
                result_value(two_facilities, total_account),
                aggregate_learning_multiplier
                * result_value(one_facility, total_account),
                delta=0.01,
            )
        for dependent_annual_account in [741, 742, 743, 747.1, 747.4]:
            self.assertAlmostEqual(
                result_value(two_facilities, dependent_annual_account),
                aggregate_learning_multiplier
                * result_value(one_facility, dependent_annual_account),
                delta=0.01,
            )
        for normalized_account in ['OCC per reactor', 'TCI per reactor']:
            self.assertAlmostEqual(
                result_value(two_facilities, normalized_account),
                aggregate_learning_multiplier
                / 2
                * result_value(one_facility, normalized_account),
                delta=0.01,
            )
        self.assertLess(
            result_value(two_facilities, 'LCOE'),
            result_value(one_facility, 'LCOE'),
        )

    def test_servicing_facility_occ_learning_stops_after_five_facilities(self):
        multipliers = servicing_facility_occ_learning_multipliers(10, 0.30, 5)
        self.assertAlmostEqual(multipliers[1], 0.70)
        self.assertAlmostEqual(multipliers[3], 0.49)
        self.assertTrue(all(value == multipliers[4] for value in multipliers[4:]))

    def test_cask_learning_adds_each_unit_and_stops_after_unit_100(self):
        exponent = np.log2(0.85)
        first_100 = sum(unit ** exponent for unit in range(1, 101))
        learned_182 = cumulative_unit_learning_multiplier(182, 0.15, 100)

        self.assertAlmostEqual(
            learned_182,
            first_100 + 82 * 100 ** exponent,
        )
        self.assertAlmostEqual(
            cumulative_unit_learning_multiplier(2, 0.15, 100),
            1 + 0.85,
        )

    def test_cask_learning_propagates_to_annual_servicing_cost(self):
        learned_params = copy.deepcopy(self.params)
        baseline_params = copy.deepcopy(self.params)
        baseline_params['Cask Learning Rate'] = 0
        learned = bottom_up_cost_estimate_servicing_campus(
            str(DATABASE), learned_params, sample_seeds=[24680]
        )
        baseline = bottom_up_cost_estimate_servicing_campus(
            str(DATABASE), baseline_params, sample_seeds=[24680]
        )

        def result_value(result, account):
            return result.loc[result['Account'] == account, 'Value'].iloc[0]

        for cask_account in [747.21, 747.22, 747.23]:
            self.assertLess(
                result_value(learned, cask_account),
                result_value(baseline, cask_account),
            )
        self.assertAlmostEqual(
            result_value(learned, 747.2),
            sum(result_value(learned, account) for account in [747.21, 747.22, 747.23]),
            delta=0.01,
        )
        self.assertAlmostEqual(
            result_value(baseline, 'Annual Cost')
            - result_value(learned, 'Annual Cost'),
            result_value(baseline, 747.2) - result_value(learned, 747.2),
            delta=0.01,
        )

    def test_radwaste_warehouse_area_is_converted_from_ft2_to_m2(self):
        servicing_scale = np.rint(self.params["Servicing Rate"] / 3)
        source_area_ft2 = 15422.94 * servicing_scale ** 0.994
        self.assertAlmostEqual(
            self.params["Radwaste Storage Warehouses Area"],
            source_area_ft2 / 3.2808 ** 2,
        )

    def test_output_uses_single_value_columns(self):
        self.assertIn("Value", self.result.columns)
        self.assertIn("Value std", self.result.columns)
        self.assertFalse(any("FOAK" in column or "NOAK" in column for column in self.result.columns))
        for account in ["OCC", "OCC per reactor", "TCI", "TCI per reactor",
                        "Annual Cost", "Annual Cost per reactor", "LCOE"]:
            self.assertTrue((self.result["Account"] == account).any())

    def test_database_ratios_populate_derived_accounts(self):
        direct_cost = self.value(20)
        staffing_cost = sum(self.value(account) for account in [711, 712, 713, 714, 715])

        self.assertAlmostEqual(self.value(30), 0.50 * direct_cost, places=4)
        self.assertAlmostEqual(self.value(716), 0.02 * staffing_cost, places=4)
        self.assertAlmostEqual(self.value(717), 0.075 * staffing_cost, places=4)
        self.assertAlmostEqual(self.value(741), 0.005 * direct_cost, places=4)
        self.assertAlmostEqual(self.value(742), 0.05 * direct_cost, places=4)
        self.assertAlmostEqual(self.value(743), 0.01 * direct_cost, places=4)
        self.assertAlmostEqual(
            self.value(747.1),
            0.05 * sum(self.value(account) for account in [741, 742, 743]),
            places=4,
        )
        self.assertAlmostEqual(self.value(747.4), 0.005 * direct_cost, places=4)

    def test_account_744_converts_average_power_to_annual_electricity_cost(self):
        unit_cost_per_mwh = 75.0
        average_power_mwe = 12.0
        expected = unit_cost_per_mwh * average_power_mwe * 365 * 24
        self.assertEqual(
            non_standard_cost_scale(
                744, unit_cost_per_mwh, average_power_mwe, 1.0, self.params
            ),
            expected,
        )

    def test_lcoe_uses_fleet_annual_electricity_production(self):
        years = self.params["Levelization Period"]
        rate = self.params["Discount Rate"]
        discount_sum = sum((1 + rate) ** -year for year in range(1, years + 1))
        annual_generation = (
            self.params["Fleet"] * self.params["Annual Electricity Production"]
        )
        expected = (
            self.value("TCI") + self.value("Annual Cost") * discount_sum
        ) / (annual_generation * discount_sum)
        self.assertTrue(np.isfinite(self.value("LCOE")))
        self.assertAlmostEqual(self.value("LCOE"), expected, places=8)

    def test_missing_annual_electricity_production_raises_clear_error(self):
        params = dict(self.params)
        del params["Annual Electricity Production"]
        with self.assertRaisesRegex(
            KeyError,
            "Annual Electricity Production.*reactor operation calculation",
        ):
            bottom_up_cost_estimate_servicing_campus(str(DATABASE), params)

    def test_fleet_mode_false_skips_servicing_estimate(self):
        params = dict(self.params)
        params["Fleet Mode"] = False
        self.assertIsNone(bottom_up_cost_estimate_servicing_campus(str(DATABASE), params))


if __name__ == "__main__":
    unittest.main()
