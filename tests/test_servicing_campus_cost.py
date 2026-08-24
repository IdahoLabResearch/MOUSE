import runpy
import unittest
from pathlib import Path

import numpy as np

from cost.cost_estimation import bottom_up_cost_estimate_servicing_campus
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
