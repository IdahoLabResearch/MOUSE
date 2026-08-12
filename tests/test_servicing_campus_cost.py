import runpy
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from cost.cost_estimation import (
    bottom_up_cost_estimate_servicing_campus,
    create_servicing_campus_cost_dictionary,
    format_servicing_campus_output,
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

    def test_output_uses_single_value_columns(self):
        self.assertIn("Value", self.result.columns)
        self.assertIn("Value std", self.result.columns)
        self.assertFalse(any("FOAK" in column or "NOAK" in column for column in self.result.columns))
        for account in ["OCC", "OCC per reactor", "TCI", "TCI per reactor",
                        "Annual Cost", "Annual Cost per reactor", "LCOE"]:
            self.assertTrue((self.result["Account"] == account).any())

    def test_excel_output_has_dynamic_integer_cost_and_lcoe_columns(self):
        output = format_servicing_campus_output(self.result, self.params)
        cost_column = "Estimated Cost ($ 2025)"
        std_column = "Estimated Cost std ($ 2025)"
        lcoe_column = "LCOE Contribution ($/MWh)"

        self.assertNotIn("Value", output.columns)
        self.assertNotIn("Value std", output.columns)
        self.assertIn(cost_column, output.columns)
        self.assertIn(std_column, output.columns)
        self.assertIn(lcoe_column, output.columns)
        self.assertEqual(str(output[cost_column].dtype), "Int64")
        self.assertEqual(str(output[std_column].dtype), "Int64")
        self.assertEqual(str(output[lcoe_column].dtype), "Int64")

        account_10 = output.loc[output["Account"] == 10, lcoe_column].iloc[0]
        discount_sum = sum(
            (1 + self.params["Discount Rate"]) ** -year
            for year in range(1, self.params["Levelization Period"] + 1)
        )
        discounted_generation = (
            self.params["Fleet"]
            * self.params["Annual Electricity Production"]
            * discount_sum
        )
        self.assertEqual(account_10, int(self.value(10) / discounted_generation))

        account_711 = output.loc[output["Account"] == 711, lcoe_column].iloc[0]
        fleet_annual_generation = (
            self.params["Fleet"] * self.params["Annual Electricity Production"]
        )
        self.assertEqual(account_711, int(self.value(711) / fleet_annual_generation))

    def test_parametric_summary_tracks_servicing_values_without_foak_noak(self):
        summary = pd.DataFrame({
            "Account": [
                "OCC", "OCC per reactor", "TCI", "TCI per reactor",
                "Annual Cost", "Annual Cost per reactor", "LCOE",
            ],
            "Value": [100, 10, 120, 12, 20, 2, 3],
            "Value std": [5, 0.5, 6, 0.6, 1, 0.1, 0.2],
        })
        tracked = create_servicing_campus_cost_dictionary(summary)

        self.assertEqual(tracked["Servicing Campus OCC"], 100)
        self.assertEqual(tracked["Servicing Campus OCC std"], 5)
        self.assertEqual(tracked["Servicing Campus LCOE"], 3)
        self.assertEqual(tracked["Servicing Campus LCOE std"], 0.2)
        self.assertFalse(any("FOAK" in key or "NOAK" in key for key in tracked))

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
