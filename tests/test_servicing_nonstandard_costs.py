import unittest
from pathlib import Path

import pandas as pd

from cost.cost_scaling import NONSTANDARD_COST_ACCOUNTS, non_standard_cost_scale


DATABASE = Path(__file__).resolve().parents[1] / "cost" / "Cost_Database.xlsx"


class ServicingNonstandardCostTest(unittest.TestCase):
    def setUp(self):
        self.params = {
            "FTEs Per Onsite Operator Per Year": 2.5,
            "FTEs Per Offsite Operator (24/7)": 5.0,
            "FTEs Per Security Staff (24/7)": 5.0,
        }

    def test_account_711_uses_old_operator_equation(self):
        self.assertEqual(
            non_standard_cost_scale(711, 100.0, 4.0, 1.0, self.params),
            2.5 * 100.0 * 4.0,
        )

    def test_account_712_preserves_old_reactor_equation(self):
        self.assertEqual(
            non_standard_cost_scale(712, 100.0, 8.0, 1.0, self.params),
            5.0 * 100.0 * (1.0 / 8.0),
        )

    def test_account_712_uses_servicing_database_fleet_variables(self):
        self.assertEqual(
            non_standard_cost_scale(
                712,
                100.0,
                1000.0,
                1.0,
                self.params,
                count_variable_value=10.0,
            ),
            100.0 * (1000.0 / 10.0) * 5.0,
        )

    def test_account_712_rejects_zero_reactors_per_operator(self):
        with self.assertRaisesRegex(
            ValueError, "Reactors Monitored Per Operator.*greater than zero"
        ):
            non_standard_cost_scale(
                712,
                100.0,
                1000.0,
                1.0,
                self.params,
                count_variable_value=0.0,
            )

    def test_account_713_uses_old_security_equation(self):
        self.assertEqual(
            non_standard_cost_scale(713, 100.0, 2.0, 1.0, self.params),
            5.0 * 100.0 * 2.0,
        )

    def test_account_714_converts_per_shift_staff_to_headcount(self):
        params = dict(self.params)
        params["Shift To Headcount"] = 5.0
        self.assertEqual(
            non_standard_cost_scale(714, 100.0, 40.0, 1.0, params),
            5.0 * 100.0 * 40.0,
        )

    def test_account_715_uses_total_engineering_headcount(self):
        self.assertEqual(
            non_standard_cost_scale(715, 100.0, 80.0, 1.0, self.params),
            100.0 * 80.0,
        )

    def test_account_744_uses_annual_electricity_equation(self):
        self.assertEqual(
            non_standard_cost_scale(744, 75.0, 12.0, 1.0, self.params),
            75.0 * 12.0 * 365 * 24,
        )

    def test_missing_nonstandard_equation_raises_clear_error(self):
        with self.assertRaisesRegex(
            NotImplementedError,
            r"Nonstandard cost equation is missing for account 999",
        ):
            non_standard_cost_scale(999, 100.0, 1.0, 1.0, self.params)

    def test_every_servicing_nonstandard_account_has_a_handler(self):
        database = pd.read_excel(DATABASE, sheet_name="Servicing Campus Database")
        accounts = set(
            database.loc[
                database["Standard Cost Equation?"] == "nonstandard", "Account"
            ]
        )
        self.assertEqual(accounts - NONSTANDARD_COST_ACCOUNTS, set())


if __name__ == "__main__":
    unittest.main()
