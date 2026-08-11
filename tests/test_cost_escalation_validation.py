import unittest
from pathlib import Path

from cost.cost_escalation import calculate_inflation_multiplier


DATABASE = Path(__file__).resolve().parents[1] / "cost" / "Cost_Database.xlsx"


class CostEscalationValidationTest(unittest.TestCase):
    def test_exact_inflation_type_is_accepted(self):
        self.assertAlmostEqual(
            calculate_inflation_multiplier(DATABASE, 2024, "General", 2025),
            1.03,
        )

    def test_missing_type_raises(self):
        with self.assertRaisesRegex(ValueError, "Inflation Type is missing"):
            calculate_inflation_multiplier(DATABASE, 2024, None, 2025)

    def test_type_must_match_inflation_sheet_exactly(self):
        with self.assertRaisesRegex(ValueError, "Inflation Type 'general' is invalid"):
            calculate_inflation_multiplier(DATABASE, 2024, "general", 2025)

    def test_non_text_type_is_invalid(self):
        with self.assertRaisesRegex(ValueError, "Inflation Type '0' is invalid"):
            calculate_inflation_multiplier(DATABASE, 2024, 0, 2025)

    def test_missing_dollar_year_raises(self):
        with self.assertRaisesRegex(ValueError, "Dollar Year is missing"):
            calculate_inflation_multiplier(DATABASE, None, "Labor", 2025)


if __name__ == "__main__":
    unittest.main()
