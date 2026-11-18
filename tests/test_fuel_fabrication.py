import math
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from cost.cost_estimation import FOAK_to_NOAK
from cost.cost_database_csv import read_cost_database_sheet
from cost.cost_escalation import escalate_cost_database
from cost.cost_scaling import scale_cost, scale_central_facility_cost
from cost.sampling import sampler


class FuelLearningTests(unittest.TestCase):
    def params(self, unit):
        return {
            'NOAK Unit Number': unit,
            'No Learning': 0.0,
            'Licensing Learning': 0.2,
            'Factory Primary Structure': 0.184432,
            'Factory Drums': 0.24138,
            'Factory Other': 0.2296284802211553,
            'Factory Be': 0.415054,
            'Factory BeO': 0.399149,
            'Non-nuclear off-the-shelf': 0.05,
            'Onsite Learning': 0.08,
        }

    def frame(self):
        return pd.DataFrame({
            'FOAK to NOAK Multiplier Type': [
                'UO2 Learning', 'TRISO Learning', 'Factory Other',
            ],
            'FOAK Estimated Cost ($2025)': [400.0, 10900.0, 1520.0],
        })

    def test_first_unit_and_hundredth_unit_cap(self):
        for unit in (1, 100, 200):
            with self.subTest(unit=unit):
                result = FOAK_to_NOAK(self.frame(), self.params(unit))
                for index, rate in enumerate((0.02, 0.13, 0.2296284802211553)):
                    expected = (1 - rate) ** math.log2(min(unit, 100))
                    self.assertAlmostEqual(result.loc[index, 'Multiplier'], expected)
                    self.assertAlmostEqual(
                        result.loc[index, 'NOAK Estimated Cost ($2025)'],
                        result.loc[index, 'FOAK Estimated Cost ($2025)'] * expected,
                    )

    def test_scenario_learning_override_wins(self):
        params = self.params(100)
        params.update({'UO2 Learning': 0.0, 'TRISO Learning': 0.05})
        result = FOAK_to_NOAK(self.frame(), params)
        self.assertEqual(result.loc[0, 'Multiplier'], 1.0)
        self.assertAlmostEqual(result.loc[1, 'Multiplier'], 0.95 ** math.log2(100))


class TriangularFabricationTests(unittest.TestCase):
    def test_dispatch_uses_low_mode_high_order(self):
        with patch('cost.sampling.np.random.triangular', return_value=410.0) as draw:
            self.assertEqual(sampler('Triangular', low_cost=230, class3_cost=400, high_cost=575), 410.0)
            draw.assert_called_once_with(230, 400, 575)

    def test_samples_obey_source_bounds_and_mean(self):
        state = np.random.get_state()
        try:
            np.random.seed(904)
            for low, mode, high in ((230, 400, 575), (3300, 10900, 29400)):
                draws = np.array([
                    sampler('Triangular', low_cost=low, class3_cost=mode, high_cost=high)
                    for _ in range(4000)
                ])
                self.assertTrue(((draws >= low) & (draws <= high)).all())
                self.assertAlmostEqual(draws.mean(), (low + mode + high) / 3, delta=(high - low) * 0.02)
        finally:
            np.random.set_state(state)

    def test_invalid_mode_is_rejected(self):
        with self.assertRaises(ValueError):
            sampler('Triangular', low_cost=230, class3_cost=600, high_cost=575)

    def frame(self, cost_kind):
        fixed = cost_kind == 'fixed'
        row = {
            'Account': 254, 'Level': 2, 'Account Title': 'Fabrication',
            'FOAK to NOAK Multiplier Type': 'UO2 Learning',
            'Fixed Cost ($)': 400.0 if fixed else np.nan,
            'Fixed Cost Low End': 230.0 if fixed else np.nan,
            'Fixed Cost High End': 575.0 if fixed else np.nan,
            'Fixed Cost Distribution': 'Triangular',
            'Unit Cost': np.nan if fixed else 400.0,
            'Unit Cost Low End': np.nan if fixed else 230.0,
            'Unit Cost High End': np.nan if fixed else 575.0,
            'Unit Cost Distribution': 'Triangular',
            'Adjusted Fixed Cost ($)': 400.0 if fixed else np.nan,
            'Adjusted Fixed Cost Low End ($)': 230.0 if fixed else np.nan,
            'Adjusted Fixed Cost High End ($)': 575.0 if fixed else np.nan,
            'Adjusted Unit Cost ($)': np.nan if fixed else 400.0,
            'Adjusted Unit Cost Low End ($)': np.nan if fixed else 230.0,
            'Adjusted Unit Cost High End ($)': np.nan if fixed else 575.0,
            'Scaling Variable': 'Uranium Mass',
            'Scaling Variable Ref Value': np.nan,
            'Count Scaling Variable': np.nan, 'Count per Variable': np.nan,
            'Exponent': 1.0, 'Exponent Min': 1.0, 'Exponent Max': 1.0,
            'Exponent std': 0.0, 'Exponent Distribution': None,
            'Standard Cost Equation?': 'standard',
        }
        return pd.DataFrame([row])

    def test_reactor_and_central_paths_sample_instead_of_using_mode(self):
        for scale in (scale_cost, scale_central_facility_cost):
            for kind in ('unit', 'fixed'):
                with self.subTest(scale=scale.__name__, kind=kind):
                    params = {'Escalation Year': 2025, 'Number of Samples': 2, 'Uranium Mass': 10}
                    with patch('cost.cost_scaling.sampler', return_value=450.0) as draw:
                        result = scale(self.frame(kind), params)
                    draw.assert_called_once_with('Triangular', low_cost=230.0, high_cost=575.0, class3_cost=400.0)
                    self.assertEqual(result.loc[0, 'FOAK Estimated Cost ($2025)'], 4500.0 if kind == 'unit' else 450.0)

    def test_single_sample_uses_mode_for_both_cost_paths(self):
        for scale in (scale_cost, scale_central_facility_cost):
            params = {'Escalation Year': 2025, 'Number of Samples': 1, 'Uranium Mass': 10}
            with patch('cost.cost_scaling.sampler') as draw:
                result = scale(self.frame('unit'), params)
            draw.assert_not_called()
            self.assertEqual(result.loc[0, 'FOAK Estimated Cost ($2025)'], 4000.0)


class FabricationWorkbookIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.directory = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.directory.cleanup)
        cls.workbook = Path(cls.directory.name) / 'Cost_Database.xlsx'
        shutil.copy2(
            Path(__file__).resolve().parents[1] / 'cost' / 'Cost_Database.xlsx',
            cls.workbook,
        )
        cls.raw = read_cost_database_sheet(cls.workbook, 'Cost Database')

    def params(self, unit):
        params = {'Escalation Year': 2025, 'Number of Samples': 1,
                  'Uranium Mass': 10.0, 'NOAK Unit Number': unit}
        # Escalation resolves every account's parameter-valued cost fields.
        # Their unrelated values do not enter the two fabrication rows below.
        for column in ('Fixed Cost ($)', 'Fixed Cost Low End', 'Fixed Cost High End',
                       'Unit Cost', 'Unit Cost Low End', 'Unit Cost High End'):
            for value in self.raw[column]:
                if isinstance(value, str):
                    params.setdefault(value, 1.0)
        return params

    def fabrication(self, params):
        escalated = escalate_cost_database(self.workbook, 2025, params)
        return escalated.loc[
            (escalated['Account'] == 254)
            & escalated['FOAK to NOAK Multiplier Type'].isin(['UO2 Learning', 'TRISO Learning'])
        ].copy().reset_index(drop=True)

    def test_workbook_inputs_escalation_and_deployment_learning(self):
        expected = {
            'UO2 Learning': (230.0, 400.0, 575.0, 0.02, 610.701200887318),
            'TRISO Learning': (3300.0, 10900.0, 29400.0, 0.13, 16641.607724179416),
        }
        for unit in (1, 100):
            params = self.params(unit)
            rows = self.fabrication(params)
            self.assertEqual(len(rows), 2)
            for _, row in rows.iterrows():
                category = row['FOAK to NOAK Multiplier Type']
                low, mode, high, rate, escalated_mode = expected[category]
                self.assertEqual(row['Dollar Year'], 2017)
                self.assertEqual(row['Unit Cost Distribution'], 'Triangular')
                self.assertEqual(row['Scaling Variable'], 'Uranium Mass')
                self.assertEqual((row['Unit Cost Low End'], row['Unit Cost'], row['Unit Cost High End']), (low, mode, high))
                self.assertEqual(params[category], rate)
                self.assertAlmostEqual(row['Adjusted Unit Cost ($)'], escalated_mode)
            result = FOAK_to_NOAK(scale_cost(rows, params), params)
            for _, row in result.iterrows():
                _, _, _, rate, escalated_mode = expected[row['FOAK to NOAK Multiplier Type']]
                multiplier = (1.0 - rate) ** math.log2(unit)
                self.assertAlmostEqual(row['FOAK Estimated Cost ($2025)'], 10.0 * escalated_mode)
                self.assertAlmostEqual(row['Multiplier'], multiplier)
                self.assertAlmostEqual(row['NOAK Estimated Cost ($2025)'], 10.0 * escalated_mode * multiplier)

    def test_zero_learning_overrides_survive_database_loading(self):
        params = self.params(100)
        params.update({'UO2 Learning': 0.0, 'TRISO Learning': 0.0})
        rows = self.fabrication(params)
        result = FOAK_to_NOAK(scale_cost(rows, params), params)
        self.assertEqual(params['UO2 Learning'], 0.0)
        self.assertEqual(params['TRISO Learning'], 0.0)
        self.assertTrue((result['Multiplier'] == 1.0).all())
        np.testing.assert_allclose(result['NOAK Estimated Cost ($2025)'], result['FOAK Estimated Cost ($2025)'])


if __name__ == '__main__':
    unittest.main()
