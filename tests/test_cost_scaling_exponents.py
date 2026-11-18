import unittest

import numpy as np
import pandas as pd

from cost.cost_scaling import scale_central_facility_cost, scale_cost


class BlankScalingExponentTests(unittest.TestCase):
    @staticmethod
    def database():
        rows = []
        for account, exponent in ((100.1, 0.5), (100.2, np.nan)):
            rows.append({
                'Account': account,
                'Level': 1,
                'Account Title': f'Account {account}',
                'FOAK to NOAK Multiplier Type': 'No Learning',
                'Fixed Cost ($)': np.nan,
                'Fixed Cost Low End': np.nan,
                'Fixed Cost High End': np.nan,
                'Fixed Cost Distribution': None,
                'Adjusted Fixed Cost ($)': np.nan,
                'Adjusted Fixed Cost Low End ($)': np.nan,
                'Adjusted Fixed Cost High End ($)': np.nan,
                'Unit Cost': 10.0,
                'Unit Cost Low End': 10.0,
                'Unit Cost High End': 10.0,
                'Unit Cost Distribution': None,
                'Adjusted Unit Cost ($)': 10.0,
                'Adjusted Unit Cost Low End ($)': 10.0,
                'Adjusted Unit Cost High End ($)': 10.0,
                'Scaling Variable': 'Quantity',
                'Scaling Variable Ref Value': 2.0,
                'Count Scaling Variable': np.nan,
                'Count per Variable': np.nan,
                'Exponent': exponent,
                'Exponent Min': exponent,
                'Exponent Max': exponent,
                'Exponent std': 0.0,
                'Exponent Distribution': None,
                'Standard Cost Equation?': 'standard',
            })
        return pd.DataFrame(rows)

    def test_blank_exponent_defaults_to_one_in_both_scaling_paths(self):
        params = {
            'Escalation Year': 2025,
            'Number of Samples': 1,
            'Quantity': 4.0,
        }
        output_column = 'FOAK Estimated Cost ($2025)'

        for scaling_function in (scale_cost, scale_central_facility_cost):
            with self.subTest(scaling_function=scaling_function.__name__):
                result = scaling_function(self.database(), params.copy())

                # The first row deliberately uses exponent 0.5. The blank
                # second row must independently use exponent 1.0: 10 * 4.
                self.assertAlmostEqual(result.loc[1, output_column], 40.0)
                self.assertNotAlmostEqual(result.loc[0, output_column], 40.0)


if __name__ == '__main__':
    unittest.main()
