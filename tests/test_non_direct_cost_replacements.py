import unittest

import pandas as pd

from cost.non_direct_cost import calculate_accounts_31_32_75_82_cost


class ReplacementAccountTests(unittest.TestCase):
    def test_replacements_and_maintenance_are_separate(self):
        foak = 'FOAK Estimated Cost ($/unit)'
        noak = 'NOAK Estimated Cost ($/unit)'
        account_costs = {
            20: 1000.0,
            21: 10.0,
            22: 20.0,
            23: 30.0,
            25: 0.0,
            31: 0.0,
            32: 0.0,
            82: 0.0,
            221.12: 100.0,
            221.13: 200.0,
            221.2: 75.0,
            221.31: 25.0,
            221.33: 50.0,
            221.34: 150.0,
            212: 100.0,
            213: 50.0,
            214: 50.0,
            215: 50.0,
            222.1: 30.0,
            222.2: 20.0,
            222.3: 20.0,
            222.6: 30.0,
            223: 10.0,
            221.11: 10.0,
            226: 10.0,
            228: 10.0,
            227: 5.0,
            23: 2.0,
            24: 2.0,
            26: 1.0,
            751: 0.0,
            752: 0.0,
            753: 0.0,
            754: 0.0,
            755: 0.0,
            757: 0.0,
            759: 0.0,
            759.1: 0.0,
            759.2: 0.0,
            759.3: 0.0,
            759.4: 0.0,
            759.5: 0.0,
            759.6: 0.0,
            759.7: 0.0,
            759.8: 0.0,
        }
        data = pd.DataFrame({
            'Account': list(account_costs),
            foak: list(account_costs.values()),
            noak: list(account_costs.values()),
        })
        params = {
            'indirect to direct field-related cost': 0.0,
            'Fuel Lifetime': 365.0,
            'Refueling Period': 0.0,
            'Startup Duration after Refueling': 0.0,
            'A75: Outer Vessel Structure Replacement Period (years)': 20.0,
            'A75: Inner Vessel Structure Replacement Period (cycles)': 10.0,
            'A75: Reflector Replacement Period (cycles)': 10.0,
            'A75: Reactor Control Devices Replacement Period (cycles)': 10.0,
            'A75: Moderator Booster Replacement Period (cycles)': 1.0,
            'Discount Rate': 0.0,
            'Maintenance to Direct Cost Ratio': 0.01,
        }

        result = calculate_accounts_31_32_75_82_cost(data, params)

        for estimated_cost_column in (foak, noak):
            booster = result.loc[
                result['Account'] == 757, estimated_cost_column
            ].iat[0]
            maintenance = result.loc[
                result['Account'].isin([
                    759.1, 759.2, 759.3, 759.4,
                    759.5, 759.6, 759.7, 759.8,
                ]),
                estimated_cost_column,
            ].tolist()
            maintenance_total = result.loc[
                result['Account'] == 759, estimated_cost_column
            ].iat[0]
            self.assertEqual(booster, 150.0)
            self.assertEqual(
                maintenance,
                [2.5, 1.0, 0.1, 0.3, 0.05, 0.02, 0.02, 0.01],
            )
            self.assertEqual(maintenance_total, 4.0)


if __name__ == '__main__':
    unittest.main()
