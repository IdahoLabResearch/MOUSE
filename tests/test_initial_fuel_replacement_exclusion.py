import unittest

import pandas as pd

from cost.non_direct_cost import calculate_accounts_31_32_75_82_cost


class InitialFuelReplacementExclusionTests(unittest.TestCase):
    def test_initial_fuel_is_excluded_from_misc_replacements(self):
        foak = 'FOAK Estimated Cost ($/unit)'
        noak = 'NOAK Estimated Cost ($/unit)'
        account_costs = {
            20: 950.0,
            21: 10.0,
            22: 20.0,
            23: 30.0,
            25: 100.0,
            31: 0.0,
            32: 0.0,
            82: 0.0,
            221.12: 100.0,
            221.13: 200.0,
            221.2: 75.0,
            221.31: 25.0,
            221.33: 50.0,
            222.1: 0.0,
            222.2: 0.0,
            222.3: 0.0,
            222.61: 0.0,
            751: 0.0,
            752: 0.0,
            753: 0.0,
            754: 0.0,
            755: 0.0,
            756: 0.0,
            759: 0.0,
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
            'A75: Vessel Replacement Period (cycles)': 10.0,
            'A75: Core Barrel Replacement Period (cycles)': 10.0,
            'A75: Reflector Replacement Period (cycles)': 10.0,
            'A75: Drum Replacement Period (cycles)': 10.0,
            'A75: Integrated HX Replacement Period (cycles)': 0.0,
            'Discount Rate': 0.0,
            'Maintenance to Direct Cost Ratio': 0.01,
        }

        result = calculate_accounts_31_32_75_82_cost(data, params)

        for estimated_cost_column in (foak, noak):
            miscellaneous = result.loc[
                result['Account'] == 759, estimated_cost_column
            ].iat[0]
            annualized_fuel = result.loc[
                result['Account'] == 82, estimated_cost_column
            ].iat[0]
            # (950 total - 450 explicit replacements - 100 fuel) * 1%.
            self.assertEqual(miscellaneous, 4.0)
            self.assertEqual(annualized_fuel, 100.0)


if __name__ == '__main__':
    unittest.main()
