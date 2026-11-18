import unittest

import pandas as pd

from cost.non_direct_cost import calculate_accounts_31_32_75_82_cost


class ReplacementAccountTests(unittest.TestCase):
    components = (
        ('Inner Vessel Structure', 752, 200.0),
        ('Moderator', 753, 50.0),
        ('Reflector', 754, 25.0),
        ('Reactor Control Devices', 755, 75.0),
        ('Moderator Booster', 757, 150.0),
    )
    foak = 'FOAK Estimated Cost ($/unit)'
    noak = 'NOAK Estimated Cost ($/unit)'

    def example(self):
        foak = 'FOAK Estimated Cost ($/unit)'
        noak = 'NOAK Estimated Cost ($/unit)'
        account_costs = {
            20: 1100.0,
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
            'A75: Moderator Replacement Period (cycles)': 2.0,
            'A75: Moderator Booster Replacement Period (cycles)': 1.0,
            'Discount Rate': 0.0,
            'Maintenance to Direct Cost Ratio': 0.01,
        }
        return data, params

    def replacement_key(self, component, units):
        return f'A75: {component} Replacement Period ({units})'

    def assert_replacement_periods(self, result, periods, rate=0.0, noak_factor=1.0):
        for component, account, capital in self.components:
            period = periods[component]
            factor = (1.0 / period if rate == 0.0 else
                      rate / (1.0 - (1.0 + rate) ** -period))
            for column, multiplier in ((self.foak, 1.0), (self.noak, noak_factor)):
                with self.subTest(component=component, column=column):
                    actual = result.loc[result['Account'] == account, column].iat[0]
                    self.assertAlmostEqual(actual, capital * multiplier * factor)

    def test_replacements_and_maintenance_are_separate(self):
        data, params = self.example()
        foak, noak = self.foak, self.noak

        result = calculate_accounts_31_32_75_82_cost(data, params)

        for estimated_cost_column in (foak, noak):
            booster = result.loc[
                result['Account'] == 757, estimated_cost_column
            ].iat[0]
            moderator = result.loc[
                result['Account'] == 753, estimated_cost_column
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
            self.assertEqual(moderator, 25.0)
            self.assertEqual(
                maintenance,
                [2.5, 1.0, 0.1, 0.3, 0.05, 0.02, 0.02, 0.01],
            )
            self.assertEqual(maintenance_total, 4.0)

    def test_year_only_inputs_work_for_all_five_components_and_both_cost_columns(self):
        data, params = self.example()
        data[self.noak] *= 1.5
        periods = {}
        for index, (component, _, _) in enumerate(self.components):
            params.pop(self.replacement_key(component, 'cycles'))
            periods[component] = 5.0 + index
            params[self.replacement_key(component, 'years')] = periods[component]
        result = calculate_accounts_31_32_75_82_cost(data, params)
        self.assert_replacement_periods(result, periods, noak_factor=1.5)

    def test_cycles_include_operation_refueling_and_restart_time(self):
        data, params = self.example()
        params.update({'Fuel Lifetime': 900.0, 'Refueling Period': 100.0,
                       'Startup Duration after Refueling': 95.0})
        periods = {}
        for component, _, _ in self.components:
            cycles = 3.0 if component == 'Moderator' else 2.0
            params[self.replacement_key(component, 'cycles')] = cycles
            periods[component] = 3.0 * cycles  # 1095 / 365 years per complete cycle
        result = calculate_accounts_31_32_75_82_cost(data, params)
        self.assert_replacement_periods(result, periods)

    def test_calendar_periods_override_cycles_and_ignore_full_cycle_duration(self):
        for lifetime, refueling, startup in ((365.0, 0.0, 0.0), (900.0, 100.0, 95.0)):
            with self.subTest(lifetime=lifetime):
                data, params = self.example()
                params.update({'Fuel Lifetime': lifetime, 'Refueling Period': refueling,
                               'Startup Duration after Refueling': startup})
                periods = {}
                for component, _, _ in self.components:
                    params[self.replacement_key(component, 'years')] = 10.0
                    params[self.replacement_key(component, 'cycles')] = 0.0
                    periods[component] = 10.0
                result = calculate_accounts_31_32_75_82_cost(data, params)
                self.assert_replacement_periods(result, periods)

    def test_units_can_be_mixed_independently_by_component(self):
        data, params = self.example()
        params[self.replacement_key('Inner Vessel Structure', 'years')] = 4.0
        params[self.replacement_key('Moderator', 'years')] = 5.0
        params[self.replacement_key('Reactor Control Devices', 'years')] = 6.0
        params[self.replacement_key('Reflector', 'cycles')] = 2.0
        params[self.replacement_key('Moderator Booster', 'cycles')] = 3.0
        result = calculate_accounts_31_32_75_82_cost(data, params)
        self.assert_replacement_periods(result, {
            'Inner Vessel Structure': 4.0, 'Moderator': 5.0,
            'Reflector': 2.0, 'Reactor Control Devices': 6.0,
            'Moderator Booster': 3.0,
        })

    def test_moderator_and_booster_retain_one_cycle_fallback(self):
        data, params = self.example()
        for component in ('Moderator', 'Moderator Booster'):
            params.pop(self.replacement_key(component, 'cycles'))
        result = calculate_accounts_31_32_75_82_cost(data, params)
        for account, capital in ((753, 50.0), (757, 150.0)):
            for column in (self.foak, self.noak):
                self.assertEqual(result.loc[result['Account'] == account, column].iat[0], capital)

    def test_calendar_periods_use_capital_recovery_factor_with_discounting(self):
        data, params = self.example()
        params['Discount Rate'] = 0.07
        periods = {}
        for component, _, _ in self.components:
            params[self.replacement_key(component, 'years')] = 8.5
            periods[component] = 8.5
        result = calculate_accounts_31_32_75_82_cost(data, params)
        self.assert_replacement_periods(result, periods, rate=0.07)

    def test_invalid_years_or_cycles_are_rejected_for_each_component(self):
        for component, _, _ in self.components:
            for units in ('years', 'cycles'):
                values = [-1.0, float('nan'), float('inf'), 'bad', None]
                if units == 'years':
                    values.append(0.0)
                for value in values:
                    with self.subTest(component=component, units=units, value=value):
                        data, params = self.example()
                        # For year inputs, the existing valid cycle input must
                        # not silently replace an explicitly invalid override.
                        params[self.replacement_key(component, units)] = value
                        with self.assertRaisesRegex(ValueError, component):
                            calculate_accounts_31_32_75_82_cost(data, params)

    def test_legacy_zero_cycles_still_disable_replacement_cost(self):
        for rate in (0.0, 0.07):
            for component, account, _ in self.components:
                with self.subTest(rate=rate, component=component):
                    data, params = self.example()
                    params['Discount Rate'] = rate
                    params[self.replacement_key(component, 'cycles')] = 0.0
                    result = calculate_accounts_31_32_75_82_cost(data, params)
                    for column in (self.foak, self.noak):
                        self.assertEqual(result.loc[result['Account'] == account, column].iat[0], 0.0)

    def test_required_component_period_does_not_acquire_a_hidden_default(self):
        data, params = self.example()
        params.pop(self.replacement_key('Inner Vessel Structure', 'cycles'))
        with self.assertRaisesRegex(KeyError, 'Inner Vessel Structure'):
            calculate_accounts_31_32_75_82_cost(data, params)


if __name__ == '__main__':
    unittest.main()
