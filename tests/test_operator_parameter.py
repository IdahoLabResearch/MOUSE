import unittest

from cost.cost_scaling import non_standard_cost_scale
from reactor_engineering_evaluation.operation import reactor_operation


class OperatorParameterTests(unittest.TestCase):
    def base_params(self, operation_mode):
        return {
            'Operation Mode': operation_mode,
            'Number of On-Site Operators per Shift': 2,
            'Levelization Period': 60,
            'Refueling Period': 7,
            'Fuel Lifetime': 365,
            'Startup Duration after Refueling': 2,
            'Startup Duration after Emergency Shutdown': 14,
            'Emergency Shutdowns Per Year': 0.2,
            'Work Hours Per Shift': 10,
            'Hours Per FTE': 1800,
            'FTEs Per Onsite Operator (24/7)': 5,
            'Power MWe': 10,
        }

    def account_711_cost(self, params):
        reactor_operation(params)
        return non_standard_cost_scale(
            account=711,
            unit_cost=178500,
            scaling_variable_value=params['Number of On-Site Operators per Shift'],
            exponent=1,
            params=params,
        )

    def test_account_711_reflects_operation_mode(self):
        onsite = self.account_711_cost(self.base_params('On-Site Staffed'))
        remote = self.account_711_cost(self.base_params('Remotely Monitored'))

        self.assertEqual(onsite, 2 * 5 * 178500)
        self.assertGreater(remote, 0)
        self.assertLess(remote, onsite)

    def test_legacy_operator_parameter_is_mapped(self):
        params = self.base_params('Remotely Monitored')
        params['Number of Operators'] = params.pop('Number of On-Site Operators per Shift')

        reactor_operation(params)

        self.assertEqual(params['Number of On-Site Operators per Shift'], 2)


if __name__ == '__main__':
    unittest.main()
