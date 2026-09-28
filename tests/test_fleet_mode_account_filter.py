import unittest

import pandas as pd

from cost.cost_estimation import remove_fleet_mode_reactor_accounts


class FleetModeAccountFilterTests(unittest.TestCase):
    def test_excluding_reactivity_control_leaves_does_not_remove_parent_hierarchy(self):
        cost_table = pd.DataFrame({
            'Account': [221.2, 221.21, 221.211, 221.212, 221.213, 221.3],
            'Level': [3, 4, 5, 5, 5, 3],
            'Account Title': [
                'Reactor Control Devices',
                'Reactivity Control System',
                'Fabrication',
                'Installation',
                'Control Drum Materials',
                'Non-Fuel Core Internals',
            ],
        })

        result = remove_fleet_mode_reactor_accounts(cost_table)

        self.assertEqual(
            result['Account'].tolist(),
            [221.2, 221.21, 221.213, 221.3],
        )


if __name__ == '__main__':
    unittest.main()
