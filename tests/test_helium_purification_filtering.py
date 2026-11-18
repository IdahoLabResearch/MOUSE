import unittest

import pandas as pd

from cost.code_of_account_processing import remove_irrelevant_account


class HeliumPurificationFilteringTests(unittest.TestCase):
    def accounts(self):
        return pd.DataFrame([
            {
                'Account': 226,
                'Account Name': 'Helium purification',
                'Account Title': 'Helium purification',
                'Optional Variable': 'reactor type',
                'Optional Value': 'GCMR',
                'Sec Optional Variable': 'Primary Loop Purification',
                'Sec Optional Value': 'True',
            },
            {
                'Account': 226,
                'Account Name': 'LTMR equipment',
                'Account Title': 'LTMR equipment',
                'Optional Variable': 'reactor type',
                'Optional Value': 'LTMR',
                'Sec Optional Variable': None,
                'Sec Optional Value': None,
            },
        ])

    def test_boolean_parameter_matches_boolean_string_from_database(self):
        selected = remove_irrelevant_account(
            self.accounts(),
            {'reactor type': 'GCMR', 'Primary Loop Purification': True},
        )
        self.assertEqual(selected['Account Name'].tolist(), ['Helium purification'])

    def test_false_purification_flag_removes_gcmr_account(self):
        selected = remove_irrelevant_account(
            self.accounts(),
            {'reactor type': 'GCMR', 'Primary Loop Purification': False},
        )
        self.assertTrue(selected.empty)

    def test_string_boolean_parameter_is_case_insensitive(self):
        selected = remove_irrelevant_account(
            self.accounts(),
            {'reactor type': 'GCMR', 'Primary Loop Purification': ' true '},
        )
        self.assertEqual(selected['Account Name'].tolist(), ['Helium purification'])

    def test_ltmr_account_does_not_require_purification_flag(self):
        selected = remove_irrelevant_account(
            self.accounts(), {'reactor type': 'LTMR'}
        )
        self.assertEqual(selected['Account Name'].tolist(), ['LTMR equipment'])


if __name__ == '__main__':
    unittest.main()
