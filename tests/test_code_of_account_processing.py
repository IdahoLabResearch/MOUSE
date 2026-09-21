import unittest

import pandas as pd

from cost.code_of_account_processing import remove_irrelevant_account


class OptionalAccountSelectionTests(unittest.TestCase):
    def _gcmr_purification_account(self):
        # CSV mirrors serialize the workbook's boolean cell as the string
        # "True", while the app supplies a Python bool in params.
        return pd.DataFrame([{
            'Account': 226,
            'Account Name': 'Other Reactor Plant Equipment (Helium purification)',
            'Account Title': 'Other Reactor Plant Equipment (Helium purification)',
            'Optional Variable': 'reactor type',
            'Optional Value': 'GCMR',
            'Sec Optional Variable': 'Primary Loop Purification',
            'Sec Optional Value': 'True',
        }])

    def test_gcmr_purification_retains_account_226(self):
        result = remove_irrelevant_account(
            self._gcmr_purification_account(),
            {'reactor type': 'GCMR', 'Primary Loop Purification': True},
        )

        self.assertEqual(result['Account'].tolist(), [226])

    def test_gcmr_without_purification_excludes_account_226(self):
        result = remove_irrelevant_account(
            self._gcmr_purification_account(),
            {'reactor type': 'GCMR', 'Primary Loop Purification': False},
        )

        self.assertTrue(result.empty)


if __name__ == '__main__':
    unittest.main()
