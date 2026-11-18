import unittest
from unittest.mock import patch

from reactor_engineering_evaluation.vessels_calcs import vessels_specs


class VesselDimensionTests(unittest.TestCase):
    @patch('reactor_engineering_evaluation.vessels_calcs.materials_densities', return_value=1.0)
    @patch('reactor_engineering_evaluation.vessels_calcs.circle_area', side_effect=lambda radius: radius ** 2)
    @patch('reactor_engineering_evaluation.vessels_calcs.ellipsoid_shell', return_value=1.0)
    def test_nested_vessel_radii_include_preceding_wall_thickness(self, *_):
        params = {
            'reactor type': 'LTMR',
            'Active Height': 100.0,
            'Axial Reflector Thickness': 10.0,
            'Vessel Lower Plenum Height': 5.0,
            'Vessel Upper Plenum Height': 6.0,
            'Vessel Upper Gas Gap': 7.0,
            'Vessel Radius': 50.0,
            'Vessel Thickness': 2.0,
            'Vessel Bottom Depth': 8.0,
            'Vessel Material': 'steel',
            'Gap Between Vessel And Guard Vessel': 3.0,
            'Guard Vessel Thickness': 4.0,
            'Guard Vessel Material': 'steel',
            'Gap Between Guard Vessel And Cooling Vessel': 5.0,
            'Cooling Vessel Thickness': 6.0,
            'Cooling Vessel Material': 'steel',
            'Gap Between Cooling Vessel And Intake Vessel': 7.0,
            'Intake Vessel Thickness': 8.0,
            'Intake Vessel Material': 'steel',
        }

        vessels_specs(params)

        self.assertEqual(params['Guard Vessel Radius'], 55.0)
        self.assertEqual(params['Cooling Vessel Radius'], 64.0)
        self.assertEqual(params['Intake Vessel Radius'], 77.0)
        self.assertEqual(params['Vessels Total Radius'], 85.0)


if __name__ == '__main__':
    unittest.main()
