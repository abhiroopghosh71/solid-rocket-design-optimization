import os
import unittest

import numpy as np

from problems.rocket_propellant_design import rocket_single_eval


class TestRocketEvaluator(unittest.TestCase):
    """Consists of unit tests for the rocket burn simulator."""

    def test_input(self):
        """Tests the output for a specific input case."""
        x_test = np.array(
            [[1, 7, 6, 7, 4, 9, 3, 10, 0, 2, 4, 8, 2, 2, 2, 5, 5, 6, 4, 1, 7, 8, 10, 8, 1, 0, 2, 10, 10, 9, 10, 0, 2, 5,
              10, 1, 2, 10, 0, 0, 7, 6, 10, 9, 0, 8, 1, 9, 10, 10, 10, 6, 10, 5, 1, 4, 4, 0, 1, 0, 9, 9, 10, 10, 7, 9,
              10, 2, 9, 0, 5, 2, 7, 6, 2, 0, 4, 0, 0, 6, 10, 9, 2, 10, 5, 8, 10, 9, 10, 10, 7, 4, 6, 0, 1, 4, 4, 0, 0,
              4, 1, 7, 6, 7, 4, 5, 9, 8, 10, 0, 2, 7, 3, 6, 2, 2, 6, 6, 0, 0, 1, 7, 6, 7, 4, 10, 1, 9, 1, 3, 9, 4, 5, 2,
              0, 1, 10, 5, 0, 0, 1, 7, 6, 7, 4, 10, 3, 1, 10, 1, 7, 3, 2, 7, 10, 9, 2, 5, 0, 0, 1, 7, 6, 7, 4, 8, 9, 5,
              9, 6, 4, 7, 0, 5, 1, 1, 1, 2, 0, 0, 1, 7, 6, 7, 4, 3, 5, 4, 10, 4, 3, 9, 6, 10, 3, 2, 1, 9, 0, 0, 1, 7, 6,
              7, 4, 9, 0, 9, 5, 0, 8, 6, 4, 2, 2, 6, 8, 5, 0, 0, 1, 7, 6, 7, 4, 5, 0, 2, 7, 5, 2, 0, 10, 6, 2, 10, 4, 3,
              0, 0, 7, 8, 10, 8, 1, 2, 0, 10, 4, 3, 6, 3, 1, 4, 1, 0, 0, 0, 0, 0, 2, 1, 3, 3, 2, 2, 14, 9, 3, 3, 2, 2,
              2, 1, 2, 6, 3, 2, 2, 2, 0, 3, 6, 9]], dtype=int)
        for i in range(2):
            simulator_output_arr = rocket_single_eval.evaluate(x_test, target_thrust_profile=None, ray_depth_flag=0)
            # KLUGE
            simulator_output = simulator_output_arr[0]

            thrust_reward = simulator_output['thrust_reward']
            simultaneity_reward = simulator_output['simultaneity_reward']
            simulation_termination_flag = simulator_output['stop_code']
            time_data_coarse = simulator_output['time_data_coarse']
            thrust_profile_coarse = simulator_output['thrust_profile_coarse']
            pressure_profile_coarse = simulator_output['pressure_profile_coarse']
            time_data_fine = simulator_output['time_data_fine']
            thrust_profile_fine = simulator_output['thrust_profile_fine']
            pressure_profile_fine = simulator_output['pressure_profile_fine']
            unburned_fuel_data = simulator_output['unburned_fuel_data']
            ray_burn_depths = simulator_output['ray_burn_depths']
            seg_star_burn_status = simulator_output['seg_star_burn_status']
            circularization_time = simulator_output['circularization_time']
            seg_residual_per_time_step_coarse = simulator_output['seg_residual_per_time_step_coarse']
            seg_residual_per_time_step_fine = simulator_output['seg_residual_per_time_step_fine']

            burn_layer_per_timestep_coarse = simulator_output['burn_layer_per_timestep_coarse']
            burn_layer_per_timestep_fine = simulator_output['burn_layer_per_timestep_fine']

            # burn_layer_per_timestep_coarse = rocket_single_eval.get_burn_layer_per_timestep(
            #     seg_residual_per_time_step_coarse)
            # burn_layer_per_timestep_fine = rocket_single_eval.get_burn_layer_per_timestep(
            #     seg_residual_per_time_step_fine)

            self.assertEqual(thrust_reward, 321114.5209178659)
            self.assertEqual(simultaneity_reward, 428436.2400948988)
            self.assertEqual(simulation_termination_flag, 0)

            self.assertTrue(
                np.allclose(thrust_profile_coarse,
                            np.array([7238.34365998, 6355.11604302, 6428.12892106, 5944.83192139, 5765.10659394,
                                     5835.07203787, 5593.34972486, 5182.88571753, 4967.62473903, 4864.53466283,
                                     4902.6610781,  4636.60217346, 4130.13557409, 4075.20247541, 3855.60918429,
                                     3647.39638416, 3461.04689244, 3364.87654633, 3114.91404598, 3134.99595683,
                                     3101.27984383, 5437.7673139, 5671.36289391, 0., 0., 0., 0., 0., 0., 0., 0.,
                                     0., 0., 0., 0., 0., 0., 0., 0., 0.])))

            self.assertTrue(
                np.allclose(pressure_profile_coarse,
                            np.array([3289372.96273456, 2888013.80011734, 2921192.54931581, 2701571.15686828,
                                      2619899.78978036, 2651693.7140755, 2541849.47542816, 2355325.08793923,
                                      2257505.49487392, 2210658.96758795, 2227984.49705745, 2107081.14539022,
                                      1876930.94818628, 1851968.07105024, 1752179.77441134, 1657563.03745279,
                                      1572881.49612103, 1529179.45506106, 1415590.683716, 1424716.37087569,
                                      1409394.9852884, 2471149.2108745, 2577300.47311134, 0., 0., 0., 0., 0.,
                                      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])))

            self.assertTrue(
                np.allclose(unburned_fuel_data,
                            np.array([2.94463669e-04, 7.80797052e-05, 2.30986229e-04, 2.31323498e-03, 2.16117296e-04,
                                      3.08452273e-04, 6.36774075e-04, 1.25373874e-03, 6.28214145e-05, 1.60950331e-04,
                                      1.22376606e-04, 2.88977000e-04, 1.05326938e-04])))


if __name__ == '__main__':
    unittest.main()
