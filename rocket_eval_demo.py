import time
import pickle

import matplotlib.pyplot as plt
import numpy as np

from problems.rocket_propellant_design import rocket_single_eval

if __name__ == '__main__':
    # TODO: Get results from original inputs for consistency with the documentation
    x_test = np.array(
             [[1, 7, 6, 7, 4, 9, 3, 10, 0, 2, 4, 8, 2, 2, 2, 5, 5, 6, 4, 1, 7, 8, 10, 8, 1, 0, 2, 10, 10, 9, 10, 0, 2, 5,
              10, 1, 2, 10, 0, 0, 7, 6, 10, 9, 0, 8, 1, 9, 10, 10, 10, 6, 10, 5, 1, 4, 4, 0, 1, 0, 9, 9, 10, 10, 7, 9,
              10, 2, 9, 0, 5, 2, 7, 6, 2, 0, 4, 0, 0, 6, 10, 9, 2, 10, 5, 8, 10, 9, 10, 10, 7, 4, 6, 0, 1, 4, 4, 0, 0,
              4, 1, 7, 6, 7, 4, 5, 9, 8, 10, 0, 2, 7, 3, 6, 2, 2, 6, 6, 0, 0, 1, 7, 6, 7, 4, 10, 1, 9, 1, 3, 9, 4, 5, 2,
              0, 1, 10, 5, 0, 0, 1, 7, 6, 7, 4, 10, 3, 1, 10, 1, 7, 3, 2, 7, 10, 9, 2, 5, 0, 0, 1, 7, 6, 7, 4, 8, 9, 5,
              9, 6, 4, 7, 0, 5, 1, 1, 1, 2, 0, 0, 1, 7, 6, 7, 4, 3, 5, 4, 10, 4, 3, 9, 6, 10, 3, 2, 1, 9, 0, 0, 1, 7, 6,
              7, 4, 9, 0, 9, 5, 0, 8, 6, 4, 2, 2, 6, 8, 5, 0, 0, 1, 7, 6, 7, 4, 5, 0, 2, 7, 5, 2, 0, 10, 6, 2, 10, 4, 3,
              0, 0, 7, 8, 10, 8, 1, 2, 0, 10, 4, 3, 6, 3, 1, 4, 1, 0, 0, 0, 0, 0, 2, 1, 3, 3, 2, 2, 14, 9, 3, 3, 2, 2,
              2, 1, 2, 6, 3, 2, 2, 2, 0, 3, 6, 9]],
              # [4, 6, 6, 3, 2, 2, 8, 0, 7, 10, 7, 1, 6, 8, 5, 8, 4, 10, 10, 0, 10, 10, 7, 2, 7, 7, 5, 8, 8, 7, 5, 3, 1, 9, 0,
              #  3, 2, 2, 7, 2, 7, 8, 5, 10, 4, 3, 5, 4, 7, 10, 8, 8, 6, 8, 5, 6, 4, 0, 1, 3, 7, 10, 9, 6, 9, 4, 2, 5, 7, 7,
              #  7, 8, 6, 4, 5, 3, 1, 3, 4, 1, 5, 9, 9, 7, 7, 10, 10, 10, 8, 8, 5, 2, 4, 2, 4, 4, 5, 4, 1, 2, 4, 6, 6, 3, 2,
              #  0, 10, 10, 1, 8, 7, 10, 10, 5, 5, 4, 4, 4, 5, 1, 4, 6, 6, 3, 2, 10, 8, 1, 9, 7, 0, 4, 8, 5, 8, 7, 9, 1, 0, 0,
              #  4, 6, 6, 3, 2, 10, 9, 9, 10, 9, 6, 6, 10, 0, 5, 9, 10, 1, 0, 0, 4, 6, 6, 3, 2, 8, 8, 6, 9, 1, 8, 0, 5, 7, 6,
              #  6, 9, 1, 0, 0, 4, 6, 6, 3, 2, 10, 9, 8, 8, 2, 10, 3, 5, 8, 9, 7, 3, 5, 0, 0, 4, 6, 6, 3, 2, 8, 7, 5, 10, 5,
              #  6, 10, 7, 7, 8, 8, 1, 0, 0, 0, 4, 6, 6, 3, 2, 8, 4, 4, 10, 1, 3, 10, 0, 9, 9, 4, 6, 4, 0, 0, 10, 10, 7, 2, 7,
              #  3, 0, 10, 3, 10, 0, 3, 1, 9, 6, 0, 0, 0, 0, 0, 2, 0, 3, 2, 2, 3, 21, 10, 3, 1, 0, 2, 1, 0, 6, 9, 3, 3, 2, 2,
              #  3, 3, 0, 6
              #  ]],
            dtype=np.int32)
    # x_test = np.array([[ 5, 10,  7, 10,  3,  7, 10,  3,  4,  9,  3,  0, 10,  8,  3,  8,  3,
    #                     0,  2,  1,  8,  9,  2,  5, 10,  1,  1,  5,  9,  5,  6, 10,  6,  6,
    #                     7,  4,  1,  1,  0,  0,  8, 10,  0,  9,  7,  7,  5,  7,  6,  4,  4,
    #                     0,  0,  6,  7,  9,  4,  9,  5,  8,  9,  5,  9,  5,  5,  7,  2, 10,
    #                     3,  3,  0,  8,  2,  9,  7,  5,  4,  7,  4,  6,  3, 10,  7,  9, 10,
    #                     7, 10,  6,  8,  3,  6,  5,  2,  7,  7,  6,  5,  9,  0,  5,  5, 10,
    #                     7, 10,  3, 10,  2,  2,  8,  5, 10,  6, 10,  0,  4,  6,  4,  9,  4,
    #                     7,  5, 10,  7, 10,  3,  6,  1,  5,  6,  1,  2,  2,  5,  8,  1,  1,
    #                     0,  9,  0,  0,  5, 10,  7, 10,  3, 10,  8,  2,  6, 10,  9,  1,  8,
    #                     4,  7,  8,  1,  1,  0,  0,  5, 10,  7, 10,  3,  8,  4,  3,  3,  3,
    #                     7,  5,  2,  3,  1,  0,  6,  4,  0,  0,  5, 10,  7, 10,  3,  5,  5,
    #                     4,  9,  8,  6,  2,  7,  9,  2, 10,  8,  0,  0,  0,  5, 10,  7, 10,
    #                     3,  3,  1,  1,  9,  3,  8,  3,  4,  9,  1,  0,  2,  5,  0,  0,  5,
    #                     10,  7, 10,  3,  8,  7,  0,  0, 10,  3,  3,  6,  8,  2,  7,  9,  3,
    #                     0,  0,  8,  9,  2,  5, 10,  3,  6,  5,  6,  2,  4,  1,  6,  8,  0,
    #                     0,  0,  0,  0,  0,  1,  2,  2,  2,  3,  1, 17,  7,  2,  2,  2,  3,
    #                     1,  2, 19,  9,  0,  3,  1,  0,  1,  0, 35,  9]])

    # x_test = np.array(
    #     [[1, 7, 6, 7, 4, 9, 3, 10, 0, 2, 4, 8, 2, 2, 2, 5, 5, 6, 4, 1, 7, 8, 10, 8, 1, 0, 2, 10, 10, 9, 10, 0, 2, 5,
    #       10, 1, 2, 10, 0, 0, 7, 6, 10, 9, 0, 8, 1, 9, 10, 10, 10, 6, 10, 5, 1, 4, 4, 0, 1, 0, 9, 9, 10, 10, 7, 9,
    #       10, 2, 9, 0, 5, 2, 7, 6, 2, 0, 4, 0, 0, 6, 10, 9, 2, 10, 5, 8, 10, 9, 10, 10, 7, 4, 6, 0, 1, 4, 4, 0, 0,
    #       4, 1, 7, 6, 7, 4, 5, 9, 8, 10, 0, 2, 7, 3, 6, 2, 2, 6, 6, 0, 0, 1, 7, 6, 7, 4, 10, 1, 9, 1, 3, 9, 4, 5, 2,
    #       0, 1, 10, 5, 0, 0, 1, 7, 6, 7, 4, 10, 3, 1, 10, 1, 7, 3, 2, 7, 10, 9, 2, 5, 0, 0, 1, 7, 6, 7, 4, 8, 9, 5,
    #       9, 6, 4, 7, 0, 5, 1, 1, 1, 2, 0, 0, 1, 7, 6, 7, 4, 3, 5, 4, 10, 4, 3, 9, 6, 10, 3, 2, 1, 9, 0, 0, 1, 7, 6,
    #       7, 4, 9, 0, 9, 5, 0, 8, 6, 4, 2, 2, 6, 8, 5, 0, 0, 1, 7, 6, 7, 4, 5, 0, 2, 7, 5, 2, 0, 10, 6, 2, 10, 4, 3,
    #       0, 0, 7, 8, 10, 8, 1, 2, 0, 10, 4, 3, 6, 3, 1, 4, 1, 0, 0, 0, 0, 0, 2, 1, 3, 3, 2, 2, 14, 9, 3, 3, 2, 2,
    #       2, 1, 2, 6, 3, 2, 2, 2, 0, 3, 6, 9]], dtype=np.int32)
    # x_test = np.loadtxt('tests/unit/x1').reshape((1, -1))
    # with open('tests/test_pymoo_res_correct.pkl', 'rb') as f:
    #     res_correct = pickle.load(f)
    # x_test = res_correct.pop.get('X')

    target_thrust_profile = np.array([[0.000, 7532.600],
                                      [0.500, 7190.600],
                                      [1.000, 6864.100],
                                      [1.500, 6552.700],
                                      [2.000, 6255.100],
                                      [2.500, 5971.300],
                                      [3.000, 5700.000],
                                      [3.500, 5441.100],
                                      [4.000, 5194.200],
                                      [4.500, 4958.400],
                                      [5.000, 4733.400],
                                      [5.500, 4518.500],
                                      [6.000, 4313.400],
                                      [6.500, 4117.300],
                                      [7.000, 3930.400],
                                      [7.500, 3752.100],
                                      [8.000, 3581.700],
                                      [8.500, 3419.300],
                                      [9.000, 3264.100],
                                      [9.500, 3115.500],
                                      [10.000, 2974.100],
                                      ])

    n_runs = 1
    for i in range(n_runs):
        t0 = time.time()
        simulator_output_arr = rocket_single_eval.evaluate(x_test, target_thrust_profile=target_thrust_profile,
                                                           ray_depth_flag=0, parallel_evaluation=False)
        print(f"Total execution time = {time.time() - t0} seconds")
        for simulator_output in simulator_output_arr:
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

            print(f"Thrust Reward = {thrust_reward}, Simultaeinity Reward = {simultaneity_reward}")
            print(f"Stop Code {simulation_termination_flag}")
            print("Thrust Profile")
            print(thrust_profile_coarse)
            print("Pressure Profile")
            print(pressure_profile_coarse)
            print("Residual fuel thickness for each segment")
            print(unburned_fuel_data)
            print("Ray burn depths for each star segment")
            print(ray_burn_depths)
            print("Segment star burning or not")
            print(seg_star_burn_status)
            print("Time of circularization for each segment")
            print(circularization_time)

            # plt.plot(time_data_fine, thrust_profile_fine)

    # plt.show()
