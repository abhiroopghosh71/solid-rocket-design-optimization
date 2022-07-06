import warnings

import numpy as np

import approcket as apr
import problems.rocket_propellant_design.rocket_single_eval as rocket_eval
from problems.rocket_propellant_design.rocket_util import RocketData


def calc_obj(x_row_indx, x, rocket_problem):
    x_out, propellant_input, layer_start_radius, ray_burn_depth_code_input, pre_cyl_burn_rate, \
        catchup_propellant_input, added_x_model_indices = rocket_problem.parse_vector_to_model_input(x)

    f, g = np.zeros(rocket_problem.n_obj), np.zeros(rocket_problem.n_constr)
    ray_depth_flag = 0  # If 1, then the model assumes actual ray depths are supplied
    if rocket_problem.star_curve_func:
        ray_depth_flag = 1
    # FIXME: Time/thrust data coarse not correct
    simulator_output = rocket_eval.evaluate(x_out.reshape([1, -1]),
                                            mode=rocket_problem.mode,
                                            ray_depth_flag=ray_depth_flag)[0]

    simulator_output['x_model'] = x_out
    simulator_output['propellant_input'] = propellant_input
    simulator_output['layer_start_radius'] = layer_start_radius
    simulator_output['ray_burn_depth_code_input'] = ray_burn_depth_code_input
    simulator_output['pre_cyl_propellant_input'] = pre_cyl_burn_rate
    simulator_output['catchup_propellant_input'] = catchup_propellant_input
    simulator_output['star_param_vector'] = x_out[-(apr.NRAYS+3)*(apr.NCYLSEGS-3):]

    if rocket_problem.calc_rewards:
        # Thrust profile error objective
        thrust_obj = rocket_problem.obj_metric_used[0]

        thrust_abs_error = 0
        thrust_squared_error = 0

        buffer_time = rocket_problem.target_thrust_profile[-1, 0] - rocket_problem.target_thrust_buffer
        for i, thrust_data in enumerate(rocket_problem.target_thrust_profile):
            current_time = thrust_data[0]

            # Don't calculate thrust error beyond the target burn time (if defined)
            # 1e-5 used as tolerance to avoid floating point errors
            if (rocket_problem.target_burn_time is not None) and (current_time > (rocket_problem.target_burn_time
                                                                                  + 1e-5)):
                break

            current_target_thrust = thrust_data[1]

            # KLUGE: thrust_profile_fine filled with garbage so comment out next 2 lines and use thrust_profile_coarse
            #  instead to calculate error. This assumes all time steps are withing 0.5 second intervals
            # thrust_indx = int(np.round(current_time / apr.DeltaTimeFine))
            # thrust_obtained = simulator_output['thrust_profile_fine'][thrust_indx]
            warnings.warn("Thrust error calculation assumes time step of 0.5. So only coarse target thrust profiles "
                          "will work")
            # FIXME: If rocket burns for lets say 9.025 sec, the thrust at 9.025 sec is written to thrust_profile_coarse
            #   Thus thrust error calc would be slightly wrong for the last time point of burn
            thrust_indx = int(np.round(current_time / apr.DeltaTimeCoarse))
            thrust_obtained = simulator_output['thrust_profile_coarse'][thrust_indx]

            # Sometimes towards the end of the burn, the thrust values change abruptly since the optimizer
            # has no incentive to continue burning once it hits the target burn time. This is undesirable
            # behavior since it means the rocket will blow up as soon as the target burn time is reached and
            # burn at one of the segments hits the shell.
            # A buffer period may be defined by the user to take care of this issue. For the buffer period at
            # the end of burn, the thrust match tolerance will be increased 3 times.
            thrust_match_tolerance = rocket_problem.thrust_match_tolerance
            if (rocket_problem.target_thrust_buffer is not None) and (current_time > (buffer_time + 1e-5)):
                thrust_match_tolerance *= 3

            if (abs(current_target_thrust - thrust_obtained) / current_target_thrust) > thrust_match_tolerance:
                thrust_abs_error += np.abs(current_target_thrust - thrust_obtained)
                thrust_squared_error += (current_target_thrust - thrust_obtained) ** 2

        thrust_metric = thrust_squared_error
        if thrust_obj == 'abs_thrust_error':
            thrust_metric = thrust_abs_error
        elif thrust_obj == 'squared_thrust_error':
            thrust_metric = thrust_squared_error
        else:
            warnings.warn('Invalid thrust metric. Defaulting to squared thrust error')

        # Residual fuel objective
        sum_residual_fuel = np.sum(np.abs(simulator_output['unburned_fuel_data']))
        mean_residual_fuel = np.mean(np.abs(simulator_output['unburned_fuel_data']))
        std_residual_fuel = np.std(np.abs(simulator_output['unburned_fuel_data']))

        longest_distance_to_mean = \
            np.max(np.abs(
                simulator_output['unburned_fuel_data'])) - mean_residual_fuel
        range_of_residual_fuel = \
            np.max(np.abs(simulator_output['unburned_fuel_data'])) \
            - np.min(np.abs(simulator_output['unburned_fuel_data']))

        # Sum of residual thicknesses used as default measure of simultaneous burnout
        simultaneity_objective = rocket_problem.obj_metric_used[1]
        simultaneity_metric = sum_residual_fuel
        if simultaneity_objective == 'mean+std':
            simultaneity_metric = mean_residual_fuel + std_residual_fuel
        elif simultaneity_objective == 'mean':
            simultaneity_metric = mean_residual_fuel
        elif simultaneity_objective == 'std':
            simultaneity_metric = std_residual_fuel / mean_residual_fuel
        elif simultaneity_objective == 'worst_dist_to_mean':
            simultaneity_metric = longest_distance_to_mean
        elif simultaneity_objective == 'min_residual_range':
            simultaneity_metric = range_of_residual_fuel
        elif simultaneity_objective == 'worst_residual':
            simultaneity_metric = np.max(np.abs(
                simulator_output['unburned_fuel_data']))
        elif simultaneity_objective == 'sum_star_residual':
            simultaneity_metric = np.sum(
                simulator_output['unburned_fuel_data'][2:apr.NCYLSEGS-1])
        elif simultaneity_objective == 'mssq-sqm':
            simultaneity_metric = np.sum(
                simulator_output['unburned_fuel_data'] ** 2) / apr.NTOTALSEGS \
                                  - (mean_residual_fuel ** 2)
        elif simultaneity_objective != 'sum':
            warnings.warn("Invalid residual metric. Reverting to sum of residuals")

        f = [thrust_metric, simultaneity_metric]
    else:
        # Use objective values obtained from rocket model
        f = [-simulator_output['thrust_reward'], -simulator_output['simultaneity_reward']]

    return x_row_indx, f, g, simulator_output
