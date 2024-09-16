import logging
import os
import warnings
from shutil import copyfile

import h5py
import numpy as np
from pymoo.core.repair import NoRepair
from pymoo.util.display.display import Display

import approcket as apr
# from optimize import env_vars
from problems import calc_obj


class OptimizationDisplay(Display):

    def _do(self, problem, evaluator, algorithm):
        super()._do(problem, evaluator, algorithm)

        f_pop = algorithm.pop.get('F')

        self.output.append("thrust_min", np.min(f_pop[:, 0]))
        self.output.append("thrust_mean", np.mean(f_pop[:, 0]))
        self.output.append("thrust_median", np.median(f_pop[:, 0]))
        self.output.append("thrust_max", np.max(f_pop[:, 0]))

        self.output.append("residue_min", np.min(f_pop[:, 1]))
        self.output.append("residue_mean", np.mean(f_pop[:, 0]))
        self.output.append("residue_median", np.median(f_pop[:, 0]))
        self.output.append("residue_max", np.max(f_pop[:, 0]))

        if algorithm.n_gen == 1:
            header_str = ""
            for attr in self.output.attrs:
                header_str += f"{attr[0]}"
                for s in range(attr[2] // 2):
                    header_str += " "
                header_str += "|"
                for s in range(attr[2] // 2):
                    header_str += " "
            logging.info(header_str)

        out_str = ""
        for attr in self.output.attrs:
            out_str += f"{attr[1]} | "
        logging.info(out_str)


env_vars_rocket = None


# A Pymoo callback function to record the current state of the optimization as a fail-safe against abrupt termination
# and for informational purposes
def record_state(algorithm):
    """Records the current progress of the optimization into multiple files to prevent loss of data due to abnormal
    termination/crash of the program.

        Args:
            algorithm (pymoo.model.algorithm.Algorithm): The pymoo algorithm object representing the optimization\
                algorithm being used to perform the optimization

    """
    # env_vars_rocket = algorithm.problem.env_vars

    # Print generation-wise statistics
    x_pop = algorithm.pop.get('X')
    f_pop = algorithm.pop.get('F')
    g_pop = algorithm.pop.get('G')
    cv_pop = algorithm.pop.get('CV')
    rank_pop = algorithm.pop.get('rank')

    interactive_interface = algorithm.problem.interactive_interface
    if interactive_interface is not None:
        rep_operator = algorithm.pop.get('repaired_by')
        # alpha = interactive_interface.alpha
        n_offsprings_survived = np.sum(rep_operator == 1) + np.sum(rep_operator == 0)

        # with open(os.path.join(env_vars_rocket['output_folder'], 'interact_log'), 'a') as f:
        #     f.write(f"{algorithm.n_gen},"
        #             f"{interactive_interface.seg234_l0108_probability},"
        #             f"{interactive_interface.seg234_l0108_range[0]},"
        #             f"{interactive_interface.seg234_l0108_range[1]},"
        #             f"{interactive_interface.seg234_l1820_probability},"
        #             f"{interactive_interface.seg234_l1820_range[0]},"
        #             f"{interactive_interface.seg234_l1820_range[1]},"
        #             f"{np.sum(rep_operator == 0)},"
        #             f"{np.sum(rep_operator == 1)},"
        #             f"{x_pop.shape[0]}\n")
        # with open(os.path.join(env_vars_rocket['output_folder'], 'interact_log_confidence'), 'a') as f:
        #     for i in range(apr.NTOTALSEGS):
        #         for j in range(apr.NLAYERS):
        #             f.write(f"{interactive_interface.confidence[i, j]}")
        #             if j < (apr.NLAYERS - 1):
        #                 f.write(",")
        #         f.write("\n")
        # with open(os.path.join(env_vars_rocket['output_folder'], 'interact_log_prop_ll'), 'a') as f:
        #     for i in range(apr.NTOTALSEGS):
        #         for j in range(apr.NLAYERS):
        #             f.write(f"{interactive_interface.prop_ll[i, j]}")
        #             if j < (apr.NLAYERS - 1):
        #                 f.write(",")
        #         f.write("\n")
        # with open(os.path.join(env_vars_rocket['output_folder'], 'interact_log_prop_ul'), 'a') as f:
        #     for i in range(apr.NTOTALSEGS):
        #         for j in range(apr.NLAYERS):
        #             f.write(f"{interactive_interface.prop_ul[i, j]}")
        #             if j < (apr.NLAYERS - 1):
        #                 f.write(",")
        #         f.write("\n")
        # TODO: replace denominator by number of offsprings as set in pymoo
        if n_offsprings_survived != 0:
            user_input_survival_rate = np.sum(rep_operator == 1) / n_offsprings_survived
        else:
            user_input_survival_rate = 0.5  # Assume 1 offspring survived for base and user input

        if algorithm.n_gen % algorithm.problem.probability_update_freq == 0:
            print(f"Updating probability, gen = {algorithm.n_gen}")
            interactive_interface.update_probability(user_input_survival_rate=user_input_survival_rate)
        # interactive_interface.seg234_l0108_probability = np.max(0.1, ((alpha * user_input_survival_rate)
        #                                                               + ((1 - alpha)
        #                                                                  * interactive_interface.seg234_l0108_probability)))
        # interactive_interface.seg234_l1820_probability = np.max(0.1, ((alpha * user_input_survival_rate)
        #                                                         + ((1 - alpha)
        #                                                            * interactive_interface.seg234_l0108_probability)))
        # else:
        #     # If no offsprings survived from base or user input, not enough information to update probabilities
        #     # Get probability of user input closer to mean
        #     interactive_interface.seg234_l0108_probability = 0.5 + (interactive_interface.seg234_l0108_probability
        #                                                             - 0.5) / 2
        #     interactive_interface.seg234_l1820_probability = 0.5 + (interactive_interface.seg234_l1820_probability
        #                                                             - 0.5) / 2
        algorithm.pop.set('repaired_by', -999*np.ones(x_pop.shape[0]))

    x_pop_rank0 = x_pop[rank_pop == 0]
    algorithm.problem.percent_rank_0 = x_pop_rank0.shape[0] / x_pop.shape[0]
    if hasattr(algorithm, 'repair') and algorithm.repair is not None and type(algorithm.repair) != NoRepair:
        algorithm.repair.learn_rules(algorithm.problem, x_pop_rank0)

    target_thrust_profile_name = algorithm.problem.target_thrust_profile_name
    if algorithm.n_gen % env_vars_rocket['logging_step'] == 0 or algorithm.n_gen == 1 \
            or algorithm.n_gen == algorithm.termination.n_max_gen:

        simulator_output_pop = [None for _ in range(x_pop.shape[0])]  # To store other rocket design data
        for indx in range(x_pop.shape[0]):
            _, _, _, simulator_output_pop[indx] = calc_obj(indx, x_pop[indx, :], algorithm.problem)

        out = {}
        for key in simulator_output_pop[0].keys():
            if type(simulator_output_pop[0][key]) == float or type(simulator_output_pop[0][key]) == int:
                out[key] = np.zeros(x_pop.shape[0])
                continue
            arr_shape = simulator_output_pop[0][key].shape
            if len(arr_shape) == 1:
                out[key] = np.zeros([x_pop.shape[0], arr_shape[0]])
            elif len(arr_shape) == 2:
                out[key] = np.zeros([x_pop.shape[0], arr_shape[0], arr_shape[1]])
            else:
                warnings.warn("3D array present in simulator output. Please check.")

        for pop_i in range(x_pop.shape[0]):
            for key in simulator_output_pop[0].keys():
                out[key][pop_i] = simulator_output_pop[pop_i][key]

        # Get rocket stats
        thrust_reward = out['thrust_reward']
        simultaneity_reward = out['simultaneity_reward']
        simulation_termination_flag = out['stop_code']
        time_data_coarse = out['time_data_coarse']
        thrust_profile_coarse = out['thrust_profile_coarse']
        pressure_profile_coarse = out['pressure_profile_coarse']
        time_data_fine = out['time_data_fine']
        thrust_profile_fine = out['thrust_profile_fine']
        pressure_profile_fine = out['pressure_profile_fine']
        unburned_fuel_data = out['unburned_fuel_data']
        ray_burn_depths = out['ray_burn_depths']
        seg_star_burn_status = out['seg_star_burn_status']
        circularization_time = out['circularization_time']
        seg_residual_per_time_step_coarse = out['seg_residual_per_time_step_coarse']
        seg_residual_per_time_step_fine = out['seg_residual_per_time_step_fine']
        seg_burn_area_per_time_step_coarse = out['seg_burn_area_per_time_step_coarse']
        seg_burn_area_per_time_step_fine = out['seg_burn_area_per_time_step_fine']
        burn_layer_per_timestep_coarse = out['burn_layer_per_timestep_coarse']
        burn_layer_per_timestep_fine = out['burn_layer_per_timestep_fine']
        x_model = out['x_model']
        propellant_input = out['propellant_input']
        layer_start_radius = out['layer_start_radius']
        ray_burn_depth_code_input = out['ray_burn_depth_code_input']
        pre_cyl_propellant_input = out['pre_cyl_propellant_input']
        catchup_propellant_input = out['catchup_propellant_input']
        star_param_vector = out['star_param_vector']

        optim_history_hdf_out = os.path.join(env_vars_rocket['output_folder'],
                                             env_vars_rocket['optimization_history_file_name'])

        with h5py.File(optim_history_hdf_out, 'a') as hf:
            # Learn rules even if online innovization not enabled for user review
            # if not (hasattr(algorithm, 'repair') and algorithm.repair is not None):
            #     ParameterlessInequalityRepair().learn_rules(
            #         algorithm.problem, x_pop, burn_layer_per_timestep_coarse=out['burn_layer_per_timestep_coarse'],
            #         seg_star_burn_status=out['seg_star_burn_status'])
            # Store data generation-wise
            hf.attrs['current_gen'] = algorithm.n_gen
            if 'generations' not in hf.attrs:
                hf.attrs['generations'] = []
            hf.attrs['generations'] = np.append(hf.attrs['generations'], algorithm.n_gen)

            g1 = hf.create_group(f'gen{algorithm.n_gen}')

            # Basic population data
            g1.create_dataset('X', data=x_pop)
            g1.create_dataset('F', data=f_pop)
            g1.create_dataset('rank', data=rank_pop)
            g1.create_dataset('G', data=g_pop)
            g1.create_dataset('CV', data=cv_pop)

            # Simulation results for every population member
            # TODO: Add stop reason
            g1.create_dataset('simulation_termination_flag', data=simulation_termination_flag)

            g1.create_dataset('thrust_reward', data=thrust_reward)
            g1.create_dataset('simultaneity_reward', data=simultaneity_reward)
            g1.create_dataset('x_model', data=x_model)

            g1.create_dataset('time_data_coarse', data=time_data_coarse)
            g1.create_dataset('thrust_profile_coarse', data=thrust_profile_coarse)
            g1.create_dataset('pressure_profile_coarse', data=pressure_profile_coarse)

            g1.create_dataset('time_data_fine', data=time_data_fine)
            g1.create_dataset('thrust_profile_fine', data=thrust_profile_fine)
            g1.create_dataset('pressure_profile_fine', data=pressure_profile_fine)

            g1.create_dataset('unburned_fuel_data', data=unburned_fuel_data)
            g1.create_dataset('ray_burn_depths', data=ray_burn_depths.astype(int))
            g1.create_dataset('seg_star_burn_status', data=seg_star_burn_status)
            g1.create_dataset('circularization_time', data=circularization_time)
            g1.create_dataset('seg_residual_per_time_step_coarse', data=seg_residual_per_time_step_coarse)
            g1.create_dataset('seg_residual_per_time_step_fine', data=seg_residual_per_time_step_fine)
            g1.create_dataset('seg_burn_area_per_time_step_coarse', data=seg_burn_area_per_time_step_coarse)
            g1.create_dataset('seg_burn_area_per_time_step_fine', data=seg_burn_area_per_time_step_fine)
            g1.create_dataset('burn_layer_per_timestep_coarse', data=burn_layer_per_timestep_coarse)
            g1.create_dataset('burn_layer_per_timestep_fine', data=burn_layer_per_timestep_fine)

            g1.create_dataset('propellant_input', data=propellant_input.astype(int))
            g1.create_dataset('layer_start_radius', data=layer_start_radius)
            g1.create_dataset('ray_burn_depth_code_input', data=ray_burn_depth_code_input.astype(int))
            g1.create_dataset('pre_cyl_propellant_input', data=pre_cyl_propellant_input.astype(int))
            g1.create_dataset('catchup_propellant_input', data=catchup_propellant_input.astype(int))
            g1.create_dataset('star_param_vector', data=star_param_vector.astype(int))

            # General run details
            if 'n_cylsegs' not in hf.keys():
                hf.create_dataset('n_cylsegs', data=apr.NCYLSEGS)
            if 'n_cornersegs' not in hf.keys():
                hf.create_dataset('n_cornersegs', data=apr.NCORNERSEGS)
            if 'n_totalsegs' not in hf.keys():
                hf.create_dataset('n_totalsegs', data=apr.NTOTALSEGS)
            if 'n_rays' not in hf.keys():
                hf.create_dataset('n_rays', data=apr.NRAYS)
            if 'n_syms' not in hf.keys():
                hf.create_dataset('n_syms', data=apr.NSYMS)

            if 'target_thrust_profile_name' not in hf.keys():
                hf.create_dataset('target_thrust_profile_name', data=algorithm.problem.target_thrust_profile_name)
            if 'target_thrust_profile' not in hf.keys():
                hf.create_dataset('target_thrust_profile', data=algorithm.problem.target_thrust_profile)
            if 'target_thrust_profile_path' not in hf.keys():
                hf.create_dataset('target_thrust_profile_path', data=algorithm.problem.target_thrust_profile_path)
            if 'target_thrust_buffer' not in hf.keys():
                hf.create_dataset('target_thrust_buffer', data=algorithm.problem.target_thrust_buffer)

            hf.attrs['obj_metric_used'] = algorithm.problem.obj_metric_used

            # Innovization data
            if hasattr(algorithm, 'repair') and algorithm.repair is not None:
                g1.create_dataset('repair', data=True)
                # g1.create_dataset('propellant_param_rule_score', data=algorithm.problem.propellant_param_rule_score)
            else:
                g1.create_dataset('repair', data=False)
            # Save "average" of design variables for good solutions
            g1.create_dataset('star_shape_param_avg', data=algorithm.problem.star_shape_param_avg)
            g1.create_dataset('star_shape_param_rule_score', data=algorithm.problem.star_shape_param_rule_score)
            g1.create_dataset('propellant_param_avg', data=algorithm.problem.propellant_param_avg)
            g1.create_dataset('propellant_param_std', data=algorithm.problem.propellant_param_std)
            # Propellant
            # Currently grouped seg-wise. Use ast to recover original dict
            hf.attrs['grouped_vars'] = str(algorithm.problem.grouped_vars)

        # Backup current HDF file as a safeguard against data corruption
        if os.path.exists(optim_history_hdf_out):
            copyfile(optim_history_hdf_out, optim_history_hdf_out + ".bak")
