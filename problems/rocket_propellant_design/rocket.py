import multiprocessing as mp
import os
import warnings
import copy

import numpy as np
from pymoo.core.problem import Problem
from pymoo.core.repair import NoRepair
# from memory_profiler import profile

import approcket as apr
import problems.rocket_propellant_design.rocket_single_eval as rocket_eval
# import star_curve
from problems.rocket_propellant_design.calc_obj import calc_obj as calc_obj
# from problems.rocket_propellant_design.calc_obj_dict import calc_obj_dict
from problems.rocket_propellant_design.rocket_util import RocketData
# from problems.rocket_propellant_design.rocket_repair_interactive import InteractiveRepair


class RocketProblem(Problem):

    def __init__(self, n_obj=2, n_constr=0, target_thrust_profile_name='baseline',
                 min_pressure=1.379e6, type_var=np.int, propellant_in='propellants.txt', resolution=1e-9,
                 calc_rewards=False, mode='const_layer_thickness', residual_metric='sum', use_parallelization=False,
                 ncores=mp.cpu_count() // 3,
                 file_io=False, obj_metric=('squared_thrust_error', 'mean+std'),
                 fixed_star=False, user_supplied_star_param=None, target_burn_time=None, star_curve_func=False,
                 env_vars=None, target_thrust_buffer=0., burn_timestep=0.5, interactive_interface=None, interact_script=None,
                 logging_freq=0, user_input_freq=None, probability_update_freq=None, **kwargs):

        # Rocket problem-specific variables mostly consisting of path to different data files
        if env_vars is None:
            problem_path = os.path.join('problems', 'rocket_propellant_design')
            self.env_vars = {'rocket_module_path': os.path.join('problems', 'rocket_propellant_design'),
                             'target_thrust_profile_path': os.path.join(problem_path, 'thrust_profiles'),
                             'propellant_path': problem_path}
        else:
            self.env_vars = env_vars
        self.logging_freq = logging_freq

        self.thrust_match_tolerance = 0.05  # The tolerance with which thrust match is to be calculated
        self.burn_timestep = burn_timestep
        self.surrogate_model = None

        # Parallelization
        self.use_parallelization = use_parallelization
        self.ncores = ncores

        # Innovization parameters
        self.star_shape_param_avg = np.zeros((apr.NRAYS + 1) * (apr.NCYLSEGS - 3))
        self.star_shape_param_std = np.zeros((apr.NRAYS + 1) * (apr.NCYLSEGS - 3))
        self.star_shape_param_rule_score = np.zeros(apr.NRAYS * (apr.NCYLSEGS - 3))
        self.x_model_avg = None
        self.x_model_std = None
        self.percent_rank_0 = 0
        self.grouped_vars = {}
        self.interactive_interface = None
        if interactive_interface is not None:
            self.interact_script = interact_script
            self.interactive_interface = interactive_interface
            self.user_input_freq = user_input_freq
            self.probability_update_freq = probability_update_freq

        # if use_parallelization is not None:
        #     elementwise_evaluation = True
        #     parallelization = use_parallelization
        # else:
        #     elementwise_evaluation = False
        #     parallelization = None
        # elementwise_evaluation = False
        # parallelization = None

        self.target_thrust_buffer = target_thrust_buffer
        self.obj_metric_used = obj_metric
        self.file_io = file_io
        self.calc_rewards = calc_rewards
        self.residual_metric = residual_metric

        self.target_burn_time = target_burn_time
        self.fixed_star = fixed_star
        self.user_supplied_star_param = user_supplied_star_param

        self.star_curve_func = star_curve_func

        self.INITIAL_BURN_RADIUS = apr.INITBURNRADIUS
        self.CHAMBER_RADIUS = apr.CHAMBERRADIUS
        self.web_thickness = self.CHAMBER_RADIUS - self.INITIAL_BURN_RADIUS
        self.zone_thickness = self.web_thickness / apr.NZONES

        self.target_thrust_profile_name = target_thrust_profile_name
        self.target_thrust_profile_path = os.path.join(self.env_vars['target_thrust_profile_path'],
                                                       self.target_thrust_profile_name + '.txt')
        self.target_thrust_profile = np.loadtxt(self.target_thrust_profile_path, delimiter=',')

        self.min_pressure = min_pressure

        # Propellant ref. burn rates are sorted in ascending order. The slowest burning propellant is at index 0.
        self.ref_burn_rates = np.sort(np.loadtxt(os.path.join(self.env_vars['propellant_path'], propellant_in)))

        self.min_resolution = resolution  # Minimum volume of propellant that can be manufactured

        self.min_layer_thickness = 0.001
        self.max_layer_thickness = self.zone_thickness - (self.min_layer_thickness * (apr.NLAYERSPERZONE - 1))

        # TODO: Add remaining optimization modes
        # 'full' - Optimize on propellant types, layer thicknesses, star burn parameters
        # 'const_layer_burn_area' - Optimize on propellant types, star burn parameters with layer thicknesses set as
        #                           constant values ensuring the surface area for all layers are equal
        # 'const_layer_thickness' - Optimize on propellant types, star burn parameters with all layer thicknesses equal
        # 'discrete_layer_thickness' - Optimize on propellant types, discretized layer thicknesses, star burn parameters
        # 'single_seg' - Optimize on propellant types, layer thicknesses, star burn parameters considering a single
        #                segment only
        self.mode = mode

        # TODO: Handle edge cases (like no. of cyl. segs = 1)
        n_var = 0
        if mode == 'const_layer_burn_area' or mode == 'const_layer_thickness':
            # No. of vars =
            #  no. of segs * no. of zones per seg * no. of layers per zone (ref. burn rate)
            #               + (no. of rays + 3)*(no. of cylindrical segments - 3)

            # Segment 0 (dome)
            n_var = apr.NLAYERS
            grp_var_no = 0
            self.grouped_vars['seg0'] = np.arange(grp_var_no, n_var)
            grp_var_no = n_var

            # Segment beside dome (Segment 5)
            n_var += (apr.NZONES - 1) * apr.NLAYERSPERZONE

            # Remaining cylindrical segments
            n_var += (apr.NCYLSEGS - 2) * apr.NLAYERS
            for cyl_seg_no in range((apr.NCYLSEGS - 2)):
                self.grouped_vars[f'seg{1 + cyl_seg_no}'] = np.arange(grp_var_no, grp_var_no + apr.NLAYERS)
                grp_var_no += apr.NLAYERS
            self.grouped_vars[f'seg{apr.NCYLSEGS - 1}'] = np.arange(grp_var_no, grp_var_no + (apr.NZONES - 1) * apr.NLAYERSPERZONE)
            grp_var_no = n_var

            # Corner segments
            n_var += ((apr.NZONES - 1) * apr.NLAYERSPERZONE * apr.NCORNERSEGS) - (2*apr.NCORNERSEGS)
            for corner_seg_no in range(apr.NCORNERSEGS):
                self.grouped_vars[f'seg{apr.NCYLSEGS + corner_seg_no}'] = np.arange(grp_var_no, grp_var_no + (apr.NZONES - 1) * apr.NLAYERSPERZONE - 2)
                grp_var_no += (apr.NZONES - 1) * apr.NLAYERSPERZONE - 2
            grp_var_no = n_var

            # Nozzle Segment
            n_var += (apr.NZONES - 2) * apr.NLAYERSPERZONE
            self.grouped_vars[f'seg{apr.NCYLSEGS + apr.NCORNERSEGS}'] = np.arange(grp_var_no, n_var)

            self.propellant_param_avg = np.zeros(n_var)
            self.propellant_param_std = np.zeros(n_var)
            self.propellant_param_rule_score = None
            # Propellant type lower limit set as -0.5 to give an equal chance of selecting propellant 0.
            # The upper limit is equal to the (no. of propellants - 0.5)
            # assuming they are numbered starting from 0
            # self.xl = np.array(-0.5 * np.ones(int(n_var)))
            # self.xu = np.array((self.ref_burn_rates.size - 0.5) * np.ones(int(n_var)))

            # updated by Zhichao, this is handled by pymoo framework now
            self.xl = np.zeros(int(n_var))
            self.xu = (self.ref_burn_rates.size - 1) * np.ones(int(n_var))

            if self.star_curve_func:
                # The star ray depths are now generated by a function with a single parameter
                if apr.EVOLVINGSTARSHAPE == 1:
                    n_var += 3 * (apr.NCYLSEGS - 3)

                    # updated by Zhichao, this is handled by pymoo framework now
                    self.xl = np.append(self.xl, np.tile(np.append(-0.2, [0.0, 0.0]),
                                                         apr.NCYLSEGS - 3))
                    self.xu = np.append(self.xu, np.tile(np.append(1.2, [35.0, 10.0]),
                                                         apr.NCYLSEGS - 3))
                else:
                    warnings.warn("Evolving star shape turned off in approcket.c. Optimization results might be wrong.")
            elif self.fixed_star is False:
                if apr.EVOLVINGSTARSHAPE == 1:
                    n_var += (apr.NRAYS + 3) * (apr.NCYLSEGS - 3)

                    # PreCylStarRDotRef (1st var) set from [0-35]. Each star section for each segment is divided
                    # in half.
                    # The rays in the inner half burn PreCylStarRDotRef / 6 + 5
                    # The rays in outer half burn PreCylStarRDotRef % 6 + 5
                    # RDotRefCatchup (2nd var) can be [0-10] but set as [5-10] to allow faster circularization
                    self.xl = np.append(self.xl, np.tile(np.append(0.0 * np.ones(apr.NRAYS + 1), [0.0, 0.0]),
                                                         apr.NCYLSEGS - 3))
                    self.xu = np.append(self.xu, np.tile(np.append(3.0 * np.ones(apr.NRAYS + 1), [35.0, 10.0]),
                                                         apr.NCYLSEGS - 3))
                else:
                    warnings.warn("Evolving star shape turned off in approcket.c. Optimization results might be wrong.")
                # else:
                #     # Warning: Code is buggy and might give unexpected results. Use with caution.
                #     n_var += (apr.NRAYS + 3) * 2
                #     self.xl = np.append(self.xl, np.tile(np.append(-0.5 * np.ones(apr.NRAYS + 1), [4.5, -0.5]), 2))
                #     self.xu = np.append(self.xu, np.tile(np.append(3.5 * np.ones(apr.NRAYS + 1), [10.5, 2.5]), 2))
            else:
                warnings.warn("Fixed star mode enabled. Star parameters will be set to user provided values and will"
                              "not be optimized")
        elif mode == 'discrete_layer_thickness':
            pass
        elif mode == 'single_seg':
            pass
        else:
            # FIXME: Due to model constraints above, the following code is invalid
            if mode == 'full':
                warnings.warn("Invalid optimization mode selected. Defaulting to 'full'")
            # No. of vars =
            #  no. of segs * no. of zones per seg * no. of layers per zone * 2 (ref. burn rate, layer thick.)
            #               + (no. of rays + 3)*(no. of cylindrical segments - 3)
            n_var = apr.NTOTALSEGS * apr.NLAYERS * 2

            # Propellant type lower limit set as -0.5 to give an equal chance of selecting propellant 0.
            # The upper limit is equal to the (no. of propellants - 0.5)
            # assuming they are numbered starting from 0
            self.xl = np.append(-0.5 * np.ones(int(n_var / 2)),
                                self.min_layer_thickness * np.ones(int(n_var / 2)))
            self.xu = np.append((self.ref_burn_rates.size - 0.5) * np.ones(int(n_var / 2)),
                                self.max_layer_thickness * np.ones(int(n_var / 2)))

            if apr.EVOLVINGSTARSHAPE == 1:
                n_var += (apr.NRAYS + 3) * (apr.NCYLSEGS - 3)
                self.xl = np.append(self.xl, np.tile(np.append(-0.5 * np.ones(apr.NRAYS + 1), [-0.5, -0.5]),
                                                     apr.NCYLSEGS - 3))
                self.xu = np.append(self.xu, np.tile(np.append(3.5 * np.ones(apr.NRAYS + 1), [35.5, 2.5]),
                                                     apr.NCYLSEGS - 3))
            else:
                # Warning: Code is buggy and might give unexpected results. Use with caution.
                n_var += (apr.NRAYS + 3) * 2
                self.xl = np.append(self.xl, np.tile(np.append(-0.5 * np.ones(apr.NRAYS + 1), [4.5, -0.5]), 2))
                self.xu = np.append(self.xu, np.tile(np.append(3.5 * np.ones(apr.NRAYS + 1), [10.5, 2.5]), 2))

        super().__init__(n_var=n_var, n_obj=n_obj, n_constr=n_constr, type_var=type_var, xl=self.xl, xu=self.xu,
                         **kwargs)

    # @profile
    def _evaluate(self, x_in, out, *args, **kwargs):
        """Evaluates the objective functions for the rocket design

        Args:
            x (np.array): The NSGA-II population matrix where each row represents one population member
            out (dict): The output dictionary which is to be filled with the objective values (F) and constraint
                        violation (G) for each population member

        Returns:
            out (dict): The output dictionary which is to be filled with the objective values (F) and constraint
                violation (G) for each population member
        """
        # Handle the case where 1D numpy array is passed to the function
        if x_in.ndim == 1:
            x = x_in.reshape([1, -1])
        else:
            x = np.copy(x_in)

        out['repaired_by'] = np.zeros(x.shape[0])
        if self.interactive_interface is not None:
            if kwargs['algorithm'].n_gen % self.user_input_freq == 0:
                print(f"Awaiting user input, gen = {kwargs['algorithm'].n_gen}")
                self.interactive_interface.get_user_feedback(kwargs['algorithm'].n_gen)
            x, repair_indx = self.interactive_interface.do(self, np.copy(x), **kwargs)
            if len(repair_indx) > 0:
                out['repaired_by'][repair_indx] = 1
        if hasattr(kwargs['algorithm'], 'repair') and kwargs['algorithm'].repair is not None and type(kwargs['algorithm'].repair) != NoRepair:
            print("Repair")
            x = kwargs['algorithm'].repair.do(self, np.copy(x), **kwargs)

        # KLUGE: Repair offspring acc to heuristic
        # prop_mean = 8
        # prop_std = 2
        # if self.interactive_interface:
        #     print("interactive")
        #     for i in range(x.shape[0]):
        #         if np.random.rand() <= 0.8:
        #             for seg_no in range(6):
        #                 indx = self.grouped_vars[f'seg{seg_no}'][:8]
        #                 x[i, indx] = np.random.randint(prop_mean - prop_std, prop_mean + prop_std + 1, len(indx))

        if self.use_parallelization:
            print("Parallel")
            pool = mp.Pool(self.ncores)
            result_objects = [pool.apply_async(calc_obj, args=(i, row, self))
                              for i, row in enumerate(x)]
            pool.close()  # Need to close the pool to prevent spawning too many processes
            pool.join()
            # Result_objects is a list of pool.ApplyResult objects
            results = [r.get() for r in result_objects]

            # apply_async() might return results in a different order
            results.sort(key=lambda r: r[0])

            if x_in.ndim == 1:
                out['X'] = x.flatten()
            else:
                out['X'] = np.copy(x)
            out['F'] = np.array([[r[1][k] for k in range(self.n_obj)] for r in results])
            out['G'] = np.array([[r[2][c] for c in range(self.n_constr)] for r in results])
            # out['rocket_data'] = np.array([r[3] for r in results])
        else:
            # print("Sequential")
            f = np.zeros([x.shape[0], self.n_obj])
            g = np.zeros([x.shape[0], self.n_constr])

            simulator_output_pop = [None for _ in range(x.shape[0])]  # To store other rocket design data
            for indx in range(x.shape[0]):
                x_row_indx, f[indx, :], g[indx, :], simulator_output_pop[indx] = calc_obj(indx, x[indx, :], self)

            if x_in.ndim == 1:
                out['X'] = x.flatten()
            else:
                out['X'] = np.copy(x)
            out['G'] = g
            out['F'] = f

            # FIXME: Its possible to induce a mem leak due to the way Pymoo stores extra data in the out array.
            #  Once a data is attached to an offspring, it might survive for multiple generations and may not be
            #  garbage collected ever. This leads to inflation in mem usage.b
            #  Any extra data attached to the out dict in one generation
            #  seems to remain in memory for the subsequent generations
            # a = {}
            # for key in simulator_output_pop[0].keys():
            #     if type(simulator_output_pop[0][key]) == float or type(simulator_output_pop[0][key]) == int:
            #         a[key] = np.zeros(x.shape[0])
            #         continue
            #     arr_shape = simulator_output_pop[0][key].shape
            #     if len(arr_shape) == 1:
            #         a[key] = np.zeros([x.shape[0], arr_shape[0]])
            #     elif len(arr_shape) == 2:
            #         a[key] = np.zeros([x.shape[0], arr_shape[0], arr_shape[1]])
            #     else:
            #         warnings.warn("3D array present in simulator output. Please check.")
            # for pop_i in range(x.shape[0]):
            #     for key in simulator_output_pop[0].keys():
            #         a[key][pop_i] = simulator_output_pop[pop_i][key]

            # out['simulation_termination_flag'] = np.zeros(x.shape[0])
            # out['pop_indx'] = np.arange(0, x_in.shape[0])
            # for key in simulator_output_pop[0].keys():
            #     if type(simulator_output_pop[0][key]) == float or type(simulator_output_pop[0][key]) == int:
            #         out[key] = np.zeros(x.shape[0])
            #         continue
            #     arr_shape = simulator_output_pop[0][key].shape
            #     if len(arr_shape) == 1:
            #         out[key] = np.zeros([x.shape[0], arr_shape[0]])
            #     elif len(arr_shape) == 2:
            #         out[key] = np.zeros([x.shape[0], arr_shape[0], arr_shape[1]])
            #     else:
            #         warnings.warn("3D array present in simulator output. Please check.")
            #
            # for pop_i in range(x.shape[0]):
            #     for key in simulator_output_pop[0].keys():
            #         out[key][pop_i] = simulator_output_pop[pop_i][key]

    def parse_vector_to_model_input(self, x):
        """Converts the pymoo decision variable into the evaluator input vector. It also splits the evaluator input
        vector into individual components like propellant inputs and layer start radii for the non-star portions of
        all segments and ray depth codes, pre-cylindrical propellants and catchup propellants for all star segments.

        Args:
            x (numpy.ndarray): The Numpy array of decision variables

        Returns:
            (tuple): Tuple containing:

                - **x_out** (*numpy.ndarray*): The model input vector obtained from the pymoo decision variable vector.
                - **propellant_input** (*numpy.ndarray*): Propellants in each layer of the non-star portions of all the segments.
                - **layer_start_radius** (*numpy.ndarray*): Start radius of each layer of the non-star portions of all the segments.
                - **ray_burn_depth_input** (*numpy.ndarray*): Defines the ray depths of for the star segments.
                - **pre_cyl_burn_rate** (*numpy.ndarray*): Defines the propellants in the pre-cylindrical sections of the star segments.
                - **catchup_propellant_input** (*numpy.ndarray*): Defines the catchup propellant to be used for circularization once burn in the star segments move into the cylindrical section.
        """
        added_var_indices = []  # To store the location of non-decision vars on the model input vector
        if x.ndim == 1:
            x_out = x.flatten()

        # Convert NSGA-II decision vector to model input vector
        # There are some additional variables which the optimization should not set. They should follow some constraints
        # dictated by the model in order to obtain accurate results. The following sections incorporate those\
        # constraints into the input vector to be sent to the rocket evaluator
        curr_indx = (apr.NCYLSEGS - 1) * apr.NLAYERS
        arr_to_insert = x_out[:apr.NLAYERSPERZONE]  # Zone 0 of Seg 0

        # Equalize Seg 0 Zone 0 propellants to Zone 0 of the adjacent segment
        x_out = np.insert(x_out, curr_indx, arr_to_insert)
        added_var_indices += [(curr_indx + k) for k in range(len(arr_to_insert))]

        # Equalize Seg 0 Zone 0 propellants to Zone 0 of the corner segments
        for i in range(apr.NCORNERSEGS):
            curr_indx = apr.NCYLSEGS * apr.NLAYERS + i * apr.NLAYERS
            # Equalize Zone 0 of corner segs to Seg 0 Zone 0
            x_out = np.insert(x_out, curr_indx, arr_to_insert)
            added_var_indices += [(curr_indx + k) for k in range(len(arr_to_insert))]
            # Insert dummy values for last two values of corner segments since they're set by the model
            x_out = np.insert(x_out, curr_indx + apr.NLAYERS - 2, np.zeros(2))
            added_var_indices += [(curr_indx + apr.NLAYERS - 2 + k) for k in range(len(np.zeros(2)))]

        # Equalize first zone of nozzle segment to Seg 1 Zone 0
        curr_indx = (apr.NCYLSEGS + apr.NCORNERSEGS) * apr.NLAYERS
        arr_to_insert = x_out[apr.NLAYERS:apr.NLAYERS + apr.NLAYERSPERZONE]
        x_out = np.insert(x_out, curr_indx, arr_to_insert)
        added_var_indices += [(curr_indx + k) for k in range(len(arr_to_insert))]

        # Set last zone of nozzle segment to dummy values since its set by the model
        arr_to_insert = np.zeros(apr.NLAYERSPERZONE)
        x_out = np.insert(x_out, apr.NTOTALSEGS * apr.NLAYERS - apr.NLAYERSPERZONE, arr_to_insert)
        added_var_indices += [(apr.NTOTALSEGS * apr.NLAYERS - apr.NLAYERSPERZONE + k) for k in range(len(arr_to_insert))]

        if self.star_curve_func:
            # Extract the star params from the decision vector.
            # x_star_params = x_out[-3 * (apr.NCYLSEGS - 3):]
            #
            # new_star_arr = np.array([])
            # for seg in range(apr.NCYLSEGS - 3):
            #     curve_param_seg = x_star_params[seg * 3]
            #     _, ray_depth_seg = star_curve.get_depths(curve_param_seg, viz=False)
            #     new_star_arr = np.concatenate((new_star_arr, ray_depth_seg[1:],
            #                                    x_star_params[seg * 3 + 1:(seg + 1) * 3]))
            #
            # x_out = x_out[:-3 * (apr.NCYLSEGS - 3)]
            # x_out = np.append(x_out, new_star_arr)
            pass
        elif self.fixed_star is True:
            x_out = np.append(x_out, self.user_supplied_star_param)

        if self.mode == 'full' or self.mode == 'discrete_layer_thickness':
            layer_thickness = x_out[apr.NTOTALSEGS*apr.NLAYERS:apr.NTOTALSEGS*apr.NLAYERS*2]
            layer_start_radius = rocket_eval.get_layer_start_radius_from_thickness(layer_thickness)
            # FIXME: Pointless statement below (??)
            x_out = np.append(x_out[:apr.NTOTALSEGS*apr.NLAYERS], x_out[apr.NTOTALSEGS*apr.NLAYERS*2:])
        else:
            layer_start_radius, layer_thickness = rocket_eval.get_layer_start_radius_and_thickness(self.mode)

        propellant_input, ray_burn_depth_input, pre_cyl_burn_rate, catchup_propellant_input\
            = rocket_eval.split_params(x_out)

        added_var_indices = np.array(added_var_indices)

        return x_out, propellant_input, layer_start_radius, ray_burn_depth_input, pre_cyl_burn_rate, \
            catchup_propellant_input, added_var_indices
