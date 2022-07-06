import numpy as np

import approcket as apr
from problems.rocket_propellant_design import rocket_single_eval


class RocketData:
    def __init__(self):
        # Termination flags:
        # 0 normal termination (burn reached shell)
        # 1 high pressure termination
        # 2 low pressure termination
        # 3 time limit reached (10 seconds)
        # 4 abnormal situation (could not complete the simulation)
        self.simulation_termination_flag = None

        self.x_model = None
        self.target_thrust_profile_name = None

        # Reward values returned by rocket model
        self.thrust_reward = None
        self.simultaneity_reward = None
        self.thrust_obj_metric = None
        self.simultaneity_obj_metric = None
        self.obj_metric_used = None

        self.time_data_coarse = None
        self.thrust_profile_coarse = None
        self.pressure_profile_coarse = None
        self.time_data_fine = None
        self.thrust_profile_fine = None
        self.pressure_profile_fine = None
        self.target_thrust_profile = None
        self.unburned_fuel_data = None  # Segment Number, Unburned Depth
        self.n_segs = 0  # Number of segments
        self.seg_data = None  # Segment no., Layer no., rDotRef, Layer start radius
        self.ray_burn_depths = None
        self.seg_star_burn_status = None
        self.time_star_burn_finish = None
        self.seg_residual_per_time_step_coarse = None
        self.seg_residual_per_time_step_fine = None

        # Rocket model inputs
        self.propellant_input = None
        self.layer_start_radius = None
        self.star_param_vector = None
        self.ray_burn_depth_code_input = None
        self.pre_cyl_propellant_input = None
        self.catchup_propellant_input = None


def parse_vector_to_model_input(x):
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
    x_out = x.flatten()

    # Convert NSGA-II decision vector to model input vector
    # There are some additional variables which the optimization should not set. They should follow some constraints
    # dictated by the model in order to obtain accurate results. The following sections incorporate those\
    # constraints into the input vector to be sent to the rocket evaluator
    curr_indx = (apr.NCYLSEGS - 1) * apr.NLAYERS
    arr_to_insert = x_out[:apr.NLAYERSPERZONE]  # Zone 0 of Seg 0

    # Equalize Seg 0 Zone 0 propellants to Zone 0 of the adjacent segment
    x_out = np.insert(x_out, curr_indx, arr_to_insert)

    # Equalize Seg 0 Zone 0 propellants to Zone 0 of the corner segments
    for i in range(apr.NCORNERSEGS):
        curr_indx = apr.NCYLSEGS * apr.NLAYERS + i * apr.NLAYERS
        # Equalize Zone 0 of corner segs to Seg 0 Zone 0
        x_out = np.insert(x_out, curr_indx, arr_to_insert)
        # Insert dummy values for last two values of corner segments since they're set by the model
        x_out = np.insert(x_out, curr_indx + apr.NLAYERS - 2, np.zeros(2))

    # Equalize first zone of nozzle segment to Seg 1 Zone 0
    curr_indx = (apr.NCYLSEGS + apr.NCORNERSEGS) * apr.NLAYERS
    arr_to_insert = x_out[apr.NLAYERS:apr.NLAYERS + apr.NLAYERSPERZONE]
    x_out = np.insert(x_out, curr_indx, arr_to_insert)

    # Set last zone of nozzle segment to dummy values since its set by the model
    arr_to_insert = np.zeros(apr.NLAYERSPERZONE)
    x_out = np.insert(x_out, apr.NTOTALSEGS * apr.NLAYERS - apr.NLAYERSPERZONE, arr_to_insert)

    layer_start_radius, layer_thickness = rocket_single_eval.get_layer_start_radius_and_thickness()

    propellant_input, ray_burn_depth_input, pre_cyl_burn_rate, catchup_propellant_input \
        = rocket_single_eval.split_params(x_out)

    return x_out, propellant_input, layer_start_radius, ray_burn_depth_input, pre_cyl_burn_rate,\
        catchup_propellant_input
