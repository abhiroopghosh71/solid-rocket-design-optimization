import math
import warnings
import importlib
import copy

import deprecation
import numpy as np
# from memory_profiler import profile

import approcket as apr


def get_layer_start_radius_and_thickness(mode='const_layer_thickness'):
    """Assign layer thicknesses based on a predefined scheme

    Args:
        mode (str): Scheme to assign layer thicknesses

    Returns:
        layer_start_radius (numpy.array): The radius at which each layer starts (in metres)
        layer_thickness (numpy.array): Thickness of each layer (in metres)
    """
    web_thickness = apr.CHAMBERRADIUS - apr.INITBURNRADIUS
    zone_thickness = web_thickness / apr.NZONES
    if mode == 'const_layer_burn_area':
        # Iterative method to calculate layer thicknesses for keeping approximately constant amount of fuel in each
        # layer
        layer_thickness_per_seg = np.zeros(apr.NLAYERS)
        layer_start_radius_per_seg = np.zeros(apr.NLAYERS)
        for zone in range(apr.NZONES):
            zone_start_indx = zone * apr.NLAYERSPERZONE
            layer_start_radius_per_seg[zone_start_indx] = apr.INITBURNRADIUS + zone * zone_thickness

            next_layer_start_radius = layer_start_radius_per_seg[zone_start_indx] + zone_thickness
            c = (next_layer_start_radius ** 2 - layer_start_radius_per_seg[
                zone_start_indx] ** 2) / apr.NLAYERSPERZONE
            for i in range(zone_start_indx + 1, zone_start_indx + apr.NLAYERSPERZONE):
                layer_start_radius_per_seg[i] = math.sqrt(c + layer_start_radius_per_seg[i - 1] ** 2)
                layer_thickness_per_seg[i - 1] = layer_start_radius_per_seg[i] - layer_start_radius_per_seg[i - 1]
                if i == zone_start_indx + apr.NLAYERSPERZONE - 1:
                    layer_thickness_per_seg[i] = next_layer_start_radius - layer_start_radius_per_seg[i]

        layer_thickness = np.tile(layer_thickness_per_seg, apr.NTOTALSEGS)
        layer_start_radius = np.tile(layer_start_radius_per_seg, apr.NTOTALSEGS)

        return layer_start_radius, layer_thickness
    elif mode == 'const_layer_thickness':
        # For all layer thicknesses equal
        layer_thickness = (web_thickness / apr.NLAYERS) * np.ones(apr.NTOTALSEGS * apr.NLAYERS)
        layer_start_radius = get_layer_start_radius_from_thickness(layer_thickness)

        return layer_start_radius, layer_thickness
    else:
        warnings.warn(f"Invalid layer thickness mode '{mode}'")
        return None


def get_layer_start_radius_from_thickness(layer_thickness):
    """Given layer thicknesses, calculate the layer start radii

    Args:
        layer_thickness (numpy.array): Thickness of each layer (in metres)

    Returns:
        layer_start_radius (numpy.array): The radius at which each layer starts (in metres)
    """
    layer_start_radius = np.zeros(apr.NTOTALSEGS * apr.NLAYERS)
    temp_indx = 0
    for count1 in range(apr.NTOTALSEGS):
        layer_start_radius[temp_indx] = apr.INITBURNRADIUS
        for count2 in range(1, apr.NLAYERS):
            layer_start_radius[temp_indx + count2] = \
                layer_start_radius[temp_indx + count2 - 1] + layer_thickness[temp_indx + count2 - 1]
        temp_indx += apr.NLAYERS

    return layer_start_radius


def split_params(x_in):
    """Split model input parameter vector into propellant inputs, ray depths, pre-cylindrical section propellants, and
    catchup propellant.

    Args:
        x (numpy.array): The model input parameter vector

    Returns:
        propellant_input (numpy.array): Propellant types for the cylindrical sections in all the segments
        ray_burn_depth_input (numpy.array): Ray depth codes for defining the geometry of the star
        pre_cyl_burn_rate (numpy.array): Propellant distribution in the non-cylindrical sections for each star segment
        catchup_propellant_input (numpy.array): Propellant used for circularization of burn surface once it comes out
            of the star portion.
    """
    if x_in.ndim == 1:
        x = x_in.reshape(1, -1)
    else:
        x = x_in
    propellant_uindx = apr.NTOTALSEGS * apr.NLAYERS  # The index of x till which propellant information is present
    propellant_input = x[:, :propellant_uindx]  # Propellant input is arranged segment-wise

    ray_burn_depth_input = np.empty((0,))
    pre_cyl_burn_rate = np.empty((0,))
    catchup_propellant_input = np.empty((0,))

    star_param_start_indx = propellant_uindx
    for count in range(apr.NCYLSEGS - 3):
        star_arr = x[:, star_param_start_indx:star_param_start_indx + apr.NRAYS + 3]

        ray_burn_depth_input = np.append(ray_burn_depth_input, star_arr[:, :apr.NRAYS + 1])
        pre_cyl_burn_rate = np.append(pre_cyl_burn_rate, star_arr[:, apr.NRAYS + 1])
        catchup_propellant_input = np.append(catchup_propellant_input, star_arr[:, apr.NRAYS + 2])

        star_param_start_indx += apr.NRAYS + 3

    return propellant_input, ray_burn_depth_input, pre_cyl_burn_rate, catchup_propellant_input


@deprecation.deprecated(  # deprecated_in="1.01", removed_in="1.1", current_version=__version__,
                        details="Use data from the simulation output tuple returned by evaluate()")
def get_burn_layer_per_timestep(seg_residual_per_time_step):
    """Returns which layer the burn is at each timestep for every segment."""
    burn_layer_per_timestep = np.zeros(seg_residual_per_time_step.shape)

    layer_thickness = (apr.CHAMBERRADIUS - apr.INITBURNRADIUS) / apr.NLAYERS * np.ones(apr.NTOTALSEGS * apr.NLAYERS)
    layer_start_radius = get_layer_start_radius_from_thickness(layer_thickness)

    # KLUGE (corner): Corner seg propellant thickness to shell not same as that of cylsegs. Using residual values to guess.
    # For the following vars, only segments 6 to 11 have valid values
    corner_layer_thickness = seg_residual_per_time_step[6:12, 0] / apr.NLAYERS
    corner_chamber_radius = np.zeros(apr.NCORNERSEGS)
    for seg_no_relative, thickness in enumerate(corner_layer_thickness):
        corner_layer_start_radius\
            = get_layer_start_radius_from_thickness(thickness * np.ones(apr.NTOTALSEGS * apr.NLAYERS))
        layer_start_radius[(6 + seg_no_relative) * apr.NLAYERS:(6 + seg_no_relative + 1) * apr.NLAYERS]\
            = corner_layer_start_radius[(6 + seg_no_relative) * apr.NLAYERS:(6 + seg_no_relative + 1) * apr.NLAYERS]
        corner_chamber_radius[seg_no_relative]\
            = layer_start_radius[(6 + seg_no_relative + 1) * apr.NLAYERS - 1] + corner_layer_thickness[seg_no_relative]

    burn_hit_shell = False
    burn_hit_shell_indx = None
    for t in range(seg_residual_per_time_step.shape[1]):
        if burn_hit_shell:
            burn_hit_shell_indx = t
            break
        # Current residual fuel for all 13 segments
        current_residue = seg_residual_per_time_step[:, t]

        for seg_no in range(apr.NTOTALSEGS):
            if current_residue[seg_no] < 1e-5:
                burn_hit_shell = True
            for layer_no in range(apr.NLAYERS):
                if layer_no == (apr.NLAYERS - 1):
                    # Checks if burn is in last layer
                    # if layer_start_radius[seg_no * layer_no] <= (apr.CHAMBERRADIUS - current_residue[seg_no]):
                    burn_layer_per_timestep[seg_no, t] = layer_no + 1
                    break

                # KLUGE (corner): For corner segs
                if 6 <= seg_no <= 11:
                    if (layer_start_radius[seg_no*apr.NLAYERS + layer_no]
                            <= np.around((corner_chamber_radius[seg_no - 6] - current_residue[seg_no]), decimals=5)
                            < layer_start_radius[seg_no*apr.NLAYERS + layer_no + 1]):
                        burn_layer_per_timestep[seg_no, t] = layer_no + 1
                        break

                if (layer_start_radius[seg_no*apr.NLAYERS + layer_no]
                        <= (apr.CHAMBERRADIUS - current_residue[seg_no])
                        < layer_start_radius[seg_no*apr.NLAYERS + layer_no + 1]):
                    burn_layer_per_timestep[seg_no, t] = layer_no + 1
                    break

    return burn_layer_per_timestep[:, :burn_hit_shell_indx]


# @profile
def evaluate(x_mat, layer_thickness=None, mode='const_layer_thickness', target_thrust_profile=None, ray_depth_flag=0,
             parallel_evaluation=False):
    """Takes an input matrix with the propellant types, layer radii (optional), and the star segment configurations of
    each design.

        Args:
            x_mat (numpy.ndarray): The input to be supplied to the rocket model.
            layer_thickness (numpy.ndarray): An array defining the layer thicknesses in the cylindrical portion of the
                solid rocket fuel. Defaults to None.
            mode (str): Defines the scheme through which layer thickness is to be set. Used only if layer_thickness
                parameter is None. Defaults to 'const_layer_thickness'.
            target_thrust_profile (numpy.ndarray): The target thrust profile to be matched. Only used for calculating
                the thrust and simultaneity rewards by the rocket model. Defaults to None
            ray_depth_flag (int): If 1, then the model assumes the actual ray depths are supplied instead of ray depth
                codes. Defaults to 0.
            env_vars (dict): Problem-specific variable dictionary. Mostly defines the paths consisting of the data
                files used by the program. Defaults to None.

        Returns:
            thrust_reward (float): A metric denoting how well the target thrust profile was matched.
            simultaneity_reward (float): A metric denoting how well the simultaneous burnout objective was achieved.
            stop_code (int): A code denoting the reason the simulation stopped.
            time_data_coarse (numpy.ndarray): Time points at which thrust is calculated. Has 0.5 second intervals.
            thrust_values_coarse (numpy.ndarray): Thrust values for the time points defined in time_data_coarse.
            pressure_values_coarse (numpy.ndarray): Pressure values for the time points defined in time_data_coarse.
            time_data_fine (numpy.ndarray): Time points at which thrust is calculated. Has 0.025 second intervals.
            thrust_values_fine (numpy.ndarray): Thrust values for the time points defined in time_data_fine.
            pressure_values_fine (numpy.ndarray): Pressure values for the time points defined in time_data_fine.
            residual_fuel (numpy.ndarray): Segment-wise residual fuel at the end of burn.
            ray_burn_depths (numpy.ndarray): The depth of the propellant in different portions of the star section.
            seg_star_burn_status (numpy.ndarray): Denotes which star segments have not circularized by marking as 1.
            time_star_burn_finish (numpy.ndarray): Time at which the star segments circularized.
            seg_residual_per_time_step_coarse (numpy.ndarray): Residual fuel for the time points defined in time_data_coarse.
            seg_residual_per_time_step_fine (numpy.ndarray): Residual fuel for the time points defined in time_data_fine.
    """

    if target_thrust_profile is None:
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

    if not layer_thickness:
        # If user does not supply layer thickness, use a predefined scheme to assign layer thicknesses
        layer_start_radius, layer_thickness = get_layer_start_radius_and_thickness(mode)
    else:
        layer_start_radius = get_layer_start_radius_from_thickness(layer_thickness)

    simulator_output_arr = []
    if parallel_evaluation is False:
        for i in range(x_mat.shape[0]):
            # importlib.reload(apr)
            propellant_input, ray_burn_depth_input, pre_cyl_burn_rate, catchup_propellant_input\
                = split_params(x_mat[i].reshape([1, -1]))

            # Call the Cython function to evaluate the current rocket model
            simulator_output = apr.evaluate_rocket_design(propellant_input.flatten().astype(np.int32),
                                                          layer_start_radius,
                                                          ray_burn_depth_input.flatten().astype(np.int32),
                                                          pre_cyl_burn_rate.flatten().astype(np.int32),
                                                          catchup_propellant_input.flatten().astype(np.int32),
                                                          np.ascontiguousarray(target_thrust_profile[:, 1]),
                                                          7532.6, ray_depth_flag)
            # simulator_output_arr.append(copy.deepcopy(simulator_output))
            simulator_output_arr.append(simulator_output)

    return simulator_output_arr
