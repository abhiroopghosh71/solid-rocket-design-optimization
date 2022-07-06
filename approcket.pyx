"""
The following piece of C code defined by Cython the function to be called in approcket.c
"""

import numpy as np
cimport numpy as np
from cython.parallel import prange


cdef extern from "approcket.h":
    # Get rocket geometry configuration parameters from approcket.h
    cdef int _NZONES "NZONES"
    cdef int _NLAYERSPERZONE "NLAYERSPERZONE"
    cdef int _NLAYERS "NLAYERS"

    cdef int _NCYLSEGS "NCYLSEGS"
    cdef int _NCORNERSEGS "NCORNERSEGS"
    cdef int _NOZZLESEG "NOZZLESEG"
    cdef int _NTOTALSEGS "NTOTALSEGS"

    cdef int _NNOZRAYS "NNOZRAYS"
    cdef int _NSEGINTERFACES "NSEGINTERFACES"

    cdef float _CHAMBERRADIUS "CHAMBERRADIUS"
    cdef float _INITBURNRADIUS "INITBURNRADIUS"

    # cdef int _STARGRAIN "STARGRAIN"

    # TODO: ifdef STARGRAIN
    cdef double _INNERSTARPOINT "INNERSTARPOINT"
    cdef int _NSYMS "NSYMS"
    cdef int _NRAYS "NRAYS"
    cdef int _EVOLVINGSTARSHAPE "EVOLVINGSTARSHAPE"

    # TODO: ifdef EVOLVINGSTARSHAPE
    cdef int _NRAYDEPTHS "NRAYDEPTHS"

    cdef int _NumTimePointsCoarse "NumTimePointsCoarse"
    cdef int _NumTimePointsFine "NumTimePointsFine"

    cdef float _DeltaTimeFine "DeltaTime"
    cdef float _DeltaTimeCoarse "DeltaTimeCoarse"

    void objectivefunction(int *reason_stopped, double *l_thrust_reward, double *l_simultaneity_reward,
                           double *l_times, double *l_thrust_profile, double *l_m_pressure_profile,
                           double *l_times_fine, double *l_thrust_profile_fine, double *l_m_pressure_profile_fine,
                           double *l_burn_dist_remaining, int *l_propellant_to_eval,
                           double *l_layer_start_radius, int *l_ray_burn_depth_input,
                           int *l_pre_cyl_rdot_ref_input, int *l_rdot_ref_catchup_input,
                           double l_throat_area_input, double *l_delta_v, double l_max_thrust, int *l_ray_burn_depths,
                           int *l_star_burning_this_seg_now, double *l_time_star_burn_finish,
                           double *l_seg_residual_per_time_step_coarse,
                           double *l_seg_residual_per_time_step_fine, double *l_target_thrusts, int l_ray_depth_flag,
                           double *l_seg_burn_area_per_time_step_coarse,
                           double *l_seg_burn_area_per_time_step_fine,
                           double *l_burn_layer_per_timestep_coarse, double *l_burn_layer_per_timestep_fine)


# Import rocket geometry parameters into Python variables
NZONES = _NZONES
NLAYERSPERZONE = _NLAYERSPERZONE
NLAYERS = _NLAYERS

NCYLSEGS = _NCYLSEGS
NCORNERSEGS = _NCORNERSEGS
NOZZLESEG = _NOZZLESEG
NTOTALSEGS = _NTOTALSEGS

NNOZRAYS = _NNOZRAYS
NSEGINTERFACES = _NSEGINTERFACES

CHAMBERRADIUS = _CHAMBERRADIUS
INITBURNRADIUS = _INITBURNRADIUS

# STARGRAIN = _STARGRAIN

INNERSTARPOINT = _INNERSTARPOINT
NSYMS = _NSYMS
NRAYS = _NRAYS
EVOLVINGSTARSHAPE = _EVOLVINGSTARSHAPE
NRAYDEPTHS = _NRAYDEPTHS

NumTimePointsCoarse = _NumTimePointsCoarse
NumTimePointsFine = _NumTimePointsFine

DeltaTimeFine = _DeltaTimeFine
DeltaTimeCoarse = _DeltaTimeCoarse


def evaluate_design_in_parallel(int [:, :] propellant_memview, double [:] layer_start_radius_memview,
                                int [:, :] ray_depth_code_memview, int [:, :] pre_cyl_memview,
                                int [:, :] rdotref_catchup_memview, double [:, :] target_thrusts_memview,
                                max_thrust, ray_depth_flag):
    cdef Py_ssize_t i

    cpdef int c_reason_stopped
    cpdef double c_thrust_reward
    cpdef double c_simultaneity_reward
    cdef double [::1] c_times = np.zeros(NumTimePointsCoarse, dtype=np.double)
    cdef double [::1] c_thrust_profile = np.zeros(NumTimePointsCoarse, dtype=np.double)
    cdef double [::1] c_m_pressure_profile = np.zeros(NumTimePointsCoarse, dtype=np.double)

    cdef double [::1] c_times_fine = np.zeros(NumTimePointsFine, dtype=np.double)
    cdef double [::1] c_thrust_profile_fine = np.zeros(NumTimePointsFine, dtype=np.double)
    cdef double [::1] c_m_pressure_profile_fine = np.zeros(NumTimePointsFine, dtype=np.double)

    cdef double [::1] c_burn_dist_remaining = np.zeros(NTOTALSEGS, dtype=np.double)
    cdef int [::1] c_ray_burn_depths = np.zeros((NCYLSEGS - 3) * (NRAYS + 1), dtype=np.int32)
    cdef int [::1] c_seg_star_burn_status = np.zeros(NTOTALSEGS, dtype=np.int32)
    cdef double [::1] c_circularization_time = np.zeros(NTOTALSEGS, dtype=np.double)

    cdef double[::1] c_seg_residual_per_time_step_coarse = np.zeros(NTOTALSEGS * NumTimePointsCoarse,
                                                                    dtype=np.double)
    cdef double[::1] c_seg_residual_per_time_step_fine = np.zeros(NTOTALSEGS * NumTimePointsFine,
                                                                  dtype=np.double)

    cdef double[::1] c_seg_burn_area_per_time_step_coarse = np.zeros(NTOTALSEGS * NumTimePointsCoarse,
                                                                     dtype=np.double)
    cdef double[::1] c_seg_burn_area_per_time_step_fine = np.zeros(NTOTALSEGS * NumTimePointsFine,
                                                                            dtype=np.double)
    cdef double[::1] c_burn_layer_per_timestep_coarse = np.zeros(NTOTALSEGS * NumTimePointsCoarse,
                                                                 dtype=np.double)
    cdef double[::1] c_burn_layer_per_timestep_fine = np.zeros(NTOTALSEGS * NumTimePointsFine,
                                                               dtype=np.double)

    cdef double c_l_throat_area_input = .001297  # Currently not used by the model
    cdef double c_l_delta_v = 0.
    cpdef double c_max_thrust = max_thrust

    cdef int c_l_ray_depth_flag = ray_depth_flag

    for i in prange(propellant_memview.shape[0], nogil=True):
        pass
        # FIXME: Cython says objectivefunction requires GIL (How??)
        # objectivefunction(&c_reason_stopped, &c_thrust_reward, &c_simultaneity_reward, &c_times[0], &c_thrust_profile[0],
        #                   &c_m_pressure_profile[0], &c_times_fine[0], &c_thrust_profile_fine[0],
        #                   &c_m_pressure_profile_fine[0], &c_burn_dist_remaining[0], &propellant_memview[0],
        #                   &layer_start_radius_memview[0], &ray_depth_code_memview[0], &pre_cyl_memview[0],
        #                   &rdotref_catchup_memview[0], c_l_throat_area_input, &c_l_delta_v, c_max_thrust, &c_ray_burn_depths[0],
        #                   &c_seg_star_burn_status[0], &c_circularization_time[0], &c_seg_residual_per_time_step_coarse[0],
        #                   &c_seg_residual_per_time_step_fine[0], &target_thrusts_memview[0], c_l_ray_depth_flag,
        #                   &c_seg_burn_area_per_time_step_coarse[0], &c_seg_burn_area_per_time_step_fine[0],
        #                   &c_burn_layer_per_timestep_coarse[0], &c_burn_layer_per_timestep_fine[0])

def evaluate_rocket_design(int [::1] propellant_memview, double [::1] layer_start_radius_memview,
                           int [::1] ray_depth_code_memview, int [::1] pre_cyl_memview,
                           int [::1] rdotref_catchup_memview, double [::1] target_thrusts_memview,
                           max_thrust, ray_depth_flag):
    """
    Evaluates a single population member and gives the thrust profile and other details
    :param
    thrust_out_file: Thrust profile output file generated by the rocket C code
    seg_out_file: Segment data output file generated by the rocket C code
    :return: True indicates output file has been generated successfully
    """

    cdef int c_reason_stopped
    cdef double c_thrust_reward
    cdef double c_simultaneity_reward
    cdef double [::1] c_times = np.zeros(NumTimePointsCoarse, dtype=np.double)
    cdef double [::1] c_thrust_profile = np.zeros(NumTimePointsCoarse, dtype=np.double)
    cdef double [::1] c_m_pressure_profile = np.zeros(NumTimePointsCoarse, dtype=np.double)

    cdef double [::1] c_times_fine = np.zeros(NumTimePointsFine, dtype=np.double)
    cdef double [::1] c_thrust_profile_fine = np.zeros(NumTimePointsFine, dtype=np.double)
    cdef double [::1] c_m_pressure_profile_fine = np.zeros(NumTimePointsFine, dtype=np.double)

    cdef double [::1] c_burn_dist_remaining = np.zeros(NTOTALSEGS, dtype=np.double)
    cdef int [::1] c_ray_burn_depths = np.zeros((NCYLSEGS - 3) * (NRAYS + 1), dtype=np.int32)
    cdef int [::1] c_seg_star_burn_status = np.zeros(NTOTALSEGS, dtype=np.int32)
    cdef double [::1] c_circularization_time = np.zeros(NTOTALSEGS, dtype=np.double)

    cdef double[::1] c_seg_residual_per_time_step_coarse = np.zeros(NTOTALSEGS * NumTimePointsCoarse, dtype=np.double)
    cdef double[::1] c_seg_residual_per_time_step_fine = np.zeros(NTOTALSEGS * NumTimePointsFine, dtype=np.double)

    cdef double[::1] c_seg_burn_area_per_time_step_coarse = np.zeros(NTOTALSEGS * NumTimePointsCoarse, dtype=np.double)
    cdef double[::1] c_seg_burn_area_per_time_step_fine = np.zeros(NTOTALSEGS * NumTimePointsFine,
                                                                            dtype=np.double)
    cdef double[::1] c_burn_layer_per_timestep_coarse = np.zeros(NTOTALSEGS * NumTimePointsCoarse, dtype=np.double)
    cdef double[::1] c_burn_layer_per_timestep_fine = np.zeros(NTOTALSEGS * NumTimePointsFine,
                                                                            dtype=np.double)

    cdef double c_l_throat_area_input = .001297  # Currently not used by the model
    cdef double c_l_delta_v = 0.
    cpdef double c_max_thrust = max_thrust

    cdef int c_l_ray_depth_flag = ray_depth_flag

    objectivefunction(&c_reason_stopped, &c_thrust_reward, &c_simultaneity_reward, &c_times[0], &c_thrust_profile[0],
                      &c_m_pressure_profile[0], &c_times_fine[0], &c_thrust_profile_fine[0],
                      &c_m_pressure_profile_fine[0], &c_burn_dist_remaining[0], &propellant_memview[0],
                      &layer_start_radius_memview[0], &ray_depth_code_memview[0], &pre_cyl_memview[0],
                      &rdotref_catchup_memview[0], c_l_throat_area_input, &c_l_delta_v, c_max_thrust, &c_ray_burn_depths[0],
                      &c_seg_star_burn_status[0], &c_circularization_time[0], &c_seg_residual_per_time_step_coarse[0],
                      &c_seg_residual_per_time_step_fine[0], &target_thrusts_memview[0], c_l_ray_depth_flag,
                      &c_seg_burn_area_per_time_step_coarse[0], &c_seg_burn_area_per_time_step_fine[0],
                      &c_burn_layer_per_timestep_coarse[0], &c_burn_layer_per_timestep_fine[0])

    time_data = np.zeros(NumTimePointsCoarse)
    thrust_profile = np.zeros(NumTimePointsCoarse)
    pressure_profile = np.zeros(NumTimePointsCoarse)

    time_data_fine = np.zeros(NumTimePointsFine)
    thrust_profile_fine = np.zeros(NumTimePointsFine)
    pressure_profile_fine = np.zeros(NumTimePointsFine)

    burn_dist_remaining = np.zeros(NTOTALSEGS)
    ray_burn_depths = np.zeros((NCYLSEGS - 3) * (NRAYS + 1), dtype=np.int32)
    seg_star_burn_status = np.zeros(NTOTALSEGS, dtype=np.int32)
    circularization_time = np.zeros(NTOTALSEGS, dtype=np.double)

    seg_residual_per_time_step_coarse = np.zeros([NTOTALSEGS, NumTimePointsCoarse], dtype=np.double)
    seg_residual_per_time_step_fine = np.zeros([NTOTALSEGS, NumTimePointsFine], dtype=np.double)

    seg_burn_area_per_time_step_coarse = np.zeros([NTOTALSEGS, NumTimePointsCoarse], dtype=np.double)
    seg_burn_area_per_time_step_fine = np.zeros([NTOTALSEGS, NumTimePointsFine], dtype=np.double)

    burn_layer_per_timestep_coarse = np.zeros([NTOTALSEGS, NumTimePointsCoarse], dtype=np.double)
    burn_layer_per_timestep_fine = np.zeros([NTOTALSEGS, NumTimePointsFine], dtype=np.double)

    cdef int i
    cdef int j
    for i in range(NumTimePointsCoarse):
        time_data[i] = c_times[i]
        # print(c_thrust_profile[i])
        thrust_profile[i] = c_thrust_profile[i]
        pressure_profile[i] = c_m_pressure_profile[i]

    for i in range(NumTimePointsFine):
        time_data_fine[i] = c_times_fine[i]
        # print(c_thrust_profile[i])
        thrust_profile_fine[i] = c_thrust_profile_fine[i]
        pressure_profile_fine[i] = c_m_pressure_profile_fine[i]

    for i in range(NTOTALSEGS):
        # KLUGE: To correct nozzle
        if i == (NTOTALSEGS - 1):
            burn_dist_remaining[i] = c_burn_dist_remaining[i] / 100
        else:
            burn_dist_remaining[i] = c_burn_dist_remaining[i]

    for i in range(len(ray_burn_depths)):
        ray_burn_depths[i] = c_ray_burn_depths[i]

    for i in range(NTOTALSEGS):
        seg_star_burn_status[i] = c_seg_star_burn_status[i]
        circularization_time[i] = c_circularization_time[i]
        for j in range(NumTimePointsCoarse):
            seg_residual_per_time_step_coarse[i, j] = c_seg_residual_per_time_step_coarse[i*NumTimePointsCoarse + j]
            seg_burn_area_per_time_step_coarse[i, j] = c_seg_burn_area_per_time_step_coarse[i*NumTimePointsCoarse + j]
            burn_layer_per_timestep_coarse[i, j] = c_burn_layer_per_timestep_coarse[i*NumTimePointsCoarse + j]

        for j in range(NumTimePointsFine):
            seg_residual_per_time_step_fine[i, j] = c_seg_residual_per_time_step_fine[i*NumTimePointsFine + j]
            seg_burn_area_per_time_step_fine[i, j] = c_seg_burn_area_per_time_step_fine[i*NumTimePointsFine + j]
            burn_layer_per_timestep_fine[i, j] = c_burn_layer_per_timestep_fine[i*NumTimePointsFine + j]

    # del propellant_memview, layer_start_radius_memview, ray_depth_code_memview, pre_cyl_memview,\
    #     rdotref_catchup_memview, c_times, c_thrust_profile, c_m_pressure_profile, c_burn_dist_remaining, c_ray_burn_depths,\
    #     c_seg_star_burn_status, c_circularization_time, c_seg_residual_per_time_step_coarse,\
    #     c_seg_residual_per_time_step_fine, c_seg_burn_area_per_time_step_coarse, c_seg_burn_area_per_time_step_fine,\
    #     c_burn_layer_per_timestep_coarse, c_burn_layer_per_timestep_fine

    simulator_output = {
        'thrust_reward': c_thrust_reward,  # Thrust reward returned by the model
        'simultaneity_reward': c_simultaneity_reward,  # Simultaneous burnout reward returned by the model
        'stop_code': c_reason_stopped,  # Simulation termination flag
        'time_data_coarse': time_data,  # Time (in seconds) for which thrusts are recorded during burn
                                        # (0.5 sec intervals)
        'thrust_profile_coarse': thrust_profile,  # Thrust Profile obtained by the input rocket design
                                                  # (0.5 sec intervals)
        'pressure_profile_coarse': pressure_profile,  # Pressure Profile obtained by the input rocket design
                                                      # (0.5 sec intervals)
        'time_data_fine': time_data_fine,  # Time (in seconds) for which thrusts are recorded during burn
                                           # (0.025 sec intervals)
        'thrust_profile_fine': thrust_profile_fine,  # Thrust Profile obtained by the input rocket design
                                                     # (0.025 sec intervals)
        'pressure_profile_fine': pressure_profile_fine,  # Pressure Profile obtained by the input rocket design
                                                         # (0.025 sec intervals)
        'unburned_fuel_data': burn_dist_remaining,  # Thickness of fuel at end of burn
        'ray_burn_depths': ray_burn_depths,  # Actual depths of the star segments
        'seg_star_burn_status': seg_star_burn_status,  # Flags denoting which segment still not circularized at the end
                                                       # of simulation
        'circularization_time': circularization_time,  # Time of circularization of each segment
        # Time-wise burn distance remaining
        'seg_residual_per_time_step_coarse': seg_residual_per_time_step_coarse,  # 0.5 seconds
        'seg_residual_per_time_step_fine': seg_residual_per_time_step_fine,  # 0.025 seconds
        # Time-wise segment burn areas (doesn't include rays)
        'seg_burn_area_per_time_step_coarse': seg_burn_area_per_time_step_coarse,  # 0.5 seconds
        'seg_burn_area_per_time_step_fine': seg_burn_area_per_time_step_fine,  # 0.025 seconds
        # Layer burning at every timestep
        'burn_layer_per_timestep_coarse': burn_layer_per_timestep_coarse,  # 0.5 seconds
        'burn_layer_per_timestep_fine': burn_layer_per_timestep_fine,  # 0.025 seconds
    }

    # gc.collect()
    return simulator_output
