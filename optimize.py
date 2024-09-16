import argparse
import logging
import multiprocessing as mp
import os
import sys
import time
import warnings
import pickle

from pymoo.operators.sampling.rnd import IntegerRandomSampling
from pymoo.operators.crossover.pntx import PointCrossover, SinglePointCrossover, TwoPointCrossover
from pymoo.operators.crossover.sbx import SimulatedBinaryCrossover
from pymoo.operators.mutation.pm import PolynomialMutation
from pymoo.operators.repair.rounding import RoundingRepair
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize

import approcket as apr
from problems import problemdef
import utils.record_data as rec

rocket_module_path = os.path.join(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                               'problems',
                                               'rocket_propellant_design'))
target_thrust_profile_path = os.path.join(rocket_module_path,
                                          'thrust_profiles')  # Location of target thrust profile input file
# Problem-specific environment variables
env_vars = {'output_folder': os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output'),
            'logging_step': None,
            'rocket_module_path': rocket_module_path,
            'target_thrust_profile_path': target_thrust_profile_path,
            'propellant_path': rocket_module_path,
            'optimization_history_file_name': 'optimization_history.hdf5'
            }
rec.env_vars_rocket = env_vars


def setup_logging(out_to_console=False, log_file=None):
    """Set up the logging system."""
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    if out_to_console:
        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(logging.DEBUG)
        handler.setFormatter(formatter)
        root.addHandler(handler)

    if log_file is not None:
        fh = logging.FileHandler(log_file)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        root.addHandler(fh)

    logging.getLogger('matplotlib').setLevel(logging.WARNING)


def run_optimization(args):
    """Runs the optimization using NSGA-II according to the user-supplied parameters.

        Args:
            args (argparse.Namespace): Stores the command line arguments supplied by the user.

    """
    global env_vars

    # Create output directory
    if args.save is None:
        # Create an output folder if none is given by the user
        args.save = f'nsga2-{args.target_thrust}-seed{args.seed}-{time.strftime("%Y%m%d-%H%M%S")}'
    env_vars['output_folder'] = os.path.join(env_vars['output_folder'], args.save)

    if not os.path.exists(env_vars['output_folder']):
        os.makedirs(env_vars['output_folder'])
    print('Experiment dir : {}'.format(args.save))

    setup_logging(log_file=os.path.join(env_vars['output_folder'], 'log.txt'))

    # Set up the log file. This will record important events during execution.
    log_format = '%(asctime)s %(message)s'
    logging.basicConfig(stream=sys.stdout, level=logging.INFO,
                        format=log_format, datefmt='%m/%d %I:%M:%S %p')
    fh = logging.FileHandler(os.path.join(env_vars['output_folder'], 'log.txt'))
    fh.setFormatter(logging.Formatter(log_format))
    logging.getLogger().addHandler(fh)

    env_vars['logging_step'] = args.report_freq
    obj_mode = ('squared_thrust_error', 'mean+std')
    args.obj_mode = obj_mode

    # save all arguments
    logging.info("args = %s", args)

    star_param = None
    if args.fixed_star is True:
        # sampling = np.loadtxt('knee_start_pop_baseline')[:, :-24]
        if args.target_thrust == 'baseline':
            # star_param = [1, 3, 3, 0, 2, 2, 21, 10, 2, 2, 2, 2, 1, 2, 0, 10, 2, 2, 2, 2, 3, 1, 4, 9]
            star_param = [2, 0, 2, 3, 3, 2, 28, 10, 3, 2, 2, 2, 1, 3, 10, 8, 3, 2, 2, 2, 2, 3, 6, 9]
        elif args.target_thrust == 'boost_sustain':
            star_param = [1, 3, 2, 2, 2, 0, 17, 8, 3, 3, 2, 2, 3, 1, 0, 2, 3, 2, 2, 2, 2, 2, 0, 2]
        elif args.target_thrust == 'boost_wait_boost':
            star_param = [2, 3, 2, 2, 1, 2, 3, 3, 2, 2, 2, 1, 2, 1, 0, 4, 3, 2, 2, 2, 0, 1, 18, 3]
        elif args.target_thrust == 'bucket':
            star_param = [3, 1, 0, 1, 1, 2, 29, 8, 1, 0, 2, 3, 3, 2, 25, 8, 3, 2, 2, 2, 2, 2, 0, 9]
        elif args.target_thrust == 'constant_thrust':
            star_param = [3, 2, 2, 3, 2, 0, 1, 5, 3, 2, 0, 2, 1, 3, 15, 10, 2, 2, 1, 1, 0, 0, 14, 8]
        elif args.target_thrust == 'hold_regress':
            star_param = [3, 2, 2, 2, 3, 2, 6, 5, 1, 1, 2, 0, 0, 2, 1, 9, 3, 2, 2, 2, 3, 3, 26, 10]
        elif args.target_thrust == 'two_step':
            star_param = [2, 3, 3, 2, 0, 2, 11, 6, 0, 1, 0, 0, 0, 0, 24, 10, 2, 2, 2, 2, 3, 1, 6, 9]
        else:
            warnings.warn("No star param defined even though fixed-star parameter given")
        # elif args.target_thrust == 'two_step':
        #     star_param = [3, 0, 3, 2, 3, 2, 32, 9, 2, 2, 2, 2, 3, 3, 3, 8, 3, 2, 2, 2, 2, 1, 12, 9]

    ncores = int(args.ncores)

    if (args.burn_timestep % apr.DeltaTimeFine) < 0.024999999999999974:
        warnings.warn(f"Burn timestep of {args.burn_timestep} seconds not a multiple of {apr.DeltaTimeFine} seconds")

    if args.interactive:
        # interactive_interface = InteractiveRepair(rocket_problem=problemdef.get_problem('rocket')())
        pass
    else:
        interactive_interface = None
    sampling = IntegerRandomSampling()
    problem = problemdef.get_problem('rocket')(target_thrust_profile_name=args.target_thrust,
                                               burn_timestep=args.burn_timestep,
                                               mode='const_layer_thickness', propellant_in='propellants.txt',
                                               calc_rewards=True, file_io=False,
                                               obj_metric=obj_mode, fixed_star=args.fixed_star,
                                               user_supplied_star_param=star_param,
                                               target_burn_time=args.target_burn_time, star_curve_func=args.star_curve,
                                               env_vars=env_vars, target_thrust_buffer=float(args.target_thrust_buffer),
                                               use_parallelization=args.parallel, ncores=ncores,
                                               logging_freq=args.report_freq,
                                               interactive_interface=interactive_interface,
                                               interact_script=args.interact_script,
                                               user_input_freq=args.user_input_freq,
                                               probability_update_freq=args.probability_update_freq)

    crossover_operator, mutation_operator = None, None

    if args.star_curve:
        if args.crossover == 'two_point':
            crossover_operator = TwoPointCrossover()
        elif args.crossover == 'one_point':
            crossover_operator = SinglePointCrossover
        elif args.crossover == 'sbx':
            crossover_operator = SimulatedBinaryCrossover()
        else:
            warnings.warn("Crossover operator not defined or invalid choice")
        
        mutation_operator = PolynomialMutation(prob=args.mutation_prob, eta=args.mutation_eta)
    else:
        if args.crossover == 'two_point':
            crossover_operator = TwoPointCrossover(repair=RoundingRepair())
        elif args.crossover == 'one_point':
            crossover_operator = SinglePointCrossover(repair=RoundingRepair())
        elif args.crossover == 'sbx':
            crossover_operator = SimulatedBinaryCrossover(repair=RoundingRepair())
        else:
            warnings.warn("Crossover operator not defined or invalid choice")
        
        mutation_operator = PolynomialMutation(prob=args.mutation_prob, eta=args.mutation_eta, repair=RoundingRepair())

    # Write optimization parameters to file
    run_details_file = open(os.path.join(env_vars['output_folder'], 'run_details.txt'), 'w')
    run_details_file.write("# NSGA-II parameters")
    run_details_file.write(f"seed={args.seed}\n")
    run_details_file.write(f"popsize={args.popsize}\n")
    run_details_file.write(f"maxgen={args.ngen}\n")
    run_details_file.write(f"crossover={args.crossover}\n")
    run_details_file.write("mutation=polynomial\n")
    run_details_file.write(f"mutation-p={args.mutation_prob}\n")
    run_details_file.write(f"mutation-eta={args.mutation_eta}\n")
    run_details_file.write(f"popsize={args.popsize}\n")

    # Write model configuration to file
    run_details_file.write("\n# Problem details\n")
    run_details_file.write(f"objectives={args.obj_mode}\n")
    run_details_file.write(f"logging frequency={int(args.report_freq)}\n")
    run_details_file.write(f"thrust profile={args.target_thrust}\n")
    run_details_file.write(f"target burn time={args.target_burn_time}\n")
    run_details_file.write(f"target thrust profile buffer={args.target_thrust_buffer}\n")
    run_details_file.write(f"fixed star={args.fixed_star}\n")
    run_details_file.write(f"parameterized star curve={args.star_curve}\n")
    run_details_file.write(f"total no. of segs={apr.NTOTALSEGS}\n")
    run_details_file.write(f"no. of star segs={apr.NCYLSEGS - 3}\n")
    run_details_file.write(f"no. of cylindrical segs={apr.NCYLSEGS}\n")
    run_details_file.write(f"no. of corner segs={apr.NCORNERSEGS}\n")
    run_details_file.write(f"no. of layers={apr.NLAYERS}\n")
    run_details_file.write(f"no. of rays={apr.NRAYS}\n")
    run_details_file.write(f"no. of ray depths allowed={apr.NRAYDEPTHS}\n")
    run_details_file.write(f"no. of symmetrical circle sectors in star segments={apr.NSYMS}\n")
    run_details_file.write(f"inner radius of cylindrical portion={apr.INITBURNRADIUS}\n")
    run_details_file.write(f"inner radius of deepest allowable star point={apr.INNERSTARPOINT}\n")
    run_details_file.close()

    data_record_function = None
    if args.report_freq > 0:
        data_record_function = rec.record_state
    method = NSGA2(pop_size=args.popsize,
                   sampling=sampling,
                   crossover=crossover_operator,
                   mutation=mutation_operator,
                   elimate_duplicates=True,
                   display=rec.OptimizationDisplay(),
                   callback=data_record_function
                   )

    args_file = open(os.path.join(env_vars['output_folder'], 'args.txt'), 'w')
    args_file.write(str(args))
    args_file.close()
    
    result = minimize(problem,
                      method,
                      verbose=True,
                      termination=('n_gen', args.ngen),
                      seed=args.seed,
                      save_history=False)

    return result


def parse_args(args):
    """Defines and parses the command line arguments that can be supplied by the user.

    Args:
        args (list): Command line arguments supplied by the user.

    """
    # Command line args accepted by the program
    parser = argparse.ArgumentParser(description='DARPA TRADES Solid Rocket Fuel Design Challenge')

    # Optimization parameters
    parser.add_argument('--seed', type=int, default=0, help='Random seed')
    parser.add_argument('--target-thrust', type=str, default='baseline', help='Thrust profile to target')
    parser.add_argument('--ngen', type=int, default=60, help='Maximum number of generations')
    parser.add_argument('--popsize', type=int, default=50, help='Population size')
    parser.add_argument('--report-freq', type=float, default=500, help='Default logging frequency in generations')
    parser.add_argument('--target-burn-time', type=float, default=None,
                        help='Minimum amount of time (in seconds) the rocket should burn')
    parser.add_argument('--burn-timestep', type=float, default=0.5,
                        help='Timesteps at which thrust error is calculated')

    # Parallelization
    parser.add_argument('--parallel', action='store_true', default=False,
                        help='Use parallel evaluation of the population every generation')
    parser.add_argument('--ncores', default=mp.cpu_count() // 3,
                        help='How many cores to use for population members to be evaluated in parallel')

    # Innovization arguments
    parser.add_argument('--repair', action='store_true', default=False, help='Apply custom repair operator')
    parser.add_argument('--shape-only', action='store_true', default=False, help='Repair only shape variables')
    parser.add_argument('--no-group', action='store_true', default=False,
                        help='Do not group propellant variables segment wise')
    parser.add_argument('--interactive', action='store_true', default=False,
                        help='Enable interactive mode. Might interfere with online innovization')
    parser.add_argument('--user-input-freq', type=float, default=100, help='Frequency with which user input is taken.')
    parser.add_argument('--probability-update-freq', type=float, default=1, help='Frequency with which probability of selection of user provided operators is updated.')
    # Following arg is deprecated
    parser.add_argument('--interact-script', default=None,
                        help='Use a script to follow a predetermined set of user actions')

    parser.add_argument('--save', type=str, help='Experiment name')
    parser.add_argument('--fixed-star', action='store_true', default=False, help='Keep star parameters fixed')
    parser.add_argument('--star-curve', action='store_true', default=False,
                        help='Generate star ray depths from a parameterized curve')
    parser.add_argument('--crossover', default='two_point', help='Choose crossover operator')
    parser.add_argument('--mutation-eta', default=3, help='Define mutation parameter eta')
    parser.add_argument('--mutation-prob', default=0.05, help='Define mutation parameter eta')
    parser.add_argument('--target-thrust-buffer', default=0., help='Time (in seconds) to keep as buffer at the end of '
                                                                   'the target thrust profile')

    return parser.parse_args(args)


if __name__ == '__main__':
    # Sample arg string:
    # arg_str = f'--ngen 60 --popsize 50 --report-freq 20 --target-thrust baseline --seed 0 --mutation-eta 3 ' \
    #           f'--mutation-prob 0.05 --crossover two_point --save nsga2-test-{time.strftime("%Y%m%d-%H%M%S")}'
    t0 = time.time()
    cmd_args = parse_args(sys.argv[1:])
    res = run_optimization(cmd_args)

    pickle.dump(res, open('test_pymoo_res_correct_coarse.pkl', 'wb'))

    print(res.F)
    print(f"Total execution time = {time.time() - t0}")
