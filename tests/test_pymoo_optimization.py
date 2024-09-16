import importlib
import os
import time
import unittest
import pickle
import multiprocessing as mp

import h5py
import numpy as np
import pandas as pd

import optimize
from problems.rocket_propellant_design import rocket_single_eval
from problems.rocket_propellant_design.rocket_util import parse_vector_to_model_input


class TestOptimization():
    """Tests the pymoo optimization framework for the rocket design problem."""

    def run_optimization(self, parallel=False):
        """Tests the optimization output for a specific set of input parameters."""

        dir_path = os.getcwd()
        # FIXME: Gives different results when optimization run twice or in parallel if thrust_profile_fine used for
        #  fitness calculation
        # Test the results returned by Pymoo (X, F, pop, etc)
        # with open('test_pymoo_res_correct.pkl', 'rb') as f:
        #     res_original = pickle.load(f)
        # x_pop_original = res_original.pop.get('X')
        # with open('test_pymoo_res_correct_bug.pkl', 'rb') as f:
        #     res_correct = pickle.load(f)

        # Uncomment following line and change calc_obj to use thrust_profile_fine for thrust error calculation to test
        # the case when target thrust profiles are not defined in steps of 0.5
        # with open('test_pymoo_res_correct.pkl', 'rb') as f:
        with open(os.path.join(dir_path, 'test_pymoo_res_correct_coarse.pkl'), 'rb') as f:
            res_correct = pickle.load(f)
        x_pop_correct = res_correct.pop.get('X')

        # FIXME: Calling rocket_single_eval multiple times before running the optimization seems to produce
        #  different results in the optimization than if it wasn't called. This does not happen if the entire
        #  module containing the optimization routine is stopped and rerun. Uncomment following lines to
        #  reproduce issue
        # for j in range(3):
        #     x_model_input_correct = []
        #     for pop_indx in range(x_pop_correct.shape[0]):
        #         # x_model_input_correct.append(parse_vector_to_model_input(x_pop_original[pop_indx, :])[0])
        #         x_model_input_correct.append(parse_vector_to_model_input(x_pop_correct[pop_indx, :])[0])
        #     x_model_input_correct = np.array(x_model_input_correct)
        #     simulator_res_arr_correct = rocket_single_eval.evaluate(x_model_input_correct)
        #     print(simulator_res_arr_correct[0]['thrust_reward'],
        #           simulator_res_arr_correct[0]['simultaneity_reward'])

        for i in range(1):
            t0 = time.time()
            arg_str = f'--ngen 60 --popsize 50 --report-freq 20 --target-thrust baseline_1sec_buffer --seed 0 ' \
                      f'--mutation-eta 3 ' \
                      f'--mutation-prob 0.05 --crossover two_point ' \
                      f'--save nsga2-test-{time.strftime("%Y%m%d-%H%M%S")} ' \
                      f'--target-thrust-buffer 0'
            if parallel:
                n_cpu = mp.cpu_count()
                if n_cpu == 1:
                    print("Parallel execution not possible on machine. Exiting.")
                    return
                print(f"Using {max(n_cpu // 4, 2)} cores.")
                arg_str += f' --parallel --ncores {max(n_cpu // 4, 2)}'
            cmd_args = optimize.parse_args(arg_str.split(' '))

            res = optimize.run_optimization(cmd_args)
            print(res.F)
            print(f"Total execution time = {time.time() - t0}")

            rocket_data_pop = res.pop.get('rocket_data')
            self.assertTrue(np.allclose(res.X, res_correct.X))
            self.assertTrue(np.allclose(res.F, res_correct.F))
            self.assertTrue(np.allclose(res.G, res_correct.G))
            self.assertTrue(np.allclose(res.pop.get('X'), x_pop_correct))
            self.assertTrue(np.allclose(res.pop.get('F'), res_correct.pop.get('F')))
            self.assertTrue(np.allclose(res.pop.get('G'), res_correct.pop.get('G')))
            self.assertTrue(np.allclose(res.pop.get('CV'), res_correct.pop.get('CV')))
            self.assertTrue(np.allclose(res.pop.get('rank'), res_correct.pop.get('rank')))

            # for pop_indx in range(x_pop_correct.shape[0]):
            #     self.assertTrue(np.allclose(simulator_res_arr_correct[pop_indx]['thrust_profile_fine'],
            #                                 rocket_data_pop[pop_indx, 0].thrust_profile_fine))

            # Test if the correct data has been written to the output file. Only final generation is checked under the
            # assumption that a correct final data will ensure the previous generation data are also correct
            with h5py.File(os.path.join(optimize.env_vars['output_folder'],
                                        optimize.env_vars['optimization_history_file_name']), 'r') as hf:
                last_gen = hf.attrs['current_gen']
                hf_x = np.array(hf[f'gen{last_gen}']['X'])
                hf_f = np.array(hf[f'gen{last_gen}']['F'])
                hf_rank = np.array(hf[f'gen{last_gen}']['rank'])
                hf_g = np.array(hf[f'gen{last_gen}']['G'])
                hf_cv = np.array(hf[f'gen{last_gen}']['CV'])

                self.assertTrue(np.allclose(hf_x, res_correct.pop.get('X')))
                self.assertTrue(np.allclose(hf_f, res_correct.pop.get('F')))
                self.assertTrue(np.allclose(hf_rank, res_correct.pop.get('rank')))
                self.assertTrue(np.allclose(hf_g, res_correct.pop.get('G')))
                self.assertTrue(np.allclose(hf_cv, res_correct.pop.get('CV')))

            # Test if the correct data has been written to the HDF5 backup file.
            with h5py.File(os.path.join(optimize.env_vars['output_folder'],
                                        optimize.env_vars['optimization_history_file_name'] + '.bak'), 'r') as hf:
                last_gen = hf.attrs['current_gen']
                hf_x = np.array(hf[f'gen{last_gen}']['X'])
                hf_f = np.array(hf[f'gen{last_gen}']['F'])
                hf_rank = np.array(hf[f'gen{last_gen}']['rank'])
                hf_g = np.array(hf[f'gen{last_gen}']['G'])
                hf_cv = np.array(hf[f'gen{last_gen}']['CV'])

                self.assertTrue(np.allclose(hf_x, res_correct.pop.get('X')))
                self.assertTrue(np.allclose(hf_f, res_correct.pop.get('F')))
                self.assertTrue(np.allclose(hf_rank, res_correct.pop.get('rank')))
                self.assertTrue(np.allclose(hf_g, res_correct.pop.get('G')))
                self.assertTrue(np.allclose(hf_cv, res_correct.pop.get('CV')))

    def test_optimization(self):
        """Tests the optimization output for a specific set of input parameters."""
        self.run_optimization(parallel=False)

    def test_optimization_parallel(self):
        """Tests the optimization output for a specific set of input parameters when run in parallel."""
        self.run_optimization(parallel=True)


if __name__ == '__main__':
    unittest.main()
