import subprocess
import sys
import os
import time

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))
sys.path.append(parent_dir)

import numpy as np

from visual import utils
from ics import generator_ics


class TestEx3extimes:
    def generate_ics(self, exp_min=1, exp_max=5, num=25):
        self.N_particles = np.logspace(exp_min, exp_max, num=num, dtype=int)
        for i in self.N_particles:
            ics = generator_ics.ic(
                f'test_{i}part', theta=1, dt=1e-2, dt_out=10, t_end=10
            )
            ics.uni_sphere(int(i / 2), 10, i, 0, (-10, 0, 0))
            ics.uni_sphere(int(i / 2), 10, i, 0, (10, 0, 0))
            ics.save_file()

    def setUp(self):
        self.input_files = np.empty(len(self.N_particles), dtype=object)
        self.output_files = np.empty(len(self.N_particles), dtype=object)
        for i, N in enumerate(self.N_particles):
            self.input_files[i] = f'ics/ic_test_{N}part.txt'
            self.output_files[i] = f'output/test_{N}part.dat'
        self.all_good = True

    def test_times_many_particles(self):
        for ic_file in self.input_files:
            try:
                print('Running: ', ic_file)
                result = subprocess.run(
                    ['mpirun',  '-np', str(num_cores), './ex3', ic_file],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                )
                print('Done!')
            except subprocess.CalledProcessError as e:
                print(f"ex3 failed with error: {e.stderr}")

        _, _, _, _, exec_time = utils.read_data(
            self.output_files[len(self.output_files) // 2]
        )
        ref_exec_time = float(exec_time)
        # if good: t(N) = A * Nlog(N) --> A = t(N_ref) / (N_ref*log(N_ref))
        A = ref_exec_time / (
            self.N_particles[len(self.output_files) // 2]
            * np.log10(self.N_particles[len(self.output_files) // 2])
        )

        for index, output_file in enumerate(self.output_files):
            _, _, _, _, exec_time = utils.read_data(output_file)
            print(f"Execution time for {self.N_particles[index]} particles: {exec_time}")
            exec_max_time = (
                100 * A * self.N_particles[index] * np.log10(self.N_particles[index])
            )
            if exec_max_time == 0.0:
                exec_max_time = 10

            if float(exec_time) > exec_max_time:
                print(
                    f"[WARNING!]Maybe something wrong, too many time to run: {output_file}\nNormal for low particles number"
                )
                self.all_good = False

        if self.all_good:
            print('Everything seems OK')

    def delete_files(self):
        for ic_file in self.input_files:
            os.remove(ic_file)
        for output_file in self.output_files:
            os.remove(output_file)


if __name__ == '__main__':
    if sys.argv[1] == 'None':
        num_cores = os.cpu_count()
    else:
        num_cores = int(sys.argv[1])
    print('-----------------------------------')
    print(f'Test execution time with {num_cores} processes')
    test = TestEx3extimes()
    print('-----------------------------------')
    print('Generating ics for test')
    test.generate_ics()
    test.setUp()
    start_time = time.time()
    print('Starting simulation tests')
    test.test_times_many_particles()
    end_time = time.time()
    print(f'Test execution time: {end_time - start_time} seconds')
    test.delete_files()
    print('-----------------------------------')
