import subprocess
import sys
import os
import time

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))
sys.path.append(parent_dir)

import numpy as np

from visual import utils


class TestEx3Output:
    def setUp(self, file):
        self.input_file = file
        self.all_good = True

    def test_2part(self, file):
        try:
            result = subprocess.run(
                ['mpirun',  '-np', str(num_cores), './ex3', self.input_file],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )
        except subprocess.CalledProcessError as e:
            print(f"ex3 failed with error: {e.stderr}")

        pos, ids, mass, t, exec_time = utils.read_data(file)
        r = np.linalg.norm(pos, axis=1)

        if len(mass) != 396:
            print("Something wrong with output time steps")
            self.all_good = False

        if r[-1] > 1.01 or r[-1] < 0.98 or r[-2] > 1e-5:
            print("Something wrong in the last position at least")
            self.all_good = False

        if np.any((r[ids[1]] > 1.01) | (r[ids[1]] < 0.98)) or np.any(r[ids[0]] > 1e-5):
            print("Something wrong with the positions")
            self.all_good = False

        if not (self.all_good):
            print('Something is wrong, check the code')
        if self.all_good:
            print('Everything seems OK')

    def test_5part(self, file):
        try:
            result = subprocess.run(
                ['mpirun',  '-np', str(num_cores), './ex3', self.input_file],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )
        except subprocess.CalledProcessError as e:
            print(f"ex3 failed with error: {e.stderr}")

        pos, ids, mass, t, exec_time = utils.read_data(file)
        r = np.linalg.norm(pos, axis=1)

        if len(mass) != 2475:
            print("Something wrong with output time steps")
            self.all_good = False

        if r[-3] > 1.01 and r[-3] < 0.98 and r[-5] > 1e-5:
            print("Something wrong in the last position at least")
            self.all_good = False

        if (
            np.any((r[ids[2]] > 1.01) | (r[ids[2]] < 0.98))
            or np.any(r[ids[0]] > 1e-5)
            or np.any((r[ids[1]] > 0.71) | (r[ids[1]] < 0.69))
            or np.any((r[ids[3]] > 1.6) | (r[ids[3]] < 1.4))
        ):
            print("Something wrong with the positions")
            self.all_good = False

        if not (self.all_good):
            print('Something is wrong, check the code')
        if self.all_good:
            print('Everything seems OK')


if __name__ == '__main__':
    if sys.argv[1] == 'None':
        num_cores = os.cpu_count()
    else:
        num_cores = int(sys.argv[1])
    print('-----------------------------------')
    print(f'Running tests for ex3 with {num_cores} processes')
    print('-----------------------------------')
    print('Test 2 particles stable orbit')
    test = TestEx3Output()
    test.setUp('./ics/ic_test_2part_stable_orbit.txt')
    start_time = time.time()
    test.test_2part('./output/test_2part_stable_orbit.dat')
    end_time = time.time()
    print(f'Test execution time: {end_time - start_time} seconds')
    print('-----------------------------------')
    print('Test 5 particles stable orbit')
    test = TestEx3Output()
    test.setUp('./ics/ic_test_5part_stable_orbit.txt')
    start_time = time.time()
    test.test_5part('./output/test_5part_stable_orbit.dat')
    end_time = time.time()
    print(f'Test execution time: {end_time - start_time} seconds')
