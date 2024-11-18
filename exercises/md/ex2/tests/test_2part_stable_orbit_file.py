import subprocess
import sys
import os
import time

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))
sys.path.append(parent_dir)

import numpy as np

from visual import utils


class TestEx2Output:
    def setUp(self):
        self.input_file = './ics/ic_test_2part_stable_orbit.txt'
        self.all_good = True

    def test_ex2_output(self):
        try:
            result = subprocess.run(
                ['./ex2', self.input_file],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )
        except subprocess.CalledProcessError as e:
            print(f"ex2 failed with error: {e.stderr}")

        pos, ids, mass, t = utils.read_data('./output/test_2part_stable_orbit.txt')
        r = np.linalg.norm(pos, axis=1)

        if len(mass) != 398:
            print("Something wrong with output time steps")
            self.all_good = False

        if r[-1] > 1.01 and r[-1] < 0.98 and r[-2] > 1e-5:
            print("Something wrong in the last position at least")
            self.all_good = False

        if not (np.all((r[ids[1]] > 1.01) & (r[ids[1]] < 0.98))) and not (
            np.all(r[ids[0]] > 1e-5)
        ):
            print("Something wrong with the positions")
            self.all_good = False

        if not (self.all_good):
            print('Something is wrong, check the code')
        if self.all_good:
            print('Everything seems OK')


if __name__ == '__main__':
    test = TestEx2Output()
    test.setUp()
    start_time = time.time()
    test.test_ex2_output()
    end_time = time.time()
    print(f"Test execution time: {end_time - start_time} seconds")
