import unittest
import numpy as np
from analysisclass.analysis import *


# Just tests for structure errors
class TestFrameAnalysis(unittest.TestCase):
    def test_dist_array0(self):
        # importing data
        unv = build_universe('testfiles/testtop.psf', 'testfiles/testtraj.dcd')
        func_act = [True, False, False, False]
        res_pairs = [[1, 2]]
        dist_results = FrameAnalysis(unv, res_pairs, func_act).run(start=0, stop=50)
        traj_time = dist_results.time_list
        array_keys = dist_results.res_keys
        array_result = dist_results.res_dists
        self.assertTrue(len(traj_time), 50)
        self.assertEqual(array_keys, ['1-2'])
        for arr in array_result['1-2']:
            self.assertEqual(np.shape(arr), (19*24,))

    def test_ca_dist0(self):
        unv = build_universe('testfiles/testtop.psf', 'testfiles/testtraj.dcd')
        func_act = [False, True, False, False]
        res_pairs = [[1, 2]]
        dist_results = FrameAnalysis(unv, res_pairs, func_act).run(start=0, stop=50)
        traj_time = dist_results.time_list
        array_keys = dist_results.ca_keys
        array_result = dist_results.ca_dists
        self.assertEqual(len(traj_time), 50)
        self.assertEqual(array_keys, ['1-2'])
        self.assertEqual(len(array_result['1-2']), 50)

    def test_cm_dist0(self):
        unv = build_universe('testfiles/testtop.psf', 'testfiles/testtraj.dcd')
        func_act = [False, False, True, False]
        res_pairs = [[1, 2]]
        dist_results = FrameAnalysis(unv, res_pairs, func_act).run(start=0, stop=50)
        traj_time = dist_results.time_list
        array_keys = dist_results.cm_keys
        array_result = dist_results.cm_dists
        self.assertEqual(len(traj_time), 50)
        self.assertEqual(array_keys, ['1-2'])
        self.assertEqual(len(array_result['1-2']), 50)

    def test_cg_dist0(self):
        unv = build_universe('testfiles/testtop.psf', 'testfiles/testtraj.dcd')
        func_act = [False, False, False, True]
        res_pairs = [[1, 2]]
        dist_results = FrameAnalysis(unv, res_pairs, func_act).run(start=0, stop=50)
        traj_time = dist_results.time_list
        array_keys = dist_results.cg_keys
        array_result = dist_results.cg_dists
        self.assertEqual(len(traj_time), 50)
        self.assertEqual(array_keys, ['1-2'])
        self.assertEqual(len(array_result['1-2']), 50)


if __name__ == '__main__':
    unittest.main()
