import unittest
from analysisclass.analysis import *


class TestFrameAnalysis(unittest.TestCase):
    def test_dist_array0(self):
        # importing data
        unv = build_universe('testtop.psf', 'testtraj.dcd')
        func_act = [True, False, False, False]
        res_pairs = [1, 2]
        dist_results = FrameAnalysis(unv, res_pairs, func_act).run(start=0, stop=50)
        traj_time = dist_results.time_list
        array_keys = dist_results.res_keys
        array_result = dist_results.res_dists
        self.assertTrue(len(traj_time), 50)
        self.assertEqual(array_keys, ['1-2'])


if __name__ == '__main__':
    unittest.main()
