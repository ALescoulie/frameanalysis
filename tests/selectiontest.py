import unittest
import numpy as np
from analysisclass.selection.selectionfuncs import *


class SelectionTest(unittest.TestCase):
    def test_build_universe0(self):
        unv = build_universe('testtop.psf', 'testtraj.dcd')
        self.assertEqual(len(unv.trajectory), 98)
        self.assertEqual(len(unv.residues), 214)

    def test_select_resnum_atoms0(self):
        unv = build_universe('testtop.psf', 'testtraj.dcd')
        res_group = select_resnum_atoms(1, unv)
        self.assertEqual(len(res_group), 19)

    def test_select_resnum_atoms1(self):
        unv = build_universe('testtop.psf', 'testtraj.dcd')
        res_group = select_resnum_atoms(2, unv)
        self.assertEqual(len(res_group), 24)

    def test_save_resnum_coords0(self):
        unv = build_universe('testtop.psf', 'testtraj.dcd')
        ca_coords = save_ca_coords(1, unv)
        self.assertEqual(np.shape(ca_coords), (1, 3))

    def test_save_resnum_coords1(self):
        unv = build_universe('testtop.psf', 'testtraj.dcd')
        resid_coords = save_resnum_coords(1, unv)
        self.assertEqual(np.shape(resid_coords), (19, 3))

    def test_build_reslist_dict0(self):
        resid_list = [[1, 2], [3, 4], [6, 7]]
        resid_dict, resid_keys = build_reslist_dict(resid_list)

        # Testing key list
        self.assertEqual(len(resid_keys), len(resid_list))
        self.assertEqual(resid_keys[0], '1-2')
        self.assertEqual(resid_keys[1], '3-4')
        self.assertEqual(resid_keys[2], '6-7')

        # Testing dictionary keys
        dict_keys = []
        for key in resid_dict:
            dict_keys.append(key)
        self.assertEqual(dict_keys, resid_keys)


if __name__ == '__main__':
    unittest.main()
