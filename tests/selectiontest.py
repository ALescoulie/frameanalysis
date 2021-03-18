import unittest
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


if __name__ == '__main__':
    unittest.main()
