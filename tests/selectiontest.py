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

<<<<<<< HEAD
    def test_save_resnum_coords(self):
        unv = build_universe('testtop.psf', 'testtraj.dcd')
        ca_coords = save_ca_coords(unv, 1)


if __name__ == '__main__':
    unittest.main()
=======

if __name__ == '__main__':
    unittest.main()
>>>>>>> origin/main
