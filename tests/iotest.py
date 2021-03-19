import unittest
from analysisclass.selection.io import *


class TestIOFunctions(unittest.TestCase):
    def test_read_input0(self):
        input0 = read_input('test0.in')
        self.assertEqual(input0['topology'], 'md_in/testtop.psf')
        self.assertEqual(input0['trajectory'], ['md_in/testtraj.dcd'])
        self.assertEqual(input0['residues'], [[1, 5], [3, 4]])
        self.assertEqual(input0['functions'], [False, True, False, False])
        self.assertEqual(input0['start'], 0)
        self.assertEqual(input0['stop'], 40)
        self.assertEqual(input0['step'], 1)
        self.assertEqual(input0['out_dir'], 'md_out')

    def test_read_inputs1(self):
        input1 = read_input('test1.in')
        self.assertEqual(input1['topology'], 'md_in/testtop.psf')
        self.assertEqual(input1['trajectory'], ['md_in/testtraj.dcd', 'md_in/testtraj.dcd'])
        self.assertEqual(input1['residues'], [[1, 2], [3, 5], [2, 7], [8, 2]])
        self.assertEqual(input1['functions'], [False, True, False, False])
        self.assertEqual(input1['start'], 0)
        self.assertEqual(input1['stop'], 50)
        self.assertEqual(input1['step'], 2)
        self.assertEqual(input1['out_dir'], 'md_out')

    def test_write_dataframe(self):
        results = {'1-2': [0, 1, 2, 3, 4]}
        time_list = [0, 1, 2, 3, 4]
        data = write_dataframe(results, time_list)
        print(data)

    def test_write_csv(self):
        pass


if __name__ == '__main__':
    unittest.main()
