import unittest
from analysisclass.analysis import *


class TestFullAnalysis(unittest.TestCase):
    def test_frame_analysis0(self):
        # Getting and testing input
        input_settings = read_input('test2.in')
        self.assertEqual(input_settings['topology'], 'testtop.psf')
        self.assertEqual(input_settings['trajectory'], ['testtraj.dcd'])
        self.assertEqual(input_settings['residues'], [[1, 5], [3, 4]])
        self.assertEqual(input_settings['functions'], [False, True, False, False])
        self.assertEqual(input_settings['start'], 0)
        self.assertEqual(input_settings['stop'], 40)
        self.assertEqual(input_settings['step'], 1)
        self.assertEqual(input_settings['out_dir'], 'md_out')

        # Defining universe and calculating output
        unv = build_universe(input_settings['topology'], input_settings['trajectory'])

        output = FrameAnalysis(unv, input_settings['residues'], input_settings['functions']) \
            .run(start=input_settings['start'], step=input_settings['step'], stop=input_settings['stop'])

        # Results Testing
        self.assertEqual(len(output.time_list), 40)
        self.assertEqual(output.ca_keys, ['1-5', '3-4'])

        # Testing outputs
        out_name = input_settings['out_dir'] + '/' + 'ca.csv'
        write_csv(output.ca_data, out_name)
        save_figure(out_name)

        with open(out_name, 'r') as outfile:
            out_data = outfile.readlines()

            for line in out_data:
                if 'Time' in line:
                    line = line.split(',')
                    self.assertEqual(line, ['', '1-5', '3-4', 'Time\n'])

    def test_frame_analysis1(self):
        # Getting and testing input
        input_settings = read_input('test3.in')
        self.assertEqual(input_settings['topology'], 'testtop.psf')
        self.assertEqual(input_settings['trajectory'], ['testtraj.dcd'])
        self.assertEqual(input_settings['residues'], [[1, 2], [3, 5], [2, 7], [8, 2]])
        self.assertEqual(input_settings['functions'], [True, False, False, False])
        self.assertEqual(input_settings['start'], 0)
        self.assertEqual(input_settings['stop'], 50)
        self.assertEqual(input_settings['step'], 1)
        self.assertEqual(input_settings['out_dir'], 'md_out')

        # Defining universe and calculating output
        unv = build_universe(input_settings['topology'], input_settings['trajectory'])

        output = FrameAnalysis(unv, input_settings['residues'], input_settings['functions']) \
            .run(start=input_settings['start'], step=input_settings['step'], stop=input_settings['stop'])

        self.assertEqual(len(output.time_list), 50)
        self.assertEqual(output.ca_keys, ['1-2', '3-5', '2-7', '8-2'])
        out_name = input_settings['out_dir'] + '/' + 'res_dist.csv'
        write_csv(output.res_data, out_name)
        save_figure(out_name)

        with open(out_name, 'r') as outfile:
            out_data = outfile.readlines()

            for line in out_data:
                if 'Time' in line:
                    line = line.split(',')
                    self.assertEqual(line, ['', '1-2', '3-5', '2-7', '8-2', 'Time\n'])


if __name__ == '__main__':
    unittest.main()
