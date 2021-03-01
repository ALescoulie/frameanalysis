import analysis
from selectionfuncs import *
import sys
import os


def read_input(path):
    # Reading input
    input_file = open(path, 'r')
    input_data = input_file.readline()

    # Input data structures
    section = 0
    topology_path = ''
    trajectory_path = []
    single_frame_funcs = []
    resid_pair = [0, 0]
    resid_pair_list = []
    start = None
    stop = None
    step = None
    out_dir = ''
    with open(path, 'r') as input_data:
        # Data pulled from list
        for line in input_data:
            print(line)
            # Establishing where loop in input file
            if '1:' in line:
                section = 1

            elif '2:' in line:
                section = 2

            elif '3:' in line:
                section = 3

            elif '4:' in line:
                section = 4

            # reading input based on position in file
            elif section == 0 and '0:':
                if 'topology_path' in line:
                    input_line = line
                    text = input_line.split('=')
                    topology_path = text[1]
                elif 'trajectory_path' in line:
                    input_line = line
                    text = input_line.split('=')
                    trajectory_path.append(text[1])

            elif section == 1 and '1:' not in line:
                input_line = line
                text = input_line.split()
                resid_pair[0] = int(text[0])
                resid_pair[1] = int(text[1])
                resid_pair_list.append(resid_pair.copy())

            elif section == 2 and '2:' not in line:
                input_line = line
                text = input_line.split('=')
                print(text)
                setting = [None]
                setting[0] = text[1]
                if setting[0] == 'False\n':
                    setting[0] = False

                elif setting[0] == 'True\n':
                    setting[0] = True

                single_frame_funcs.append(setting[0])

            elif section == 3 and '3:' not in line:
                if 'start' in line:
                    input_line = line
                    text = input_line.split('=')
                    start = int(text[1])

                if 'stop' in line:
                    input_line = line
                    text = input_line.split('=')
                    stop = int(text[1])

                if 'step' in line:
                    input_line = line
                    text = input_line.split('=')
                    step = int(text[1])

            elif section == 4 and '4:' not in line:
                input_line = line
                text = input_line.split('=')
                print(text)
                out_dir = text[1]

    # Output as data
    input_vars = {'topology': topology_path,
                  'trajectory': trajectory_path,
                  'residues': resid_pair_list,
                  'functions': single_frame_funcs,
                  'start': start,
                  'stop': stop,
                  'step': step,
                  'out_dir': out_dir}
    return input_vars


if __name__ == '__main__':
    # Accepting arguments
    #input_path = str(sys.argv[0])
    #output_path = str(sys.argv[1])
    input_test = 'template_input'

    # Getting inputs
    input_variables = read_input(input_test)

    # Defining universe
    #unv = build_universe(input_variables['topology'], input_variables['trajectory'])

    # Running analysis
    #analysis_output = analysis.FrameAnalysis(unv, input_variables['residues'], input_variables['functions'])
