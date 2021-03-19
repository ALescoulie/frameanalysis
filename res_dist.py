import MDAnalysis as mda
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.analysis.distances import dist


def remove_return(string):
    ret_line = string.split('\n')
    return ret_line[0]


def read_input(path):
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
                    topology_path = remove_return(text[1])
                elif 'trajectory_path' in line:
                    input_line = line
                    text = input_line.split('=')
                    trajectory_path.append(remove_return(text[1]))

            elif section == 1 and '1:' not in line:
                input_line = line
                text = input_line.split()
                resid_pair[0] = int(text[0])
                resid_pair[1] = int(text[1])
                resid_pair_list.append(resid_pair.copy())

            elif section == 2 and '2:' not in line:
                input_line = line
                text = input_line.split('=')
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


def write_dataframe(result_dict, time_list):
    column_key = []
    result_dict['Time'] = time_list
    for key in result_dict:
        column_key.append(key)

    data = np.zeros((len(time_list), len(column_key)))

    for header in range(len(column_key)):
        key = column_key[header]
        for i in range(data.shape[0]):
            data[i, header] = result_dict[key][i]

    result_frame = pd.DataFrame(data, columns=column_key)
    return result_frame


def write_csv(dataframe, out_dir):
    dataframe.to_csv(out_dir, header=True)
    return


def save_figure(csv_path):
    def read_csv(path):
        file_headers = []

        with open(path, 'r') as datafile:
            data = datafile.readlines()
            header_line = data[0]
            if '\n' in header_line:
                header_line = remove_return(header_line)

            items = header_line.split(',')
            for item in items:
                file_headers.append(item)

        file_headers = file_headers[1:]
        csv_data = np.genfromtxt(path, delimiter=',', skip_header=1)
        return csv_data, file_headers

    def build_figure(data, figure_headers, png_name):
        for col in range(0, len(headers) - 1):
            fig = plt.plot(data[:, col + 1], label=figure_headers[col])
            plt.legend()

        plt.xlabel('Time (ps)')
        plt.ylabel('Pair Distance (angstrom)')
        plt.savefig(png_name, dpi=600)
        return

    dist_data, headers = read_csv(csv_path)
    csv_path = csv_path.split('.')
    file_name = csv_path[0] + '.png'
    build_figure(dist_data, headers, file_name)
    return


def build_universe(topology_path, *args):
    universe = mda.Universe(topology_path, args)
    return universe


# Selects residue atoms from an its integer residue number
def select_resnum_atoms(res_num, universe):
    resid_num = 'resid ' + str(res_num)
    atoms = universe.select_atoms(resid_num)
    return atoms


# Selects residue alpha carbon from its integer residue number
def select_resid_ca(res_num, universe):
    selection = 'resid ' + str(res_num) + ' and name CA'
    res_ca = universe.select_atoms(selection)
    return res_ca


def save_resid_cm(res_num, universe):
    res_atoms = select_resnum_atoms(res_num, universe)
    coords = res_atoms.center_of_mass()
    return coords


def save_resid_cg(res_num, universe):
    res_atoms = select_resnum_atoms(res_num, universe)
    coords = res_atoms.center_of_geometry()
    return coords


# Writes atom group coordinates into numpy array
def save_atom_coords(atom_group):
    coords = atom_group.positions
    return coords


# Selects residue atoms and writes coordinates into numpy array
def save_resnum_coords(res_num, universe):
    res_atoms = select_resnum_atoms(res_num, universe)
    res_coords = save_atom_coords(res_atoms)
    return res_coords


def save_ca_coords(res_num, universe):
    ca = select_resid_ca(res_num, universe)
    ca_coords = save_atom_coords(ca)
    return ca_coords


def build_reslist_dict(res_pair_list):
    res_pairs = {}
    for respair in res_pair_list:
        key_name = str(respair[0]) + '-' + str(respair[1])
        res_pairs[key_name] = []
    res_keys = list(res_pairs.keys())
    return res_pairs, res_keys


def dist_calc(coordinates_a, coordinates_b):
    group_dist = distance_array(coordinates_a, coordinates_b)
    return group_dist


def atom_dist(cooridnates_a, cooridnates_b):
    distance = dist(cooridnates_a, cooridnates_b)
    return distance


class FrameAnalysis(AnalysisBase):
    def __init__(self, universe, res_pairs, func_set):
        super(FrameAnalysis, self).__init__(universe.trajectory)
        self._unv = universe        # Universe being analyzed
        self._rpl = res_pairs       # First list of int residue numbers
        self._fxn = func_set        # Bool values enabling functions

    def _prepare(self):
        """ Builds dictionary to store function outputs computed for inputted residue pairs.
            Keys are string of the respective residue numbers punctuated by a hyphen ex: '121-145'
            dist_calc outputs 2d numpy array with dimensions equal to atoms of each residue.
            Each frame will output an array of the distances of calculated for the given reside pairs"""
        # Time list
        self.time_list = []
        # Distance array
        if self._fxn[0] is True:
            self.res_dists, self.res_keys = build_reslist_dict(self._rpl)

        # Distance between alpha carbons
        if self._fxn[1] is True:
            self.ca_dists, self.ca_keys = build_reslist_dict(self._rpl)

        # Distance between resid center of mass
        if self._fxn[2] is True:
            self.cm_dists, self.cm_keys = build_reslist_dict(self._rpl)

        # Distance between resid center of geometry
        if self._fxn[3] is True:
            self.cg_dists, self.cg_keys = build_reslist_dict(self._rpl)

    def _single_frame(self):
        # Saving time values
        self.time_list.append(self._unv.trajectory.time)
        # Only running functions specified in input file
        if self._fxn[0] is True:
            # Iterating through index pairs and returning distance array to dictionary
            for key in range(len(self.res_keys)):
                self.cord1 = save_resnum_coords(self._rpl[key][0], self._unv)
                self.cord2 = save_resnum_coords(self._rpl[key][1], self._unv)
                self.res_dists[self.res_keys[key]].append(dist_calc(self.cord1, self.cord2))

        if self._fxn[1] is True:
            # Iterating through index pairs and returning alpha carbon distance
            for key in range(len(self.ca_keys)):
                self.cord1 = save_ca_coords(self._rpl[key][0], self._unv)
                self.cord2 = save_ca_coords(self._rpl[key][1], self._unv)
                self.ca_dists[self.ca_keys[key]].append(dist_calc(self.cord1, self.cord2))

        if self._fxn[2] is True:
            # Iterating through index pairs and returning residue center of mass distance
            for key in range(len(self.cm_keys)):
                self.cord1 = save_resid_cm(self._rpl[key][0], self._unv)
                self.cord2 = save_resid_cm(self._rpl[key][1], self._unv)
                self.cm_dists[self.cm_keys[key]].append(dist_calc(self.cord1, self.cord2))

        if self._fxn[3] is True:
            # Iterating through index pairs and returning residue center of geometry distance
            for key in range(len(self.cg_keys)):
                self.cord1 = save_resid_cg(self._rpl[key][0], self._unv)
                self.cord2 = save_resid_cg(self._rpl[key][1], self._unv)
                self.cg_dists[self.cg_keys[key]].append(dist_calc(self.cord1, self.cord2))

    def _conclude(self):
        if self._fxn[0] is True:
            pass

        if self._fxn[1] is True:
            self.ca_data = write_dataframe(self.ca_dists, self.time_list)

        if self._fxn[2] is True:
            self.cm_data = write_dataframe(self.cm_dists, self.time_list)

        if self._fxn[3] is True:
            self.cg_data = write_dataframe(self.cg_dists, self.time_list)


if __name__ == '__main__':
    # Accepting arguments
    input_path = str(sys.argv[1])
    input_test = 'template_input'

    # Getting inputs
    input_variables = read_input(input_test)

    # Defining universe
    unv = build_universe(input_variables['topology'], input_variables['trajectory'])

    # Running analysis
    analysis_output = FrameAnalysis(unv, input_variables['residues'], input_variables['functions']) \
        .run(start=input_variables['start'], step=input_variables['step'], stop=input_variables['stop'])

    # Writing output
    if input_variables['functions'][0]:
        pass

    elif input_variables['functions'][1]:
        out_name = input_variables['out_dir'] + '/' + 'ca.csv'
        write_csv(analysis_output.ca_keys, out_name)
        save_figure(out_name)

    elif input_variables['functions'][2]:
        out_name = input_variables['out_dir'] + '/' + 'cm.csv'
        write_csv(analysis_output.cm_keys, out_name)
        save_figure(out_name)

    elif input_variables['functions'][3]:
        out_name = input_variables['out_dir'] + '/' + 'cg.csv'
        write_csv(analysis_output.cg_keys, out_name)
        save_figure(out_name)
