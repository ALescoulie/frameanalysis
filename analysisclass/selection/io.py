import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def remove_return(string: str) -> str:
    ret_line = string.split('\n')
    return ret_line[0]


def read_input(path: str) -> dict:
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

#TODO finish time array function
def write_time_array(result_dict: dict, time_list: list) -> np.array():
    pass


def write_csv(dataframe: pd.DataFrame, out_dir: str) -> None:
    dataframe.to_csv(out_dir, header=True)
    return


def save_figure(csv_path: str) -> None:
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

    def build_figure(data: np.array(), figure_headers: list, png_name: str) -> None:
        for col in range(0, len(headers) - 1):
            fig = plt.plot(data[:, col + 1], label=figure_headers[col])
            plt.legend()

        plt.xlabel('Time (ps)')
        plt.ylabel('Pair Distance (angstrom)')
        plt.savefig(png_name, dpi=1200)
        return

    dist_data, headers = read_csv(csv_path)
    csv_path = csv_path.split('.')
    file_name = csv_path[0] + '.png'
    build_figure(dist_data, headers, file_name)
    return


def write_spat_in(residue_coords, out_dir):
    pass
