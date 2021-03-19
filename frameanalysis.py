import sys
from analysisclass.analysis import *


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

    # TODO fix writing outputs
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
