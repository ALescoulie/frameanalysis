import os
import sys
from analysisclass.analysis import *


if __name__ == '__main__':
    # Accepting arguments
    input_path = str(sys.argv[1])
    output_path = str(sys.argv[2])
    input_test = 'template_input'

    # Getting inputs
    input_variables = read_input(input_test)

    # Defining universe
    unv = build_universe(input_variables['topology'], input_variables['trajectory'])

    # Running analysis
    analysis_output = FrameAnalysis(unv, input_variables['residues'], input_variables['functions']) \
        .run(start=input_variables['start'], step=input_variables['step'], stop=input_variables['stop'])
