import analysisclass/analysis
import analysisclass/selections/selection.py
import sys
import os


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
