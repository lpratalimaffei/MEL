"""
Main script to run MEL. 
executes the jobs of input_lumping.txt and generates folders with the output
"""

import os
import sys

from . import license_message
from . import A_read_input as readinp


def main():
    # print license message
    print(license_message.message)
    # read input file
    try:
        inputfile = sys.argv[1]
    except IndexError:
        print('input file required')
        exit()

    # read the input
    cwd = os.path.abspath(os.getcwd())
    inp_instructions = readinp.READ_INPUT(cwd,inputfile)
    # extract input + job_list and corresponding subdictionaries
    input_parameters,job_list = inp_instructions.read_file_lines()
    print(input_parameters,job_list)
    # jobs list and corresponding subdictionaries

    # read the mechanism and do first checks on input compatibility
    

    # loop over the job list
    for jobid,subdict in job_list.items():
        print('\nstarting task: ' + jobid)
        # create subfolder if it does not exist

if __name__ == "__main__":
    main()
