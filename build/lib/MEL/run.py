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
    inp_instr = readinp.READ_INPUT(cwd,inputfile)
    # extract input + job_list and corresponding subdictionaries
    input_par,job_list = inp_instr.read_file_lines()
    print(input_par,job_list)

    # first check of the input: exit in case of exceptions
    try:
        inp_instr.CHECK_INPUT()
    except Exception as e:
        print(str(e))
        exit()
        
    # read the mechanism and do first checks on input compatibility
    

    # loop over the job list
    for key, value in job_list.items():
        print('\nstarting task: ' + key)
        
        # for single simulation: re-set job type
        if key == 'single_simulation':
            jobtype = value['simul_type']
        elif key == 'prescreening_equilibrium' or key == 'prescreening_allreactive' or key == 'composition_selection' or key == 'lumping' or key == 'validation':
            jobtype = key
        else:
            print('\nValue Error: dictionary type not recognized')
        # call subdictionaries
        subdict = value
        # check the input based on job type
        try:
            inp_instr.CHECK_INPUT_JOBTYPE(jobtype,key,subdict)
        except Exception as e:
            print(str(e))
            exit()

        # set rest of input parameters (except reac/prod) based on simulation type
        input_par_jobtype = inp_instr.set_inputparam_job(jobtype)

        # create subfolder if it does not exist

        # for single simulation: check input compatibility

        # for other simulations: set iterative operations

        

if __name__ == "__main__":
    main()
