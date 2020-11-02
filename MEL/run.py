"""
Main script to run MEL. 
executes the jobs of input_lumping.txt and generates folders with the output
"""

import os
import sys
import pandas as pd
from time import perf_counter as clock
from . import license_message
from . import A_read_inputGIT as readinp
from . import B_extract_rates as extr_rates

def main():
    # print license message
    print(license_message.message)

    #################### input reading and checking #######################
    tic_inpread = clock()

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

    toc_inpread = clock()  

    ################ mechanism reading and checking ##########################
    tic_mechread = clock()
    mech_dict = extr_rates.readmechanism(input_par['mech_type'],cwd)
    toc_mechread = clock()

    # check the consistency of the input: parameters section
    try:
        inp_instr.COMPARE_INPUT_PT(input_par['mech_type'],mech_dict)
    except Exception as e:
        print(str(e))
        exit()

    # check the input based on job type
    try:
        inp_instr.CHECK_INPUT_JOBTYPE(job_list,mech_dict)
    except Exception as e:
        print(str(e))
        exit()

    # create dataframe where to store times for each simulation
    tictoc = pd.Series(index=job_list.keys())
    tictoc['input reading'] = toc_inpread-tic_inpread
    tictoc['mech reading'] = toc_mechread-tic_mechread

    ################ loops over each job list ##########################
    for key, value in job_list.items():
        
        # for single simulation: re-set job type
        if key == 'single_simulation':
            jobtype = value['simul_type']
        elif key == 'prescreening_equilibrium' or key == 'prescreening_allreactive' or key == 'composition_selection' or key == 'lumping' or key == 'validation':
            jobtype = key

        # call subdictionaries
        job_subdict = value

        print('\nstarting task: ' + key)

        # set rest of input parameters (except reac/prod) based on simulation type
        input_par_jobtype = inp_instr.set_inputparam_job(jobtype)
        print(input_par_jobtype)

        # set iterative operations for each type of simulations

            # iterate over the selected set of species
            # check if files exist and if not run the simulation / generate them

        

if __name__ == "__main__":
    main()
