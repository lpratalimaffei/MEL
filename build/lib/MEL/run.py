"""
Main script to run MEL. 
executes the jobs of input_lumping.txt and generates folders with the output
"""

import os
import sys
import shutil
import pandas as pd
from time import perf_counter as clock
from . import license_message
from . import A_read_input as readinp
from . import B_extract_rates as extr_rates
from . import extract_RATES_V0 as sim
from . import H_SET_LOOPS_FLD as set_sim
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
        # set iterative operations for each type of simulations
        sim_DF = set_sim.set_simul_loop(cwd,jobtype,job_subdict,mech_dict)
        print(sim_DF)

        # iterate over the selected set of species
        for i in sim_DF.index:
            sim_series = sim_DF.loc[i]
            # check if folder exists, otherwise create it
            YE_NO = set_sim.setfolder(sim_series['fld'])
            # create the mechanism folder and copy the input preprocessor
            set_mechfld = set_sim.setfolder(os.path.join(cwd,'mech_tocompile'))
            shutil.copy(os.path.join(cwd,'inp','input_preproc.dic'),os.path.join(cwd,'mech_tocompile','input_preproc.dic'))
            if YE_NO == 0:
                # perform the simulation
                sim.main_simul(cwd,jobtype,input_par,input_par_jobtype,mech_dict,sim_series)

            # delete folders to avoid confusion and do cleaning
            set_sim.rmfolder(os.path.join(cwd,'mech_tocompile'))
            set_sim.rmfolder(os.path.join(cwd,'Output'))
            # delete simulation files
            os.remove(os.path.join(cwd,'OS_output.txt'))
            os.remove(os.path.join(cwd,'input_OS.dic'))

        # for lumping: derive the full lumped mechanism from the submechs in each subfolder

        


            

        

if __name__ == "__main__":
    main()
