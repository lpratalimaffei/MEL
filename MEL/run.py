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
from . import C_preprocessing as preproc
from . import main_flow as sim
from . import H_SET_LOOPS_FLD as set_sim
from . import preproc_irreversible as preproc_irr


def main():
    # print license message
    print(license_message.message)

    #################### input reading and checking #######################
    tic_inpread = clock()

    try:
        inputfile = sys.argv[1]
    except IndexError:
        print('input file required')
        sys.exit()

    # read the input
    cwd = os.path.abspath(os.getcwd())
    try:
        inp_instr = readinp.READ_INPUT(cwd, inputfile)
    except RuntimeError as e:
        print(str(e))
        sys.exit()

    # extract input + job_list and corresponding subdictionaries
    input_par, job_list = inp_instr.read_file_lines()

    # first check of the input: exit in case of exceptions
    try:
        inp_instr.CHECK_INPUT()
    except Exception as e:
        print(str(e))
        sys.exit()

    toc_inpread = clock()

    ################ mechanism reading and checking ##########################
    tic_mechread = clock()
    mech_dict = extr_rates.readmechanism(input_par['mech_type'], cwd)
    toc_mechread = clock()

    # check the consistency of the input: parameters section
    try:
        inp_instr.COMPARE_INPUT_PT(input_par['mech_type'], mech_dict)
    except Exception as e:
        print(str(e))
        sys.exit()

    # check the input based on job type
    try:
        inp_instr.CHECK_INPUT_JOBTYPE(job_list, mech_dict)
    except Exception as e:
        print(str(e))
        sys.exit()

    # create dataframe where to store times for each simulation
    tictoc = pd.Series(index=job_list.keys())
    tictoc['input reading'] = toc_inpread-tic_inpread
    tictoc['mech reading'] = toc_mechread-tic_mechread

    ################ loops over each job list ##########################
    for key, value in job_list.items():

        # for single simulation: re-set job type
        if key == 'single_simulation':
            jobtype = value['simul_type']
        elif key in ['prescreening_equilibrium', 'prescreening_allreactive', 'composition_selection',
                    'lumping', 'validation', 'preproc_irreversible']:
            jobtype = key

        # call subdictionaries
        job_subdict = value

        # set optional parameters for composition_selection
        if jobtype == 'composition_selection':
            BF_tol = job_subdict['BF_tol']
            maxiter = job_subdict['maxiter']
            opts = [BF_tol, maxiter]
        else:
            opts = ['', '']

        print('\nstarting task: ' + key + '\n')

        # set rest of input parameters (except reac/prod) based on simulation type
        input_par_jobtype = inp_instr.set_inputparam_job(jobtype)
        # set iterative operations for each type of simulations
        sim_DF = set_sim.set_simul_loop(cwd, jobtype, job_subdict, mech_dict)
        print(sim_DF)
        # iterate over the selected set of species
        for i in sim_DF.index:
            sim_series = sim_DF.loc[i]

            # check if folder exists, otherwise create it
            YE_NO = set_sim.setfolder(sim_series['fld'])
            # create the mechanism folder and copy the input preprocessor
            set_mechfld = set_sim.setfolder(
                os.path.join(cwd, 'mech_tocompile'))
            shutil.copy(os.path.join(cwd, 'inp', 'input_preproc.dic'),
                        os.path.join(cwd, 'mech_tocompile', 'input_preproc.dic'))

            if YE_NO == 0 and jobtype != 'preproc_irreversible':

                print('\nStart with set of reactants [{}]'.format(
                    sim_series['REAC']))
                # perform the simulation
                sim.main_simul(cwd, jobtype, input_par,
                            input_par_jobtype, mech_dict, sim_series, opts)

                # delete folders to avoid confusion and do cleaning
                set_sim.rmfolder(os.path.join(cwd, 'mech_tocompile'))
                set_sim.rmfolder(os.path.join(cwd, 'Output'))
                # delete simulation files
                os.remove(os.path.join(cwd, 'OS_output.txt'))
                os.remove(os.path.join(cwd, 'input_OS.dic'))

            elif YE_NO == 0 and jobtype == 'preproc_irreversible':
                preproc_irr.run_preproc(cwd, input_par['opensmoke_folder'])

        # for lumping: derive the full lumped mechanism from the submechs in each subfolder
        if jobtype == 'lumping' and key != 'single_simulation':
            # list of folders
            fld_list = sim_DF['fld']
            lumpedmech_fld = os.path.join(cwd, 'lumpedmech')
            YE_NO = set_sim.setfolder(lumpedmech_fld)
            set_sim.renamefiles(cwd, lumpedmech_fld)

            # write the new lumped mech
            shutil.copy(os.path.join(fld_list[0], 'therm.txt'), os.path.join(
                lumpedmech_fld, 'therm.txt'))
            # write kinetics
            preproc.COMBINE_CKI(lumpedmech_fld, fld_list)





if __name__ == "__main__":
    main()
