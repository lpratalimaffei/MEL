# import modules of python
import os
import sys
import numpy as np
import pandas as pd
import scipy.linalg as la
import scipy.optimize as op
import scipy.interpolate as interp
from scipy import integrate
from scipy.integrate import ode
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import shutil
import copy

# subprocess.call('Run.bat',stdout=subprocess.DEVNULL)

# import custom modules
from . import A_read_input as readinp
from . import B_extract_rates as extr_rates
from . import C_preprocessing as preproc
# I call the function without the need of pre-pending anything
from . import D_ODESYSTEM as odesys
from . import E_PLOTTING as mine_plt
from . import F_FITS as fitall
from . import G_PROFILES as prof_CKImech
from time import perf_counter as clock

def get_OS():
    if os.name == 'nt':
        OS_exe = 'exe'
    elif os.name == 'posix':
        OS_exe = 'sh'
    return OS_exe

def get_libpath():
    checkpath = os.path.expanduser('~')+'/.bash_profile'
    if sys.platform == 'darwin' and os.path.exists(checkpath):
        # mac: add libraries
        my_env = subprocess.check_output('source ' +  checkpath + ' && echo $DYLD_LIBRARY_PATH', shell=True)
        my_env = str(my_env)[2:-3]
        exec0 = 'export DYLD_LIBRARY_PATH=' + my_env + '; '
    else:
        exec0 = ''

    return exec0

def main_simul(cwd, jobtype, input_par, input_par_jobtype, mech_dict, sim_series, opts):
    '''
    Perform simulations for 1 reacting pseudospecies at all T,P provided
    '''
    
    # Save times as dataframes for every pressure; the first will be a dictionary
    Dt_names_Pi = ['ode solving', 'time Pi', 'plotting and saving figs']

    ############### DERIVE VARIABLES REQUIRED FOR CODE FLOW #######################################

    OS_folder = input_par['opensmoke_folder']
    input_type = input_par['mech_type']
    P_VECT = input_par['P_vect']
    T_VECT = input_par['T_vect']
    T_VECT_SKIP = input_par['T_skip']
    UNITS_BIMOL = input_par['units_bimol']
    STOICH = input_par['stoich']
    CUTOFF = input_par['cutoff']

    ISOM_EQUIL = input_par_jobtype['isom_equil']
    PRODSINKS = input_par_jobtype['Prods_sinks']
    PLOT_CMP = input_par_jobtype['plot_cmp']

    ################## DERIVE VARIABLES OF THE CORRESPONDING MECHANISM  ################################

    if input_type == 'MESS':
        P_VECT_MESS = mech_dict['P_VECT_MESS']
        T_VECT_MESS = mech_dict['T_VECT_MESS']
        SPECIES = mech_dict['SPECIES']
        SPECIES_BIMOL = mech_dict['SPECIES_BIMOL']
        rates = mech_dict['rates']

    elif input_type == 'CKI':
        SPECIES = mech_dict['SPECIES']
        SPECIES_BIMOL = mech_dict['SPECIES_BIMOL']

    ################ DERIVE REACTANTS AND PRODUCTS FROM THE INPUT DICTIONARY #####################
    fld = sim_series['fld']
    REAC = sim_series['REAC']
    REACLUMPED = sim_series['REACLUMPED']
    PRODS = sim_series['PRODS']
    PRODSLUMPED = sim_series['PRODSLUMPED']

    ################# OS EXECUTION FOLDERS AND COMMANDS ########################################
    OS_exe = get_OS()
    exec0 = get_libpath()
    preproc_exe = os.path.join('"' + OS_folder, "OpenSMOKEpp_CHEMKIN_PreProcessor." + OS_exe + '"')
    osbatch_exe = os.path.join('"' + OS_folder, "OpenSMOKEpp_BatchReactor." + OS_exe + '"')
    input_preproc = os.path.join(os.path.join(".", "mech_tocompile", "input_preproc.dic"))    
    output_preproc = os.path.join(".", "mech_tocompile", "preproc_output.txt")
    
    ################## READ OPTIONS CONCERNING MAXIT AND BF TOLERANCE (ONLY FOR COMPOSITION_SELECTIO; OTHERWISE EMPTY) ######
    BF_tol = opts[0]
    maxiter = opts[1]

    # PROCESSING OF SPECIES NAMES AND REACTANT:

    SPECIES_SERIES = pd.Series(np.arange(0, len(SPECIES)), index=SPECIES)
    # species values are the names of the bimol species corresponding to the "primary" one
    SPECIES_BIMOL_SERIES = pd.Series(SPECIES_BIMOL, index=SPECIES)

    i_PRODS = SPECIES_SERIES[PRODS].values  # array
    if isinstance(REAC, np.ndarray):
        # value (single reactant) or array (lumped reactant)
        i_REAC = SPECIES_SERIES[REAC].values
    else:
        i_REAC = SPECIES_SERIES[REAC]

    # PREPROCESSING: WRITE THE THERMODYNAMIC FILE
    if input_type == 'MESS':
        print('writing therm.txt for OS preprocessor ...'), preproc.WRITE_THERM(
            os.path.join(cwd, 'mech_tocompile'), STOICH, SPECIES_SERIES, SPECIES_BIMOL_SERIES)
    elif input_type == 'CKI':
        print('copy therm.txt to preproc folder ...')
        shutil.copy(os.path.join(cwd, 'inp', 'therm.txt'),
                    os.path.join(cwd, 'mech_tocompile', 'therm.txt'))
    # ASSIGN THE INITIAL NUMBER OF MOLES
    N_INIT_REAC = preproc.INITIALMOLES(REAC, SPECIES_BIMOL_SERIES, 1)

    # reduce T range if necessary
    # skip the temperature selected in the input
    for TSKIP in T_VECT_SKIP:
        mask = np.where(TSKIP != T_VECT)
        T_VECT = T_VECT[mask]

    ######################## PREPROCESSING OF CKI MECHANISM #########################
    if input_type == 'CKI':
        #
        extr_rates.copy_CKI_processed(os.path.join(cwd, 'inp'), os.path.join(
            cwd, 'mech_tocompile'), PRODSINKS, ISOM_EQUIL, REAC, PRODS)

        toexecute = exec0 + preproc_exe + " --input " + input_preproc + ">" + output_preproc

        print('compiling mech ...'),  subprocess.run(toexecute, shell=True)

    ####################################################################################

    ######################## LOOP OVER THE SELECTED PRESSURES ########################
    # call the class of the output processing
    postproc = odesys.ODE_POSTPROC(cwd)

    # PREALLOCATIONS
    Dt_Pi = np.zeros((len(P_VECT), len(Dt_names_Pi)))
    Pi = 0
    # dictionary for each pressure
    profiles_P = dict.fromkeys(list(P_VECT))
    arrfit_P = dict.fromkeys(list(P_VECT))
    # reactant profiles for lumped reactant
    if isinstance(REAC, np.ndarray):
        profiles_P_reac = dict.fromkeys(list(P_VECT))

    T_VECT_0 = T_VECT

    for P in P_VECT:

        print('processing: P = ' + str(P) + ' atm ... ')
        ticP = clock()

        # CHECK REACTIVITY UNTIL THE SELECTED TEMPERATURE
        # extract the rates in the selected pressure
        if input_type == 'MESS':
            # SET T_VECT=T_VECT_0 TO MAKE TEMPERATURE REDUCTION INDEPENDENT OF THE PREVIOUS PRESSURE
            T_VECT = T_VECT_0
            if isinstance(REAC, np.ndarray):
                # IN CASE OF A LUMPED REACTANT: DERIVE A TOTAL RATE CONSTANT WITH A LOOP
                kR_j = np.zeros((len(T_VECT_MESS), len(SPECIES)-1))
                # STORE REACTIVITY OF EACH REACTANT TO SELECT BF APPROPRIATELY
                T_VECT_rr_len = pd.Series(index=REAC, dtype=int)

                for rr in REAC:
                    kR_j_rr = extr_rates.REAC_P(
                        P, rr, P_VECT_MESS, T_VECT_MESS, SPECIES, rates)
                    flag_rr, Err_rr, T_VECT_rr = preproc.CHECK_REACTIVITY(
                        kR_j_rr, T_VECT, T_VECT_MESS)
                    T_VECT_rr_len[rr] = len(T_VECT_rr)
                    kR_j = kR_j + kR_j_rr

            elif isinstance(REAC, str):
                # NON LUMPED REACTANT: DO A SINGLE DERIVATION OF THE REACTIVITY MATRIX
                kR_j = extr_rates.REAC_P(
                    P, REAC, P_VECT_MESS, T_VECT_MESS, SPECIES, rates)

            flag_exit, Err, T_VECT = preproc.CHECK_REACTIVITY(
                kR_j, T_VECT, T_VECT_MESS)
            print(Err)
            if flag_exit == 1:
                exit()

        elif input_type == 'CKI' and isinstance(REAC, np.ndarray):
            T_VECT_rr_len = pd.Series(
                len(T_VECT)*np.ones(REAC.shape), index=REAC, dtype=int)

        # PRE PROCESSING FOR LUMPED REACTANTS AND PRODUCTS   ##################à

        # FOR LUMPED REACTANT: EXTRACT IN A DATAFRAME THE COMPOSITION OF THE REACTANT AT EVERY PRESSURE
        # COMPARE THE BRANCHINGS WITH THE T RANGE AND PRINT A WARNING IF YOU REDUCE IT

        if isinstance(REAC, np.ndarray):
            try:
                BR_L_REAC, T_VECT = preproc.BRANCHING_LUMPEDREAC(os.path.join(
                    cwd, jobtype), REACLUMPED, REAC, T_VECT, T_VECT_rr_len, P, P_VECT)
            except ValueError as e:
                print('error while extracting the branching fractions of ' +
                      REACLUMPED.index[0] + ': ' + str(e))
                exit()

        elif isinstance(REAC, str):
            # create empty variable (so that the input will exist anyways)
            BR_L_REAC = ''

        # FOR LUMPED PRODUCTS or reactants: GENERATE DATAFRAMES FOR BRANCHINGS
        if len(PRODS) != len(PRODSLUMPED) or len(REAC) != len(REACLUMPED):
            try:
                postproc.MAKE_BRANCHINGFOLDERS(
                    jobtype, REACLUMPED, PRODSLUMPED, T_VECT)
            except ValueError as e:
                print(e)

        # inside the profile_P dictionary: allocate a profiles_T dictionary to include the profiles at every temperature
        profiles_P[P] = dict.fromkeys(list(T_VECT))
        if isinstance(REAC, np.ndarray):
            profiles_P_reac[P] = dict.fromkeys(list(T_VECT))

        # todo: for the composition_seleciton: iterate until BFs converge or until you reach max number of iterations
        # at the beginning of each iteration: set BR_L_REAC = BF_OUTPUT
        # at the end: copy the selected BF to lumping/BF_INPUT and validation/BF_INPUT
        flag_BFiter = 0
        it = 0
        while flag_BFiter == 0:
            ########################### LOOP OVER ALL THE TEMPERATURES OF INTEREST ####################

            Ti = 0

            for T in T_VECT:
                print('processing: T = ' + str(T) + ' K ... ')
                Ti += 1

                # make the folder of the output to optiSMOKE++
                postproc.MAKE_FOLDERS(fld, P, T, REACLUMPED)

                ################# PREPROCESSING - MESS INPUT ONLY #######################################################
                if input_type == 'MESS':
                    # EXTRACT THE KINETIC MATRIX
                    k_ij_TP = extr_rates.MATRIX_TP(
                        T, P, T_VECT_MESS, P_VECT_MESS, SPECIES, rates)

                    # PREPROCESSING: SET THE PRODUCTS AS IRREVERSIBLE SINKS
                    k_ij_TP_proc = preproc.PREPROCESSING(
                        k_ij_TP, SPECIES_SERIES, PRODS)
                    k_ij_TP_prodsinks = k_ij_TP_proc.SET_PRODS_SINKS(PRODSINKS)
                    k_ij_TP_OK = k_ij_TP_proc.MULTIPLY_BIMOL(
                        SPECIES_BIMOL, UNITS_BIMOL)

                    ################# WRITE FICTITIOUS MECH IN CKI FORMAT ##################################

                    k_to_CKI = preproc.WRITE_MECH_CKI(k_ij_TP_OK, np.zeros(k_ij_TP_OK.shape), np.zeros(
                        k_ij_TP_OK.shape), np.zeros(k_ij_TP_OK.shape, dtype=str), SPECIES_SERIES, i_REAC, i_PRODS, SPECIES_BIMOL)

                    print('writing kin.txt for OS preprocessor ...')
                    # ISOM_EQUIL is put here: if you want to study the isomer equilibrium, all the other reactions are not copied
                    CKI_lines = k_to_CKI.MAKE_CKI(PRODSINKS, ISOM_EQUIL)
                    try:
                        k_to_CKI.WRITE_CKI(os.path.join(
                            cwd, 'mech_tocompile'), CKI_lines.values)
                    except RuntimeError as e:
                        print(str(e))
                        exit()

                    toexecute = exec0 + preproc_exe + " --input " + input_preproc + ">" + output_preproc
                    print('compiling mech ...'), subprocess.run(toexecute, shell=True)

                ################### WRITE OPENSMOKE INPUT ##################################################
                try:
                    OS_write = preproc.WRITE_OS_INPUT(
                        cwd, T, P, SPECIES, REACLUMPED, 1, SPECIES_BIMOL_SERIES, BR_L_REAC)
                except RuntimeError as e:
                    print(str(e))
                    exit()
                print('writing new input ...'), OS_write.w_input_OS()

                ################## SOLUTION OF THE ODE SYSTEM ##########################################
                # CALL OPENSMOKE
                tic = clock()

                toexecute = exec0 + osbatch_exe + " --input input_OS.dic > OS_output.txt"
                print('solving OS Batch Reactor ...'), subprocess.run(toexecute, shell=True)

                toc = clock()
                Dt_Pi[Pi, 0] = toc-tic

                ##################### EXTRACT THE OUTPUT ###############################################
                try:
                    tW_DF, PV = postproc.EXTRACT_PROFILES(
                        SPECIES, i_REAC, N_INIT_REAC, SPECIES_BIMOL_SERIES, ISOM_EQUIL, CUTOFF)
                except ValueError as e:
                    print(str(e))

                # process the output: rewrite tW_DF if there are lumped species
                print('Rewriting profiles for lumped species ...')
                tW_DF, REAC_L, i_REAC_L, SPECIES_L, SPECIES_SERIES_L, SPECIES_BIMOL_SERIES_L, PRODS_L = postproc.LUMP_PROFILES(
                    PRODS, PRODSLUMPED)
                # SAVE PROFILES IN A DICTIONARY FOR LATER POSTPROCESSING AND PLOTTING
                profiles_P[P][T] = tW_DF
                if isinstance(REAC, np.ndarray):
                    tW_DF_reac = postproc.PROFILES_REAC_COMPOSITION()
                    profiles_P_reac[P][T] = tW_DF_reac
                # IF YOU HAVE LUMPING (I.E.: PRODSINKS=1)
                # LINEAR FIT OF THE PROFILES OBTAINED: GUESS FOR SUCCESSIVE optiSMOKE++
                if PRODSINKS == 1:
                    # if it is the first temperature: call the class
                    # generate dataframe at every T,P to allocate the fits of k_ij
                    if Ti == 1:
                        kfit_P = fitall.FITTING(T_VECT, REAC_L, PRODS_L)
                    # then do the fits
                    k_prods_T, fit_error_T = kfit_P.fit_profiles(
                        tW_DF, i_REAC_L, SPECIES_SERIES_L, T, PV, SPECIES_BIMOL_SERIES_L, 1)
                # I think I need to do something with this fit_error_T, save it somewhere
                # maybe I can do a 3D plot with the fitting error at a certain T,P or with the difference between the sum(ki) and k of the reactant
                # and finally the plot with the branching among different channels with the lumped rates

                # WRITE OUTPUT PROFILES TO OPTISMOKE
                print('Writing profiles for Optismoke ...'), postproc.WRITE_PROFILES(
                    PRODS_L)

                # GENERATE THE NEW OS INPUT AND WRITE IT TO THE FOLDERS
                print(
                    'Generate OS input for the lumped mech and copy to the folder'), postproc.WRITE_NEW_OSINPUT(1)

            # BF processing before getting to the following cycle; fitting if needed

            # Return the dataframes of the fits of the rate constants
            if PRODSINKS == 1:
                print('Fitting the arrhenius profiles ...')
                # write the original rate constants in the corresponding folder
                out_fld = os.path.join(fld, str(P) + 'atm')
                kfit_P.write_originalk(out_fld)
                # fit
                rates_P_CKI = kfit_P.fits_lumped_k(
                    cwd, P, SPECIES_BIMOL_SERIES_L)
                arrfit_P[P] = rates_P_CKI

            if len(REAC) != len(REACLUMPED):
                print(
                    'writing branchings of lumped reactant at {} atm ...'.format(str(P)))
                BF_OUTPUT = postproc.WRITE_BRANCHINGS_REACS()

            # FOR COMPOSITION_SELECTION: COMPARE THE BF OBTAINED

            if jobtype == 'composition_selection':
                # the reactant is lumped; compare the BFs obtained
                max_deltaBF = preproc.COMPARE_BRANCHINGS(BR_L_REAC, BF_OUTPUT)
                print('iteration: {} ; max_deltaBF: {} '.format(it, max_deltaBF))

                if max_deltaBF < BF_tol or it > maxiter:

                    # write the products BF
                    if len(PRODS) != len(PRODSLUMPED):
                        print(
                            'writing branchings of lumped products at {} atm ...'.format(str(P)))
                        postproc.WRITE_BRANCHINGS_PRODS(PRODS)

                    flag_BFiter = 1
                    # COPY THE COMPOSITION TO validation and lumping directories
                    BF_validation_path = os.path.join(
                        cwd, 'validation', 'BF_INPUT', REAC_L)
                    BF_lumping_path = os.path.join(
                        cwd, 'lumping', 'BF_INPUT', REAC_L)
                    file_tocopy = os.path.join(
                        cwd, 'composition_selection', 'BF_OUTPUT', REAC_L, REAC_L, str(P)+'atm.txt')

                    if os.path.isdir(BF_validation_path) == False:
                        os.makedirs(BF_validation_path)
                    if os.path.isdir(BF_lumping_path) == False:
                        os.makedirs(BF_lumping_path)

                    shutil.copy(file_tocopy, os.path.join(
                        BF_validation_path, str(P)+'atm.txt'))
                    shutil.copy(file_tocopy, os.path.join(
                        BF_lumping_path, str(P)+'atm.txt'))

                else:
                    # continue with the while loop
                    BR_L_REAC = BF_OUTPUT.iloc[:, 1:]
                    it += 1

            else:
                # ALL OTHER SIMULATION TYPES: ONLY 1 CYCLE REQUIRED; also write products BF
                # Write branchings to the lumped products
                if len(PRODS) != len(PRODSLUMPED):
                    print(
                        'writing branchings of lumped products at {} atm ...'.format(str(P)))
                    postproc.WRITE_BRANCHINGS_PRODS(PRODS)

                    # FOR PRESCREENING_ALLREACTIVE: IF YOU HAVE LUMPED PRODUCTS, COPY THE BF TO THE INPUT OF COMPOSITION SELECTION
                    if jobtype == 'prescreening_allreactive':
                        for PRi_L in PRODS_L:
                            BF_comp_sel = os.path.join(
                                cwd, 'composition_selection', 'BF_INPUT', PRODS_L)
                            file_tocopy = os.path.join(
                                cwd, 'prescreening_allreactive', 'BF_OUTPUT', REAC_L, PRODS_L, str(P)+'atm.txt')

                            if os.path.isdir(BF_comp_sel) == False:
                                os.makedirs(BF_comp_sel)

                            shutil.copy(file_tocopy, os.path.join(
                                BF_comp_sel, str(P)+'atm.txt'))

                # FOR PRESCREENING_EQUILIBRIUM : COPY BF TO composition_selection-prescreening_allreactive/BF_INPUT folder?
                if jobtype == 'prescreening_equilibrium':
                    # COPY THE COMPOSITION TO validation and lumping directories
                    BF_prescreen_allreac = os.path.join(
                        cwd, 'prescreening_allreactive', 'BF_INPUT', REAC_L)
                    BF_comp_sel = os.path.join(
                        cwd, 'composition_selection', 'BF_INPUT', REAC_L)
                    file_tocopy = os.path.join(
                        cwd, 'prescreening_equilibrium', 'BF_OUTPUT', REAC_L, REAC_L, str(P)+'atm.txt')

                    if os.path.isdir(BF_prescreen_allreac) == False:
                        os.makedirs(BF_prescreen_allreac)
                    if os.path.isdir(BF_comp_sel) == False:
                        os.makedirs(BF_comp_sel)

                    shutil.copy(file_tocopy, os.path.join(
                        BF_prescreen_allreac, str(P)+'atm.txt'))
                    shutil.copy(file_tocopy, os.path.join(
                        BF_comp_sel, str(P)+'atm.txt'))

                flag_BFiter = 1

        tocP = clock()

        Dt_Pi[Pi, 1] = tocP-ticP
        # update Pi
        Pi += 1
    # write the pathways to the experimental datasets, the final mechanism, and the thermodynamic file
    print('Writing the list of pathways to experimental data ...'), postproc.WRITE_FINAL_PATHS()
    if PRODSINKS == 1:
        print('Writing Arrhenius fits in PLOG form ...'), kfit_P.WRITE_PLOG_FITS(
            arrfit_P, P_VECT, fld)
        # write the thermodynamic file
        preproc.WRITE_THERM(fld, STOICH, SPECIES_SERIES_L,
                            SPECIES_BIMOL_SERIES_L)

    ############### IN CASE YOU HAVE PLOT_CMP : PERFORM SIMULATIONS WITH THE LUMPED or OPTIMIZED MECHANISM #############
    profiles_P_all = {'detailed': profiles_P}
    i_reac_all = {'detailed': i_REAC_L}
    # derive the profiles
    if PLOT_CMP == 'YES' and jobtype == 'lumping':
        plt_cmp = prof_CKImech.PROFILES_FROM_CKI(cwd, fld, OS_folder)
        try:
            plt_cmp.COMPILE_MECH()
        except ValueError as e:
            print('error while compiling mech: ' + str(e))
        profiles_P_guess = plt_cmp.DERIVE_PROFILES(
            fld, P_VECT, T_VECT, SPECIES_L, PRODS_L, i_REAC_L, N_INIT_REAC, SPECIES_BIMOL_SERIES_L, ISOM_EQUIL, CUTOFF)
        profiles_P_all.update({'lumped_guess': profiles_P_guess})

    # derive the profiles with the optimized mech, if indicated
    # remove the variable "OPT_MECH" and replace it with checks on the simulation type
    if PLOT_CMP == 'YES' and jobtype == 'validation':
        OPT_MECH = os.path.join(cwd, 'lumpedmech')
        plt_cmp = prof_CKImech.PROFILES_FROM_CKI(cwd, OPT_MECH, OS_folder)
        try:
            plt_cmp.COMPILE_MECH()
        except ValueError as e:
            print('error while compiling mech: ' + str(e))
        profiles_P_opt = plt_cmp.DERIVE_PROFILES(
            fld, P_VECT, T_VECT, SPECIES_L, PRODS_L, i_REAC_L, N_INIT_REAC, SPECIES_BIMOL_SERIES_L, ISOM_EQUIL, CUTOFF)
        profiles_P_all.update({'lumped_opt': profiles_P_opt})

    # if you find the optimized lumped mechanism: plot also with that
    if os.path.isdir(os.path.join(cwd, 'lumpedmech_opt')):
        OPT_MECH2 = os.path.join(cwd, 'lumpedmech_opt')
        plt_cmp2 = prof_CKImech.PROFILES_FROM_CKI(cwd, OPT_MECH2, OS_folder)
        plt_cmp2.COMPILE_MECH()
        profiles_P_opt2 = plt_cmp2.DERIVE_PROFILES(
            fld, P_VECT, T_VECT, SPECIES_L, PRODS_L, i_REAC_L, N_INIT_REAC, SPECIES_BIMOL_SERIES_L, ISOM_EQUIL, CUTOFF)
        profiles_P_all.update({'opt2': profiles_P_opt2})

    ###################### PLOTS #######################################################
    print('Plotting ...')
    # generate the class and the folder "PLOTS"
    plots_callclass = mine_plt.PLOTTING(os.path.join(cwd, jobtype))
    Pi = 0
    profiles_Pi = dict.fromkeys(list(profiles_P_all.keys()))

    for P in P_VECT:
        tic = clock()
        # set the data to plot
        for key in list(profiles_P_all.keys()):
            profiles_Pi[key] = profiles_P_all[key][P]
        # generate the common frame for all the plots
        plots_callclass.plot_data(profiles_Pi, REAC_L, PRODS_L, SPECIES_L, P)
        # generate the different figures and save them
        plots_callclass.plot_exp_profiles(N_INIT_REAC)
        # for lumped reactant: plot also the reactant composition
        if isinstance(REAC, np.ndarray):
            plots_callclass.plot_data_reac(
                profiles_P_reac[P], REAC, P, N_INIT_REAC)
        toc = clock()
        Dt_Pi[Pi, 2] = toc-tic
        Pi += 1

    ##################### SAVE THE TIMES REQUIRED TO PERFORM EVERY STEP ##################

    # create series with the times and the description and print it to excel
    times_Pi = pd.DataFrame(Dt_Pi, index=list(P_VECT), columns=Dt_names_Pi)
    # times_Pi.to_excel('clock_Pi.xlsx',sheet_name='times_Pi')
