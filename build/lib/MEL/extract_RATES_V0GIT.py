# import modules of python
import os
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

#subprocess.call('Run.bat',stdout=subprocess.DEVNULL)

# import custom modules
from . import A_read_input as readinp
from . import B_extract_rates as extr
from . import C_preprocessing as preproc 
from . import D_ODESYSTEM as odesys # I call the function without the need of pre-pending anything
from . import E_PLOTTING as mine_plt
from . import F_FITS as fitall
from . import G_PROFILES as prof_CKImech
from time import perf_counter as clock

def set_simul_loop(cwd,jobtype,job_subdict,mech_dict):
    '''
    Set the list of operations to perform for the selected jobtype
    Output: dataframe
    rows: N of simulation
    columns: main subfolder where output is stored; reactant; lumpedreactant; product; lumpedproduct
    '''
    # read species list from mech_dict
    SPECIES = mech_dict['SPECIES']
    SPECIES_BIMOL = mech_dict['SPECIES_BIMOL']

    # initialize the array of simulations
    simul_array = np.array([])
    # SINGLE SIMULATION
    if 'simul_type' in job_subdict.keys():
        N_simul = 1
        # for single simulations: for lumping and validation, reactants and products are already indicated
        if jobtype == 'lumping' or jobtype == 'validation':
            # look for reactants and products in the subdictionaries
            REAC = job_subdict['REAC']
            REACLUMPED = job_subdict['REACLUMPED']
            PROD = job_subdict['PROD']
            PRODSLUMPED = job_subdict['PRODSLUMPED']
            # set output folder
            fld = os.path.join(cwd,jobtype,REACLUMPED)
            simul_array = np.append(simul_array,[fld,REAC,REACLUMPED,PROD,PRODSLUMPED])

        elif jobtype == 'prescreening_equilibrium' or jobtype == 'prescreening_allreactive' or jobtype == 'composition_selection':

            # select the reactant
            if isinstance(job_subdict['pseudospecies'].iloc[0],np.ndarray):
                REAC = job_subdict['pseudospecies'].iloc[0]
                REACLUMPED = job_subdict['pseudospecies'] # series

            elif jobtype == 'prescreening_equilibrium' and job_subdict['pseudospecies'].iloc[0] == 'all':
                # only accepted case where you consider all pseudospecies together: all isomers as one
                REAC = SPECIES # all species together
                REACLUMPED = pd.Series([SPECIES],index=[SPECIES[0]+'_L'])

            # PRODUCTS: irrelevant to the reactivity; set them as the "other" species
            PRODSUMPED = pd.Series(SPECIES,index=[SPECIES])
            PRODSLUMPED = PRODSLUMPED.drop(index=REAC)
            PRODS = PRODSLUMPED.values

    else:             
    ############ ALL THE OTHER SIMULATIONS
    # PRESCREENING_EQUILIBRIUM
    '''
        if jobtype == 'prescreening_equilibrium':


            for SP in subdict['pseudospecies'].index :
                # CHECK THAT THE PSEUDOSPECIES ARE PRESENT IN THE LIST
                if isinstance(subdict['pseudospecies'].loc[SP],np.ndarray):
                    if np.array([[sp == SPECIES for sp in subdict['pseudospecies'].loc[SP]][ii].any() for ii in range(0,len(subdict['pseudospecies'].loc[SP]))]).all() != True:
                        error_list = error_list + '\nNot all species indicated in the set of pseudospecies match the species list. \n \t please select among {SPECIES} '.format(SPECIES = str(SPECIES))
                    # WITHIN EACH GROUP: CHECK THAT PSEUDOSPECIES BELONG TO THE SAME TYPE (I.E. HAVE SAME BIMOL FRAGMENT)
                    i_sp_group = np.where([sp_group == SPECIES for sp_group in subdict['pseudospecies'].loc[SP]])[1]
                    if np.array([i_sp_bim == SPECIES_BIMOL[i_sp_group] for i_sp_bim in SPECIES_BIMOL[i_sp_group]]).all == False:
                        error_list = error_list + '\nPseudospecies do not belong to the same type: mix of unimol/bimol or bimol with different fragments '

                elif subdict['pseudospecies'].loc[SP] == 'all':
                    # all pseudospecies together
                    species = subdict['pseudospecies'].loc[SP]
                    if np.array([species == SP_list for SP_list in SPECIES]).any() != True:
                        error_list = error_list + '\nNot all species indicated in the set of pseudospecies match the species list. \n \t please select among {SPECIES} '.format(SPECIES = str(SPECIES))
                 REAC = SPECIES # all species together
                REACLUMPED = pd.Series([SPECIES],index=[SPECIES[0]+'_L'])

            # PRODUCTS: irrelevant to the reactivity; set them as the "other" species
            PRODSUMPED = pd.Series(SPECIES,index=[SPECIES])
            PRODSLUMPED = PRODSLUMPED.drop(index=REAC)
            PRODS = PRODSLUMPED.values 
    '''              
    # PRESCREENING_ALLREACTIVE

    # COMPOSITION_SELECTION

    # LUMPING

    # VALIDATION

    output_DF = pd.DataFrame(simul_array,index=np.arange(0,N_simul),columns=['fld','REAC','REACLUMPED','PROD','PRODLUMPED'])
    return output_DF


def main_simul(cwd,jobtype,input_par,input_par_jobtype,mech_dict):
    '''
    Perform simulations for 1 reacting pseudospecies at all T,P provided
    '''

    # Save times as dataframes for every pressure; the first will be a dictionary
    Dt_names_Pi = ['ode solving','time Pi','plotting and saving figs']

    ############### DERIVE VARIABLES REQUIRED FOR CODE FLOW #######################################

    OS_folder = input_par['opensmoke_folder']
    input_type = input_par['mech_type']
    P_VECT = input_par['P_vect']
    T_VECT = input_par['T_vect']
    T_VECT_SKIP = input_par['T_skip']
    UNITS_BIMOL = input_par['units_bimol']
    STOICH = input_par['stoich']

    ISOM_EQUIL = input_par_jobtype['isom_equil']
    PRODSINKS = input_par_jobtype['Prods_sinks']
    PLOT_CMP = input_par_jobtype['isom_equil']    


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

    REAC
    REACLUMPED
    PRODS
    PRODSLUMPED


    # PROCESSING OF SPECIES NAMES AND REACTANT:
    print(SPECIES,SPECIES_BIMOL)
    SPECIES_SERIES = pd.Series(np.arange(0,len(SPECIES)),index = SPECIES)
    SPECIES_BIMOL_SERIES = pd.Series(SPECIES_BIMOL,index = SPECIES) #species values are the names of the bimol species corresponding to the "primary" one

    i_PRODS = SPECIES_SERIES[PRODS].values  # array
    if isinstance(REAC,np.ndarray):
        i_REAC = SPECIES_SERIES[REAC].values          # value (single reactant) or array (lumped reactant)
    else:
        i_REAC = SPECIES_SERIES[REAC]

    # PREPROCESSING: WRITE THE THERMODYNAMIC FILE
    print('writing therm.txt for OS preprocessor ...'),preproc.WRITE_THERM(cwd + '/mech_tocompile',STOICH,SPECIES_SERIES,SPECIES_BIMOL_SERIES)

    # ASSIGN THE INITIAL NUMBER OF MOLES
    N_INIT_REAC = preproc.INITIALMOLES(REAC,SPECIES_BIMOL_SERIES,1)

    # reduce T range if necessary
    # skip the temperature selected in the input
    for TSKIP in T_VECT_SKIP:
        mask = np.where(TSKIP != T_VECT)
        T_VECT = T_VECT[mask]
    ######################## PREPROCESSING OF CKI MECHANISM #########################
    if input_type == 'CKI':
        # 
        case.copy_CKI_processed(cwd + '/mech_tocompile',PRODSINKS,ISOM_EQUIL,REAC,PRODS)
        toexecute = '"' + OS_folder + "\OpenSMOKEpp_CHEMKIN_PreProcessor.exe" + '"' + " --input .\mech_tocompile\input_preproc.dic >.\mech_tocompile\preproc_output.txt"
        #toexecute = r'"C:\Users\Luna Pratali Maffei\OpenSMOKE++Suite\bin\OpenSMOKEpp_CHEMKIN_PreProcessor.exe" --input .\mech_tocompile\input_preproc.dic >.\mech_tocompile\preproc_output.txt'
        #toexecute = r'"%OPENSMOKEPP_EXE_FOLDER%\OpenSMOKEpp_CHEMKIN_PreProcessor.exe" --input .\mech_tocompile\input_preproc.dic >.\mech_tocompile\preproc_output.txt'
        print(OS_folder)
        print(toexecute)
        print('compiling mech ...'),os.system(toexecute)

    ######################## LOOP OVER THE SELECTED PRESSURES ########################
    # call the class of the output processing
    postproc = odesys.ODE_POSTPROC(cwd)

    # PREALLOCATIONS
    Dt_Pi = np.zeros((len(P_VECT),len(Dt_names_Pi)))
    Pi = 0
    # dictionary for each pressure  
    profiles_P = dict.fromkeys(list(P_VECT))
    arrfit_P = dict.fromkeys(list(P_VECT))
    # reactant profiles for lumped reactant
    if isinstance(REAC,np.ndarray):
        profiles_P_reac = dict.fromkeys(list(P_VECT))

    T_VECT_0=T_VECT
    for P in P_VECT:

        
        print('processing: P = ' + str(P) + ' atm ... ')
        ticP = clock()

        ###### CHECK REACTIVITY UNTIL THE SELECTED TEMPERATURE
        # extract the rates in the selected pressure
        if input_type == 'MESS':
            # SET T_VECT=T_VECT_0 TO MAKE TEMPERATURE REDUCTION INDEPENDENT OF THE PREVIOUS PRESSURE
            T_VECT=T_VECT_0
            if isinstance(REAC,np.ndarray):
                # IN CASE OF A LUMPED REACTANT: DERIVE A TOTAL RATE CONSTANT WITH A LOOP
                kR_j = np.zeros((len(T_VECT_MESS),len(SPECIES)-1))
                for rr in REAC:
                    kR_j = kR_j + case.REAC_P(P,rr)

            elif isinstance(REAC,str):
                # NON LUMPED REACTANT: DO A SINGLE DERIVATION OF THE REACTIVITY MATRIX
                kR_j = case.REAC_P(P,REAC)

            #print(kR_j)
            flag_exit,Err,T_VECT = preproc.CHECK_REACTIVITY(kR_j,T_VECT,T_VECT_MESS)
            print(Err)
            if flag_exit == 1:
                exit()

        ############### PRE PROCESSING FOR LUMPED REACTANTS AND PRODUCTS   ##################à

        # FOR LUMPED REACTANT: EXTRACT IN A DATAFRAME THE COMPOSITION OF THE REACTANT AT EVERY PRESSURE
        # COMPARE THE BRANCHINGS WITH THE T RANGE AND PRINT A WARNING IF YOU REDUCE IT
        if isinstance(REAC,np.ndarray):
            try:
                BR_L_REAC,T_VECT = preproc.BRANCHING_LUMPEDREAC(cwd,REACLUMPED.index[0],REAC,T_VECT,P,STOICH,ISOM_EQUIL)
            except ValueError as e:
                print('error while extracting the branching fractions of ' + REACLUMPED.index[0] + ': ' + str(e))
                exit()

        elif isinstance(REAC,str):
            # create empty variable (so that the input will exist anyways)
            BR_L_REAC = ''

        # FOR LUMPED PRODUCTS or reactants: GENERATE DATAFRAMES FOR BRANCHINGS
        if len(PRODS) != len(PRODSLUMPED) or len(REAC) != len(REACLUMPED):
            try:
                postproc.MAKE_BRANCHINGFOLDERS(REACLUMPED,PRODSLUMPED,T_VECT,STOICH)
            except ValueError as e:
                print(e)
                
        # inside the profile_P dictionary: allocate a profiles_T dictionary to include the profiles at every temperature
        profiles_P[P] = dict.fromkeys(list(T_VECT))
        if isinstance(REAC,np.ndarray):
            profiles_P_reac[P] = dict.fromkeys(list(T_VECT))

        ########################### LOOP OVER ALL THE TEMPERATURES OF INTEREST ####################

        Ti = 0

        for T in T_VECT:
            print('processing: T = ' + str(T) + ' K ... ')
            Ti += 1
            # make the folder of the output to optiSMOKE++
            postproc.MAKE_FOLDERS(STOICH,P,T,REACLUMPED)

            
            ################# PREPROCESSING - MESS INPUT ONLY #######################################################
            if input_type == 'MESS':
                # EXTRACT THE KINETIC MATRIX
                k_ij_TP = case.MATRIX_TP(T,P)
                
                # PREPROCESSING: SET THE PRODUCTS AS IRREVERSIBLE SINKS
                k_ij_TP_proc = preproc.PREPROCESSING(k_ij_TP,SPECIES_SERIES,PRODS)
                k_ij_TP_prodsinks = k_ij_TP_proc.SET_PRODS_SINKS(PRODSINKS)
                k_ij_TP_OK = k_ij_TP_proc.MULTIPLY_BIMOL(SPECIES_BIMOL,UNITS_BIMOL)

                ################# WRITE FICTITIOUS MECH IN CKI FORMAT ##################################
                
                k_to_CKI = preproc.WRITE_MECH_CKI(k_ij_TP_OK,np.zeros(k_ij_TP_OK.shape),np.zeros(k_ij_TP_OK.shape),np.zeros(k_ij_TP_OK.shape,dtype=str),SPECIES_SERIES,i_REAC,i_PRODS,SPECIES_BIMOL)
                
                print('writing kin.txt for OS preprocessor ...')
                # ISOM_EQUIL is put here: if you want to study the isomer equilibrium, all the other reactions are not copied
                CKI_lines = k_to_CKI.MAKE_CKI(PRODSINKS,ISOM_EQUIL)
                try:
                    k_to_CKI.WRITE_CKI(cwd + '/mech_tocompile',CKI_lines.values)
                except RuntimeError as e:
                    print(str(e))
                    exit()
                
                #toexecute = r'"C:\Users\Luna Pratali Maffei\OpenSMOKE++Suite\bin\OpenSMOKEpp_CHEMKIN_PreProcessor.exe" --input .\mech_tocompile\input_preproc.dic >.\mech_tocompile\preproc_output.txt'
                toexecute = '"' + OS_folder + "\OpenSMOKEpp_CHEMKIN_PreProcessor.exe" + '"' + " --input .\mech_tocompile\input_preproc.dic >.\mech_tocompile\preproc_output.txt"
                #toexecute = r'"%OPENSMOKEPP_EXE_FOLDER%\OpenSMOKEpp_CHEMKIN_PreProcessor.exe" --input .\mech_tocompile\input_preproc.dic >.\mech_tocompile\preproc_output.txt'
                print(OS_folder)
                print(toexecute)
                print('compiling mech ...'),os.system(toexecute)
                #subprocess.run([r'"%OPENSMOKEPP_EXE_FOLDER%\OpenSMOKEpp_CHEMKIN_PreProcessor.exe"'],input=(cwd + '/mech_tocompile/input_preproc.dic'),stdout=(cwd+ '/mech_tocompile/preproc_output.txt'))
                #subprocess.run([r'%OPENSMOKEPP_EXE_FOLDER%\OpenSMOKEpp_CHEMKIN_PreProcessor.exe --input' + cwd + r'\mech_tocompile\input_preproc.dic' +  '>' + cwd + r'\mech_tocompile\preproc_output.txt'])
            
            ################### WRITE OPENSMOKE INPUT ##################################################
            try:
                OS_write = preproc.WRITE_OS_INPUT(cwd,T,P,SPECIES,REACLUMPED,1,SPECIES_BIMOL_SERIES,BR_L_REAC)
            except RuntimeError as e:
                print(str(e))
                exit()
            print('writing new input ...'),OS_write.w_input_OS()
            
            ################## SOLUTION OF THE ODE SYSTEM ##########################################
            # CALL OPENSMOKE
            tic = clock()
            #toexecute = r'"C:\Users\Luna Pratali Maffei\OpenSMOKE++Suite\bin\OpenSMOKEpp_BatchReactor.exe" --input input_OS.dic > OS_output.txt'
            toexecute = '"' + OS_folder + "\OpenSMOKEpp_BatchReactor.exe" + '"' + " --input input_OS.dic > OS_output.txt"
            # toexecute = r'"%OPENSMOKEPP_EXE_FOLDER%\OpenSMOKEpp_BatchReactor.exe" --input input_OS.dic > OS_output.txt'
            print('solving OS Batch Reactor ...'),os.system(toexecute)
            toc = clock()
            Dt_Pi[Pi,0]= toc-tic

            ##################### EXTRACT THE OUTPUT ###############################################
            try:
                tW_DF,PV = postproc.EXTRACT_PROFILES(SPECIES,i_REAC,N_INIT_REAC,SPECIES_BIMOL_SERIES,ISOM_EQUIL)
            except ValueError as e:
                print(str(e))
            print(tW_DF)

            # process the output: rewrite tW_DF if there are lumped species
            print('Rewriting profiles for lumped species ...')
            tW_DF,REAC_L,i_REAC_L,SPECIES_L,SPECIES_SERIES_L,SPECIES_BIMOL_SERIES_L,PRODS_L = postproc.LUMP_PROFILES(PRODS,PRODSLUMPED)
            print(tW_DF)

            # SAVE PROFILES IN A DICTIONARY FOR LATER POSTPROCESSING AND PLOTTING
            profiles_P[P][T] = tW_DF
            if isinstance(REAC,np.ndarray):
                tW_DF_reac = postproc.PROFILES_REAC_COMPOSITION()
                profiles_P_reac[P][T] = tW_DF_reac
            # IF YOU HAVE LUMPING (I.E.: PRODSINKS=1)
            # LINEAR FIT OF THE PROFILES OBTAINED: GUESS FOR SUCCESSIVE optiSMOKE++
            if PRODSINKS == 1:
                # if it is the first temperature: call the class
                # generate dataframe at every T,P to allocate the fits of k_ij
                if Ti == 1:
                    kfit_P = fitall.FITTING(T_VECT,REAC_L,PRODS_L)
                # then do the fits
                k_prods_T,fit_error_T = kfit_P.fit_profiles(tW_DF,i_REAC_L,SPECIES_SERIES_L,T,PV,SPECIES_BIMOL_SERIES_L,1)
            # I think I need to do something with this fit_error_T, save it somewhere
            # maybe I can do a 3D plot with the fitting error at a certain T,P or with the difference between the sum(ki) and k of the reactant
            # and finally the plot with the branching among different channels with the lumped rates

            # WRITE OUTPUT PROFILES TO OPTISMOKE
            print('Writing profiles for Optismoke ...'),postproc.WRITE_PROFILES(PRODS_L)

            # GENERATE THE NEW OS INPUT AND WRITE IT TO THE FOLDERS
            print('Generate OS input for the lumped mech and copy to the folder'),postproc.WRITE_NEW_OSINPUT(1)

            ############################### GENERATE RANDOM PROFILES AND SAVE THEM IN OUTPUT_TO_OPTISMOKE_RANDOM

        tocP = clock()

        Dt_Pi[Pi,1]=tocP-ticP

        # Return the dataframes of the fits of the rate constants
        if PRODSINKS == 1:
            print('Fitting the arrhenius profiles ...')
            # write the original rate constants in the corresponding folder
            out_fld = cwd+'/Outputs_to_optiSMOKE_' + ''.join(STOICH) +  '/' + REAC_L + '/' + str(P) + 'atm'
            kfit_P.write_originalk(out_fld)
            # fit
            rates_P_CKI = kfit_P.fits_lumped_k(cwd,P,SPECIES_BIMOL_SERIES_L)
            arrfit_P[P] = rates_P_CKI

        # Write branchings to the lumped products
        if len(PRODS) != len(PRODSLUMPED) or len(REAC) != len(REACLUMPED):
            print('writing branchings of lumped reactant/products at the given pressure ...')
            postproc.WRITE_BRANCHINGS(PRODS)
        # update Pi
        Pi += 1
    # write the pathways to the experimental datasets, the final mechanism, and the thermodynamic file
    print('Writing the list of pathways to experimental data ...'),postproc.WRITE_FINAL_PATHS()
    if PRODSINKS == 1:
        out_fld = cwd+'/Outputs_to_optiSMOKE_' + ''.join(STOICH) +  '/' + REAC_L
        print('Writing Arrhenius fits in PLOG form ...'),kfit_P.WRITE_PLOG_FITS(arrfit_P,P_VECT,out_fld)
        #write the thermodynamic file
        preproc.WRITE_THERM(out_fld,STOICH,SPECIES_SERIES_L,SPECIES_BIMOL_SERIES_L)


    ############### IN CASE YOU HAVE PLOT_CMP : PERFORM SIMULATIONS WITH THE OPTIMIZED MECHANISM #############
    profiles_P_all = {'detailed':profiles_P}
    i_reac_all = {'detailed':i_REAC_L}
    # derive the profiles 
    if PLOT_CMP == 'YES':
        path_guessmech = cwd+'/Outputs_to_optiSMOKE_' + ''.join(STOICH) +  '/' + REAC_L
        plt_cmp = prof_CKImech.PROFILES_FROM_CKI(cwd,path_guessmech,OS_folder)
        try:
            plt_cmp.COMPILE_MECH()
        except ValueError as e:
            print('error while compiling mech: ' + str(e))
        profiles_P_guess = plt_cmp.DERIVE_PROFILES(path_guessmech,P_VECT,T_VECT,SPECIES_L,PRODS_L,i_REAC_L,N_INIT_REAC,SPECIES_BIMOL_SERIES_L,ISOM_EQUIL)
        profiles_P_all.update({'lumped_guess':profiles_P_guess})
        i_reac_all.update({'lumped_guess':0})

    # derive the profiles with the optimized mech, if indicated
    # remove the variable "OPT_MECH" and replace it with checks on the simulation type
    if OPT_MECH and PLOT_CMP == 'YES':
        path_exp = cwd+'/Outputs_to_optiSMOKE_' + ''.join(STOICH) +  '/' + REAC_L
        plt_cmp = prof_CKImech.PROFILES_FROM_CKI(cwd,OPT_MECH,OS_folder)
        try:
            plt_cmp.COMPILE_MECH()
        except ValueError as e:
            print('error while compiling mech: ' + str(e))
        profiles_P_opt = plt_cmp.DERIVE_PROFILES(path_exp,P_VECT,T_VECT,SPECIES_L,PRODS_L,i_REAC_L,N_INIT_REAC,SPECIES_BIMOL_SERIES_L,ISOM_EQUIL)
        profiles_P_all.update({'lumped_opt':profiles_P_opt})
        i_reac_all.update({'lumped_opt':0})

    ###################### PLOTS #######################################################
    print('Plotting ...')
    plots_callclass = mine_plt.PLOTTING(cwd,STOICH) # generate the class and the folder "PLOTS"
    Pi = 0
    profiles_Pi = dict.fromkeys(list(profiles_P_all.keys()))

    for P in P_VECT:
        tic = clock()
        # set the data to plot
        for key in list(profiles_P_all.keys()):
            profiles_Pi[key] = profiles_P_all[key][P]
        # generate the common frame for all the plots
        plots_callclass.plot_data(profiles_Pi,REAC_L,PRODS_L,SPECIES_L,P)
        # generate the different figures and save them
        plots_callclass.plot_exp_profiles(N_INIT_REAC)
        # for lumped reactant: plot also the reactant composition
        if isinstance(REAC,np.ndarray):
            plots_callclass.plot_data_reac(profiles_P_reac[P],REAC,P,N_INIT_REAC)
        toc = clock()
        Dt_Pi[Pi,2]=toc-tic
        Pi += 1

    ##################### SAVE THE TIMES REQUIRED TO PERFORM EVERY STEP ##################

    # create series with the times and the description and print it to excel
    times_Pi = pd.DataFrame(Dt_Pi,index=list(P_VECT),columns=Dt_names_Pi)
    times_Pi.to_excel('clock_Pi.xlsx',sheet_name='times_Pi')