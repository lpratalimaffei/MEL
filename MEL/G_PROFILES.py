import shutil
import os
import subprocess
from . import main_flow
from . import D_ODESYSTEM as odesys
import numpy as np

class PROFILES_FROM_CKI:
    def __init__(self,cwd,path,OS_folder):
        self.cwd = cwd      # execution folder
        self.path = path    # path where the kin.txt mech is
        self.OS_folder = OS_folder # folder containing OS executables
        self.OS_exe = main_flow.get_OS()
        self.exec0 = main_flow.get_libpath()
        self.preproc_exe = os.path.join('"' + OS_folder, "OpenSMOKEpp_CHEMKIN_PreProcessor." + self.OS_exe + '"')
        self.osbatch_exe = os.path.join('"' + OS_folder, "OpenSMOKEpp_BatchReactor." + self.OS_exe + '"')
        self.input_preproc = os.path.join(os.path.join(".", "mech_tocompile", "input_preproc.dic"))    
        self.output_preproc = os.path.join(".", "mech_tocompile", "preproc_output.txt")

    # derive the profiles in the selected range of T,P with the optimized mech
    # 1 compile the mechanism
    def COMPILE_MECH(self):
        '''
        Compile the mechanism of the "path" folder after copying it into the mech_tocompile folder
        '''
        # if kin and therm files exist in the selected path: remove them in the destination folder, so as to overwrite them
        if os.path.isfile(os.path.join(self.path, 'kin.txt')):
            os.remove(os.path.join(self.cwd,'mech_tocompile','kin.txt'))
            shutil.copyfile(os.path.join(self.path, 'kin.txt'), os.path.join(self.cwd, 'mech_tocompile', 'kin.txt'))
        else:
            raise ValueError('kin.txt not found: file not substituted in the mech_tocompile folder')

        if os.path.isfile(os.path.join(self.path, 'therm.txt')):
            os.remove(os.path.join(self.cwd, 'mech_tocompile', 'therm.txt'))
            shutil.copyfile(os.path.join(self.path, 'therm.txt'), os.path.join(self.cwd, 'mech_tocompile', 'therm.txt'))
        
        # compile the optimized mechanism
        toexecute = self.exec0 + self.preproc_exe + " --input " + self.input_preproc + ">" + self.output_preproc
        print('compiling mech ...'), subprocess.run(toexecute, shell=True)

    # 2 derive again all the profiles with the compiled mechanism
    def DERIVE_PROFILES(self,path_exp,P_VECT,T_VECT,SPECIES,PRODS,i_REAC,N_INIT_REAC,SPECIES_BIMOL_SERIES,ISOM_EQUIL,CUTOFF):
        '''
        Derive the output profiles of the OS simulations with the appropriate compiled mechanism
        INPUT:
            path_exp: path to the experimental profiles, where you find the input_OS.dic to use for each simulation
            P_VECT,T_VECT: vectors with the temperatures and pressures of interest
        OUTPUT: profiles_P => dictionary with all the pressure; in each dictionary, another dictionary with the temperatures of the simulations is contained.
                in each T dictionary, a DataFrame with the profiles of the selected output species is contained
        '''
        # SPECIES NAMES: REAC_PRODS AS WRITTEN IN OPENSMOKE
        PRODS = np.array(PRODS, dtype='<U16')
        SPECIES_PR = np.insert(PRODS,0,SPECIES[i_REAC])
        SPECIES_BIMOL_SERIES_NEW = SPECIES_BIMOL_SERIES[SPECIES_PR]
        # initialize postprocessor to extract profiles
        postproc = odesys.ODE_POSTPROC(self.cwd)
        # dictionary where to store values
        profiles_P = dict.fromkeys(list(P_VECT))
        # loops
        Pi = 0
        for P in P_VECT:
            # inside the profile_P dictionary: allocate a profiles_T dictionary to include the profiles at every temperature
            profiles_P[P] = dict.fromkeys(list(T_VECT))
            print('processing lumped mech: P = ' + str(P) + ' atm ... ')

            ########################### LOOP OVER ALL THE TEMPERATURES OF INTEREST ####################

            Ti = 0

            for T in T_VECT:
                Ti += 1
                ################### retrieve OPENSMOKE INPUT ##################################################
                
                # if the pathway exists: perform the simulation
                inputpath = os.path.join(path_exp,str(P)+'atm',str(T)+'K','input_OS.dic')
                if os.path.isfile(inputpath):
                    os.remove(os.path.join(self.cwd,'input_OS.dic'))
                    shutil.copyfile(inputpath,os.path.join(self.cwd,'input_OS.dic'))

                    ################## SOLUTION OF THE ODE SYSTEM ##########################################
                    # CALL OPENSMOKE
                    toexecute = self.exec0 + self.osbatch_exe + " --input input_OS.dic > OS_output.txt"
                    subprocess.run(toexecute, shell=True)


                    ##################### EXTRACT THE OUTPUT ###############################################
                    try:
                        tW_DF,PV = postproc.EXTRACT_PROFILES(SPECIES_PR,0,N_INIT_REAC,SPECIES_BIMOL_SERIES_NEW,ISOM_EQUIL,CUTOFF)
                    except ValueError as e:
                        print(str(e))
                    # SAVE PROFILES IN A DICTIONARY FOR LATER POSTPROCESSING AND PLOTTING
                    profiles_P[P][T] = tW_DF

            # update Pi
            Pi += 1

        return profiles_P

