import shutil
import os
import D_ODESYSTEM as odesys
import numpy as np

class PROFILES_FROM_CKI:
    def __init__(self,cwd,path,OS_folder):
        self.cwd = cwd      # execution folder
        self.path = path    # path where the CKI mech is
        self.OS_folder = OS_folder # folder containing OS executables
    # derive the profiles in the selected range of T,P with the optimized mech
    # 1 compile the mechanism
    def COMPILE_MECH(self):
        '''
        Compile the mechanism of the "path" folder after copying it into the /mech_tocompile folder
        '''
        # if kin and therm files exist in the selected path: remove them in the destination folder, so as to overwrite them
        if os.path.isfile(self.path + '/kin.txt'):
            os.remove(self.cwd +'/mech_tocompile/kin.txt')
            shutil.copyfile(self.path + '/kin.txt',self.cwd + '/mech_tocompile/kin.txt')
        else:
            raise ValueError('kin.txt not found: file not substituted in the mech_tocompile folder')

        if os.path.isfile(self.path + '/therm.txt'):
            os.remove(self.cwd + '/mech_tocompile/therm.txt')
            shutil.copyfile(self.path + '/therm.txt',self.cwd + '/mech_tocompile/therm.txt')
        
        # compile the optimized mechanism
        #toexecute = r'"C:\Users\Luna Pratali Maffei\OpenSMOKE++Suite\bin\OpenSMOKEpp_CHEMKIN_PreProcessor.exe" --input .\mech_tocompile\input_preproc.dic >.\mech_tocompile\preproc_output.txt'
        #toexecute = r'"%OPENSMOKEPP_EXE_FOLDER%\OpenSMOKEpp_CHEMKIN_PreProcessor.exe" --input .\mech_tocompile\input_preproc.dic >.\mech_tocompile\preproc_output.txt'
        toexecute = '"' + self.OS_folder + "\OpenSMOKEpp_CHEMKIN_PreProcessor.exe" + '"' + " --input .\mech_tocompile\input_preproc.dic >.\mech_tocompile\preproc_output.txt"
    
        #print(toexecute)
        print('compiling mech ...'),os.system(toexecute)

    # 2 derive again all the profiles with the compiled mechanism
    def DERIVE_PROFILES(self,path_exp,P_VECT,T_VECT,SPECIES,PRODS,i_REAC,N_INIT_REAC,SPECIES_BIMOL_SERIES,ISOM_EQUIL):
        '''
        Derive the output profiles of the OS simulations with the appropriate compiled mechanism
        INPUT:
            path_exp: path to the experimental profiles, where you find the input_OS.dic to use for each simulation
            P_VECT,T_VECT: vectors with the temperatures and pressures of interest
        OUTPUT: profiles_P => dictionary with all the pressure; in each dictionary, another dictionary with the temperatures of the simulations is contained.
                in each T dictionary, a DataFrame with the profiles of the selected output species is contained
        '''
        # SPECIES NAMES: REAC_PRODS AS WRITTEN IN OPENSMOKE
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
                inputpath = path_exp + '/' + str(P) + 'atm/' + str(T) + 'K/input_OS.dic'
                if os.path.isfile(inputpath):
                    os.remove(self.cwd + '/input_OS.dic')
                    shutil.copyfile(inputpath,self.cwd + '/input_OS.dic')

                    ################## SOLUTION OF THE ODE SYSTEM ##########################################
                    # CALL OPENSMOKE
                    #toexecute = r'"C:\Users\Luna Pratali Maffei\OpenSMOKE++Suite\bin\OpenSMOKEpp_BatchReactor.exe" --input input_OS.dic > OS_output.txt'
                    #toexecute = r'"%OPENSMOKEPP_EXE_FOLDER%\OpenSMOKEpp_BatchReactor.exe" --input input_OS.dic > OS_output.txt'
                    toexecute = '"' + self.OS_folder + "\OpenSMOKEpp_BatchReactor.exe" + '"' + " --input input_OS.dic > OS_output.txt"
                    os.system(toexecute)


                    ##################### EXTRACT THE OUTPUT ###############################################
                    try:
                        tW_DF,PV = postproc.EXTRACT_PROFILES(SPECIES_PR,0,N_INIT_REAC,SPECIES_BIMOL_SERIES_NEW,ISOM_EQUIL)
                    except ValueError as e:
                        print(str(e))
                    # SAVE PROFILES IN A DICTIONARY FOR LATER POSTPROCESSING AND PLOTTING
                    profiles_P[P][T] = tW_DF

            # update Pi
            Pi += 1

        return profiles_P

