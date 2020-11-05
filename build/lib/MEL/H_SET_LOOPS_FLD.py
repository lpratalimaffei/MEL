import sys
import os
import numpy as np
import pandas as pd
import shutil

def setfolder(fld):
    '''
    Create folder if it does not exist, and return 0
    If folder exists and is not empty: return 1
    '''
    fld_ex = int(os.path.isdir(fld))
    if fld_ex == 0:
        # create folder
        os.makedirs(fld)
        YE_NO = 0
    elif fld_ex == 1:
        # check if it has subfolders
        sub_fld = len(os.listdir(fld))   
        if sub_fld == 0:
            YE_NO = 0
        else:
            YE_NO = 1

    return YE_NO

def rmfolder(fld):
    '''
    remove folder if exists
    '''
    if os.path.isdir(fld):
        shutil.rmtree(fld)


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
    SPECIES_UNIMOL = SPECIES[np.where(''==SPECIES_BIMOL)]

    # SINGLE SIMULATION
    if 'simul_type' in job_subdict.keys():
        N_simul = 1
        # for single simulations: for lumping and validation, reactants and products are already indicated
        if jobtype == 'lumping' or jobtype == 'validation':
            # look for reactants and products in the subdictionaries
            REAC = job_subdict['REAC']
            REACLUMPED = job_subdict['REACLUMPED']
            PRODS = job_subdict['PRODS']
            PRODSLUMPED = job_subdict['PRODSLUMPED']
            # set output folder
            

        elif jobtype == 'prescreening_equilibrium' or jobtype == 'prescreening_allreactive' or jobtype == 'composition_selection':

            # select the reactant
            if isinstance(job_subdict['pseudospecies'].iloc[0],np.ndarray):
                REAC = job_subdict['pseudospecies'].iloc[0]
                REACLUMPED = job_subdict['pseudospecies'] # series

            elif jobtype == 'prescreening_equilibrium' and job_subdict['pseudospecies'].iloc[0] == 'all':
                # only accepted case where you consider all pseudospecies together: all isomers as one
                REAC = SPECIES_UNIMOL # all species together
                REACLUMPED = pd.Series([SPECIES_UNIMOL],index=[SPECIES_UNIMOL[0]+'_L'])

            # PRODUCTS: irrelevant to the reactivity; set them as the "other" species
            PRODSLUMPED = pd.Series(SPECIES,index=[SPECIES])
            PRODSLUMPED = PRODSLUMPED.drop(index=REAC)
            PRODS = PRODSLUMPED.values

        # set output folder and simul_array
        fld = os.path.join(cwd,jobtype,REACLUMPED.index[0])
        simul_array = np.array([fld,REAC,REACLUMPED,PRODS,PRODSLUMPED])
        simul_array = simul_array[np.newaxis,:]
        output_DF = pd.DataFrame(simul_array,index=np.arange(0,N_simul),columns=['fld','REAC','REACLUMPED','PRODS','PRODSLUMPED'])

    else:             
    ############ ALL THE OTHER SIMULATIONS
    # PRESCREENING_EQUILIBRIUM
        if jobtype == 'prescreening_equilibrium' or jobtype == 'prescreening_allreactive' :

            if job_subdict['pseudospecies'].index[0] == 'all':
                ######### equilibrium: all isomers together
                if jobtype == 'prescreening_equilibrium':
                    N_simul = 1
                    # all pseudospecies together
                    REAC = SPECIES_UNIMOL # all species together
                    REACLUMPED = pd.Series([SPECIES_UNIMOL],index=[SPECIES_UNIMOL[0]+'_L'])
                    # PRODUCTS: irrelevant to the reactivity; set them as the "other" species
                    PRODSLUMPED = pd.Series(SPECIES,index=[SPECIES])
                    PRODSLUMPED = PRODSLUMPED.drop(index=REAC)
                    PRODS = PRODSLUMPED.values 
                    # set output folder and simul_array
                    fld = os.path.join(cwd,jobtype,REACLUMPED.index[0])
                    simul_array = np.array([fld,REAC,REACLUMPED,PRODS,PRODSLUMPED]) 
                    simul_array = simul_array[np.newaxis,:] 
                    output_DF = pd.DataFrame(simul_array,index=np.arange(0,N_simul),columns=['fld','REAC','REACLUMPED','PRODS','PRODSLUMPED'])

                if jobtype == 'prescreening_allreactive' :
                    # create one set of simulations for each species
                    N_simul = len(SPECIES)
                    output_DF = pd.DataFrame(index=np.arange(0,N_simul),columns=['fld','REAC','REACLUMPED','PRODS','PRODSLUMPED'])
                    ind = 0
                    for SP in SPECIES:
                        REAC = SP
                        REACLUMPED = pd.Series([SP],index=[SP])
                        # PRODUCTS: ALL THW OTHERS
                        PRODSLUMPED = pd.Series(SPECIES,index=[SPECIES])
                        PRODSLUMPED = PRODSLUMPED.drop(index=REAC)
                        PRODS = PRODSLUMPED.values 
                        # set output folder and simul_array
                        fld = os.path.join(cwd,jobtype,REACLUMPED.index[0])
                        output_DF.loc[ind] = np.array([fld,REAC,REACLUMPED,PRODS,PRODSLUMPED])
                        ind += 1 

            else:
                # same for both prescreenings
                N_simul = job_subdict['pseudospecies'].index.size
                output_DF = pd.DataFrame(index=np.arange(0,N_simul),columns=['fld','REAC','REACLUMPED','PRODS','PRODSLUMPED'])
                ind = 0
                for SP in job_subdict['pseudospecies'].index :
                    # for each set of pseudospecies: generate dictionaries
                    if isinstance(job_subdict['pseudospecies'].loc[SP],np.ndarray):
                        # generate reactant arrays
                        REAC = job_subdict['pseudospecies'].loc[SP]
                        REACLUMPED = job_subdict['pseudospecies'] # series                   
                        # PRODUCTS: irrelevant to the reactivity; set them as the "other" species
                        PRODSLUMPED = pd.Series(SPECIES,index=[SPECIES])
                        PRODSLUMPED = PRODSLUMPED.drop(index=REAC)
                        PRODS = PRODSLUMPED.values 
                        # select output folder and allocate to dataframe
                        fld = os.path.join(cwd,jobtype,SP)
                        output_DF.loc[ind] = np.array([fld,REAC,REACLUMPED,PRODS,PRODSLUMPED])
                        ind += 1  

        elif jobtype == 'composition_selection' or 'lumping' or 'validation':
            # read the pseudospecies file
            pseudospecies_file = os.path.join(cwd,'inp','pseudospecies.txt')
            pseudospecies_array = np.genfromtxt(pseudospecies_file,delimiter='',dtype=str)
            # set an array of all pseudospecies
            stable_species = [] 
            # generate pseudospecies series
            pseudospecies_series = pd.Series(pseudospecies_array[:,1],index=pseudospecies_array[:,0])
            for SP in pseudospecies_series.index:
                if pseudospecies_series.loc[SP].split('+') != [pseudospecies_series.loc[SP]]:
                    for i in pseudospecies_series.loc[SP].split('+'):
                        stable_species.append(i)
                    pseudospecies_series.loc[SP] = np.array(pseudospecies_series.loc[SP].split('+'),dtype=str)
                else :
                    stable_species.append(pseudospecies_series.loc[SP])
                    # single species: composition_selection does not require them
                    if jobtype == 'composition_selection':
                         pseudospecies_series = pseudospecies_series.drop(index=[SP])
            stable_species = np.array(stable_species,dtype=str)

            # organize dataframe
            N_simul = pseudospecies_series.index.size
            output_DF = pd.DataFrame(index=np.arange(0,N_simul),columns=['fld','REAC','REACLUMPED','PRODS','PRODSLUMPED'])
            ind = 0
            
            # now do the loop on the right indices
            for SP in pseudospecies_series.index:
                # for each set: setup simulations
                REACLUMPED = pd.Series([pseudospecies_series.loc[SP]],index=[SP])
                REAC = pseudospecies_series.loc[SP]
                PRODSLUMPED = pseudospecies_series.drop(index=[SP])
                if isinstance(REAC,np.ndarray):
                    PRODS = np.delete(stable_species,np.where([rr==stable_species for rr in REAC])[1])
                else:
                    PRODS = np.delete(stable_species,np.where([ss==REAC for ss in stable_species])[0])
                    
                # select output folder and allocate to dataframe
                fld = os.path.join(cwd,jobtype,SP)
                output_DF.loc[ind] = np.array([fld,REAC,REACLUMPED,PRODS,PRODSLUMPED])
                ind += 1  


    return output_DF