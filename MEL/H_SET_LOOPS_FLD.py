import sys
import os
import numpy as np
import pandas as pd
import shutil
from . import A_read_input as readinp


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


def renamefiles(cwd, fld):
    '''
    Renames the files in the given folder with name_old
    '''
    files = os.listdir(fld)
    N = str(round(len(files)/2))
    for F in files:
        oldname = os.path.join(fld, F)
        newname = os.path.join(fld, F+'_old'+N)
        os.rename(oldname, newname)


def set_simul_loop(cwd, jobtype, job_subdict, mech_dict):
    '''
    Set the list of operations to perform for the selected jobtype
    Output: dataframe
    rows: N of simulation
    columns: main subfolder where output is stored; reactant; lumpedreactant; product; lumpedproduct
    '''
    # read species list from mech_dict
    SPECIES = mech_dict['SPECIES']
    SPECIES_BIMOL = mech_dict['SPECIES_BIMOL']
    SPECIES_UNIMOL = SPECIES[np.where('' == SPECIES_BIMOL)]

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
            if (isinstance(job_subdict['pseudospecies'].iloc[0], np.ndarray)
                    or (jobtype == 'prescreening_allreactive' and isinstance(job_subdict['pseudospecies'].iloc[0], str))):
                REAC = job_subdict['pseudospecies'].iloc[0]
                REACLUMPED = job_subdict['pseudospecies']  # series

            elif jobtype == 'prescreening_equilibrium' and job_subdict['pseudospecies'].iloc[0] == 'all':
                # only accepted case where you consider all pseudospecies together: all isomers as one
                REAC = SPECIES_UNIMOL  # all species together
                REACLUMPED = pd.Series([SPECIES_UNIMOL], index=[
                                       SPECIES_UNIMOL[0]+'_L'])

            # PRODUCTS: irrelevant to the reactivity; set them as the "other" species
            PRODSLUMPED = pd.Series(SPECIES, index=[SPECIES])
            PRODSLUMPED = PRODSLUMPED.drop(index=REAC)
            PRODS = PRODSLUMPED.values

        # set output folder and simul_array
        fld = os.path.join(cwd, jobtype, REACLUMPED.index[0])
        simul_array = np.array(
            [fld, REAC, REACLUMPED, PRODS, PRODSLUMPED], dtype=object)
        simul_array = simul_array[np.newaxis, :]
        output_DF = pd.DataFrame(simul_array, index=np.arange(0, N_simul), columns=[
                                 'fld', 'REAC', 'REACLUMPED', 'PRODS', 'PRODSLUMPED'])

    else:
        # ALL THE OTHER SIMULATIONS
        if jobtype == 'preproc_irreversible':
            # empty dataframe; it must have 1 index for consistency with the rest of the simulation set
            N_simul = 1
            output_DF = pd.DataFrame(
                index=np.arange(0, N_simul), columns=['fld'])
            output_DF.loc[0]['fld'] = os.path.join(cwd, jobtype)

        elif jobtype == 'prescreening_equilibrium' or jobtype == 'prescreening_allreactive':

            if job_subdict['pseudospecies'].index[0] == 'all':
                # equilibrium: all isomers together
                if jobtype == 'prescreening_equilibrium':
                    N_simul = 1
                    # all pseudospecies together
                    REAC = SPECIES_UNIMOL  # all species together
                    REACLUMPED = pd.Series([SPECIES_UNIMOL], index=[
                                           SPECIES_UNIMOL[0]+'_L'])
                    # PRODUCTS: irrelevant to the reactivity; set them as the "other" species
                    PRODSLUMPED = pd.Series(SPECIES, index=[SPECIES])
                    PRODSLUMPED = PRODSLUMPED.drop(index=REAC)
                    PRODS = PRODSLUMPED.values
                    # set output folder and simul_array
                    fld = os.path.join(cwd, jobtype, REACLUMPED.index[0])
                    simul_array = np.array(
                        [fld, REAC, REACLUMPED, PRODS, PRODSLUMPED], dtype=object)
                    simul_array = simul_array[np.newaxis, :]
                    output_DF = pd.DataFrame(simul_array, index=np.arange(0, N_simul), columns=[
                                             'fld', 'REAC', 'REACLUMPED', 'PRODS', 'PRODSLUMPED'])

                if jobtype == 'prescreening_allreactive':
                    # create one set of simulations for each species
                    N_simul = len(SPECIES)
                    output_DF = pd.DataFrame(index=np.arange(0, N_simul), columns=[
                                             'fld', 'REAC', 'REACLUMPED', 'PRODS', 'PRODSLUMPED'])
                    ind = 0
                    for SP in SPECIES:
                        REAC = SP
                        REACLUMPED = pd.Series([SP], index=[SP])
                        # PRODUCTS: ALL THW OTHERS
                        PRODSLUMPED = pd.Series(SPECIES, index=[SPECIES])
                        PRODSLUMPED = PRODSLUMPED.drop(index=REAC)
                        PRODS = PRODSLUMPED.values
                        # set output folder and simul_array
                        fld = os.path.join(cwd, jobtype, REACLUMPED.index[0])
                        output_DF.loc[ind] = np.array(
                            [fld, REAC, REACLUMPED, PRODS, PRODSLUMPED], dtype=object)
                        ind += 1

            else:
                # same for both prescreenings
                N_simul = job_subdict['pseudospecies'].index.size
                output_DF = pd.DataFrame(index=np.arange(0, N_simul), columns=[
                                         'fld', 'REAC', 'REACLUMPED', 'PRODS', 'PRODSLUMPED'])
                ind = 0
                for SP in job_subdict['pseudospecies'].index:
                    # for each set of pseudospecies: generate dictionaries
                    # if isinstance(job_subdict['pseudospecies'].loc[SP],np.ndarray):
                    # also valid for single reactants in case of prescreening_allreactive
                    # generate reactant arrays
                    REAC = job_subdict['pseudospecies'].loc[SP]
                    REACLUMPED = pd.Series([REAC], index=[SP])  # series
                    # PRODUCTS: irrelevant to the reactivity; set them as the "other" species
                    PRODSLUMPED = pd.Series(SPECIES, index=[SPECIES])
                    PRODSLUMPED = PRODSLUMPED.drop(index=REAC)
                    PRODS = PRODSLUMPED.values
                    # select output folder and allocate to dataframe
                    fld = os.path.join(cwd, jobtype, SP)
                    output_DF.loc[ind] = np.array(
                        [fld, REAC, REACLUMPED, PRODS, PRODSLUMPED], dtype=object)
                    ind += 1

        elif jobtype == 'composition_selection' or 'lumping' or 'validation':
            # read the pseudospecies file
            pseudospecies_file = os.path.join(cwd, 'inp', 'pseudospecies.txt')
            pseudospecies_series, stable_species = readinp.read_pseudospecies(
                pseudospecies_file)

            # remove single species for composition_selection jobtype
            if jobtype == 'composition_selection':
                for SP in pseudospecies_series.index:
                    if isinstance(pseudospecies_series.loc[SP], np.ndarray) == False:
                        pseudospecies_series = pseudospecies_series.drop(index=[
                                                                         SP])

            # organize dataframe
            N_simul = pseudospecies_series.index.size
            output_DF = pd.DataFrame(index=np.arange(0, N_simul), columns=[
                                     'fld', 'REAC', 'REACLUMPED', 'PRODS', 'PRODSLUMPED'])
            ind = 0

            # now do the loop on the right indices
            for SP in pseudospecies_series.index:
                # for each set: setup simulations
                REACLUMPED = pd.Series(
                    [pseudospecies_series.loc[SP]], index=[SP])
                REAC = pseudospecies_series.loc[SP]
                PRODSLUMPED = pseudospecies_series.drop(index=[SP])
                if isinstance(REAC, np.ndarray):
                    PRODS = np.delete(stable_species, np.where(
                        [rr == stable_species for rr in REAC])[1])
                else:
                    PRODS = np.delete(stable_species, np.where(
                        [ss == REAC for ss in stable_species])[0])

                # select output folder and allocate to dataframe
                fld = os.path.join(cwd, jobtype, SP)
                output_DF.loc[ind] = np.array(
                    [fld, REAC, REACLUMPED, PRODS, PRODSLUMPED], dtype=object)
                ind += 1

    return output_DF
