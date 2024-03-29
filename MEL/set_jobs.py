import os
import numpy as np
import pandas as pd
import shutil
from . import read_input as readinp


def checkNsubfolders(fld, fldperlevel):
    '''
    Check in fld that you have all levels of subfolders required
    level1: fldperlevel[0]; level2: fldperlevel[1] ...
    '''
    flds = [fld]
    CHECK = 1
    for nflds in fldperlevel:
        fldlist = []
        for fld in flds:
            # check N of subdirs
            CHECK *= (len(os.listdir(fld)) == nflds)
            # add subdirs
            for subfld in os.listdir(fld):
                fldlist.append(os.path.join(fld, subfld))
        # now the new level to check is the list of subfolders
        flds = fldlist

    return CHECK


def setfolder(fld, filescheck=[]):
    '''
    Create folder if it does not exist
    and if it does not have a certain list of files or folders in it
    if folder was generated from scratch, return 1
    '''

    if not os.path.isdir(fld):
        # create folder
        os.makedirs(fld)
        CHECK = 1

    else:
        if filescheck:
            allfiles = all([file in os.listdir(fld) for file in filescheck])

            if allfiles:
                CHECK = 0  # do nothing
            else:
                rmfolder(fld)
                os.makedirs(fld)
                CHECK = 1
        else:
            # do nothing
            CHECK = 0

    return CHECK


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
        F_name, F_format = F.split('.')
        oldname = os.path.join(fld, F)
        newname = os.path.join(fld, F_name + '_old' + N + '.' + F_format)
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
            # check that all pseudospecies are available
            for SP_stable in stable_species:
                if SP_stable not in SPECIES:
                    raise(ValueError(
                        '*Species {} of pseudospecies file is not among the mechanism species - check and retry'.format(SP_stable)))

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
