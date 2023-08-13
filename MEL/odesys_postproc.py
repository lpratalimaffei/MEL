import os
import numpy as np
import pandas as pd
import shutil
from . import preprocessing as preproc
from .set_jobs import  setfolder
import warnings
# warnings.simplefilter(action='ignore', category=RuntimeWarning) # ignore float division by 0 /other producing nan/inf, then removed
warnings.filterwarnings('ignore')

class ODE_POSTPROC:
    '''
    In this class the output of the ODE is post-processed and the output is written as required by optiSMOKE++
    a. Read the Output\Output.OUT file and extract the columns => generate t,W
    b. Once extracted: generate the files corresponding to the profile of every species and write them
    c. Write "path_to_Exp_Datasets.txt" indicating the experimental dataset considered
    d. Copy input_OS.dic in the new folder
    '''

    def __init__(self, cwd):
        '''
        Indicates the path where to store the output. 
        It can be called either at the beginning of the main or at the first iteration
        '''
        self.cwd = cwd
        # preallocate the list of pathways to the experimental datasets you will generate.
        # the paths to the OS input correspond to the datasets.
        self.path_to_Exp_Datasets = []
        self.path_to_OS_inputs = []

    def SET_BRANCHINGFOLDERS(self, jobtype, REACLUMPED, PRODSLUMPED, T_VECT):
        '''
        Set the subfolder "BF_OUTPUT" to store the data of the branching of the lumped products
        It also allocates the dictionary where to store the branchings, which will be written at the end of each single-P simulation
        NB this function is to be called only if PRODSLUMPED contains arrays. if not: raise exception
        '''
        self.lumped_branching = {}

        # allocate the composition of the lumped reactant
        if isinstance(REACLUMPED[0], np.ndarray):
            # generate empty dataframe:
            self.lumped_branching_reac = pd.DataFrame(
                index=T_VECT, columns=REACLUMPED[0])

        # allocate only the lumped products with arrays
        for PRi_L in PRODSLUMPED.index:
            # the indexes are the names of the lumped species
            if isinstance(PRODSLUMPED[PRi_L], np.ndarray):
                # generate empty dataframe
                df_PRi_L = pd.DataFrame(
                    index=T_VECT, columns=PRODSLUMPED[PRi_L])
                # allocate it in dictionary
                self.lumped_branching[PRi_L] = df_PRi_L

        # set branching path
        self.branchingpath = os.path.join(
            self.cwd, jobtype, 'BF_OUTPUT', REACLUMPED.index[0])

    def MAKE_FOLDERS(self, fld, P, T, REACLUMPED):
        '''
        Makes the subfolders at the selected conditions of reactant, temperature and pressure
        '''
        # extract the name of the reactant
        REAC = REACLUMPED.index[0]
        self.fld = fld
        self.dir_PT = os.path.join(fld, str(P) + 'atm', str(T) + 'K')
        _ = setfolder(self.dir_PT)

        # Allocate T,P,REAC to "self" for successive use
        self.T = T
        self.P = P
        self.REACNAME = REAC
        self.REAC = REACLUMPED[0]
        # NB self.REAC is the array of reactants, self.REACNAME is the name (single reactant or R_L)

    def EXTRACT_PROFILES(self, SPECIES, i_REAC, N_INIT_REAC, SPECIES_BIMOL_SERIES, ISOM_EQUIL, CUTOFF):
        '''
        Read the profiles obtained by OpenSMOKE++ and store them into arrays
        '''

        # if the reaction is bimolecular: read also the xi of the second species.
        # for bimolecular lumped species: set i_REAC as the first species
        if isinstance(i_REAC, np.ndarray):
            i_REAC_0 = i_REAC[0]
        else:
            i_REAC_0 = i_REAC

        if SPECIES_BIMOL_SERIES.iloc[i_REAC_0] != '' and SPECIES_BIMOL_SERIES.iloc[i_REAC_0] != SPECIES[i_REAC_0]:
            # if the second reactant has all same fragment: 1 extracolumn
            extracol = 1
        else:
            extracol = 0

        # read the output file
        filename = os.path.join(self.cwd, 'Output', 'Output.out')
        if os.path.isfile(filename):
            cols_species = np.arange(9, 9+len(SPECIES)+extracol)
            # 0 IS THE INDEX, AND [0] IS THE VALUE TO INSERT
            n_cols = np.insert(cols_species, 0, [0])
            data = np.genfromtxt(filename, dtype=float,
                                 skip_header=1, usecols=(n_cols))

            # extract also PV (needed for the total number of moles)
            data_PV = np.genfromtxt(
                filename, dtype=float, skip_header=1, usecols=(np.array([5, 6], dtype=int)))

            # check if array is 1D - something went wrong with the simulations, give a warning and make it so that the function works anyways
            if len(data.shape) == 1:
                print('*Warning: shape of array is 1D - something wrong with the simulation - check convergence / tolerance params')
                data = np.array([data]*3)
                data_PV = np.array([data_PV]*3)
                
            # lumped reactant: sum the profiles, redefine a new reactant index and name, new species and species bimol series
            if isinstance(i_REAC, np.ndarray):
                # 0: for the calculation of the derivatives: compute the absolute variation of each isomer
                dreac_i_dt = abs((-data[2:, i_REAC+1]+data[1:-1, i_REAC+1])/(
                    data[2:, 0]-data[1:-1, 0]).reshape(len(data[2:, 0]), 1))
                dreac_dt = np.sum(dreac_i_dt, axis=1)
                # 1: sum the profiles of all the reactants; save the reactant composition for later
                Wreac_new = np.sum(data[:, i_REAC+1], axis=1)
                Wreac_new = Wreac_new[np.newaxis, :]
                self.Wreac_composition = data[:, i_REAC+1]
                # 2: delete all the columns corresponding to the reactant and add the reactant in the first column
                #    do the same for the names of the species and the bimolecular species
                data = np.delete(data, i_REAC+1, axis=1)
                data = np.insert(data, 1, Wreac_new, axis=1)
                SPECIES = np.delete(SPECIES, i_REAC)
                bim_reac = pd.Series(
                    SPECIES_BIMOL_SERIES[i_REAC_0], index=[self.REACNAME])
                SPECIES_BIMOL_SERIES = SPECIES_BIMOL_SERIES[SPECIES]
                SPECIES = np.insert(SPECIES, 0, self.REACNAME)
                SPECIES_BIMOL_SERIES = pd.concat(
                    [bim_reac, SPECIES_BIMOL_SERIES])
                # 3: assign new indices
                i_REAC = 0  # now the reactant is in the first position
                reaclumped = 'YES'
            else:
                # compute the derivative of the reactant consumption
                dreac_dt = (-data[2:, i_REAC+1]+data[1:-1,
                                                    i_REAC+1])/(data[2:, 0]-data[1:-1, 0])
                # RUNTIMEWARNING
                # save variables forlater
                reaclumped = 'NO'

            # cut the profiles where needed

            i_in = np.where(data[:, i_REAC+1] <= (1-CUTOFF[0])*N_INIT_REAC)
            i_fin = np.where(data[:, i_REAC+1] <= (1-CUTOFF[1])*N_INIT_REAC)
            # if the reactant does not reach the minimum consumption (possible for lumped reactants): set the initial value as 0
            if len(i_in[0]) == 0:
                i_in = 0
            else:
                i_in = i_in[0][0]
            # impose to cut when the DERIVATIVE of the reactant (consumption) reaches a small value
            if len(i_fin[0]) == 0:
                #i_fin = np.where(dreac_dt[3:] < dreac_dt[3]*1e-4)
                # remove infinite values and keep only positive values
                dreac2_dt = dreac_dt[(np.isinf(dreac_dt) == False) & (
                    np.isnan(dreac_dt) == False) & (dreac_dt > 0)]

                if len(dreac2_dt) <= 1:
                    # include also length = 1 otherwise it means that i_in and i_fin will be the same
                    i_in = 0
                    i_fin = len(dreac_dt)
                else:
                    maxderiv = max(dreac2_dt)
                    minderiv = min(dreac2_dt)
                    i_in = np.where(dreac_dt == maxderiv)[0][0]
                    if minderiv <= maxderiv*1e-4:
                        cutoff_deriv = dreac2_dt[dreac2_dt < maxderiv*1e-4][0]
                        i_fin = np.where(dreac_dt >= cutoff_deriv)[0][-1]
                    elif minderiv > maxderiv*1e-4:
                        i_fin = np.where(dreac_dt >= minderiv)[0][-1]

            else:
                i_fin = i_fin[0][0]
            # check that i_fin > i_in, otherwise set i_in to 0
            if i_fin < i_in:
                i_in = 0
            # save data in the appropriate range of consumption of the reactant
            data = data[i_in:i_fin, :]
            data_PV = data_PV[i_in:i_fin, :]
            t = data[:, 0]
            t = t[:, np.newaxis]
            # W will have the mole fractions multiplied by PV/RT (CtotV = Ntot)
            # IF YOU HAVE AN EXTRA SPECIES (NAMELY: BIMOLECULAR REACTANT), MULTIPLY THE COLUMN OF THE REACTANT BY THAT
            W = data[:, 1:]
            if extracol == 1:
                # multiply by the number of the second fragments s.t. you reach the total fraction of N_ABU
                W_reac2 = W[:, -1]
                W = W[:, :-1]
                # self.W will be used to write the profiles to the optimizer, so no modification should be done
                self.W = pd.DataFrame(W, columns=SPECIES)
                # after this, multiply W[reac] by the bimolecular reactant to obtain xi^2
                W[:, i_REAC] = W[:, i_REAC]*W_reac2
            else:
                self.W = pd.DataFrame(W, columns=SPECIES)

            self.t = t
            tW_DF = np.concatenate((t, W), axis=1)
            tW_DF = pd.DataFrame(tW_DF, columns=np.insert(SPECIES, 0, 't'))

            # FOR LUMPED REACTANTS: SAVE THE BRANCHING FRACTIONS FOR LATER
            if reaclumped == 'YES':
                # reactant composition
                # print(self.Wreac_composition,i_in,i_fin)
                self.Wreac_composition = self.Wreac_composition[i_in:i_fin, :]
                # if ISOM_EQUIL is active: take only the last BF
                if ISOM_EQUIL == 1:
                    self.lumped_branching_reac.loc[self.T, self.REAC] = self.Wreac_composition[-1, :]/np.sum(
                        self.Wreac_composition[-1, :])
                elif ISOM_EQUIL == 0:
                    Wreac_tot = np.sum(self.Wreac_composition[1:], axis=1)
                    Wreac_tot = Wreac_tot[:, np.newaxis]
                    if len(self.t) > 0:
                        dtweight = ((self.t[1:]-self.t[:-1])/self.t[-1])
                        br_weighted = self.Wreac_composition[1:,
                                                            :]/Wreac_tot*dtweight
                        self.lumped_branching_reac.loc[self.T, self.REAC] = np.sum(
                            br_weighted, axis=0)
                    else:
                        # keep initial composition
                        self.lumped_branching_reac.loc[self.T, self.REAC] = self.Wreac_composition[0, :]
                # save the reactant composition separately for plotting

        else:
            raise ValueError('OS output file not found')
        
        self.tW_DF = tW_DF

        # save species and bimol series for the following steps
        self.i_REAC = i_REAC
        self.SPECIES = SPECIES
        self.SPECIES_BIMOL_SERIES = SPECIES_BIMOL_SERIES

        return tW_DF, data_PV

    def PROFILES_REAC_COMPOSITION(self):
        '''
        Method to be called only in presence of a lumped reactant.
        It returns the profile of lumped set of reactants in case it is needed for later use
        '''
        try:
            tW_DF_reac = np.concatenate(
                (self.t, self.Wreac_composition), axis=1)
            tW_DF_reac = pd.DataFrame(
                tW_DF_reac, columns=np.insert(self.REAC, 0, 't'))
            return tW_DF_reac

        except ValueError as e:
            print(str(e))

    def LUMP_PROFILES(self, PRODS, PRODSLUMPED):
        '''
        IN THIS METHOD:
        - TAKE THE PROFILES OF self.W and sum those of the lumped products
        - redefine self.PRODS as the names of the lumped products and rewrite self.W and tW_DF
          NB also include the species that are not part of PRODSLUMPED or reac.
        - Allocate the branchings within each lumped products to the appropriate dictionary
        - for later processing: define new
               REAC_L
               i_REAC_L
               SPECIES_SERIES_L
               SPECIES_BIMOL_SERIES_L
               PRODS_L
               these have the same format as the initial ones.
        '''
        self.PRODS = PRODS
        self.PRODSLUMPED = PRODSLUMPED
        # non-lumped products: return all the values the same as before
        if len(self.PRODS) == len(PRODSLUMPED):
            tW_DF = self.tW_DF
            i_REAC_L = self.i_REAC
            SPECIES_L = self.SPECIES
            SPECIES_SERIES_L = pd.Series(
                np.arange(0, len(self.SPECIES)), index=self.SPECIES)
            SPECIES_BIMOL_SERIES_L = self.SPECIES_BIMOL_SERIES
            PRODS_L = self.PRODS
        else:
            # redefine the product names
            self.PRODS = np.array(PRODSLUMPED.index, dtype='<U16')
            # empty dataframe for the products
            W_prods_L = pd.DataFrame(columns=self.PRODS, dtype=float)
            # empty series for bimolecular species
            PRODS_L_BIMOL = pd.Series(index=self.PRODS, dtype=str)
            # lumped products: go over the lumped products and generate W_prods
            # delete the corresponding columns in the dataframe
            for PRi_L in self.PRODS:
                # PRODSLUMPED[PRi_L] will be the value (either string or array, depending on whether the product is lumped or not)
                PRi_L_value = PRODSLUMPED[PRi_L]
                # if the product is just a string (single product): just save the corresponding line and delete it
                if isinstance(PRi_L_value, str):
                    # NB here we have PRi_L_value==PRi_L, so using one or the other makes no difference
                    # delete the column from self.W and move it to W_prods_L
                    W_prods_L[PRi_L] = self.W[PRi_L]
                    # now delete the column from self.W
                    self.W = self.W.drop(columns=PRi_L)
                    # reconstruct the series of bimolecular species
                    PRODS_L_BIMOL[PRi_L] = self.SPECIES_BIMOL_SERIES[PRi_L]

                elif isinstance(PRi_L_value, np.ndarray):
                    # for each of the products: sum the columns
                    W_prods_L[PRi_L] = np.sum(self.W[PRi_L_value], axis=1)

                    # compute the weighted average branching within each product and save them to dataframe
                    Wtot = W_prods_L[PRi_L].values[1:, np.newaxis]
                    dtweight = ((self.t[1:]-self.t[:-1])/self.t[-1])
                    BR_PRi_L = self.W.loc[1:, PRi_L_value]/Wtot*dtweight
                    self.lumped_branching[PRi_L].loc[self.T,
                                                     PRi_L_value] = np.sum(BR_PRi_L, axis=0)

                    # delete the corresponding columns
                    self.W = self.W.drop(columns=PRi_L_value)
                    # reconstruct the series of bimolecular species: take the corresponding species of the first lumped  species
                    PRODS_L_BIMOL[PRi_L] = self.SPECIES_BIMOL_SERIES[PRi_L_value[0]]

            # now you deleted from self.W all the product profiles. Put them back at the end by concatenating the dataframes
            W_noprods = self.W  # save it for later
            self.W = pd.concat([self.W, W_prods_L], axis=1)
            # new dataframe
            timeseries = pd.DataFrame(self.t, columns=['t'])
            tW_DF = pd.concat([timeseries, self.W], axis=1)
            # new names
            SPECIES_L = self.W.columns
            SPECIES_SERIES_L = pd.Series(
                np.arange(0, len(SPECIES_L)), index=SPECIES_L)
            i_REAC_L = SPECIES_SERIES_L[self.REACNAME]
            PRODS_L = self.PRODS
            # bimolecular species: # first select the non-product species
            SPECIES_BIMOL_L = self.SPECIES_BIMOL_SERIES[W_noprods.columns]
            SPECIES_BIMOL_SERIES_L = pd.Series(
                SPECIES_BIMOL_L, index=W_noprods.columns)
            # now concatenate this series with that of the products
            SPECIES_BIMOL_SERIES_L = pd.concat(
                [SPECIES_BIMOL_SERIES_L, PRODS_L_BIMOL])
            # UPDATE THE VALUE FOR THE FOLLOWING STEPS
            self.SPECIES_BIMOL_SERIES = SPECIES_BIMOL_SERIES_L

        return tW_DF, self.REACNAME, i_REAC_L, SPECIES_L, SPECIES_SERIES_L, SPECIES_BIMOL_SERIES_L, PRODS_L

    def WRITE_PROFILES(self, PRODS):
        '''
        Writes single files in Output_to_optiSMOKE/P_reac/T
        With the profile of each species
        Finally write the file Path_to_Exp_Datasets.txt with the names of the files
        '''
        PRODS = np.array(PRODS, dtype='<U20')
        self.PRODS = PRODS
        # indices of the reactant and of the products
        indices_R_prods = np.insert(self.PRODS, 0, self.REACNAME)
        # empty matrix for the results
        exp_dataset = np.zeros((np.shape(self.t)[0], 3*len(indices_R_prods)))
        # first column with the time
        exp_dataset[:, 0::3] = self.t
        # if you have only 1 species to write:
        W_reduced = self.W[indices_R_prods]
        exp_dataset[:, 1::3] = W_reduced
        # third column with the error
        exp_dataset[:, 2::3] = 0.1*np.ones(self.t.shape)
        # Write the profiles ONLY FOR THE REACTANT AND THE LUMPED PRODUCTS
        header = 'Batch m_SP {} '.format(len(self.PRODS)+1)
        header += '\t\t\t\t\t\t\t\t\t'.join(indices_R_prods)

        np.savetxt(os.path.join(self.dir_PT, str(self.T) + '.txt'),
                   exp_dataset, header=header, delimiter='\t', fmt='%.2e', comments='')
        self.path_to_Exp_Datasets.append(os.path.join(
            str(self.P) + 'atm', str(self.T) + 'K', str(self.T) + '.txt'))
        self.path_to_OS_inputs.append(os.path.join(
            str(self.P) + 'atm', str(self.T) + 'K', 'input_OS.dic'))

    def CHECK_PROD_SELECTIVITY(self, profiles_P, min_sel = 0.001):
        '''
        Check the products which accumulate above min_sel% for at least one of the T,P investigated
        '''
        PRODS_SEL = pd.Series(np.NaN, index=self.PRODS)
        # max BF will be the starting one of the reactant
        for P in profiles_P.keys():
            profiles_T = profiles_P[P]
            for T in profiles_T.keys():
                profiles = profiles_T[T]
                # get the max profile of the reactant
                max_reac = max(profiles[self.REACNAME])
                min_reac = min(profiles[self.REACNAME])
                denom = (max_reac-min_reac)
                if denom > 1e-30: 
                    # otherwise it means the reactant was not consumed at all and it does not make sense to check the selectivity
                    # get product selectivity
                    for pr in self.PRODS:
                        max_bf_pr = max(profiles[pr])
                        sel_pr = max_bf_pr/denom
                        if all(sel_pr >= np.array([min_sel, PRODS_SEL[pr]])):
                            PRODS_SEL[pr] = sel_pr
                        elif sel_pr >= min_sel:
                            PRODS_SEL[pr] = sel_pr

        PRODS_SEL = PRODS_SEL.dropna() # drop the NaN values
        PRODS_SEL = PRODS_SEL.sort_values(ascending=False) # order by selectivity
        # write file
        prod_sel = np.concatenate(
            (np.array(PRODS_SEL.index)[:, np.newaxis], np.array(PRODS_SEL.values)[:, np.newaxis]), axis=1, dtype=object)
        # save the accumulating species in the corresponding folder
        np.savetxt(os.path.join(self.fld, 'prods_selectivity.txt'),
                   prod_sel, header= ' '*12 + 'prod' + ' '*3 + 'max_sel', fmt=['%16s','%10.2e'], comments='')

    def WRITE_BRANCHINGS_PRODS(self):
        '''
        This method writes the profiles of the lumped products in the folder "Branchings"
        '''
        # products:
        # self.lumped_branching knows already if prods are lumped
        # #if len(self.PRODSLUMPED) != len(PRODS):
        for PRi_L in self.lumped_branching:
            # make corresponding subfolder if it does not exist
            _ = setfolder(os.path.join(self.branchingpath, PRi_L))
            BRfile = os.path.join(self.branchingpath, PRi_L,
                                str(self.P)+'atm.txt')
            # concatenate 2 dataframes: write also the values of the temperature
            # NB for concatenation along rows, the same index is needed!
            T_DF = pd.DataFrame(
                self.lumped_branching[PRi_L].index, index=self.lumped_branching[PRi_L].index, columns=['T[K]'])
            BRall = pd.concat([T_DF, self.lumped_branching[PRi_L]], axis=1)
            formats = pd.Series(index=BRall.columns, dtype=str)
            formats[self.lumped_branching[PRi_L].columns] = '%1.5f'
            formats['T[K]'] = '%d'
            formats_list = list(formats.values)
            head = '\t'.join(BRall.columns)
            np.savetxt(BRfile, BRall, fmt=formats_list,
                        header=head, comments='\t')

    def WRITE_BRANCHINGS_REACS(self):
        # lumped reactant:
        if isinstance(self.REAC, np.ndarray):
            # folder
            _ = setfolder(os.path.join(self.branchingpath, self.REACNAME))
            BRfile = os.path.join(self.branchingpath,
                               self.REACNAME, str(self.P)+'atm.txt')
            #fld = self.branchingpath + '/' + self.REACNAME + '_from' + self.REACNAME + '_' + str(self.P) + 'atm.txt'
            # concatenation
            T_DF = pd.DataFrame(self.lumped_branching_reac.index,
                                index=self.lumped_branching_reac.index, columns=['T[K]'])
            BRall = pd.concat([T_DF, self.lumped_branching_reac], axis=1)
            formats = pd.Series(index=BRall.columns, dtype=str)
            formats[self.REAC] = '%1.5f'
            formats['T[K]'] = '%d'
            formats_list = list(formats.values)
            head = '\t'.join(BRall.columns)
            np.savetxt(BRfile, BRall, fmt=formats_list,
                       header=head, comments='\t')
            return BRall

    def WRITE_NEW_OSINPUT(self, N_INIT):
        '''
        This method writes the new OS input in the selected subfolders
        '''
        # with the new indices, the reactants and products are lumped
        new_indices = np.insert(self.PRODS, 0, self.REACNAME)
        _ = setfolder(os.path.join(self.dir_PT, 'inp'))
        # Copy the input to OS simulations and substitute the values of interest
        shutil.copyfile(os.path.join(self.cwd, 'inp', 'input_OS_template.dic'), os.path.join(
            self.dir_PT, 'inp', 'input_OS_template.dic'))
        # write new input in the selected folder
        write_OS_new = preproc.WRITE_OS_INPUT(self.dir_PT, self.T, self.P, new_indices, pd.Series(
            self.REACNAME, index=[self.REACNAME]), N_INIT, self.SPECIES_BIMOL_SERIES, '')
        write_OS_new.w_input_OS()
        # Remove the template from the subfolders
        shutil.rmtree(os.path.join(self.dir_PT, 'inp'))

    def WRITE_FINAL_PATHS(self):
        '''
        Write the paths to the experimental datasets in the folder named after the reactant
        '''
        np.savetxt(os.path.join(self.fld, 'Path_to_Exp_Datasets.txt'),
                   self.path_to_Exp_Datasets, fmt='%s')
        np.savetxt(os.path.join(self.fld, 'Path_to_OS_inputs.txt'),
                   self.path_to_OS_inputs, fmt='%s')


def OVERALL_SELECTIVITY(jobfld, fld_list, threshold):
    tokeep = []
    for fld in fld_list:
        with open(os.path.join(fld, 'prods_selectivity.txt')) as selfile:
            for line in selfile:
                prod, sel = line.split()
                if prod != 'prod' and prod not in tokeep and sel >= threshold:
                    tokeep.append(prod)
        selfile.close()
        
    tokeep.sort()
    newfile = open(os.path.join(jobfld, 'speciestokeep.txt'), 'w')
    newfile.write('\n'.join(tokeep))
    newfile.close()
                
        
        