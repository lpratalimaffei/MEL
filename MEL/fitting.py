import os
import sys
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit
from . import preprocessing as preproc
import warnings
warnings.filterwarnings('ignore')

def arrhenius_fit(T_VECT, k_VECT):
    '''
    This function provides an Arrhenius fit of the k(T) provided as
    k = k0 * T^alpha * exp(-EA/R/T) with non linear least squares
    The function also returns the quality of the fit R2
    '''
    def ln_k(T, ln_k0, alpha, EA):
        return ln_k0+alpha*np.log(T)-EA/1.987/T

    ln_k_VECT = np.log(k_VECT)
    popt, _ = curve_fit(ln_k, T_VECT, ln_k_VECT)
    # get the model parameters
    k0 = np.exp(popt[0])
    alpha = popt[1]
    EA = popt[2]
    # get the adjusted R2
    R2 = 1 - np.sum((ln_k_VECT-ln_k(T_VECT,
                                    popt[0], alpha, EA))**2)/np.sum((ln_k_VECT-np.mean(ln_k_VECT))**2)
    # 2 is the number of parameters in the model excluding the constant, and len(T_VECT) is the number of observations
    R2adj = 1-(1-R2)*(len(T_VECT)-1)/(len(T_VECT)-1-2)

    # check on floating precision for k0 value
    if k0 <= sys.float_info.min or k0 == 0:
        k0 = np.inf
        
    return k0, alpha, EA, R2adj


class FITTING:
    '''
    class providing a dataframe with the fits and the fit qualities of the profiles obtained from solving the ode
    as linear fits.
    '''

    def __init__(self, T_VECT, REAC, PRODS):
        '''
        Generates the dataframes to allocate the rate constants and the fitting error
        '''
        self.REAC = REAC
        self.PRODS = PRODS
        self.T_VECT = T_VECT
        prods_col = np.copy(PRODS)
        prods_col = list(prods_col)
        prods_col.append(REAC)
        # frame to allocate the R2 in the fits of the profiles
        self.data_P_R2 = pd.DataFrame(np.zeros((len(T_VECT), len(
            PRODS)+1), dtype=np.dtype('<U50')), index=list(T_VECT), columns=prods_col)
        # frame to allocate the fitted rate constants and the corresponding error of the linear fit
        self.data_P_fits_wERR = pd.DataFrame(np.zeros((len(T_VECT), len(
            PRODS)+1), dtype=np.dtype('<U50')), index=list(T_VECT), columns=prods_col)
        # append another column to the products
        prods_col.append('sum_k_Pr_i')
        # frame to allocate the fitted rate constants at every temperature
        self.data_P_fits = pd.DataFrame(
            np.zeros((len(T_VECT), len(PRODS)+2)), index=list(T_VECT), columns=prods_col)

    def fit_profiles(self, tW_DF, i_REAC, SPECIES_SERIES, T, PV, SPECIES_BIMOL_SERIES, N_INIT):
        '''
        This function takes the profiles dWi/dt VS W_reac and fits it with a constant k for every Product.
        Also the fit of the profile of the reactant is performed and compared with the sum of the k_Pr_i obtained.
        It returns the values of the rate constants k_Pr_i at the defined T.
        The rates are still normalized because they should be used to simulate the profiles with the lumped model.

        next: compare the sum of the k_Pr_i obtained with the reactant decay and print a warning if the ks are different
        '''
        # reconstruct the derivatives in time _ Central Difference Scheme: dWi/dt = (W(i+1/2))-W(i-1/2))/(t(i+1/2)-t(i-1/2))
        t = tW_DF['t'].values
        t = t[:, np.newaxis]
        W = tW_DF[SPECIES_SERIES.index].values
        # turn the W (xi) into moles by multiplication by Ntot
        Pprofile = PV[:, 0]
        Pprofile = Pprofile[:, np.newaxis]
        Vprofile = PV[:, 1]
        Vprofile = Vprofile[:, np.newaxis]
        W = W*Pprofile*Vprofile/8.314/T
        dWi = (W[2:, :]-W[0:-2, :])/(t[2:]-t[0:-2])
        # bimol reaction: multiply by x_ABU*P/RT (x_ABU approximated with N_INIT)
        if SPECIES_BIMOL_SERIES.iloc[i_REAC] != '' and SPECIES_BIMOL_SERIES.iloc[i_REAC] != SPECIES_SERIES.index[i_REAC]:
            # the W[i_REAC] is already xi^2*Ntot, so I need to multiply it again by Ntot
            Wreac_fordWi = W[1:-1, i_REAC]*PV[1:-1, 0]*PV[1:-1, 1]/8.314/T
            # self reaction: the W has to be multiplied by itself (not done at the beginning)
        elif SPECIES_BIMOL_SERIES.iloc[i_REAC] == SPECIES_SERIES.index[i_REAC]:
            # take the square (already in moles)
            Wreac_fordWi = W[1:-1, i_REAC]**2
        else:
            Wreac_fordWi = W[1:-1, i_REAC]
        Wreac_fordWi = Wreac_fordWi[:, np.newaxis]
        # delete nan or inf values of dWi
        index_nan = np.where(np.isnan(dWi))[0]  # rows corresponding to nan
        index_inf = np.where(np.isinf(dWi))[0]  # rows corresponding to inf
        index_remove = np.concatenate(
            (index_nan, index_inf))  # all the rows to delete
        dWi = np.delete(dWi, index_remove, axis=0)
        Wreac_fordWi = np.delete(Wreac_fordWi, index_remove, axis=0)
        # fits
        for Pr_i in self.PRODS:
            Pr_i_index = SPECIES_SERIES[Pr_i]
            # IF THE PRODUCT IS OF THE KIND 2A (SELF-PRODUCT): SET THE CONSTANT TO 2, OTHERWISE KEEP 1
            if Pr_i == SPECIES_BIMOL_SERIES[Pr_i]:
                c = 2
            else:
                c = 1
            # alternative:
            model_np = np.linalg.lstsq(
                c*Wreac_fordWi, dWi[:, Pr_i_index], rcond=None)
            slope_np = model_np[0]
            model = LinearRegression(fit_intercept=False).fit(
                c*Wreac_fordWi, dWi[:, Pr_i_index])  # set the intercept to 0
            fiterr = model.score(c*Wreac_fordWi, dWi[:, Pr_i_index])
            # save directly the dimensional coefficient
            self.data_P_fits[Pr_i][T] = model.coef_[0]
            self.data_P_R2[Pr_i][T] = 'R2 = {:.2f}'.format(fiterr)
            self.data_P_fits_wERR[Pr_i][T] = '{:.4e}({:.2f})'.format(
                model.coef_[0], fiterr)
            # save directly the dimensional coefficient
            self.data_P_fits['sum_k_Pr_i'][T] += model.coef_[0]

        # fit for the reactant: multiply by 2 in case of self-reaction
        if SPECIES_BIMOL_SERIES.iloc[i_REAC] == SPECIES_SERIES.index[i_REAC]:
            model_reac = LinearRegression(fit_intercept=False).fit(
                2*Wreac_fordWi, -dWi[:, i_REAC])
            fiterr = model_reac.score(2*Wreac_fordWi, -dWi[:, i_REAC])
        else:
            model_reac = LinearRegression(fit_intercept=False).fit(
                Wreac_fordWi, -dWi[:, i_REAC])
            fiterr = model_reac.score(Wreac_fordWi, -dWi[:, i_REAC])

        self.data_P_fits[self.REAC][T] = model_reac.coef_[0]
        
        # save directly the dimensional coefficient
        self.data_P_R2[self.REAC][T] = 'R2 = {:.2f}'.format(fiterr)
        self.data_P_fits_wERR[self.REAC][T] = '{:.4e}({:.2f})'.format(
            model_reac.coef_[0], fiterr)

    def prod_bfs(self):
        # get product branching fractions at a given P
        BFS = self.data_P_fits[self.PRODS].values/np.array(self.data_P_fits['sum_k_Pr_i'])[:, np.newaxis]
        self.BF_P_fits =  pd.DataFrame(BFS, columns = self.PRODS, index = self.data_P_fits.index)        
        
    def write_originalk(self, out_fld, verbose = False):

        if verbose:
            print(self.data_P_fits)
            print(self.data_P_fits_wERR)
        # Save the profiles in the appropriate folder
        head = np.append('T[K]', self.data_P_fits.columns)
        # head = head[np.newaxis,:]
        headbf = ''.join([' '*(14-len(el)) + el for el in head[:-2]])
        head = ''.join([' '*(14-len(el)) + el for el in head])
        head_err = ''.join([' '*(14-len(el)) + el for el in head[:-1]])
        
        BF_towrite = np.concatenate(
            (self.T_VECT[:, np.newaxis], self.BF_P_fits.values), axis=1)
        profiles_towrite = np.concatenate(
            (self.T_VECT[:, np.newaxis], self.data_P_fits.values), axis=1)
        profiles_wERR_towrite = np.concatenate(
            (self.T_VECT[:, np.newaxis], self.data_P_fits_wERR.values), axis=1)
        #profiles_final = np.concatenate((header,profiles_towrite),axis=0)
        # save the rates (before fitting) and the corresponding fitting error
        np.savetxt(os.path.join(out_fld, 'bf.txt'),
                   BF_towrite, header=headbf, delimiter='\t', fmt='%.2e')
        np.savetxt(os.path.join(out_fld, 'rates.txt'),
                   profiles_towrite, header=head, delimiter='\t', fmt='%.2e')
        np.savetxt(os.path.join(out_fld, 'rates_wERR.txt'),
                   profiles_wERR_towrite, header=head_err, delimiter='\t', fmt='%s')

    def fits_lumped_k(self, SPECIES_BIMOL_SERIES, bfthreshold = 0.):
        '''
        Return the dataframes in the current status
        Generate matrices with arrhenius fits and comments
        BFTHRESHOLD: minimum product branching fraction to keep the rate constant
        '''

        # Allocate the matrices of the rate constants
        k0_MAT = np.ones((len(self.PRODS)+1, len(self.PRODS)+1))
        alpha_MAT = np.ones((len(self.PRODS)+1, len(self.PRODS)+1))
        EA_MAT = np.ones((len(self.PRODS)+1, len(self.PRODS)+1))
        comments_MAT = np.zeros(
            (len(self.PRODS)+1, len(self.PRODS)+1), dtype=np.dtype('<U100'))

        # Allocate new indices
        new_indices = np.append(np.array(self.REAC), np.array(self.PRODS))
        new_SPECIES_SERIES = pd.Series(
            np.arange(0, len(self.PRODS)+1), index=new_indices)
        new_SPECIES_BIMOL_SERIES = SPECIES_BIMOL_SERIES[new_indices].values
        new_i_REAC = 0
        new_i_PRODS = np.arange(1, len(self.PRODS)+1)

        # Arrhenius fits for every product, ant put them in the matrices
        for Pr_i in self.PRODS:
            i_prod = new_SPECIES_SERIES[Pr_i]
            k_tofit = self.data_P_fits.loc[:][Pr_i].values
            # if overflow is encountered in exp, it means that there is a discontinuity in the fitting range: reduce it
            # rows corresponding to small values of k
            index_ovf = np.where(k_tofit < 1e-90)[0]
            # fit and store the values in the matrices

            if len(index_ovf) == 0:
                k0_MAT[0, i_prod], alpha_MAT[0, i_prod], EA_MAT[0,
                                                                i_prod], R2adj = arrhenius_fit(self.T_VECT, k_tofit)
                comments_MAT[0, i_prod] = 'R2_adj = {:.2f} , T RANGE = {T_IN} - {T_FIN} K'.format(
                    R2adj, T_IN=self.T_VECT[0], T_FIN=self.T_VECT[-1])
            else:
                k_tofit = np.delete(k_tofit, index_ovf, axis=0)
                T_VECT_NEW = np.delete(self.T_VECT, index_ovf, axis=0)
                try:
                    k0_MAT[0, i_prod], alpha_MAT[0, i_prod], EA_MAT[0,
                                                                    i_prod], R2adj = arrhenius_fit(T_VECT_NEW, k_tofit)
                    comments_MAT[0, i_prod] = 'R2_adj = {:.2f} , T RANGE = {T_IN} - {T_FIN} K'.format(
                        R2adj, T_IN=T_VECT_NEW[0], T_FIN=T_VECT_NEW[-1])
                except TypeError:
                    print('not enough values for the fit of %s -> %s: setting val to inf' % (self.REAC, Pr_i))
                    k0_MAT[0, i_prod], alpha_MAT[0, i_prod], EA_MAT[0, i_prod] = [np.inf]*3
                    comments_MAT[0, i_prod] = 'T RANGE = {T_IN} - {T_FIN} K: too small for fit'.format(
                        T_IN=T_VECT_NEW[0], T_FIN=T_VECT_NEW[-1])
                except ValueError:
                    print('fitting unavailable - probably species are not connected')
                    k0_MAT[0, i_prod], alpha_MAT[0, i_prod], EA_MAT[0, i_prod] = [np.inf]*3
                    comments_MAT[0, i_prod] = 'probably species not connected'
                    
            # now check max prod BF
            maxbf = max(self.BF_P_fits[Pr_i])
            maxbf_idx = np.argmax(self.BF_P_fits[Pr_i])
            if maxbf < bfthreshold:
                comments_MAT[0, i_prod] += ' ! LOW BF: MAX IS {:.2e} AT {} K'.format(maxbf, self.T_VECT[maxbf_idx])
                
   # at this point, call the class WRITE_MECH_CKI from the C_preprocessing module and write the new mech
        # generate the DataFrame with the lines of the mechanism

        self.new_k_to_CKI = preproc.WRITE_MECH_CKI(
            k0_MAT, alpha_MAT, EA_MAT, comments_MAT, new_SPECIES_SERIES, new_i_REAC, new_i_PRODS, new_SPECIES_BIMOL_SERIES)
        # the products are always sinks and there is no isomer equilibrium when fitting
        CKI_lines = self.new_k_to_CKI.MAKE_CKI(1, 0)
        return CKI_lines

    def WRITE_FITS_SINGLEP(self, cwd, P, CKI_lines):
        '''
        This method writes the arrhenius fits at a SINGLE given pressure in the /Arr_fits_CKI folder
        '''

        if os.path.isdir(os.path.join(cwd, 'Arr_fits_CKI')) == False:
            os.mkdir(os.path.join(cwd, 'Arr_fits_CKI'))
        # Now check if the file already exists and remove it
        if os.path.isfile(os.path.join(cwd, 'Arr_fits_CKI', str(P) + '_' + self.REAC + '.CKI')):
            os.remove(os.path.join(cwd, 'Arr_fits_CKI',
                      str(P) + '_' + self.REAC + '.CKI'))
        # generate the file kin.txt and rename it accordingly
        self.new_k_to_CKI.WRITE_CKI(os.path.join(cwd, 'Arr_fits_CKI'), CKI_lines.values)
        os.rename(os.path.join(cwd, 'Arr_fits_CKI', 'kin.txt'), os.path.join(cwd,
                  'Arr_fits_CKI', str(P) + '_' + self.REAC + '.CKI'))

        # check if the selected folder exists
        if os.path.isdir(cwd) == False:
            raise ValueError('The selected folder ' + cwd + ' does not exist')

    def WRITE_PLOG_FITS(self, FITS_DICT, P_VECT, cwd):
        '''
        This method turns the dataframes contained in FITS_DICT into a larger dataframe with the arrhenius fits in the plog form
        The fits are then written in the appropriate folder
        '''
        # generate new dataframes based on the reaction name
        reactions = FITS_DICT[P_VECT[0]]['reac_name'].index

        for reac in reactions:
            # generate the dataframe
            DF_reac = pd.DataFrame(
                index=np.arange(-1, len(P_VECT) -1*int(len(P_VECT)==0)), columns=FITS_DICT[P_VECT[0]].columns)

            # Allocate the first row to the lowest pressure
            DF_reac.loc[-1] = FITS_DICT[P_VECT[0]].loc[reac]
            DF_reac.loc['empty'] = ' '       # add empty line between reactions
            print(DF_reac.loc[-1])
            # if fits are nans or infs: comment those lines
            k0_all = np.concatenate(
                [FITS_DICT[P].loc[reac][['k0']].values for P in P_VECT])
            comms = [FITS_DICT[P].loc[reac]['comments'] for P in P_VECT]
            listinfsnans = [(par == 'nan' or par == 'inf') for par in k0_all]
            if (all(listinfsnans)
                or all(['LOW BF' in comm for comm in comms])):
                flag = '!'
                DF_reac.loc[-1]['reac_name'] = '!' + DF_reac.loc[-1]['reac_name']
            elif (len(P_VECT) > 1 and any(listinfsnans) and not all(listinfsnans)):
                DF_reac.loc[-1][['k0', 'alpha', 'EA']] = [1.0, 1.0, 1.0]
            else:
                flag = ''
            print(DF_reac.loc[-1])
            # Scan the pressures and write in PLOG form
            Pi = 0
            # if the P_VECT has only one pressure: don't write in PLOG form
            if len(P_VECT) > 1:
                for P in P_VECT:
                    # assign the row to the selected pressure
                    DF_reac.loc[Pi] = FITS_DICT[P].loc[reac]
                    # rewrite the line of the reac_name and that of the comments
                    DF_reac.loc[Pi]['reac_name'] = flag + 'PLOG' + '/' + str(P)
                    DF_reac.loc[Pi]['comments'] = '/ ' + \
                        DF_reac.loc[Pi]['comments']
                    # if only some P should be commented, the flag will be ''. update the rxn name
                    if flag == '' and any((par == 'nan' or par == 'inf') for par in DF_reac.loc[Pi][['k0', 'alpha', 'EA']]):
                        DF_reac.loc[Pi]['reac_name'] = '!' + DF_reac.loc[Pi]['reac_name']
                    # update the value of the pressure
                    Pi += 1

            # concatenate values
            if reac == 0:
                plog_fits_array = DF_reac.values
            else:
                plog_fits_array = np.concatenate(
                    (plog_fits_array, DF_reac.values))

        # write the output
        self.new_k_to_CKI.WRITE_CKI(cwd, plog_fits_array)
