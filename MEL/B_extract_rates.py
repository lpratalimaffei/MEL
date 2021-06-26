import sys
import numpy as np
import os
import pandas as pd
import re
import copy
# extract names, temperature, pressures and make folder
# function for mech reading


def readmechanism(input_type, cwd):
    '''
    method to read the mechanism file depending on the input type
    '''
    if input_type == 'MESS':
        # extract parameters
        P_LIST, T_LIST, species_names, species_names_bimol_frag2 = data_names_mess(
            os.path.join(cwd, 'inp'))
        # extract matrix of rate constants
        rates = MATRIX(os.path.join(cwd, 'inp'), P_LIST, T_LIST, species_names)
        mech_dict = dict(zip(['P_VECT_MESS', 'T_VECT_MESS', 'SPECIES', 'SPECIES_BIMOL', 'rates'], [
                         P_LIST, T_LIST, species_names, species_names_bimol_frag2, rates]))

    elif input_type == 'CKI':
        species_names, species_names_bimol_frag2 = data_names_CKI(
            os.path.join(cwd, 'inp'))
        mech_dict = dict(zip(['SPECIES', 'SPECIES_BIMOL'], [
                         species_names, species_names_bimol_frag2]))

    return mech_dict

########################### extract and process MESS type mechanism #####################


def data_names_mess(cwd):
    """
    Extract from the input file me_ktp.inp useful data:
    - list of pressures
    - list of temperatures
    - list of species (names)
    """
    species_names_unimol = np.array([], dtype='<U16')
    species_names_bimol = np.array([], dtype='<U16')
    species_names_unimol_frag2 = np.array([], dtype='<U16')
    species_names_bimol_frag2 = np.array([], dtype='<U16')
    look_for_species = 0
    look_for_bimol_fragment = 0
    bad_wellwrds = ['WellDepth', 'WellCutoff', 'WellExtension',
                    'WellReductionThreshold', 'WellPartitionMethod', 'WellProjectionThreshold']
    with open(os.path.join(cwd, 'me_ktp.inp')) as myfile:
        for line in myfile:
            #  do not read comments
            if line[0] == '!':
                line = ''  # empty line
            elif len(line.split('!')) > 1:  # remove commented part
                line = line.split('!')[0]
            # if you enter the model section: start looking for species
            if line.find('Model') != -1 and line.find('ModelEnergyLimit') == -1:
                look_for_species = 1

            if line.find('PressureList') != -1:
                # verifica di non poter usare semplicemente line.split
                pressures = [x.strip() for x in line.split()]
                del pressures[0]
                # print(pressures)
                # print(len(pressures))
            if line.find('TemperatureList') != -1:
                temperatures = [x.strip() for x in line.split()]
                del temperatures[0]

            if (line.find('Well') != -1 and all(line.find(bad) == -1 for bad in bad_wellwrds)) and look_for_species == 1:
                full_line = [x.strip() for x in line.split()]
                species_names_unimol = np.append(
                    species_names_unimol, full_line[1])
                species_names_unimol_frag2 = np.append(
                    species_names_unimol_frag2, '')

            if (line.find('Bimolecular')) != -1 and look_for_species == 1:
                look_for_bimol_fragment = 1
                full_line = [x.strip() for x in line.split()]
                species_names_bimol = np.append(
                    species_names_bimol, full_line[1])

            if (line.find('Fragment')) != -1 and look_for_bimol_fragment > 0:
                look_for_bimol_fragment += 1
                if look_for_bimol_fragment == 3:
                    full_line = [x.strip() for x in line.split()]
                    species_names_bimol_frag2 = np.append(
                        species_names_bimol_frag2, full_line[1])
                    look_for_bimol_fragment = 0

    myfile.close()

    # write files
    P_LIST = np.unique(np.array(pressures, dtype=np.float32))
    T_LIST = np.unique(np.array(temperatures, dtype=np.int16))

    species_names = np.append(species_names_unimol, species_names_bimol)
    species_names_frag2 = np.append(
        species_names_unimol_frag2, species_names_bimol_frag2)
    print(species_names, species_names_frag2, '\n')
    # check that bimol fragments have different names
    if len(list(set(species_names_bimol_frag2))) != len(list(species_names_bimol_frag2)):
        print('*Warning: some bimol fragments share the same names. check that they are isomers')

    return P_LIST, T_LIST, species_names, species_names_frag2


def MATRIX(cwd, P_LIST, T_LIST, species_names):
    """
    Extract from rate.out all the rate constants in the form of a list
    """
    # pre-allocation
    matrix_list = []
    capture_list = []
    # define checks for temperature and pressure and read the file
    # len works on b
    check_P = len(T_LIST)*len(species_names)*len(P_LIST)
    # define current checks
    check_P_curr = 0
    check_list = 0
    with open(os.path.join(cwd, 'rate.out')) as myfile:
        for line in myfile:
            # find the 'Temperature-Species rate tables to extract the rates
            if line.find('Temperature-Species Rate Tables:') != -1:
                check_list = 1
            if ((check_list == 1) and (check_P_curr < check_P)):
                # add the check on 'Pressure' in case the values of temperature and pressure are accidentally the same
                if any(line.find(T) != -1 for T in np.array(T_LIST, dtype=str)) and (line.find('Pressure') == -1):
                    # print(line)
                    check_P_curr += 1
                    rates = [x.strip() for x in line.split()]
                    # replace '***' values with 0
                    rates = [x.replace('***', '0') for x in rates]
                    matrix_list.append(rates[1:-2])  # -2 excluded
                    capture_list.append(float(rates[-1]))
            if line.find('Temperature-Pressure Rate Tables:') != -1:
                check_list = 0  # don't read the file anylonger
    myfile.close()

    matrix_float = np.array(matrix_list, dtype=np.float64)

    # remove negative values from the matrix
    n_T = len(T_LIST)
    n_P = len(P_LIST)
    warnings_neg = ''  # generate list of warnings for negative values
    for ii, row in enumerate(matrix_float):
        mask_neg = np.where(row < 0)
        mask_toohigh = np.where(row > capture_list[ii])
        row[mask_toohigh] = 0
        for mask_neg_i in mask_neg[0]:
            row[mask_neg_i] = 0
            R = int(ii/n_T/n_P)
            P = int((ii-R*n_T*n_P)/n_T)
            T = ii-R*n_T*n_P-P*n_T
            warnings_neg = warnings_neg + \
                'removed negative k_{R} at {T} K and {P} atm, be careful \n'.format(
                    R=species_names[R], T=T_LIST[T], P=P_LIST[P])

    np.savetxt('warnings_negval_messrates.txt', [warnings_neg], fmt='%s')

    return matrix_float


def MATRIX_TP(T, P, T_LIST, P_LIST, species_names, matrix_float):
    """
    Extract square matrix of k_ij at the selected temperature and pressure
    T,P are expected to be numbers, either floating or integers, to be compared with T,P in the lists
    """
    # transform input types
    if not isinstance(T, int):
        try:
            T = int(T)
        except:
            print('input T not convertible to number')
            return None

    if not isinstance(P, float):
        try:
            P = float(P)
        except:
            print('input P not convertible to number')
            return None

    n_T = len(T_LIST)
    n_P = len(P_LIST)
    # preallocate the matrix
    n_species = len(species_names)
    mat_TP = np.zeros((n_species, n_species))

    # reconstruct indices
    # set as option the exception case: if you put a value not present
    # index [0][0] to save just the index itself
    P_index = np.where((P == P_LIST))[0][0]
    T_index = np.where((T == T_LIST))[0][0]

    for ii in range(0, n_species):
        # identify the string to place in the row: ki->prods
        rates_row = matrix_float[ii*(n_P)*(n_T)+P_index*(n_T)+T_index, :]
        # use list comprehension: kij = rate from reactant i to product j
        col_indices = np.array([jj != ii for jj in range(0, n_species)])
        mat_TP[ii, col_indices] = rates_row

    # this returns the matrix at a certain temperature and pressure
    return mat_TP


def REAC_P(P, reac, P_LIST, T_LIST, species_names, matrix_float):
    """
    Method: only needed for MESS input. it extracts from the matrix of rate constants the reactivity of the selected reactant at a certain pressure
    output: matrix [n_T*n_species-1]|P
    useful if: you need to check the reactivity of the reactant in the full range of temperature
    """
    n_T = len(T_LIST)
    n_P = len(P_LIST)
    # select the reactant
    if not isinstance(reac, str):
        print('input is not a string. please insert "reac" as string')
        return None
        # if the input is not a string: show an error

    elif sum(np.array([reac == r_list for r_list in species_names])) == 0:
        print('selected reactant not found in the list of species. select another reactant')
        return None
        # if the reactant is not in the list: show an error

    else:
        try:
            float(P)
        except ValueError as e:
            print('P not convertible to number, error: ' + str(e))
            return None

        # derive the reactivity matrix at all the different temperatures
        P_index = np.where((P == P_LIST))[0][0]
        reac_index = np.array(
            [reac == r_list for r_list in species_names], dtype=int)
        ii_reac = np.where(reac_index == 1)[0][0]  # index of the reactant
        ii_in = ii_reac*(n_P)*(n_T)+P_index*(n_T)
        rates_reac = matrix_float[ii_in:ii_in+n_T, :]

        return rates_reac

# extract and process CHEMKIN type mechanism ########################à


def data_names_CKI(cwd):
    '''
    extract species names of CKI input mechanism
    '''

    # names of all the primary species
    species_names = np.array([], dtype='<U32')
    # names of the "second" reactant found in bimolecular reaction channels
    species_names_bimol = np.array([], dtype='<U32')
    # check in the file if you reading the part with all the reactions
    check_reactions = 0

    with open(os.path.join(cwd, 'kin.CKI')) as myfile:
        for line in myfile:
            if line.find('REACTIONS') != -1:
                check_reactions += 1
            # read the lines with the reaction rates;
            # only irreversible reactions are considered ('=>')
            # lines starting with a comment ('!') are not considered
            if check_reactions == 1 and line.find('=>') != -1 and line.strip()[0] != '!':
                line = line.split('!')[0]  # remove the comments
                REACS = [x.strip() for x in line.split('=>')][0]
                REACS = [x.strip() for x in REACS.split('+')]
                rest_ofline = [x.strip() for x in line.split('=>')][1]
                # the arrhenius parameters will be the last three elements in the split
                # ARR_PAR = rest_ofline.split()[-3:]
                PRODS = ''.join(rest_ofline.split()[:-3])
                PRODS = PRODS.split('+')
                reacting_species = [REACS, PRODS]
                # ALLOCATE THE SPECIES INTO THE ARRAYS
                # if species are not into the array of species: append them
                for x in reacting_species:
                    # IF LEN IS 1 AND REACTION IS LIKE 2A=>PRODS OR REACS=>2B: RECOGNIZE SPECIES B
                    # and add another product
                    if x[0][0] == '2':
                        x[0] = x[0][1:]
                        x.append(x[0])
                    if np.array([x[0] == SP for SP in species_names]).any() != True:
                        species_names = np.append(species_names, x[0])
                        # if the species is bimolecular, also append the second species
                        if len(x) == 1:
                            species_names_bimol = np.append(
                                species_names_bimol, '')
                        elif len(x) == 2:
                            species_names_bimol = np.append(
                                species_names_bimol, x[1])
                        else:
                            raise ValueError(
                                'Wrong number of elements when reading the reaction: line ' + line)

    return species_names, species_names_bimol


def copy_CKI_processed(oldpath, newpath, PRODSINKS, ISOM_EQUIL, REAC, PRODS):
    '''
    IN THIS METHOD, THE MECHANISM IS COPIED TO NEWPATH FOLDER AFTER PREPROCESSING
    PRODSINKS = 1: THE PRODUCTS ARE SET AS IRREVERSIBLE SINKS, SO THE LINES ARE COMMENTED WITH !
    ISOM_EQUIL = 1: ALL BIMOLECULAR REACTIONS ARE DELETED,SO THAT EQUILIBRIUM WITHIN SETS OF ISOMERS IS REACHED.
                  => ALL BIMOLECULAR REACTIONS COMMENTED WITH !
    '''

    if os.path.isdir(newpath) == False:
        raise RuntimeError(
            'The destination folder for the new mech does not exist ')
    else:
        with open(os.path.join(oldpath, 'kin.CKI'), mode='r') as mech_orig_file:
            mech_orig = mech_orig_file.readlines()

    newfile = copy.deepcopy(mech_orig)

    for idx, row in enumerate(newfile):

        if ISOM_EQUIL == 1 and row.find('=>') != -1 and row.strip()[0] != '!':
            # check first if you only want isomer equilibrium: so you delete all bimolecular reactions
            reactant = [x.strip() for x in row.split('=>')][0]
            reactant = [x.strip() for x in reactant.split('+')][0]
            product = [x.strip() for x in row.split('=>')][1]
            # only bimolecular products will be meaningful
            product = [x.strip() for x in product.split('+')][0]
            product = product.split()[0]
            # delete all the reactions involving species other than the reactants
            if np.array([reactant == REAC]).any() and np.array([product == REAC]).any():
                delete = 'NO'
            else:
                delete = 'YES'

        elif PRODSINKS == 1 and row.find('=>') != -1 and row.strip()[0] != '!':
            # if there is no isom_equil, you did not delete all bimol. reactions:
            # so if you set product as sinks, then you have to delete those lines.
            # NB product sinks are incompatible with isom_equil
            # extract the first reactant of each row
            reactant = [x.strip() for x in row.split('=>')][0]
            reactant = [x.strip() for x in reactant.split('+')][0]
            # if the reactant of the line is in the list of products, comment the line
            if np.array([reactant == np.array(PRODS)]).any():
                delete = 'YES'
            else:
                delete = 'NO'

        else:  # any other case (read a reaction or an empty line with no need of deleting anything)
            delete = 'NO'

        # ON THE BASIS OF THE CONDITIONS ABOVE: DELETE THE REACTION OR NOT
        if delete == 'YES':
            # comment the line
            newfile[idx] = '!' + row
            # if duplicate reaction: comment also that
            # the "lower" notation is to trasfer all to lower cases so that you make the search independent of the font
            if newfile[idx+1].lower().find('DUPLICATE'.lower()) != -1:
                newfile[idx+1] = '!' + newfile[idx+1]
            # if PLOG: comment also all the lines below until you don't find PLOG anymore
            check_plog = 0
            iline = 1
            while check_plog == 0:
                if newfile[idx+iline].find('PLOG') != -1:
                    newfile[idx+iline] = '!' + newfile[idx+iline]
                    iline += 1
                else:
                    check_plog = 1
        elif delete == 'NO':
            # copy the line as it is
            newfile[idx] = row

    # remove the file if it exists and write the new one
    if os.path.isfile(os.path.join(newpath, 'kin.txt')):
        os.remove(os.path.join(newpath, 'kin.txt'))

    with open(os.path.join(newpath, 'kin.txt'), mode='x') as inp:
        inp.writelines(newfile)
