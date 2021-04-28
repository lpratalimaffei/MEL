import numpy as np
import os
import pandas as pd
import copy
import re
import shutil
import subprocess
from . import main_flow
from . import C_preprocessing as preproc


def run_preproc(cwd, OS_folder):
    '''
    preprocesses CKI input mechanism and re-writes it as mech with only fw reactions
    writes the bw reactions directly below the fw
    '''
    # copy the thermo to the compilation folder
    thermfile = os.path.join(cwd, 'inp', 'therm.txt')
    thermfile_compile = os.path.join(cwd, 'mech_tocompile', 'therm.txt')
    shutil.copy(thermfile, thermfile_compile)
    # reads kinetics mechanism:
    # preserve the headers (to be used later)
    element_block, species_block, reaction_block = read_CKI(
        os.path.join(cwd, 'inp'), 'kin.CKI')

    # sort the reactions
    rxn_irrev, rxn_rev_HPLIM, rxn_rev_PDEP = sort_CKI(reaction_block)
    rxn_irrev_df = block_to_df(rxn_irrev)
    rxn_irrev = df_to_block(rxn_irrev_df)
    # preserve the only-fw reactions as a separate string block

    # define 2 groups:
    # pressure independent reactions

    # 1. preproc for pressure independent reactions:
    filename = os.path.join(cwd, 'mech_tocompile', 'kin.txt')
    # write the mechanism from the blocks
    # convert the reaction block to dataframe
    rxn_rev_HPLIM_df = block_to_df(rxn_rev_HPLIM)
    write_CKI_blocks(filename, element_block, species_block, rxn_rev_HPLIM)
    rxn_bw_HPLIM_df = run_extract_fittedkin(cwd, OS_folder)

    rxn_irrev_HPLIM = convert_df_FWBW_to_irrev(
        rxn_rev_HPLIM_df, rxn_bw_HPLIM_df)

    # 2. preproc for pressure dependent reactions: at each P, set the rate as the high P limit
    #    Then read the backward reaction from the processed mechanism
    rxn_rev_PDEP_df = block_to_df(rxn_rev_PDEP)
    rxn_bw_PDEP_df = pdep_torev(
        cwd, OS_folder, filename, element_block, species_block, rxn_rev_PDEP_df)

    rxn_irrev_PDEP = convert_df_FWBW_to_irrev(rxn_rev_PDEP_df, rxn_bw_PDEP_df)
    # join all the blocks together
    print(type(rxn_irrev), type(rxn_irrev_HPLIM), type(rxn_irrev_PDEP))
    blocks_all_irrev = rxn_irrev + rxn_irrev_HPLIM + rxn_irrev_PDEP

    # write the final mechanism: preproc_irreversible/kin_irr.CKI
    filename = os.path.join(cwd, 'preproc_irreversible', 'kin.CKI')
    write_CKI_blocks(filename, element_block, species_block, blocks_all_irrev)


def read_CKI(path, name):
    '''
    reads ELEMENTS, SPECIES AND REACTIONS block of CKI mech
    '''

    if os.path.isfile(os.path.join(path, name)) == False:
        print('*Error: the mechanism file does not exist')
    else:
        with open(os.path.join(path, name), mode='r') as mech_block_file:
            mech_block = mech_block_file.readlines()

    mech_block_clean = rem_comments(mech_block, '!')

    start_el = np.where(np.array(
        ['ELEMENTS' in line for line in mech_block_clean], dtype=int) == 1)[0][0]
    start_sp = np.where(
        np.array(['SPEC' in line for line in mech_block_clean], dtype=int) == 1)[0][0]
    start_rxn = np.where(
        np.array(['REAC' in line for line in mech_block_clean], dtype=int) == 1)[0][0]
    # ends = np.where(np.array(['END' in line for line in mech_block_clean], dtype=int) == 1)[0]

    element_block = mech_block_clean[start_el:start_sp]
    species_block = mech_block_clean[start_sp:start_rxn]
    reaction_block = mech_block_clean[start_rxn:]

    return element_block, species_block, reaction_block


def rem_comments(block, cmt_sym):
    '''
    block: string block (list of strings)
    cmt_sym: string with comment sumbol
           - removes lines starting with '!' from the list
           - deletes the portion after the '!'
    '''
    new_block = []

    for line in block:
        if line[0] != cmt_sym and line != '\n':
            new_line = line.split(cmt_sym)[0]
            if line != new_line:
                new_line = new_line + '\n'
            new_block.append(new_line)

    return new_block


def sort_CKI(rxn_block):
    '''
    rxn_block: line block of CKI mechanism
    extracts 3 blocks:
           1. irreversible reactions
           2. HP lim reversible reactions: extracts this block
           3. PLOG reversible reactions: extracts this block
    '''

    rxn_irrev = []
    rxn_rev_HPLIM = []
    rxn_rev_PDEP = []
    iline = 1

    for idx, rxn in enumerate(rxn_block):

        if '=' in rxn:
            if '=>' in rxn and '<=>' not in rxn:
                obj = rxn_irrev
                # irreversible reaction
                obj.append(rxn)

            else:
                # reversible reaction
                if 'PLOG' in rxn_block[idx+iline] or 'LOW' in rxn_block[idx+iline]:
                    obj = rxn_rev_PDEP
                    obj.append(rxn)
                else:
                    # HP lim reaction
                    obj = rxn_rev_HPLIM
                    obj.append(rxn)

        if 'dup' in rxn.lower() or 'LOW' in rxn or 'TROE' in rxn or 'PLOG' in rxn:
            # append that line
            obj.append(rxn)

    return rxn_irrev, rxn_rev_HPLIM, rxn_rev_PDEP


def write_CKI_blocks(filename, element_block, species_block, reaction_block):
    '''
    Write chemkin file from element, species and reactions block (lists of strings)
    '''

    element_str = ''.join(element_block)
    species_str = ''.join(species_block)
    reaction_str = ''.join(reaction_block)
    # check that the reaction string contains the REACTION and END, otherwise add it
    if 'REACTIONS' not in reaction_str:
        reaction_str = 'REACTIONS \n' + reaction_str
    if 'END' not in reaction_str:
        reaction_str += 'END \n'

    total_str = element_str + '\n' + species_str + '\n' + reaction_str

    mechfile = open(filename, "w")
    mechfile.writelines(total_str)
    mechfile.close()


def extract_fittedkin(filename):
    '''
    Extract all the values of the reversible reactions from the fittedkinetics
    '''
    # dataframe
    # index = rxn_name
    # columns = 'HP', 'DUPLICATE'
    # open file

    if os.path.isfile(filename) == False:
        print('*Error: the mechanism file does not exist')
    else:
        with open(filename, mode='r') as fittedkin_file:
            fittedkin = fittedkin_file.readlines()

    # extract reaction block
    start = np.where(
        np.array(['index' in line for line in fittedkin], dtype=int) == 1)[0][1]+2

    rxn_block = fittedkin[start:]

    # generate rxn list and parameter list
    rxn_idx_list = []
    rxn_list = []
    param_list = []

    # order reactions by index so that duplicate reactions are read in the correct order
    for line in rxn_block:
        rxn_idx = int(line.split()[0])
        rxn_name = rxn_name = line.split()[4]
        # convert parameters in appropriate format and append them
        params = np.array(line.split()[1:4], dtype=float)
        # change units to the first parameter if required
        prods = rxn_name.split('=')[1]
        if len(prods.split('+')) >= 2:
            params[0] *= 1.e3

        str_params = '{:.2e}\t{:.2f}\t{:.2f}'.format(
            params[0], params[1], params[2])

        rxn_idx_list.append(rxn_idx)
        rxn_list.append(rxn_name)
        param_list.append(str_params)

    rxn_bw_df_all = pd.DataFrame(np.array(
        [rxn_list, param_list], dtype=object).T,
        index=rxn_idx_list, columns=['RXN_NAME', 'PARAMS'], dtype=object)

    # sort by index
    rxn_bw_df_all = rxn_bw_df_all.sort_index()

    # generate new rxn list and parameter list without duplicates
    rxn_list_nodup = []
    param_list_nodup = []
    dup_param_list = []

    for idx in rxn_bw_df_all.index:
        rxn_name = rxn_bw_df_all.loc[idx]['RXN_NAME']
        str_params = rxn_bw_df_all.loc[idx]['PARAMS']

        # check if duplicate
        if rxn_name not in rxn_list_nodup:
            rxn_list_nodup.append(rxn_name)
            param_list_nodup.append(str_params)
            dup_param_list.append([])
        else:
            # find the corresponding reaction index
            idx_dup = np.where(
                np.array([rxn_name == rxn for rxn in rxn_list_nodup], dtype=int) == 1)[0][0]
            dup_param_list[idx_dup].append(str_params)

    # generate dataframe
    rxn_bw_df = pd.DataFrame(np.array(
        [param_list_nodup, dup_param_list], dtype=object).T,
        index=rxn_list_nodup, columns=['HP', 'DUPLICATE'], dtype=object)

    # add troe and plog to dataframe
    rxn_bw_df['TROE'] = None
    rxn_bw_df['PLOG'] = None

    return rxn_bw_df


def block_to_df(rxn_block):
    '''
    rxn_block: list of strings WITHOUT HEADERS
    NB: MUST HAVE - COMMENTS REMOVED; NO SPACES IN RXN STRINGS
    Convert a given reaction block to dataframe
    df[rxn_name]['HP','PLOG','TROE','DUPLICATE']
    'HP': string with HP params
    'PLOG': dct{HP: string_params, P: string_params}
    'TROE': dct{HP: string_params, LOW: string_params, TROE: string_params}
    'DUPLICATE': same format as the non-empty column
    '''
    # lists for DF generation
    rxn_list = []
    PARAMS = []  # LIST OF ARRAYS: [HP - STR, PLOG -DCT, TROE - DCT, DUP - STR/DCT]

    # remove '\n'

    for idx, line in enumerate(rxn_block):
        # remove '\n'
        iline = 1*(int(idx < len(rxn_block)-1))
        line = line.replace('\n', '')
        if '=' in line:
            # preallocate the array
            PARAM_ARRAY = np.array([None, None, None, []], dtype=object)
            line_split = line.split()
            rxn_name = line_split[0]
            # convert parameters in appropriate format and append them
            params = '\t'.join(line_split[1: 4])
            # PDEP
            if 'PLOG' in rxn_block[idx+iline]:
                plog_dct = {}
                PARAM_ARRAY[1] = plog_dct
                plog_dct['HP'] = params
                while 'PLOG' in rxn_block[idx+iline]:
                    line2 = rxn_block[idx+iline]
                    line2 = line2.replace('\n', '')
                    all_params = line2.split()
                    plog_dct[all_params[2]] = '\t'.join(all_params[3:6])
                    iline += 1
                    if idx+iline == len(rxn_block):
                        break

            elif 'LOW' in rxn_block[idx+iline]:
                troe_dct = {}
                PARAM_ARRAY[2] = troe_dct
                troe_dct['HP'] = params
                line2 = rxn_block[idx+iline]
                line2 = line2.replace('\n', '')
                iline += 1
                line3 = rxn_block[idx+iline]
                line3 = line3.replace('\n', '')
                troe_dct['LOW'] = line2.split('/')[1]
                troe_dct['TROE'] = line3

            else:
                # HP lim reaction
                PARAM_ARRAY[0] = params

            # DUPLICATE CHECK
            if rxn_name not in rxn_list:
                rxn_list.append(rxn_name)
                PARAMS.append(PARAM_ARRAY)

            else:
                # find the corresponding reaction index
                idx = np.where(
                    np.array([rxn_name == rxn for rxn in rxn_list], dtype=int) == 1)[0][0]
                # select the nonempty param array
                check_notempty = np.array([bool(i)
                                           for i in PARAM_ARRAY], dtype=int)
                notempty = np.where(check_notempty == 1)[0][0]
                PARAMS[idx][3].append(PARAM_ARRAY[notempty])

    # generate dataframe
    rxn_df = pd.DataFrame(PARAMS, index=rxn_list, columns=[
                          'HP', 'PLOG', 'TROE', 'DUPLICATE'])

    return rxn_df


def df_to_block(rxn_df):
    '''
    rxn_df[rxn_name]['HP','PLOG','TROE','DUPLICATE']
    'HP': string with HP params
    'PLOG': dct{HP: string_params, P: string_params}
    'TROE': dct{HP: string_params, LOW: string_params, TROE: string_params}
    'DUPLICATE': same format as the non-empty column
    rxn_block:  strings WITHOUT HEADERS
    '''

    rxn_block = ''

    for rxn in rxn_df.index:

        # derive params
        check_notempty = np.array(
            [bool(i) for i in rxn_df.loc[rxn].values[0:-1]], dtype=int)
        notempty = np.where(check_notempty == 1)[0][0]
        params = rxn_df.loc[rxn].values[notempty]

        rxn_block += rxn_params_to_str(rxn, params)

        if rxn_df.loc[rxn]['DUPLICATE']:
            rxn_block += 'DUPLICATE\n'
            for params in rxn_df.loc[rxn]['DUPLICATE']:
                rxn_block += rxn_params_to_str(rxn, params)
                rxn_block += 'DUPLICATE\n'

    return rxn_block


def rxn_params_to_str(rxn, params):
    '''
    params = string: HP
    params = {HP: , 'P1': , 'P2': , ..}
    params = {'HP', 'LOW', 'TROE'}
    returns string to write the reaction
    '''
    rxn_str = rxn + '\t'

    if isinstance(params, str):
        # HP limit
        rxn_str += params + '\n'

    elif isinstance(params, dict):
        rxn_str += params['HP'] + '\n'
        if 'TROE' in list(params.keys()):
            # troe params
            rxn_str += 'LOW/\t' + params['LOW'] + '\t/\n'
            rxn_str += params['TROE'] + '\n'
        else:
            # plog params
            for P in list(params.keys())[1:]:
                rxn_str += 'PLOG/\t' + P + '\t' + params[P] + '\t/\n'

    return rxn_str


def run_extract_fittedkin(cwd, OS_folder):
    '''
    run preprocessor and extract bw HP kinetics dataframe
    '''
    OS_exe = main_flow.get_OS()
    exec0 = main_flow.get_libpath()
    preproc_exe = os.path.join('"' + OS_folder, "OpenSMOKEpp_CHEMKIN_PreProcessor." + OS_exe + '"')
    input_preproc = os.path.join(os.path.join(".", "mech_tocompile", "input_preproc.dic"))    
    output_preproc = os.path.join(".", "mech_tocompile", "preproc_output.txt")
    toexecute = exec0 + preproc_exe + " --input " + input_preproc + ">" + output_preproc
    print('compiling mech ...'), subprocess.run(toexecute, shell=True)

    # derive the dataframe for the corresponding backward reactions
    rxn_bw_HPLIM_df = extract_fittedkin(os.path.join(
        cwd, 'mech_tocompile', 'kinetics', 'Reaction_FittedKinetics.out'))

    return rxn_bw_HPLIM_df


def pdep_torev(cwd, OS_folder, filename, element_block, species_block, rxn_rev_PDEP_df):
    '''
    Takes a dataframe of pressure dependent reactions
    returns a dataframe with the backward rate constants
    reaction names are UNCHANGED
    '''

    # preallocate dataframe
    rxn_bw_PDEP_df = pd.DataFrame(index=rxn_rev_PDEP_df.index, columns=rxn_rev_PDEP_df.columns, dtype=object)
    rxn_bw_PDEP_df[rxn_rev_PDEP_df.columns] = np.array([None, None, None, []], dtype=object)
    # derive the reaction block
    rxn_block_PDEP = df_to_block(rxn_rev_PDEP_df)
    write_CKI_blocks(filename, element_block, species_block, rxn_block_PDEP)
    # run preprocessor

    rxn_bw_HP_df = run_extract_fittedkin(cwd, OS_folder)

    # for each reaction: derive pressure dependent backward rate constants
    for rxn in rxn_bw_PDEP_df.index:

        # different runs for PLOG or TROE: the only differences are
        # the column you look at : PLOG or TROE
        # the N of loops: for TROE, you just have HP and LOW, not the last one

        if rxn_rev_PDEP_df.loc[rxn]['PLOG']:
            col = 'PLOG'
            rxn_bw_PDEP_df.loc[rxn][col] = {}
        elif rxn_rev_PDEP_df.loc[rxn]['TROE']:
            col = 'TROE'
            rxn_bw_PDEP_df.loc[rxn][col] = {}
            rxn_bw_PDEP_df.loc[rxn][col]['TROE'] = rxn_rev_PDEP_df.loc[rxn][col]['TROE']

        # rewrite the HP limit value
        hp_params = rxn_bw_HP_df.loc[rxn]['HP']
        rxn_bw_PDEP_df.loc[rxn][col]['HP'] = hp_params
        # run the preprocessor for each P: rewrite the parameters
        params = rxn_rev_PDEP_df.loc[rxn][col]
        if rxn_rev_PDEP_df.loc[rxn]['PLOG']:
            P_list = list(params.keys())[1:]
        elif rxn_rev_PDEP_df.loc[rxn]['TROE']:
            P_list = ['LOW']

        for P in P_list:
            # assign params for writing
            params_write = copy.deepcopy(params)
            for key in params.keys():
                if key != 'TROE':
                    params_write[key] = params[P]
            rxn_block_single = rxn_params_to_str(rxn, params_write)
            write_CKI_blocks(filename, element_block,
                            species_block, rxn_block_single)
            # run preprocessor and extract kinetics
            rxn_bw_single_HP_df = run_extract_fittedkin(cwd, OS_folder)
            # rewrite the params
            rxn_bw_PDEP_df.loc[rxn][col][P] = rxn_bw_single_HP_df.loc[rxn]['HP']

        print(rxn, rxn_rev_PDEP_df.loc[rxn][col], rxn_bw_PDEP_df.loc[rxn][col])
        # for duplicate reactions:

        if rxn_rev_PDEP_df.loc[rxn]['DUPLICATE']:
            rxn_bw_PDEP_df.loc[rxn]['DUPLICATE'] = []
            # run over all the i:
            for i in rxn_rev_PDEP_df.loc[rxn]['DUPLICATE']:
                param_bw = {}
                hp_params = rxn_bw_HP_df.loc[rxn]['DUPLICATE'][i]
                param_bw['HP'] = hp_params
                params = rxn_rev_PDEP_df.loc[rxn]['DUPLICATE'][i]

                if 'TROE' in list(params.keys()):
                    P_list = ['LOW']
                    param_bw['TROE'] = rxn_rev_PDEP_df.loc[rxn][col]['TROE']
                else:
                    P_list = list(params.keys())[1:]

                for P in P_list:
                    params_write = copy.deepcopy(params)
                    for key in params_write.keys():
                        if key != 'TROE':
                            params_write[key] = params[P]
                    rxn_block_single = rxn_params_to_str(rxn, params_write)
                    write_CKI_blocks(filename, element_block,
                                     species_block, rxn_block_single)
                    # run preprocessor and extract kinetics
                    rxn_bw_single_HP_df = run_extract_fittedkin(cwd, OS_folder)
                    # rewrite the params
                    param_bw[P] = rxn_bw_single_HP_df.loc[rxn]['HP']

                rxn_bw_PDEP_df.loc[rxn]['DUPLICATE'].append(param_bw)

    return rxn_bw_PDEP_df


def convert_df_FWBW_to_irrev(rxn_fw_df, rxn_bw_df):
    """
    rxn_fw_df: dataframe with forward rate constants
    rxn_bw_df: dataframe with backward rate constants
    a. give proper indices to reactions
    b. rename reactions as irreversible
    c. write new dataframe
    d. write ordered reaction block
    """
    rxn_fw = []
    rxn_bw = []
    # replace the '=' sign with '=>'
    for rxn in list(rxn_fw_df.index):
        if ('=' in rxn and '=>' not in rxn):
            rxn = rxn.replace('=', '=>')
        elif ('<=>' in rxn):
            rxn = rxn.replace('<=>', '=>')
        rcts, prds = rxn.split('=>')
        rxn_fw.append(rxn)
        rxn_bw.append(prds + '=>' + rcts)
    rxn_all = rxn_fw + rxn_bw
    # order reactions in the dataframes to write them consistently
    rxn_fw_df['order'] = np.arange(0, 2*len(rxn_fw), 2)
    rxn_bw_df['order'] = np.arange(1, 2*len(rxn_fw)+1, 2)
    # new_dataframe
    rxn_irrev_df = pd.DataFrame(np.concatenate((
        rxn_fw_df[rxn_fw_df.columns], rxn_bw_df[rxn_fw_df.columns])), index=rxn_all, columns=rxn_fw_df.columns)
    # sort by 'order'
    rxn_irrev_df = rxn_irrev_df.sort_values(by='order')
    # drop the order
    rxn_irrev_df = rxn_irrev_df.drop('order', axis=1)
    print(rxn_irrev_df)

    # turn to block
    rxn_irrev_block = df_to_block(rxn_irrev_df)

    return rxn_irrev_block
