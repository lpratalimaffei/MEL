import numpy as np
import os
import pandas as pd
import copy
import re
import shutil
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
       element_block, species_block, reaction_block = read_CKI(os.path.join(cwd, 'inp'), 'kin.CKI')
       
       # sort the reactions
       rxn_irrev, rxn_rev_HPLIM, rxn_rev_PLOG = sort_CKI(reaction_block)
       
       # preserve the only-fw reactions as a separate string block

       # define 2 groups:
       # pressure independent reactions

       # 1. preproc for pressure independent reactions: 
       # write the mechanism from the blocks
       filename = os.path.join(cwd, 'mech_tocompile', 'kin.txt')
       write_CKI_blocks(filename, element_block, species_block, rxn_rev_HPLIM)
       toexecute = '"' + OS_folder + "\OpenSMOKEpp_CHEMKIN_PreProcessor.exe" + '"' + " --input .\mech_tocompile\input_preproc.dic >.\mech_tocompile\preproc_output.txt"
       print('compiling mech ...'),os.system(toexecute)

       # derive the dataframe for the corresponding backward reactions
       rxn_bw_df = extract_fittedkin(os.path.join(cwd, 'mech_tocompile', 'kinetics', 'Reaction_FittedKinetics.out'))
       print(rxn_bw_df)
       print(rxn_bw_df['DUPLICATE']['C5H10OOH2-4O2=C5H91-2,4OOH'])

       # convert the reaction block to dataframe
       rxn_rev_HPLIM_df = block_to_df(rxn_rev_HPLIM)

       # 2. preproc for pressure dependent reactions: at each P, set the rate as the high P limit
       #    Then read the backward reaction from the processed mechanism


       # join all the blocks together

       # remove the mechanism after execution
       if os.path.isfile(filename):
              os.remove(filename)

       # write the final mechanism: preproc_irreversible/kin_irr.CKI


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

       start_el = np.where(np.array(['ELEMENTS' in line for line in mech_block_clean], dtype=int) == 1)[0][0]
       start_sp = np.where(np.array(['SPEC' in line for line in mech_block_clean], dtype=int) == 1)[0][0]
       start_rxn = np.where(np.array(['REAC' in line for line in mech_block_clean], dtype=int) == 1)[0][0]
       #ends = np.where(np.array(['END' in line for line in mech_block_clean], dtype=int) == 1)[0]

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
       start = np.where(np.array(['index' in line for line in fittedkin], dtype=int) == 1)[0][1]+2

       rxn_block = fittedkin[start:]

       # generate rxn list and parameter list
       rxn_list = []
       param_list = []
       dup_param_list = []
       
       for line in (rxn_block):
              line_split = line.split()
              rxn_name = line_split[4]
              # convert parameters in appropriate format and append them
              params = np.array(line_split[1:4], dtype=float)
              # change units to the first parameter
              params[0] *= 1.e6
              str_params = '{:.2e}     {:.2f}     {:.2f}'.format(params[0], params[1], params[2])

              # check if duplicate
              if rxn_name not in rxn_list:
                     rxn_list.append(rxn_name)
                     dup_param_list.append(None)
                     # append to list
                     param_list.append(str_params)
              else:
                     # find the corresponding reaction index
                     idx = np.where(np.array([rxn_name == rxn for rxn in rxn_list], dtype=int) == 1)[0][0]
                     dup_param_list[idx] = str_params

       # generate dataframe
       rxn_bw_df = pd.DataFrame(np.array([param_list, dup_param_list]).T, index=rxn_list, columns=['HP', 'DUPLICATE'])

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
       PARAMS = [] # LIST OF ARRAYS: [HP - STR, PLOG -DCT, TROE - DCT, DUP - STR/DCT]

       # remove '\n'

       for idx, line in enumerate(rxn_block):
              # remove '\n'
              iline = 1
              line = line.replace('\n','')
              if '=' in rxn:
                     # preallocate the array
                     PARAM_ARRAY = np.array([None, None, None, None], dtype=object)
                     line_split = line.split()
                     rxn_name = line_split[0]
                     # convert parameters in appropriate format and append them
                     params = line_split[1:4]
                     # PDEP
                     flag_plog = 0
                     while flag_plog = 0:
                            plog_dct = {}
                            if 'PLOG' in rxn_block[idx+iline]:
                                   line2 = rxn_block[idx+iline]
                                   line2 = line2.replace('\n','')
                                   all_params = line2.split('/')
                                   plog_dct[all_params[0]] = all_params[1:]
                            
                                   iline += 1
                            else:
                                   flag_plog = 1
                     troe_dct = {}
                     if 'LOW' in rxn_block[idx+1]:
                            line2 = rxn_block[idx+1]
                            line2 = line2.replace('\n','')


                     # DUPLICATE CHECK
                     if rxn_name not in rxn_list:
                            rxn_list.append(rxn_name)
                            dup_param_list.append(None)
                            # append to list
                            param_list.append(str_params)
                     else:
                            # find the corresponding reaction index
                            idx = np.where(np.array([rxn_name == rxn for rxn in rxn_list], dtype=int) == 1)[0][0]
                            dup_param_list[idx] = str_params

 


       # generate dataframe
       rxn_bw_df = pd.DataFrame(np.array([param_list, dup_param_list]).T, index=rxn_list, columns=['HP', 'DUPLICATE'])
