import numpy as np
import os
import pandas as pd
import copy
import re
from . import C_preprocessing as preproc

def run_preproc(cwd):
       '''
       preprocesses CKI input mechanism and re-writes it as mech with only fw reactions
       writes the bw reactions directly below the fw
       '''

       # reads kinetics mechanism:
       # preserve the headers (to be used later)
       element_block, species_block, reaction_block = read_CKI(os.path.join(cwd, 'inp'), 'kin.CKI')
       
       # sort the reactions
       rxn_irrev, rxn_rev_HPLIM, rxn_rev_PLOG = sort_CKI(reaction_block)
       print(rxn_irrev, rxn_rev_HPLIM, rxn_rev_PLOG)
       # preserve the only-fw reactions as a separate string block

       # define 2 groups:
       # pressure independent reactions
       # pressure dependent reactions: preproc to be done at each pressure for each simulation

       # 1. preproc for pressure independent reactions: read them from the fitted kinetics

       # 2. preproc for pressure dependent reactions: at each P, set the rate as the high P limit
       #    Then read the backward reaction from the processed mechanism


       # join all the blocks together

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
       rxn_rev_PLOG = []


       for idx, rxn in enumerate(rxn_block):

              iline = 1

              if '=' in rxn:
                     if '=>' in rxn and '<=>' not in rxn:
                            obj = rxn_irrev
                            # irreversible reaction
                            obj.append(rxn)

                     else:
                            # reversible reaction
                            if 'PLOG' in rxn_block[idx+iline]:
                                   obj = rxn_rev_PLOG
                                   obj.append(rxn)
                            else:
                                   # HP lim reaction
                                   obj = rxn_rev_HPLIM
                                   obj.append(rxn)
 
                     # checks: PLOG
                     check_plog = 0

                     while check_plog == 0:
                            if 'PLOG' in rxn_block[idx+iline]:
                                   obj.append(rxn_block[idx+iline])
                                   iline += 1
                            else:
                                   check_plog = 1


              if 'dup' in rxn.lower():
                     # append that line
                     obj.append(rxn)

       return rxn_irrev, rxn_rev_HPLIM, rxn_rev_PLOG