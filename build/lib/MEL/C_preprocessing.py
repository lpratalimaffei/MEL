import numpy as np
import os
import pandas as pd
import copy
import re

def CHECK_REACTIVITY(kR_j,T_VECT,T_VECT_MESS):
       '''
       CHECK IF THE SELECTED REACTANT IS ACTIVE UNTIL THE SELECTED TEMPERATURE
       '''
       #print(kR_j,np.sum(kR_j,axis=1))
       REACTIVITY = np.sum(kR_j,axis=1) # sum by row (each temperature)
       mask_REAC = np.where(sum(np.array([Ti==T_VECT_MESS for Ti in T_VECT])))
       REACTIVITY = REACTIVITY[mask_REAC] # now you have the reactivity in the selected range of temperature
       mask_REAC_zero = np.where(REACTIVITY==0)
       mask_REAC_nzero = np.where(REACTIVITY!=0)
       if len(mask_REAC_zero[0])==len(REACTIVITY):
              Err = ('Error: reactant not active at any of the selected temperatures')
              flag_exit = 1
       elif len(mask_REAC_zero[0]>0):
              # reduce the range of T accordingly
              flag_exit = 0
              Err = (('Warning: temperature range reduced; reactant inactive at {T} K').format(T=T_VECT[mask_REAC_zero]))
              T_VECT = T_VECT[mask_REAC_nzero]
       else:
              flag_exit = 0
              Err = ''

       return flag_exit,Err,T_VECT

def INITIALMOLES(REAC,SPECIES_BIMOL_SERIES,N_INIT):
       '''
       this method returns the number of initial moles depending on whether the reaction is bimolecular (init moles = 1e-5),
       unimolecular or self-reaction (2A->B) (ninit = 1)
       '''
       # for lumped reactants: the first species is chosen to be representative of the set. print a warning
       if isinstance(REAC,np.ndarray):
              REAC = REAC[0]
              print('warning: the first lumped reactant is chosen to be representative of the full set')
       # CHECK IF THE REACTION IS BIMOLECULAR AND SET A FLAG
       if SPECIES_BIMOL_SERIES[REAC] != '' and SPECIES_BIMOL_SERIES[REAC] != REAC:
              N_INIT_REAC =  N_INIT*0.00001
       else:
              N_INIT_REAC =  N_INIT

       return N_INIT_REAC

def BRANCHING_LUMPEDREAC(cwd,REACLUMPED,REAC,T_VECT,P,STOICH,ISOM_EQUIL):
       '''
       This method extracts at every pressure the branching fractions of the inlet reactant species
       Warnings may be given by:
       - Comparison between T and T_vect: if T_VECT is not present in the branching fractions read, the range of T is reduced and a warning is given
       - Comparison between the reactants and the branchings found: if the species don't correspond, the function is stopped
       - If the file is not found: the program is stopped
       '''
       warnings = ''
       fld = cwd + '/Branchings_' + ''.join(STOICH) +'/'+ REACLUMPED + '_' + str(P) + 'atm.txt'
       #if os.path.isfile(fld) == False and ISOM_EQUIL != 1:
       #              raise ValueError('The file ' + fld + 'does not exist: impossible to retrieve the fractions of the reactant pool ')
       #elif os.path.isfile(fld) == False and ISOM_EQUIL == 1:
       if os.path.isfile(fld) == False :
              warnings = warnings + ' warning: Branching fractions not found: initial BFs randomly generated from normal distribution \n'
              # you need the equilibrium composition: start from a random composition of isomers
              rand_comp = abs(np.random.randn(len(T_VECT),len(REAC)))
              rand_comp = rand_comp/np.sum(rand_comp,axis=1).reshape(len(T_VECT),1)
              # divide by the sum of each column to normalize
              BR_L_REAC = pd.DataFrame(rand_comp,index=np.array(T_VECT),columns=np.array(REAC))
              
       else:
              BR_L_Pi = np.genfromtxt(fld,skip_header=1)
              labels = np.genfromtxt(fld,dtype=str,max_rows=1)
              # generate a dataframe with the info you read
              BR_L_REAC = pd.DataFrame(BR_L_Pi[:,1:],index=BR_L_Pi[:,0],columns=labels[1:])
              print(BR_L_Pi[:,0],print(labels[1:]))
              # 1. check that the names in the labels and REAC correspond, independent of the order
              if np.array([[RR == REAC for RR in BR_L_REAC.columns][ii].any() for ii in range(0,len(BR_L_REAC.columns))]).all() != True:
                     raise ValueError('The species read in the file with the branchings and the selected reactant do not match, please check ')

              # 2. check that the temperature range is larger or equal to T_VECT. otherwise, reduce the range
              for T in T_VECT:
                     if  (T==BR_L_REAC.index).any() != True:
                            T_VECT = np.delete(T_VECT,np.where(T_VECT==T))
                            warnings = warnings + ' warning: lumped reactivity at {TK} K not available. T skipped. \n'.format(TK=T)
       
       print(warnings)
       return BR_L_REAC,T_VECT


class PREPROCESSING:
       '''
       In this class, the matrix of rate constants is processed such that:
       a) the products are set as irreversible sinks
       b) the units of bimolecular reactions are modified if necessary
       '''

       def __init__(self,MATRIX_TP,SPECIES_SERIES,PRODS):
              self.MAT = MATRIX_TP                 #ARRAY
              self.SPECIES_SERIES = SPECIES_SERIES #SERIES
              self.PRODS = PRODS                   #LIST

       def SET_PRODS_SINKS(self,PRODSINKS):
              '''
              In this method, the reactivity of the products PRODS is set to 0 in the reactivity matrix
              This means that the columns kiP of the matrix are set to 0 
              Input: PRODSINKS number 0 or 1; if it is 1, it means the reactivity of the products must be set to 0
              '''
              MAT_PRODSINKS = self.MAT
              if PRODSINKS == 1:
                     # identify the indices of the products
                     i_PRODS = self.SPECIES_SERIES[self.PRODS].values
                     MAT_PRODSINKS[i_PRODS,:] = 0  # also self.mat is updated
              
              return MAT_PRODSINKS


       def MULTIPLY_BIMOL(self,SPECIES_BIMOL,UNITS_BIMOL):
              '''
              In this method, the rate constants of the kinetic matrix are multiplied by the avogadro number if the units are in cm3/molec/s
              input: isbimol: 0 if not, 1 if reaction bimolecular
                     UNITS_BIMOL: 'mol' or 'molec' for cm3/mol/s or cm3/molec/s, respectively
              '''
              MAT_BIMOL = self.MAT
              if UNITS_BIMOL == 'molec':
                     for SPECIES_i in self.SPECIES_SERIES.values:
                            if SPECIES_BIMOL[SPECIES_i] != '':
                                   # the species is bimolecular: multiply the corresponding rate constants                    
                                   MAT_BIMOL[SPECIES_i,:] *= 6.022E+23
                                   # convert to cm3/mol/s

              return MAT_BIMOL



class WRITE_MECH_CKI:
       '''
       In this class, the matrix of the rate constants at a selected T,P is used to write a mechanism in CKI format with fake thermodynamics
       This is required in order to be able to solve the system with opensmoke
       Input: MATRIX_TP_NORM: matrix of k0, either normalized or not (for generic mechanism)
              MATRIX_ALPHA = matrix of alpha in the expression k0*T^alpha*exp(-EA/R/T)
              MATRIX_EA = matrix of activation energies in cal/mol
              MATRIX_COMMENTS = matrix of strings with the comments for each reaction
              SPECIES_SERIES = series of species: names as indices, numbers as values
              i_REAC, i_PRODS = indices of reactant and product
       '''
       def __init__(self,MATRIX_TP_NORM,MATRIX_ALPHA,MATRIX_EA,MATRIX_COMMENTS,SPECIES_SERIES,i_REAC,i_PRODS,SPECIES_BIMOL):
              self.k0 = MATRIX_TP_NORM
              self.alpha = MATRIX_ALPHA
              self.EA = MATRIX_EA
              self.comments = MATRIX_COMMENTS
              self.SPECIES_SERIES = SPECIES_SERIES
              self.i_REAC = i_REAC  # array if lumped reactant
              self.i_PRODS = i_PRODS
              self.SPECIES_BIMOL = SPECIES_BIMOL


       def MAKE_CKI(self,PRODSINKS,ISOM_EQUIL):
              '''
              In this function:
              Organize the rate constants of self.MAT in chemkin format, at a defined temperature and pressure
              input: PRODSINKS => if it is 1, the products are irreversible sinks, therefore they are not considered in the reactivity
              if isom_equil=1 and i_reac is an array: exclude all the reactions not involving i_reac
              '''
              # organize the dataframe reactant => product k0 alpha EA !comments
              SPECIES_i = self.SPECIES_SERIES.values
              SPECIES_NAMES = np.copy(self.SPECIES_SERIES.index)
              #print(SPECIES_i,SPECIES_NAMES,self.SPECIES_BIMOL)
              # if the reaction is bimolecular: change the species name in S+S_ABU
              for S_i in SPECIES_i:
                     if self.SPECIES_BIMOL[S_i] != '':
                            SPECIES_NAMES[S_i] = SPECIES_NAMES[S_i] + '+' + self.SPECIES_BIMOL[S_i]

              if PRODSINKS == 1:
                     # EXCLUDE THE PRODUCTS FROM THE LOOP, BECAUSE THEY ARE NEVER REACTANTS
                     SPECIES_EFFECTIVE = np.delete(SPECIES_i,self.i_PRODS)
              elif PRODSINKS == 0:
                     SPECIES_EFFECTIVE = SPECIES_i

              if ISOM_EQUIL == 1 and isinstance(self.i_REAC,np.ndarray):
                     # EXCLUDE ALL THE NON-REACTIVE SPECIES FROM SPECIES_I, SPECIES_NAMES, AND SPECIES_EFFECTIVE
                     SPECIES_i = self.i_REAC
                     SPECIES_EFFECTIVE = SPECIES_i
                     SPECIES_NAMES = SPECIES_NAMES[self.i_REAC]

              CKI_lines = pd.DataFrame(index=list(np.arange(0,len(SPECIES_EFFECTIVE)*(len(SPECIES_i)-1))),columns=['reac_name','k0','alpha','EA','comments'])
              ii_row = 0       

              for S_i in SPECIES_EFFECTIVE:
                     # loop over the products with the exception of the  reactant selected
                     for Pr_i in np.delete(SPECIES_i,S_i): # S_i deleted excludes self reactions
                            CKI_lines['reac_name'][ii_row] = SPECIES_NAMES[S_i] + ' => ' + SPECIES_NAMES[Pr_i]
                            CKI_lines['k0'][ii_row] = '{:.2e}'.format(self.k0[S_i,Pr_i])
                            CKI_lines['alpha'][ii_row] = '{:.2f}'.format(self.alpha[S_i,Pr_i])
                            CKI_lines['EA'][ii_row] = '{:.2f}'.format(self.EA[S_i,Pr_i])
                            CKI_lines['comments'][ii_row] = '! ' + self.comments[S_i,Pr_i]
                            ii_row += 1

              return CKI_lines
              
              #np.savetxt(r'try.txt',CKI_lines.values,delimiter='\t',fmt='%s')

       def WRITE_CKI(self,cwd,CKI_lines_values):
              '''
              In this function:
              Write the rate constants in the file kin.txt
              '''
              self.path = cwd
              # check if the path exists
              if os.path.exists(self.path) == False:
                     raise RuntimeError('The folder ' + self.path + 'does not exist: impossible to compile the mech ')
              else:
                     kin = self.path + '/kin.txt'
                     # write arrays to concatenate 
                     elements = np.array([['ELEMENTS','','','',''],['C','H','O','N',''],['END','','','',''],['!','*','*','*','*'],['SPECIES','','','','']],dtype=str)
                     from_sp_to_reac = np.array([['END','','','',''],['!','*','*','*','*'],['REACTIONS','','','','']],dtype=str)
                     end = np.array(['END','','','',''],dtype=str)
                     end = end[np.newaxis,:]
                     AR_line = np.array(['AR','','','',''],dtype=str)
                     AR_line = AR_line[np.newaxis,:]
                     # write array with species: max 5 species per row
                     SPECIES_NAMES = np.copy(self.SPECIES_SERIES.index)
                     # for bimolecular reactants/products: append the name of the second reactant
                     SPECIES_i = self.SPECIES_SERIES.values
                     for S_i in SPECIES_i:
                            if self.SPECIES_BIMOL[S_i] != '' and (SPECIES_NAMES==self.SPECIES_BIMOL[S_i]).any() != True:
                                   # check also that the species is not already present in the list (possible for bimolecular reactions with the same second product)
                                   SPECIES_NAMES = np.append(SPECIES_NAMES,self.SPECIES_BIMOL[S_i])
                     N_rows_species = int(len(SPECIES_NAMES)/5)+1*(len(SPECIES_NAMES)/5>0)
                     empty_val = np.zeros(5-len(SPECIES_NAMES)%5,dtype=str)       # generates empty array of the selected size
                     SPECIES_TOWRITE = np.append(SPECIES_NAMES,empty_val)
                     SPECIES_TOWRITE = np.reshape(SPECIES_TOWRITE,(N_rows_species,5))           # reshape array appropriately
                     # concatenate arrays
                     data_final = np.concatenate((elements,AR_line,SPECIES_TOWRITE,from_sp_to_reac,CKI_lines_values,end))
                     # save mechanism
                     np.savetxt(kin,data_final,delimiter='\t',fmt='%s')

def WRITE_THERM(cwd,STOICH,SPECIES_SERIES,SPECIES_BIMOL_SERIES):
       
       '''
       In this function: write the fictitious thermodynamics in therm.txt
       input: STOICH provides the stoichiometry of the inlet species with carbon, hydrogen and oxygen
       THIS FUNCTION IS NOT IN THE CLASS WRITE_MECH_CKI because it does not require the rates
       '''
       # READ "STOICH" AND DERIVE THE STOICHIOMETRY TO SUBSTITUTE
       xxx = ' '*(3-len(STOICH[0][1:])) + STOICH[0][1:]
       yyy = ' '*(3-len(STOICH[1][1:])) + STOICH[1][1:]
       zzz = ' '*(3-len(STOICH[2][1:])) + STOICH[2][1:]

       # FIXED LINES TO WRITE
       header = 'THERMO \n'
       sub_header = '   300.000  1500.000  5000.000 \n'
       L1_part2 = '      C xxxH yyyO zzz     G    300.00   4000.00 1000.00      1\n'
       L2 = ' 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2\n'
       L3 = ' 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    3\n'
       L4 = ' 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00                   4\n'
       end = 'END'
       L1_part2_S_i = L1_part2.replace('xxx',xxx)
       L1_part2_S_i = L1_part2_S_i.replace('yyy',yyy)
       L1_part2_S_i = L1_part2_S_i.replace('zzz',zzz)
       # SPECIES NAMES
       SPECIES_NAMES = np.copy(SPECIES_SERIES.index)
       # WRITE FILE
       thermo = open(cwd + '/therm.txt','w')
       thermo.write(header)
       thermo.write(sub_header)
       # WRITE THERM FOR ARGON - INERT
       AR = 'AR'
       len_S_i = len(AR)   # characters in the species
       empty_spaces = 18-len_S_i # right number of spaces to reach L1_part2
       L1 = AR + ' '*empty_spaces + 'ATcT3EAr  1    0    0    0G    200.00   6000.00 1000.00      1\n'
       thermo.write(L1)
       thermo.write(L2)
       thermo.write(L3)
       thermo.write(L4)              
       # WRITE THERM FOR EVERY SPECIES
       for S_i in SPECIES_NAMES:
              len_S_i = len(S_i)   # characters in the species
              empty_spaces = 18-len_S_i # right number of spaces to reach L1_part2
              # check if the species is bimolecular reactant and remove 1 hydrogen
              # if the reaction is self-reaction: divide the stoichiometry by 2
              if SPECIES_BIMOL_SERIES[S_i] != '' and SPECIES_BIMOL_SERIES[S_i] != S_i :
                     stoich_reac = str(int(STOICH[1][1:])-1)
                     yyy = ' '*(3-len(stoich_reac)) + stoich_reac
                     L1_part2_reac = L1_part2.replace('xxx',xxx)
                     L1_part2_reac = L1_part2_reac.replace('yyy',yyy)
                     L1_part2_reac = L1_part2_reac.replace('zzz',zzz)
                     L1 = S_i + ' '*empty_spaces + L1_part2_reac
              # for self reactions: divide the stoichiometry by 2
              elif SPECIES_BIMOL_SERIES[S_i] != '' and SPECIES_BIMOL_SERIES[S_i] == S_i :
                     stoich_reacxB = str(int(int(STOICH[0][1:])/2))
                     stoich_reacyB = str(int(int(STOICH[1][1:])/2))
                     stoich_reaczB = str(int(int(STOICH[2][1:])/2))
                     xxxB = ' '*(3-len(stoich_reacxB)) + stoich_reacxB
                     yyyB = ' '*(3-len(stoich_reacyB)) + stoich_reacyB
                     zzzB = ' '*(3-len(stoich_reaczB)) + stoich_reaczB
                     L1_part2_reac = L1_part2.replace('xxx',xxxB)
                     L1_part2_reac = L1_part2_reac.replace('yyy',yyyB)
                     L1_part2_reac = L1_part2_reac.replace('zzz',zzzB)
                     L1 = S_i + ' '*empty_spaces + L1_part2_reac
              else:
                     L1 = S_i + ' '*empty_spaces + L1_part2_S_i
              thermo.write(L1)
              thermo.write(L2)
              thermo.write(L3)
              thermo.write(L4)

       # WRITE THE THERMODYNAMICS FOR THE SECOND REACTANT
       # NB this does not apply to self reactions where the second reactant is the same as the first one
       for S_i in SPECIES_NAMES:
              if SPECIES_BIMOL_SERIES[S_i] != ''  and SPECIES_BIMOL_SERIES[S_i] != S_i :
                     # write lines for the last species
                     S_i_last = SPECIES_BIMOL_SERIES[S_i]
                     len_S_i = len(S_i_last)   # characters in the species
                     empty_spaces = 18-len_S_i # right number of spaces to reach L1_part2
                     L1 = S_i_last + ' '*empty_spaces + '      C   0H   1O   0     G    300.00   4000.00 1000.00      1\n'
                     thermo.write(L1)
                     thermo.write(L2)
                     thermo.write(L3)
                     thermo.write(L4)

       thermo.write(end)
       thermo.close()


class WRITE_OS_INPUT:
       '''
       In this class:
       1. Extract the info you want to write in input_OS.dic
       2. Copy input_OS.dic and write
       '''
       def __init__(self,cwd,T,P,SPECIES,REACLUMPED,N_INIT,SPECIES_BIMOL_SERIES,BR_L_REAC):
              '''
              Check if the template exist,
              Generate its copy to edit,
              Generate the variables to overwrite
              '''
              self.cwd = cwd
              self.path = cwd + '/input_OS_template.dic'
              if os.path.isfile(self.path) == False:
                     raise RuntimeError('The file input_OS_template.dic does not exist: impossible to generate OS input ')
              else:
                     with open(self.path,mode='r') as OS_templ:
                            OS_file = OS_templ.readlines()
                     self.OS_file = OS_file
                     self.T = T
                     self.P = P
                     # at the end of the species output: also put the bimolecular reactant if present and if not self reaction
                     # distinguish between lumped and single reactant
                     if isinstance(REACLUMPED[0],np.ndarray):
                            # extract the mole fractions at the corresponding temperature
                            REAC_BRANCHING =  pd.Series(N_INIT*BR_L_REAC.loc[T,:].values,index=BR_L_REAC.columns)
                            # check that the sum of the BFs is <1, otherwise re-normalize it (else: AR fraction might be negative)
                            if np.sum(REAC_BRANCHING.values) > 1:
                                   REAC_BRANCHING[BR_L_REAC.columns] = REAC_BRANCHING[BR_L_REAC.columns]/np.sum(REAC_BRANCHING.values) 
                            # LUMPED REACTANTS
                            # take the first value and check if a second species should be added
                            if SPECIES_BIMOL_SERIES[REACLUMPED[0][0]] != '' and SPECIES_BIMOL_SERIES[REACLUMPED[0][0]] != REACLUMPED[0][0] :
                                   REAC_BRANCHING =REAC_BRANCHING*0.00001
                                   # check if all names of bimol fragments are identical; if so, assign N_REAC_ABU; otherwise, divide it among the fragments
                                   BIMOL_FRAGMENTS = SPECIES_BIMOL_SERIES[REACLUMPED[0]]
                                   N_AR = 1-np.sum(REAC_BRANCHING.values)-N_INIT*0.99998
                                   if (BIMOL_FRAGMENTS==BIMOL_FRAGMENTS[0]).all():
                                          SPECIES = np.append(SPECIES,BIMOL_FRAGMENTS[0])
                                          # you have just 1 "abundant" fragment (generally happens with CKI inputs)
                                          REAC_ABU = BIMOL_FRAGMENTS[0]
                                          N_REAC_ABU = N_INIT*0.99998
                                          REAC_ABU_WRITE = REAC_ABU + ' ' + str(N_REAC_ABU)  + ' '
                                   else:
                                          raise ValueError('Bimolecular II fragment is different among isomers: check lumping ')

                            else:
                                   REAC_ABU_WRITE = ''
                                   N_REAC_ABU = ''
                                   N_AR = 1-np.sum(REAC_BRANCHING.values)
                            
                            REAC_BRANCHING_WRITE = ''
                            for RR in BR_L_REAC.columns:
                                   REAC_BRANCHING_WRITE = REAC_BRANCHING_WRITE + RR + ' ' + str(REAC_BRANCHING[RR]) + ' '   

                            self.MOLE_FRACTIONS = REAC_BRANCHING_WRITE + REAC_ABU_WRITE + ' AR {:.5f}'.format(N_AR) 
                     else:
                            REAC = REACLUMPED[0]
                            if SPECIES_BIMOL_SERIES[REAC] != '' and SPECIES_BIMOL_SERIES[REAC] != REAC:
                                   SPECIES = np.append(SPECIES,SPECIES_BIMOL_SERIES[REAC])
                                   REAC_ABU = SPECIES_BIMOL_SERIES[REAC]
                                   N_REAC = N_INIT*0.00001
                                   N_REAC_ABU = N_INIT*0.99999
                                   self.MOLE_FRACTIONS = REAC + ' ' + str(N_REAC) + ' ' + REAC_ABU + ' ' + str(N_REAC_ABU)  + ' AR ' + str(1-N_INIT)
                            # differentiate the input for bimolecular reactions with a second "diffrent" reactant
                            else:
                                   self.MOLE_FRACTIONS = REAC + ' ' + str(N_INIT) + ' AR ' + str(1-N_INIT)

                     self.SPECIES = ' '.join(SPECIES)

       def w_input_OS(self):
              '''
              1. IDENTIFY KEYWORDS TO SUBSTITUTE
              2. OVERWRITE
              '''
              newfile = copy.deepcopy(self.OS_file)
              keywords = ['TEMP','PR','MF','SPECIES']
              values = [str(self.T),str(self.P),self.MOLE_FRACTIONS,self.SPECIES]
              # associate indices to the keywords:
              key_val = pd.Series(values,index=keywords)

              for kk in keywords:
                     newval = key_val[kk]
                     for idx,row in enumerate(newfile):
                            newfile[idx] = re.sub(kk,newval,row)

              if os.path.isfile(self.cwd + '/input_OS.dic'):
                     os.remove(self.cwd + '/input_OS.dic')

              with open(self.cwd + '/input_OS.dic',mode='x') as inp:
                     inp.writelines(newfile)




       



              


              

                                          
