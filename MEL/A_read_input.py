
import numpy as np
import pandas as pd
import re
import os

class READ_INPUT:

    def __init__(self,cwd):
        self.filename = cwd + '/input_lumping.txt'
        self.cwd = cwd
        # the main variable is the name of the file


    def read_file_lines(self):
        """
        read the file line by line and sort it into the variables you need
        """
        self.plot_cmp = 'NO' # default
        self.opt_mech = ''   # default
        
        with open(self.filename) as myfile:
            for line in myfile:
                # find OS folder
                if line.find('opensmoke_folder') != -1:
                    self.OS_folder = '"'+line.split()[2] + '"'

                # select the type of input
                if line.find('input_type') != -1:
                    self.inp_type = line.split()[2]

                if line.find('N_init') != -1:
                    try:
                        self.moles_init = float(line.split()[2])
                        # check if the mole fraction 
                        if self.moles_init > 1:
                            print('Initial mole fraction of the reactant cannot exceed 1')
                            return None
                    except:
                        print('Initial number of moles is not a number. Please put a number instead')
                        return None

                if line.find('P_vect') != -1:
                    # select things in between square brackets
                    line_Pvect = re.split('\[|\]',line)[1]
                    #get the list of pressures and turn it into an array
                    list_Pvect = line_Pvect.split()
                    self.Pvect = np.array(list_Pvect,dtype=np.float32)

                if line.find('T_range') != -1:
                    line_Tvect = re.split('\[|\]',line)[1]
                    list_Tvect = line_Tvect.split()
                    T_range = np.array(list_Tvect,dtype=np.int16)
                    self.Tvect = np.arange(T_range[0],T_range[1]+100,100)

                if line.find('T_skip') != -1:
                    line_Tvect_skip = re.split('\[|\]',line)[1]
                    list_Tvect_skip = line_Tvect_skip.split()
                    self.Tvect_skip = np.array(list_Tvect_skip,dtype=np.int16)

                if line.find('Prods_sinks') != -1:
                    try:
                        self.Prods_sinks = int(line.split()[2])
                        if self.Prods_sinks != 1 and self.Prods_sinks != 0:
                            print('The Prods_sinks option must be 0 or 1')
                            return None
                    except:
                        print('The Prods_sinks option must be 0 or 1')
                        return None

                if line.find('Bimol_activate') != -1:
                    try:
                        self.bimol_activate = int(line.split()[2])
                        if self.bimol_activate != 1 and self.bimol_activate != 0:
                            print('The Bimol_activate option must be 0 or 1')
                            return None
                    except:
                        print('The Bimol_activate option must be 0 or 1')
                        return None

                if line.find('Stoichiometry ') != -1:
                    line_stoich = re.split('\[|\]',line)[1]
                    self.stoich = line_stoich.split() 

                if line.find('Reac ') != -1:
                    Reac = re.split('\[|\]',line)[1]
                    # if the reactant is lumped: generate an array
                    if len(Reac.split('+')) > 1: 
                        self.Reac = np.array(Reac.split('+'),dtype=str)
                        reaclumped = self.Reac[0] + '_L'
                    else:
                        # the reactant is just 1 so you don't need to do anything
                        self.Reac = Reac
                        reaclumped = Reac

                if line.find('Prod ') != -1:
                    line_Prod = re.split('\[|\]',line)[1]
                    Prod = line_Prod.split()
                    # generate the list of products: if you have lumped products, allocate them too
                    self.Prod = []              # list of single species
                    Prods_Lumped = []           # names of the lumped products, including the single products
                    Prods_Lumped_array = []     # list of arrays with the set of products contained in each species of Prods_Lumped
                    # if you have lumped products: generate a series and allocate the products to a new list
                    for Pr in Prod:
                        if len(Pr.split('+')) == 1:
                            # single product: append it directly to the list
                            self.Prod.append(Pr)
                            # append it also in the lumped products
                            Prods_Lumped.append(Pr)
                            Prods_Lumped_array.append(Pr)

                        elif len(Pr.split('+')) > 1:
                            # append the first name, modified, to the list of lumped species
                            Prods_Lumped.append(Pr.split('+')[0] + '_L')
                            # append the array of products to Prods_Lumped_array:
                            Prods_Lumped_array.append(np.array(Pr.split('+'),dtype=str))
                            for Pr_L in Pr.split('+'):
                                # append the products to the array of products
                                self.Prod.append(Pr_L)

                if line.find('Isomer_equilibrium') != -1:
                    self.isom_equil = int(line.split()[2])                             

                if line.find('units_bimol') != -1:
                    self.units_bimol = line.split()[2]

                # read plotting options
                if line.find('plot_compare') != -1:
                    self.plot_cmp = line.split()[2]

                if line.find('opt_mech') != -1:
                    self.opt_mech = line.split()[2]

        myfile.close()
        # optimized mech: substitute with '' if no path was found 
        if self.opt_mech == '#':
            self.opt_mech = ''
        # units bimol: if the input is MESS, set molec by default
        if self.inp_type=='MESS' and (self.units_bimol != 'molec' or self.units_bimol != 'mol'):
            self.units_bimol = 'molec'

        # series of prodslumped
        prodslumped = pd.Series(Prods_Lumped_array,index=Prods_Lumped)
        reaclumped = pd.Series([self.Reac],index=[reaclumped])
        print(self.OS_folder)
        return self.OS_folder,self.inp_type,self.moles_init,self.Pvect,self.Tvect,self.Tvect_skip,self.stoich,self.Reac,reaclumped,self.Prod,prodslumped,self.isom_equil,self.units_bimol,self.bimol_activate,self.Prods_sinks,self.plot_cmp,self.opt_mech

    
    def CHECK_INPUT(self):
        """
        Now check that the variables saved are present when necessary, and that reactants and products are different species
        """
        error_list = ''
        # check if present - always necessary
        # type of input

        if self.inp_type != 'MESS' and self.inp_type != 'CKI':
            error_list = error_list + ' Wrong keyword for input: set MESS or CKI \n'

        if not list(self.Reac) or not self.Prod:
            # if not defined, the output of the readline above will be #
            error_list = error_list + ' No reactants or products defined \n'

        if isinstance(self.Reac,np.ndarray):
            for Rr in self.Reac:
                if np.array([Rr == Pr for Pr in self.Prod]).any():
                    error_list = error_list + ' The reactant matches one of the product species. Inconsistent! \n'

        elif isinstance(self.Reac,str):
            if np.array([self.Reac == Pr for Pr in self.Prod]).any():
                error_list = error_list + ' The reactant matches one of the product species. Inconsistent! \n'
            if self.isom_equil == 1:
                error_list = error_list + ' The reactant is not a set of isomers, but Isomer_equilibrium keyword is activated. Inconsistent! \n'

        if self.inp_type != 'MESS':
            # extra checks
            if self.Tvect.size == 0:
                error_list = error_list + ' Temperature range not defined \n'
            
            if self.units_bimol != 'molec' and self.units_bimol != 'mol':
                error_list = error_list + ' Units for bimolecular reactions wrong or not defined \n'

        # check that the stoichiometry is written correctly
        if self.stoich[0][0] != 'C' or self.stoich[1][0] != 'H' or self.stoich[2][0] != 'O':
            error_list = error_list + ' The order of the stoichiometry is incorrect: please define C,H,O \n'
        try:
            int(self.stoich[0][1:]) + int(self.stoich[1][1:]) + int(self.stoich[2][1:])
        except:
            error_list = error_list + ' The stoichiometry coefficients cannot be converted to integers \n'

        # check that isomer_equilibrium and prod_sinks are compatible
        if self.isom_equil ==1 and self.Prods_sinks == 1:
            error_list = error_list + ' The products are set as infinite sinks, but Isomer_equilibium keyword is active: incompatible \n'
        if self.isom_equil ==1 and self.bimol_activate == 0:
            error_list = error_list + ' Isom_equil requires Bimol_activate = 1 otherwise bimolecular reactions are not recognized \n'
        # check if the path to optimized mech. exists
        if os.path.exists(self.cwd + '/' + self.opt_mech) == False and self.opt_mech != '':
            error_list = error_list + ' The path to the optimized mechanism does not exist \n'

        if self.plot_cmp != 'YES' and self.plot_cmp != 'NO':
            error_list = error_list + ' Error in variable of plot comparison: must be "YES" or "NO" \n'

        if self.plot_cmp == 'YES' and self.Prods_sinks == 0:
            error_list = error_list + ' Prods_sinks == 0 is incompatible with plotting the lumped mech. set it to 1 or set the plot comparison to "NO" \n'

        if not error_list:
            return None
        else:
            raise RuntimeError('Errors detected when reading the input: \n' + str(error_list))


        # no need to check P_vect now, because if it is empty you can read it elsewhere or assume everything is high P limit
        
    def COMPARE_INPUT_PVECT(self,input_type,SPECIES,P_VECT_MESS,T_VECT_MESS):
        """
        CHECK that the input provided externally is compatible with the one given by the user
        It can be used also for CKI inputs, however the P_VECT_MESS and T_VECT_MESS will not be provided
        """
        # useful to do with pandas if you want to restructure, and then you use &
        error_list = ''
        # for MESS input types: check the pressure and temperature vectors
        if input_type == 'MESS':
            if np.array([[p == P_VECT_MESS for p in self.Pvect][ii].any() for ii in range(0,len(self.Pvect))]).all() != True:
                # YOU ARE CHECKING IF, INDEPENDENT OF THE ORDER, ALL THE ELEMENTS OF self.Pvect are contained in P_VECT_MESS
                error_list = error_list + ' Some of the selected pressures are not available in MESS pressure list \n \t please select among {PRESSURES} \n'.format(PRESSURES=P_VECT_MESS)

            if (T_VECT_MESS > self.Tvect[0]).all() or (T_VECT_MESS < self.Tvect[-1]).all():
                # CHECK IF THE RANGE OF TEMPERATURE SELECTED IS INCLUDED IN THAT AVAILABLE IN THE MESS OUTPUT
                error_list = error_list + ' The temperature range in the input is outside of the mess range available: \n \t please select between {minT} and {maxT} K \n'.format(minT=min(T_VECT_MESS),maxT=max(T_VECT_MESS))

        # the other things will always be checked
        # check reactant: the check will be different if the reactant is lumped (array of reactants)
        if isinstance(self.Reac,np.ndarray):
            if np.array([[Rr == SPECIES for Rr in self.Reac][ii].any() for ii in range(0,len(self.Reac))]).all() != True:
                # CHECK IF THE PRODUCTS SELECTED ARE AVAILABLE IN THE SPECIES LIST OF MESS
                error_list = error_list + ' The Reactants selected are not all available in the list of species! \n \t please select among {SPECIES} \n'.format(SPECIES = str(SPECIES))                
        elif isinstance(self.Reac,str):
            if np.array([self.Reac == SP for SP in SPECIES]).any() != True:
                # CHECK IF THE REACTANT SELECTED IS AVAILABLE IN THE RANGE OF SPECIES
                error_list = error_list + ' The Reactant selected is not available in the list of species! \n \t please select among {SPECIES} \n'.format(SPECIES = str(SPECIES))

        if np.array([[Pr == SPECIES for Pr in self.Prod][ii].any() for ii in range(0,len(self.Prod))]).all() != True:
            # CHECK IF THE PRODUCTS SELECTED ARE AVAILABLE IN THE SPECIES LIST OF MESS
            error_list = error_list + ' The Products selected are not all available in the list of species! \n \t please select among {SPECIES} \n'.format(SPECIES = str(SPECIES))


        if not error_list:
            return None
        else:
            raise RuntimeError('Errors detected when comparing the input with MESS file: \n' + str(error_list))



        
        
            