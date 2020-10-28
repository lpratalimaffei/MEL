
import numpy as np
import pandas as pd
import re
import os

class READ_INPUT:

    def __init__(self,cwd,inputfile):
        self.filename = os.path.join(cwd,inputfile)
        self.cwd = cwd
        # the main variable is the name of the file

    def read_file_lines(self):


        readinput = 0
        readjobs = 0
        readsubdict = 0
        jobid = []  # empty list
        read_prescreen_equil = 0
        read_prescreen_allreactive = 0
        read_comp_sel = 0
        read_single_simul = 0

        with open(self.filename) as myfile:
            for line in myfile:
                # skip empty lines
                if len(line.split()) > 0:
                    # if it starts with "#": do not read
                    if line.split()[0][0] != '#':

                        if line.find('input') != -1:
                            readinput = 1
                        if line.find('jobs') != -1:
                            readjobs = 1

                        # read the input section
                        if readinput == 1:
                            # find OS folder
                            if line.find('opensmoke_folder') != -1:
                                OS_folder = '"'+line.split()[2] + '"'

                            # select the type of input
                            if line.find('mech_type') != -1:
                                self.inp_type = line.split()[2]

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
                                Tvect_skip = np.array(list_Tvect_skip,dtype=np.int16)

                            if line.find('units_bimol') != -1:
                                self.units_bimol = line.split()[2]

                            if line.find('Stoichiometry ') != -1:
                                line_stoich = re.split('\[|\]',line)[1]
                                self.stoich = line_stoich.split() 

                            if line.find('cutoff') != -1:
                                line_cutoff = line.split()[2]
                                lower_cutoff = re.split('-',line_cutoff)[0] 
                                upper_cutoff = re.split('-',line_cutoff)[1] 

                            if line.find('end') != -1:
                                readinput = 0

                        if readjobs == 1:
                            # read the list of jobs
                            if line.find('single_simulation') != -1:
                                jobid.append('single_simulation')

                            if line.find('prescreening_equilibrium') != -1:
                                jobid.append('prescreening_equilibrium')

                            if line.find('prescreening_allreactive') != -1:
                                jobid.append('prescreening_allreactive')

                            if line.find('composition_selection') != -1:
                                jobid.append('composition_selection')

                            if line.find('lumping') != -1:
                                jobid.append('lumping')

                            if line.find('validation') != -1:
                                jobid.append('validation')

                            if line.find('end') != -1:
                                readjobs = 0
                                readsubdict = 1
                                # generate the job_list dictionary
                                job_list = dict.fromkeys(jobid)

                        if readsubdict == 1:
                            # read the subdictionaries of the tasks you perform

                            ############### subdictionaries for the prescreening ################
                            if line.find('prescreening_equilibrium') != -1 and any('prescreening_equilibrium'==x for x in jobid) :
                                read_prescreen_equil = 1

                            if line.find('prescreening_allreactive') != -1 and any('prescreening_allreactive'==x for x in jobid) :
                                read_prescreen_allreactive = 1

                            if line.find('pseudospecies') != -1 and (read_prescreen_equil == 1 or read_prescreen_allreactive == 1):
                                print(line)
                                print(re.split('\[|\]',line))
                                line_pseudospecies = re.split('\[|\]',line)[1]
                                pseudospecies = line_pseudospecies.split()
                                # generate the list of psudospecies
                                pseudospecies_names = []           # names of the lumped products, including the single products
                                pseudospecies_components = []     # list of arrays with the set of products contained in each species of Prods_Lumped
                                
                                for PS in pseudospecies:
                                    if len(PS.split('+')) == 1 and read_prescreen_equil == 1:
                                        # error: you cannot have single pseudospecies for isomer equilibrium simulations
                                        print(ValueError('subdictionary prescreening_equilibrium: cannot have single isomers'))
                                    elif len(PS.split('+')) == 1 and read_prescreen_allreactive == 1:
                                        # append it also in the lumped products
                                        pseudospecies_names.append(PS)
                                        pseudospecies_components.append(PS)

                                    elif len(PS.split('+')) > 1:
                                        # append the first name, modified, to the list of pseudospecies
                                        pseudospecies_names.append(PS.split('+')[0] + '_L')
                                        # append the array of products to Prods_Lumped_array:
                                        pseudospecies_components.append(np.array(PS.split('+'),dtype=str))

                                # dataframe with pseudospecies
                                pseudospecies_df = pd.Series(pseudospecies_components,index=pseudospecies_names)

                            if line.find('end') != -1 and read_prescreen_equil == 1:
                                # store values in dictionary
                                job_list['prescreening_equilibrium'] = pseudospecies_df
                                # exit from reading this portion
                                read_prescreen_equil = 0

                            if line.find('end') != -1 and read_prescreen_allreactive == 1:
                                # store values in dictionary
                                job_list['prescreening_allreactive'] = pseudospecies_df
                                # exit from reading this portion
                                read_prescreen_allreactive = 0

                            ############# subdictionaries for composition selection ################
                            if line.find('composition_selection') != -1 and any('composition_selection'==x for x in jobid) :
                                read_comp_sel = 1

                            if read_comp_sel == 1:
                                if line.find('BF_tolerance') != -1 :
                                    BF_tol = float(line.split()[2])
                                    
                                    if BF_tol >= 1:
                                        print('error in composition_selection subdictionary: BF tolerance should be below 1')
                                        exit()

                                if line.find('maxiter') != -1 :
                                    maxiter = line.split()[2]

                                if line.find('end') != -1:
                                    # save variabile in series
                                    try:
                                        subdict_comp_sel = dict(zip(['BF_tol','maxiter'],[BF_tol,maxiter]))
                                        job_list['composition_selection'] = subdict_comp_sel
                                    except NameError as e:
                                        print('\nmissing subdictionary specification in composition_selection: \n {}'.format(e))
                                        exit()
                                    # exit from reading this portion
                                    read_comp_sel = 0

                            if line.find('single_simulation') != -1 and any('single_simulation'==x for x in jobid) :
                                read_single_simul = 1

                            if read_single_simul == 1:

                                simul_type= ''
                                if line.find('simul_type') != -1 :
                                    simul_type = line.split()[2]

                                if line.find('maxiter') != -1 and simul_type == 'composition_selection':
                                    maxiter = line.split()[2]

                                if line.find('BF_tolerance') != -1 and simul_type == 'composition_selection':
                                    BF_tol = line.split()[2]
                                    if BF_tol >= 1:
                                        print('error in composition_selection subdictionary: BF tolerance should be below 1')
                                        exit()

                                # read reactants and products
                                if line.find('Reac ') != -1:
                                    Reac = re.split('\[|\]',line)[1]
                                    # if the reactant is lumped: generate an array
                                    if len(Reac.split('+')) > 1: 
                                        self.Reac = np.array(Reac.split('+'),dtype=str)
                                        reaclumped = Reac[0] + '_L'
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

                                if line.find('end') != -1:
                                    # subdictionaries
                                    job_list['single_simulation'] = dict.fromkeys(['main_spec','additional_spec'])
                                    if simul_type == 'composition_selection':
                                        # build subdictionary if needed
                                        try:
                                            subdict_comp_sel = dict(zip(['BF_tol','maxiter'],[BF_tol,maxiter]))
                                            job_list['single_simulation']['additional_spec'] = subdict_comp_sel
                                        except NameError as e:
                                            print('\nmissing subdictionary specification in composition_selection: \n {}'.format(e))
                                            exit()

                                    # subdictionary storing reactant and products
                                    # series of reactant/products lumped
                                    try:
                                        prodslumped = pd.Series(Prods_Lumped_array,index=Prods_Lumped)
                                        reaclumped = pd.Series([self.Reac],index=[reaclumped])
                                        job_list['single_simulation']['main_spec'] = dict(zip(['REAC','REACLUMPED','PRODS','PRODSLUMPED'],[self.Reac,reaclumped,self.Prod,prodslumped]))
                                    except NameError as e:
                                        print('\nmissing reactant/product specifications in single_simulation subdictionary: \n {}'.format(e))
                                        exit()

                                    # exit from reading this portion
                                    read_single_simul = 0
        myfile.close()

        # units bimol: if the input is MESS, set molec by default
        if self.inp_type=='MESS' and (self.units_bimol != 'molec' or self.units_bimol != 'mol'):
            self.units_bimol = 'molec'

        # input_parameters dictionary
        keys_inputpar = ['opensmoke_folder','mech_type','P_vect','T_vect','T_skip','units_bimol','stoich','cutoff'] 
        values_inputpar = [OS_folder,self.inp_type,self.Pvect,self.Tvect,Tvect_skip,self.units_bimol,self.stoich]
        input_parameters = dict(zip(keys_inputpar,values_inputpar))


        return input_parameters,job_list

    def read_file_lines_old(self):
        """
        read the file line by line and sort it into the variables you need
        """
        self.plot_cmp = 'NO' # default
        self.opt_mech = ''   # default
        
        with open(self.filename) as myfile:
            for line in myfile:
                # find OS folder

                if line.find('Prods_sinks') != -1:
                    try:
                        self.Prods_sinks = int(line.split()[2])
                        if self.Prods_sinks != 1 and self.Prods_sinks != 0:
                            print('The Prods_sinks option must be 0 or 1')
                            return None
                    except:
                        print('The Prods_sinks option must be 0 or 1')
                        return None


                if line.find('Isomer_equilibrium') != -1:
                    self.isom_equil = int(line.split()[2])                             

                # read plotting options
                if line.find('plot_compare') != -1:
                    self.plot_cmp = line.split()[2]

                if line.find('opt_mech') != -1:
                    self.opt_mech = line.split()[2]

        myfile.close()
        # optimized mech: substitute with '' if no path was found 
        if self.opt_mech == '#':
            self.opt_mech = ''



        print(self.OS_folder)
        return self.OS_folder,self.inp_type,self.Pvect,self.Tvect,self.Tvect_skip,self.stoich,self.Reac,reaclumped,self.Prod,prodslumped,self.isom_equil,self.units_bimol,self.Prods_sinks,self.plot_cmp,self.opt_mech

    
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


# functions to create the input of each type of simulation
        
        
            