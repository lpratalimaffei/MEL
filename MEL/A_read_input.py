import numpy as np
import pandas as pd
import re
import os
import sys

def read_pseudospecies(pseudospecies_file):
    '''
    Read pseudospecies file and generate series and array of pseudospecies
    '''
    pseudospecies_array = np.genfromtxt(
        pseudospecies_file, delimiter='', dtype=str)
    pseudospecies_array = np.array(pseudospecies_array, dtype='<U160')
    # set an array of all pseudospecies
    stable_species = []
    # generate pseudospecies series
    try:
        pseudospecies_series = pd.Series(
            pseudospecies_array[:, 1], index=pseudospecies_array[:, 0])
        pseudospecies_series.index.astype('<U160')
    except IndexError:
        print('error: only 1 pseudospecies provided. Reactivity will be meaningless')
        sys.exit()

    for SP in pseudospecies_series.index:
        if pseudospecies_series.loc[SP].split('+') != [pseudospecies_series.loc[SP]]:
            for i in pseudospecies_series.loc[SP].split('+'):
                stable_species.append(i)
            pseudospecies_series.loc[SP] = np.array(
                pseudospecies_series.loc[SP].split('+'), dtype='<U160')
        else:
            stable_species.append(pseudospecies_series.loc[SP])

    stable_species = np.array(stable_species, dtype='<U160')
    return pseudospecies_series, stable_species


class READ_INPUT:

    def __init__(self, cwd, inputfile):
        self.filename = os.path.join(cwd, inputfile)
        self.cwd = cwd
        # the main variable is the name of the file

    def read_file_lines(self):

        errorlist = ''
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
                                self.OS_folder = line.split()[2]
                                # if env variable: convert it to absolute
                                if self.OS_folder[0] == '%':
                                    self.OS_folder = os.environ[self.OS_folder[1:-1]]

                            # select the type of input
                            if line.find('mech_type') != -1:
                                self.inp_type = line.split()[2]

                            if line.find('P_vect') != -1:
                                # select things in between square brackets
                                line_Pvect = re.split('\[|\]', line)[1]
                                # get the list of pressures and turn it into an array
                                list_Pvect = line_Pvect.split()
                                self.Pvect = np.array(
                                    list_Pvect, dtype=np.float32)

                            if line.find('T_range') != -1:
                                line_Tvect = re.split('\[|\]', line)[1]
                                list_Tvect = line_Tvect.split()
                                T_range = np.array(list_Tvect, dtype=np.int16)
                                self.Tvect = np.arange(
                                    T_range[0], T_range[1]+100, 100)

                            if line.find('T_skip') != -1:
                                line_Tvect_skip = re.split('\[|\]', line)[1]
                                list_Tvect_skip = line_Tvect_skip.split()
                                Tvect_skip = np.array(
                                    list_Tvect_skip, dtype=np.int16)

                            if line.find('units_bimol') != -1:
                                self.units_bimol = line.split()[2]

                            if line.find('cutoff') != -1:
                                line_cutoff = line.split()[2]
                                lower_cutoff = re.split('-', line_cutoff)[0]
                                upper_cutoff = re.split('-', line_cutoff)[1]
                                self.cutoff = np.array(
                                    [lower_cutoff, upper_cutoff], dtype=np.float32)

                            if line.find('verboseoutput') != -1:
                                verbose = line.split()[2]
                                if verbose == 'True':
                                    self.verbose = True
                                else:
                                    self.verbose = None

                            if line.find('end') != -1:
                                readinput = 0

                        if readjobs == 1:
                            # read the list of jobs
                            job_ids = ['preproc_irreversible', 'single_simulation', 'prescreening_equilibrium', 'prescreening_allreactive',
                                       'composition_selection', 'lumping', 'validation']

                            for job in job_ids:
                                if line.find(job) != -1:
                                    jobid.append(job)

                            if line.find('end') != -1:
                                readjobs = 0
                                readsubdict = 1
                                # generate the job_list dictionary
                                job_list = dict.fromkeys(jobid)
                                # empty subdictionaries for lumping and validation - if they are present in the list
                                job_ids_empty = ['lumping', 'validation', 'preproc_irreversible']
                                for job_empty in job_ids_empty:
                                    if job_empty in job_list:
                                        job_list[job_empty] = dict.fromkeys(
                                        [''], '')

                        if readsubdict == 1:
                            # read the subdictionaries of the tasks you perform

                            ############### subdictionaries for the prescreening ################
                            if line.find('prescreening_equilibrium') != -1 and any('prescreening_equilibrium' == x for x in jobid) and line.find('simul_type') == -1:
                                read_prescreen_equil = 1

                            if line.find('prescreening_allreactive') != -1 and any('prescreening_allreactive' == x for x in jobid) and line.find('simul_type') == -1:
                                read_prescreen_allreactive = 1

                            if line.find('pseudospecies') != -1 and (read_prescreen_equil == 1 or read_prescreen_allreactive == 1) and read_single_simul == 0:
                                line_pseudospecies = re.split('\[|\]', line)[1]
                                pseudospecies = line_pseudospecies.split()
                                # generate the list of psudospecies
                                # names of the lumped products, including the single products
                                pseudospecies_names = []
                                # list of arrays with the set of products contained in each species of Prods_Lumped
                                pseudospecies_components = []

                                for PS in pseudospecies:
                                    if len(PS.split('+')) == 1 and read_prescreen_equil == 1 and PS != 'all':
                                        # error: you cannot have single pseudospecies for isomer equilibrium simulations
                                        errorlist = errorlist + \
                                            '\nsubdictionary prescreening_equilibrium: cannot have single isomers'

                                    elif len(PS.split('+')) == 1 and read_prescreen_equil == 1 and PS == 'all':
                                        # you will consider the full isomer pool
                                        pseudospecies_names.append(PS)
                                        pseudospecies_components.append(PS)

                                    elif len(PS.split('+')) == 1 and read_prescreen_allreactive == 1:
                                        # append it also in the lumped products
                                        pseudospecies_names.append(PS)
                                        pseudospecies_components.append(PS)

                                    elif len(PS.split('+')) > 1:
                                        # append the first name, modified, to the list of pseudospecies
                                        pseudospecies_names.append(
                                            PS.split('+')[0] + '_L')
                                        # append the array of products to Prods_Lumped_array:
                                        pseudospecies_components.append(
                                            np.array(PS.split('+'), dtype=str))

                                # dataframe with pseudospecies
                                pseudospecies_df = pd.Series(
                                    pseudospecies_components, index=pseudospecies_names)

                            if line.find('end') != -1 and read_prescreen_equil == 1:
                                # store values in dictionary
                                job_list['prescreening_equilibrium'] = dict.fromkeys(
                                    ['pseudospecies'], pseudospecies_df)
                                # exit from reading this portion
                                read_prescreen_equil = 0

                            if line.find('end') != -1 and read_prescreen_allreactive == 1:
                                # store values in dictionary
                                job_list['prescreening_allreactive'] = dict.fromkeys(
                                    ['pseudospecies'], pseudospecies_df)
                                # exit from reading this portion
                                read_prescreen_allreactive = 0

                            ############# subdictionaries for composition selection ################
                            if line.find('composition_selection') != -1 and any('composition_selection' == x for x in jobid) and line.find('simul_type') == -1:
                                read_comp_sel = 1

                            if read_comp_sel == 1:

                                if line.find('maxiter') != -1:
                                    try:
                                        maxiter = int(line.split()[2])
                                    except ValueError as e:
                                        errorlist = errorlist + \
                                            '\ninvalid input in composition_selection subdictionary: number of iterations should be an integer number'

                                if line.find('BF_tolerance') != -1:
                                    try:
                                        BF_tol = float(line.split()[2])
                                        if BF_tol >= 1 or BF_tol <= 0:
                                            errorlist = errorlist + \
                                                '\nerror in composition_selection subdictionary: BF tolerance should be betweeen 0 and 1'
                                    except ValueError as e:
                                        errorlist = errorlist + \
                                            '\ninvalid input in composition_selection subdictionary: BF tolerance should be a number between 0 and 1'

                                if line.find('end') != -1:
                                    # save variabile in series
                                    try:
                                        subdict_comp_sel = dict(
                                            zip(['BF_tol', 'maxiter'], [BF_tol, maxiter]))
                                        job_list['composition_selection'] = subdict_comp_sel
                                    except NameError as e:
                                        errorlist = errorlist + \
                                            '\nmissing subdictionary specification in composition_selection: \n {}'.format(
                                                e)
                                    # exit from reading this portion
                                    read_comp_sel = 0

                            ########## subdictionaries for single simulation ####################

                            if line.find('single_simulation') != -1 and any('single_simulation' == x for x in jobid):
                                read_single_simul = 1
                                simul_type = ''

                            if read_single_simul == 1:

                                if line.find('simul_type') != -1:
                                    simul_type = line.split()[2]

                                if line.find('pseudospecies') != -1 and (simul_type == 'prescreening_allreactive' or simul_type == 'prescreening_equilibrium' or simul_type == 'composition_selection'):
                                    line_pseudospecies = re.split(
                                        '\[|\]', line)[1]
                                    pseudospecies = line_pseudospecies.split()

                                    # generate the list of psudospecies
                                    # names of the lumped products, including the single products
                                    pseudospecies_names = []
                                    # list of arrays with the set of products contained in each species of Prods_Lumped
                                    pseudospecies_components = []

                                    for PS in pseudospecies:
                                        if len(PS.split('+')) == 1 and (simul_type == 'prescreening_equilibrium' or simul_type == 'composition_selection') and PS != 'all':
                                            # error: you cannot have single pseudospecies for isomer equilibrium simulations
                                            errorlist = errorlist + \
                                                '\nsubdictionary prescreening_equilibrium/composition_selection: cannot have single isomers'

                                        elif len(PS.split('+')) == 1 and (simul_type == 'prescreening_equilibrium' or simul_type == 'composition_selection') and PS == 'all':
                                            pseudospecies_names.append(PS)
                                            pseudospecies_components.append(PS)

                                        elif len(PS.split('+')) == 1 and simul_type == 'prescreening_allreactive':
                                            # append it also in the lumped products
                                            pseudospecies_names.append(PS)
                                            pseudospecies_components.append(PS)

                                        elif len(PS.split('+')) > 1:
                                            # append the first name, modified, to the list of pseudospecies
                                            pseudospecies_names.append(
                                                PS.split('+')[0] + '_L')
                                            # append the array of products to Prods_Lumped_array:
                                            pseudospecies_components.append(
                                                np.array(PS.split('+'), dtype=str))

                                    # dataframe with pseudospecies
                                    pseudospecies_df = pd.Series(
                                        pseudospecies_components, index=pseudospecies_names)

                                if line.find('maxiter') != -1 and simul_type == 'composition_selection':
                                    try:
                                        maxiter = int(line.split()[2])
                                    except ValueError as e:
                                        errorlist = errorlist + \
                                            '\ninvalid input in single_simulation subdictionary: number of iterations should be an integer number'

                                if line.find('BF_tolerance') != -1 and simul_type == 'composition_selection':
                                    try:
                                        BF_tol = float(line.split()[2])
                                        if BF_tol >= 1 or BF_tol <= 0:
                                            errorlist = errorlist + \
                                                '\nerror in single_simulation subdictionary: BF tolerance should be betweeen 0 and 1'
                                    except ValueError as e:
                                        errorlist = errorlist + \
                                            '\ninvalid input in single_simulation subdictionary: BF tolerance should be a number between 0 and 1'

                                # read reactants and products
                                if line.find('Reac ') != -1:
                                    Reac = re.split('\[|\]', line)[1]
                                    # if the reactant is lumped: generate an array
                                    if len(Reac.split('+')) > 1:
                                        self.Reac = np.array(
                                            Reac.split('+'), dtype=str)
                                        reaclumped = self.Reac[0] + '_L'
                                    else:
                                        # the reactant is just 1 so you don't need to do anything
                                        self.Reac = Reac
                                        reaclumped = Reac

                                if line.find('Prod ') != -1:
                                    line_Prod = re.split('\[|\]', line)[1]
                                    Prod = line_Prod.split()
                                    # generate the list of products: if you have lumped products, allocate them too
                                    self.Prod = []              # list of single species
                                    Prods_Lumped = []           # names of the lumped products, including the single products
                                    # list of arrays with the set of products contained in each species of Prods_Lumped
                                    Prods_Lumped_array = []
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
                                            Prods_Lumped.append(
                                                Pr.split('+')[0] + '_L')
                                            # append the array of products to Prods_Lumped_array:
                                            Prods_Lumped_array.append(
                                                np.array(Pr.split('+'), dtype=str))
                                            for Pr_L in Pr.split('+'):
                                                # append the products to the array of products
                                                self.Prod.append(Pr_L)

                                if line.find('end') != -1:
                                    # subdictionaries

                                    job_list['single_simulation'] = dict.fromkeys(
                                        ['simul_type'], simul_type)

                                    if simul_type == 'composition_selection':
                                        # build subdictionary if needed
                                        try:
                                            job_list['single_simulation']['BF_tol'] = BF_tol
                                            job_list['single_simulation']['maxiter'] = maxiter
                                            job_list['single_simulation']['pseudospecies'] = pseudospecies_df
                                        except NameError as e:
                                            errorlist = errorlist + \
                                                '\nmissing subdictionary specification in composition_selection: \n {}'.format(
                                                    e)

                                    elif simul_type == 'prescreening_equilibrium' or simul_type == 'prescreening_allreactive':
                                        # check for pseudospecies keyword
                                        try:
                                            job_list['single_simulation']['pseudospecies'] = pseudospecies_df
                                        except NameError as e:
                                            errorlist = errorlist + \
                                                '\nmissing subdictionary specification in single_simulation (prescreening): \n {}'.format(
                                                    e)

                                    elif simul_type == 'lumping' or simul_type == 'validation':
                                        # subdictionary storing reactant and products
                                        # series of reactant/products lumped
                                        try:
                                            prodslumped = pd.Series(
                                                Prods_Lumped_array, index=Prods_Lumped)
                                            reaclumped = pd.Series(
                                                [self.Reac], index=[reaclumped])
                                            job_list['single_simulation'].update(dict(zip(['REAC', 'REACLUMPED', 'PRODS', 'PRODSLUMPED'], [
                                                                                 self.Reac, reaclumped, self.Prod, prodslumped])))
                                        except NameError as e:
                                            errorlist = errorlist + \
                                                '\nmissing reactant/product specifications in single_simulation subdictionary: \n {}'.format(
                                                    e)

                                    else:
                                        # wrong simul type
                                        errorlist = errorlist+'\nNameError: simulation type unavailable '

                                    # exit from reading this portion
                                    read_single_simul = 0
        myfile.close()

        if not errorlist:
            # np.savetxt('pseudospecies_prova.txt',pseudospecies_df,delimiter='\t',fmt='%s')
            # units bimol: if the input is MESS, set molec by default
            if self.inp_type == 'MESS' and (self.units_bimol != 'molec' or self.units_bimol != 'mol'):
                self.units_bimol = 'molec'

            # if verboseoutput is missing: set to None
            if not hasattr(self, 'verbose'):
                self.verbose = None
            # input_parameters dictionary
            keys_inputpar = ['opensmoke_folder', 'mech_type', 'P_vect',
                             'T_vect', 'T_skip', 'units_bimol', 'cutoff', 'verbose']
            values_inputpar = [self.OS_folder, self.inp_type, self.Pvect,
                               self.Tvect, Tvect_skip, self.units_bimol, self.cutoff, self.verbose]
            input_parameters = dict(zip(keys_inputpar, values_inputpar))

            return input_parameters, job_list

        else:
            raise RuntimeError(
                'Errors detected when reading the input: ' + str(errorlist))

    def CHECK_INPUT(self):
        """
        Now check that the variables saved are present when necessary, and that reactants and products are different species
        """
        error_list = ''
        # check that the OS folder exists

        if os.path.exists(self.OS_folder) == False:
            error_list = error_list + '\nThe opensmoke bin folder indicated does not exist'

        # check that the input_OS_template exists
        OS_inpfile = os.path.join(self.cwd, 'inp', 'input_OS_template.dic')
        if os.path.isfile(OS_inpfile) == False:
            error_list = error_list + \
                '\nMissing OS template input file: {} required '.format(
                    OS_inpfile)

        # check that the input_preproc exists
        preproc_inpfile = os.path.join(self.cwd, 'inp', 'input_preproc.dic')
        if os.path.isfile(preproc_inpfile) == False:
            error_list = error_list + \
                '\nMissing OS preprocessor input file: {} required '.format(
                    preproc_inpfile)

        # check type of input
        if self.inp_type != 'MESS' and self.inp_type != 'CKI':
            error_list = error_list + '\nWrong keyword for input: set MESS or CKI '

        elif self.inp_type == 'MESS':
            # check that mech file and output file exist in the input folder
            mechinp = os.path.join(self.cwd, 'inp', 'me_ktp.inp')
            mechout = os.path.join(self.cwd, 'inp', 'rate.out')

            if os.path.isfile(mechinp) == False:
                error_list = error_list + \
                    '\nMissing MESS input file: {} required '.format(mechinp)
            if os.path.isfile(mechout) == False:
                error_list = error_list + \
                    '\nMissing MESS output file: {} required '.format(mechout)

        # check that MESS input parameters are defined appropriately
        elif self.inp_type == 'CKI':

            # chech that mech files exist in the input folder
            mechfile = os.path.join(self.cwd, 'inp', 'kin.CKI')
            if os.path.isfile(mechfile) == False:
                error_list = error_list + \
                    '\nMissing CKI input file: {} required '.format(mechfile)

            # extra checks
            if self.Tvect.size == 0:
                error_list = error_list + '\nTemperature range not defined '

            if self.units_bimol != 'molec' and self.units_bimol != 'mol':
                error_list = error_list + '\nUnits for bimolecular reactions wrong or not defined '

        # check that the cutoff is between 0 and 1 and that the lower cutoff is higher
        if self.cutoff[0] > self.cutoff[1]:
            error_list = error_list + '\nLower cutoff limit must be smaller than upper cutoff limit'
        if (self.cutoff[0] < 0 or self.cutoff[0] > 1) or (self.cutoff[1] < 0 or self.cutoff[1] > 1):
            error_list = error_list + '\nCutoff limits must be between 0 and 1'

        if not error_list:
            return None
        else:
            raise RuntimeError(
                'Errors detected when reading the input: ' + str(error_list))

        # no need to check P_vect now, because if it is empty you can read it elsewhere or assume everything is high P limit

    def COMPARE_INPUT_PT(self, input_type, mech_dict):
        """
        CHECK that the input provided externally is compatible with the one given by the user
        It can be used also for CKI inputs, however the P_VECT_MESS and T_VECT_MESS will not be provided
        """
        # useful to do with pandas if you want to restructure, and then you use &
        error_list = ''

        if input_type == 'MESS':
            P_VECT_MESS = mech_dict['P_VECT_MESS']
            T_VECT_MESS = mech_dict['T_VECT_MESS']
            if np.array([[p == P_VECT_MESS for p in self.Pvect][ii].any() for ii in range(0, len(self.Pvect))]).all() != True:
                # YOU ARE CHECKING IF, INDEPENDENT OF THE ORDER, ALL THE ELEMENTS OF self.Pvect are contained in P_VECT_MESS
                error_list = error_list + \
                    '\nSome of the selected pressures are not available in MESS pressure list \n \t please select among {PRESSURES}'.format(
                        PRESSURES=P_VECT_MESS)

            if (T_VECT_MESS > self.Tvect[0]).all() or (T_VECT_MESS < self.Tvect[-1]).all():
                # CHECK IF THE RANGE OF TEMPERATURE SELECTED IS INCLUDED IN THAT AVAILABLE IN THE MESS OUTPUT
                error_list = error_list + '\nThe temperature range in the input is outside of the mess range available: \n \t please select between {minT} and {maxT} K'.format(
                    minT=min(T_VECT_MESS), maxT=max(T_VECT_MESS))

        if not error_list:
            return None
        else:
            raise RuntimeError(
                'Errors detected when comparing the input with NESS mechanism file: \n' + str(error_list))

    def CHECK_INPUT_JOBTYPE(self, job_list, mech_dict):
        '''
        Perform checks for each type of job
        a. input values provided by the user
        b. compatibility with the mechanism input parameters
        '''
        # extract variables from mech_dict
        SPECIES = mech_dict['SPECIES']
        SPECIES_BIMOL = mech_dict['SPECIES_BIMOL']
        error_list = ''

        # loop over the dictionaries
        for key, value in job_list.items():

            # for single simulation: re-set job type
            if key == 'single_simulation':
                jobtype = value['simul_type']
            elif key in ['prescreening_equilibrium', 'prescreening_allreactive', 'composition_selection',
                    'lumping', 'validation', 'preproc_irreversible']:
                jobtype = key
            else:
                print('\nValue Error: dictionary type not recognized')
            # call subdictionaries
            subdict = value

            if key != 'single_simulation' and (jobtype == 'composition_selection' or jobtype == 'lumping' or jobtype == 'validation'):
                # check that file pseudospecies.txt exists
                pseudospecies_file = os.path.join(
                    self.cwd, 'inp', 'pseudospecies.txt')
                if os.path.isfile(pseudospecies_file) == False:
                    error_list = error_list + \
                        '\n composition_selection/lumping/validation require file {}, not found!'.format(
                            pseudospecies_file)

            # for preprocessing: throw error if mech is not chemkin
            if jobtype == 'preproc_irreversible' and self.inp_type != 'CKI':
                error_list = error_list + \
                    '\n preprocessing to irreversible mech is only available for chemkin inputs'

            if jobtype == 'composition_selection':
                # check if values are numbers
                try:
                    # MAXITERATIONS CHECK
                    if subdict['maxiter'] > 1000:
                        error_list = error_list + \
                            '\nIterations for species composition convergence set above 1000: unreasonable '
                    # CHECK ON BRANCHING FRACTION TOLERANCE
                    if subdict['BF_tol'] <= 0 or subdict['BF_tol'] > 1:
                        error_list = error_list + \
                            '\nBranching fraction tolerance for composition convergence is out of valid 0<BF<1 range '

                except ValueError as e:
                    error_list = error_list + \
                        '\nInconsistent values for composition_selection subdictionary: {e}'.format(
                            e)

            # VALIDATION: DOES THE LUMPED MECH EXIST?
            if jobtype == 'validation':
                lumpedmech_kinfile = os.path.join(
                    self.cwd, 'lumpedmech', 'kin.txt')
                lumpedmech_thermfile = os.path.join(
                    self.cwd, 'lumpedmech', 'therm.txt')
                if not 'lumping' in job_list.keys():
                    if os.path.isfile(lumpedmech_kinfile) == False:
                        error_list = error_list + \
                            '\n validation requires lumped mech {}, not found!'.format(
                                lumpedmech_kinfile)

                    if os.path.isfile(lumpedmech_thermfile) == False:
                        error_list = error_list + \
                            '\n validation requires lumped thermo {}, not found!'.format(
                                lumpedmech_thermfile)

            # PSEUDOSPECIES: CHECK THAT THEY ARE NOT REPEATED IN EACH GROUP
            if jobtype == 'prescreening_equilibrium' or jobtype == 'prescreening_allreactive':

                # check if you have more than one species/pseudospecies group
                if subdict['pseudospecies'].index.size > 1:
                    for SP in subdict['pseudospecies'].index:
                        for SP_other in subdict['pseudospecies'].index:
                            if SP_other != SP:
                                if isinstance(subdict['pseudospecies'].loc[SP], np.ndarray) and isinstance(subdict['pseudospecies'].loc[SP_other], np.ndarray):
                                    for species in subdict['pseudospecies'].loc[SP]:
                                        if np.array([species == sp_other for sp_other in subdict['pseudospecies'].loc[SP_other]]).any():
                                            error_list = error_list + '\nMatching species in different pseudospecies group '

                                elif isinstance(subdict['pseudospecies'].loc[SP], str) and isinstance(subdict['pseudospecies'].loc[SP_other], np.ndarray):
                                    species = subdict['pseudospecies'].loc[SP]
                                    if np.array([species == sp_other for sp_other in subdict['pseudospecies'].loc[SP_other]]).any():
                                        error_list = error_list + '\nMatching species in different pseudospecies group '

                # check that the pseudospecies indicated match one of the species in the list
                for SP in subdict['pseudospecies'].index:
                    # CHECK THAT THE PSEUDOSPECIES ARE PRESENT IN THE LIST
                    if isinstance(subdict['pseudospecies'].loc[SP], np.ndarray):
                        if np.array([[sp == SPECIES for sp in subdict['pseudospecies'].loc[SP]][ii].any() for ii in range(0, len(subdict['pseudospecies'].loc[SP]))]).all() != True:
                            error_list = error_list + \
                                '\nNot all species indicated in the set of pseudospecies match the species list. \n \t please select among {SPECIES} '.format(
                                    SPECIES=str(SPECIES))
                        # WITHIN EACH GROUP: CHECK THAT PSEUDOSPECIES BELONG TO THE SAME TYPE (I.E. HAVE SAME BIMOL FRAGMENT)
                        i_sp_group = np.where(
                            [sp_group == SPECIES for sp_group in subdict['pseudospecies'].loc[SP]])[1]
                        if np.array([i_sp_bim == SPECIES_BIMOL[i_sp_group] for i_sp_bim in SPECIES_BIMOL[i_sp_group]]).all == False:
                            error_list = error_list + \
                                '\nPseudospecies do not belong to the same type: mix of unimol/bimol or bimol with different fragments '

                    elif isinstance(subdict['pseudospecies'].loc[SP], str) and subdict['pseudospecies'].loc[SP] != 'all':
                        species = subdict['pseudospecies'].loc[SP]
                        if np.array([species == SP_list for SP_list in SPECIES]).any() != True:
                            error_list = error_list + \
                                '\nNot all species indicated in the set of pseudospecies match the species list. \n \t please select among {SPECIES} '.format(
                                    SPECIES=str(SPECIES))

            # SINGLE_SIMULATIONS
            # check that you don't have products or reactants matching
            if key == 'single_simulation':
                if jobtype == 'lumping' or jobtype == 'validation':
                    if not list(self.Reac) or not self.Prod:
                        # if not defined, the output of the readline above will be #
                        error_list = error_list + '\nNo reactants or products defined '

                    if np.array([[Pr == SPECIES for Pr in self.Prod][ii].any() for ii in range(0, len(self.Prod))]).all() != True:
                        # CHECK IF THE PRODUCTS SELECTED ARE AVAILABLE IN THE SPECIES LIST OF MESS
                        error_list = error_list + \
                            '\nThe Products selected are not all available in the list of species! \n \t please select among {SPECIES} '.format(
                                SPECIES=str(SPECIES))

                    if isinstance(self.Reac, np.ndarray):

                        if np.array([[Rr == SPECIES for Rr in self.Reac][ii].any() for ii in range(0, len(self.Reac))]).all() != True:
                            # CHECK IF THE PRODUCTS SELECTED ARE AVAILABLE IN THE SPECIES LIST OF MESS
                            error_list = error_list + \
                                '\nThe Reactants selected are not all available in the list of species! \n \t please select among {SPECIES} '.format(
                                    SPECIES=str(SPECIES))

                        for Rr in self.Reac:
                            # check that reactants and products do not match
                            if np.array([Rr == Pr for Pr in self.Prod]).any():
                                error_list = error_list + \
                                    '\nThe reactant matches one of the product species. Inconsistent! '

                    elif isinstance(self.Reac, str):
                        if np.array([self.Reac == Pr for Pr in self.Prod]).any():
                            error_list = error_list + \
                                '\nThe reactant matches one of the product species. Inconsistent! '
                        if np.array([self.Reac == SP for SP in SPECIES]).any() != True:
                            # CHECK IF THE REACTANT SELECTED IS AVAILABLE IN THE RANGE OF SPECIES
                            error_list = error_list + \
                                '\nThe Reactant selected is not available in the list of species! \n \t please select among {SPECIES} '.format(
                                    SPECIES=str(SPECIES))

                elif jobtype == 'prescreening_equilibrium' or jobtype == 'prescreening_allreactive' or jobtype == 'composition_selection':
                    if subdict['pseudospecies'].index.size > 1:
                        error_list = error_list + \
                            '\nMore than 1 set of pseudospecies: incompatible with single simulations'
                    elif subdict['pseudospecies'].index.size == 1:
                        # CHECK THAT THE PSEUDOSPECIES ARE PRESENT IN THE LIST
                        if isinstance(subdict['pseudospecies'].iloc[0], np.ndarray):
                            if np.array([[sp == SPECIES for sp in subdict['pseudospecies'].iloc[0]][ii].any() for ii in range(0, len(subdict['pseudospecies'].iloc[0]))]).all() != True:
                                error_list = error_list + \
                                    '\nNot all species indicated in the set of pseudospecies match the species list. \n \t please select among {SPECIES} '.format(
                                        SPECIES=str(SPECIES))
                            # WITHIN EACH GROUP: CHECK THAT PSEUDOSPECIES BELONG TO THE SAME TYPE (I.E. HAVE SAME BIMOL FRAGMENT)
                            i_sp_group = np.where(
                                [sp_group == SPECIES for sp_group in subdict['pseudospecies'].iloc[0]])[1]
                            if np.array([i_sp_bim == SPECIES_BIMOL[i_sp_group] for i_sp_bim in SPECIES_BIMOL[i_sp_group]]).all == False:
                                error_list = error_list + \
                                    '\nPseudospecies do not belong to the same type: mix of unimol/bimol or bimol with different fragments '

                        elif isinstance(subdict['pseudospecies'].iloc[0], str) and subdict['pseudospecies'].iloc[0] != 'all':
                            species = subdict['pseudospecies'].loc[SP]
                            if np.array([species == SP_list for SP_list in SPECIES]).any() != True:
                                error_list = error_list + \
                                    '\nNot all species indicated in the set of pseudospecies match the species list. \n \t please select among {SPECIES} '.format(
                                        SPECIES=str(SPECIES))

            # if you find inp/pseudospecies.txt: check that all species indicated:
            # are present in the species list
            # have consistent names (1st species + '_L' for lumped groups)

        if not error_list:
            return None
        else:
            raise RuntimeError(
                'Errors detected when checking input_lumping.txt: \n' + str(error_list))

    def set_inputparam_job(self, jobtype):
        """
        Based on job type: 
        set the rest of input parameters:
        - plot comparison : only in case of lumping
        - prodsinks (in case of lumping)
        - isom_equil (in case of prescreening_equil)
        """

        keys_inputpar = ['plot_cmp', 'Prods_sinks', 'isom_equil']

        if jobtype == 'lumping' or jobtype == 'validation':
            plot_cmp = 'YES'
        else:
            plot_cmp = 'NO'

        if jobtype == 'lumping':
            Prods_sinks = 1
        else:
            Prods_sinks = 0

        if jobtype == 'prescreening_equilibrium':
            isom_equil = 1
        else:
            isom_equil = 0

        # set the dictionary
        input_par_jobtype = dict(
            zip(keys_inputpar, [plot_cmp, Prods_sinks, isom_equil]))

        return input_par_jobtype
