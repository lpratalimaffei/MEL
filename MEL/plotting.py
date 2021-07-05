import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class PLOTTING:
       '''
       This function produces subplots of the profiles of the products at a given pressure
       '''
       def __init__(self,cwd):
              '''
              Checks if the output folder "PLOTS" exists and if not it creates it.
              '''
              self.path = cwd
              self.path_plots = os.path.join(cwd,'PLOTS')
              if os.path.exists(self.path_plots) == False:
                     os.mkdir(self.path_plots)


       def plot_data(self,profiles_Pi,REAC,PRODS,SPECIES,P):
              '''
              Extract the data common to each kind of plot
              To be called at every new pressure
              '''
              # extract the T_VECT from the dictionary keys of profiles_Pi
              T_LIST = list(profiles_Pi['detailed'].keys())
              self.T_VECT = np.array(T_LIST,dtype=int)
              # Define the number of rows and columns
              if int(len(self.T_VECT)**0.5) == (len(self.T_VECT)**0.5):
                     self.n_cols = int(len(self.T_VECT)**0.5)
                     self.n_rows = self.n_cols
              else:
                     self.n_cols = int(len(self.T_VECT)**0.5)+1
                     self.n_rows = int(len(self.T_VECT)**0.5)+int((len(self.T_VECT)/self.n_cols)>int(len(self.T_VECT)**0.5))
              self.profiles = profiles_Pi #dictionaries

              self.REAC = REAC
              self.SPECIES = SPECIES
              self.PRODS = PRODS
              self.P = P

              # linestyles for plotting
              self.lstyles = ['-','--',':']
              palette = plt.rcParams['axes.prop_cycle'].by_key()['color']
              if len(palette) < len(SPECIES):
                     ntot = round(len(SPECIES)/len(palette)) +1
                     for i in np.arange(0,ntot):
                            palette = palette + palette
              self.palette_series = pd.Series(palette[0:len(SPECIES)],index=SPECIES)
              self.palette_series[REAC] = 'black'

              # fontsizes
              if len(self.T_VECT) > 9:
                     self.fsize_vect = [6,7,8]
              else:
                     self.fsize_vect = [8,9,10]
              # legend label size
              if len(self.SPECIES) <= 10 and len(self.T_VECT) <= 9:
                     self.legsize = 6
              elif len(self.SPECIES) <= 10 and len(self.T_VECT) > 9:
                     self.legsize = 5
              elif len(self.SPECIES) > 10 and len(self.T_VECT) <= 9:
                     self.legsize = 4
              elif len(self.SPECIES) > 10 and len(self.T_VECT) > 9:
                     self.legsize = 3


       def plot_exp_profiles(self,N_INIT):
              '''
              Produces profiles of the form Wi VS time
              '''
              # generate the figure
              fig,axes = plt.subplots(self.n_rows,self.n_cols)
              # set axis
              row_plot = 0
              col_plot = 0
              Ti = 0

              for T in self.T_VECT:
                     Ti += 1
                     # loop on the dictionaries in self.profiles
                     kk = 0
                     # set vect for maxt
                     maxt = []
                     for key in list(self.profiles.keys()):
                            profiles_i = self.profiles[key]
                            # extract the matrix with the profiles
                            tW_DF = profiles_i[T]
                            t = tW_DF['t'].values
                            maxt.append(t[-1])
                            t = t[:,np.newaxis]
                            # plot the species only if present in the list
                            for Si in self.SPECIES:
                                   if np.where(tW_DF.columns==Si)[0].size > 0:
                                          axes[row_plot,col_plot].plot(t,tW_DF[Si],color=self.palette_series[Si],ls=self.lstyles[kk])      
                            kk += 1

                     # labels and formatting
                     axes[self.n_rows-1,col_plot].set_xlabel(r'$t$ [s]',fontsize=self.fsize_vect[1])
                     # if you don't have any plot below in the column: put the label
                     #miss_plots = (self.n_rows*self.n_cols)-len(self.T_VECT)
                     #miss_col = np.arange(self.n_cols-miss_plots,self.n_cols)
                     if T == self.T_VECT[-1] and col_plot != self.n_cols-1:
                            for miss_col in np.arange(col_plot+1,self.n_cols):
                                   axes[self.n_rows-2,miss_col].set_xlabel(r'$t$ [s]',fontsize=self.fsize_vect[1])

                     axes[row_plot,0].set_ylabel(r'$X_{i}$ [-]',fontsize=self.fsize_vect[1])
                     axes[row_plot,col_plot].set_title((str(T) + ' K'),fontsize=self.fsize_vect[2],fontweight='bold')
                     axes[row_plot,col_plot].set_xlim([0,min(maxt)])
                     axes[row_plot,col_plot].set_ylim([0,N_INIT])
                     axes[row_plot,col_plot].axes.ticklabel_format(axis='both', style='sci', scilimits=(0,0),useLocale=True,useMathText=True)
                     axes[row_plot,col_plot].yaxis.offsetText.set_fontsize(self.fsize_vect[0])
                     axes[row_plot,col_plot].xaxis.set_tick_params(labelsize=self.fsize_vect[0])
                     axes[row_plot,col_plot].xaxis.offsetText.set_fontsize(self.fsize_vect[0])
                     axes[row_plot,col_plot].yaxis.set_tick_params(labelsize=self.fsize_vect[0])
                     
                     
                     # update of the axis
                     row_plot = row_plot + int((Ti%self.n_cols)==0) 
                     col_plot = col_plot + 1 - self.n_cols*int((Ti%self.n_cols)==0) 

              if row_plot == self.n_rows-1:
                     for cols in np.arange(col_plot,self.n_cols):
                            fig.delaxes(axes[-1,cols])        
              axes[0,0].legend(list(self.SPECIES),loc=4,fontsize=self.legsize)
              # set tight layout
              fig.tight_layout(pad=0.1,w_pad=0.1)
              # save the figure
              fig_path = (self.path_plots + '/Prof_EXP_' + self.REAC + '_' + str(self.P) + '_atm.png')
              if os.path.isfile(fig_path):
                     os.remove(fig_path)
              fig.savefig(fig_path,dpi=200)
              
              plt.close()


       def plot_data_reac(self,profiles_Pi,REAC,P,N_INIT_REAC):
              '''
              Plot the profiles of the reactant composition
              '''

              # generate the figure
              fig,axes = plt.subplots(self.n_rows,self.n_cols)
              # set axis
              row_plot = 0
              col_plot = 0
              Ti = 0

              for T in self.T_VECT:
                     Ti += 1
                     # extract the profiles
                     tW_DF = profiles_Pi[T]
                     t = tW_DF['t'].values
                     t = t[:,np.newaxis]
                     # redefine tW_DF based on the sum of the reactants mole fraction
                     tW_DF_REAC_tot = np.sum(tW_DF[REAC].values,axis=1)
                     tW_DF_REAC_tot = tW_DF_REAC_tot[:,np.newaxis]
                     tW_DF[REAC] = tW_DF[REAC]/tW_DF_REAC_tot
                     # plot the species
                     axes[row_plot,col_plot].plot(t,tW_DF[REAC])      
                     # labels and formatting
                     axes[self.n_rows-1,col_plot].set_xlabel(r'$t$ [s]',fontsize=self.fsize_vect[1])
                     axes[row_plot,0].set_ylabel(r'$X_{i}$ [-]',fontsize=self.fsize_vect[1])
                     axes[row_plot,col_plot].set_title((str(T) + ' K'),fontsize=self.fsize_vect[2],fontweight='bold')
                     axes[row_plot,col_plot].set_ylim([0,np.max(np.max(tW_DF[REAC]))])
                     axes[row_plot,col_plot].axes.ticklabel_format(axis='both', style='sci', scilimits=(0,0),useLocale=True,useMathText=True)
                     axes[row_plot,col_plot].yaxis.offsetText.set_fontsize(self.fsize_vect[0])
                     axes[row_plot,col_plot].xaxis.set_tick_params(labelsize=self.fsize_vect[0])
                     axes[row_plot,col_plot].xaxis.offsetText.set_fontsize(self.fsize_vect[0])
                     axes[row_plot,col_plot].yaxis.set_tick_params(labelsize=self.fsize_vect[0])
                        
                     # update of the axis
                     row_plot = row_plot + int((Ti%self.n_cols)==0) 
                     col_plot = col_plot + 1 - self.n_cols*int((Ti%self.n_cols)==0) 

              if row_plot == self.n_rows-1:
                     for cols in np.arange(col_plot,self.n_cols):
                            fig.delaxes(axes[-1,cols])        
              axes[0,0].legend(list(REAC),loc=4,fontsize=self.legsize)
              # set tight layout
              fig.tight_layout(pad=0.1,w_pad=0.1)
              # save the figure
              fig_path = (os.path.join(self.path_plots, 'Prof_EXP_' + self.REAC + 'composition_' + str(self.P) + '_atm.png'))
              if os.path.isfile(fig_path):
                     os.remove(fig_path)
              fig.savefig(fig_path,dpi=200)
              
              plt.close()