#!/usr/bin/env pythonw
import tkinter as Tk
from tkinter import ttk
import matplotlib
import numpy as np
import numpy.ma as ma
import new_cmaps
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as PathEffects
from scipy.special import kn # Modified Bessel function
from scipy.stats import linregress
from scipy.integrate import simps, cumtrapz
class SpectralPanel:
    # A dictionary of all of the parameters for this plot with the default parameters

    plot_param_dict = {'spectral_type': 0, #0 dn/dp, 1 = dn/dE
                       'twoD': 0, # Plot is not 2D
                       'show_ions': 1, # Include ion spectrum?
                       'show_electrons': 1, # Include electron Spectrum?
                       'rest_frame': False, # Show the data in the rest frame?
                       'set_ylim': True,
                       'set_xlim': True,
                       'x_min': 0.05,
                       'x_max': 200,
                       'spatial_x': False, # Important info for the main program
                       'spatial_y': False,
                       'show_legend': True, # Include legend?
                       'normalize_spectra': True, # To normalze the spectrum so it's integral = 1
                       'y_min': -6, # in logspace
                       'y_max': 0, # in logspace
                       'MeasureEpsP': False, #
                       'IonLeft': -10000.0,
                       'PowerLawIonMax': 10.0,
                       'GammaIonInjection': 1.0,
                       'DoPowerLawFitIon': False,
                       'DoPowerLawFitElectron': False,
                       'DelGami': 0.06,
                       'DelGame': 0.03,
                       'SetTi': False,
                       'PowerLawIonMin': 1.0,
                       'SetTe': False,
                       'GammaElectronInjection': 30.0,
                       'PowerLawElectronMin': 1.0,
                       'ElectronRight': 0.0,
                       'PrtlIntegrationRelative': True,
                       'PowerLawElectronMax': 10.0,
                       'ElectronLeft': -10000.0,
                       'IonRight': 0.0,
                       'MeasureEpsE': False,
                       'eNormalizer' : 0.,
                       'iNormalizer' : 0.,
                       'BoostedIons' : False, #Show a energy boosted ion plot to compare to electrons
                       'T_legend_loc':  'N/A', # location of temperature legend.
                       'PL_legend_loc': 'N/A'} # lcation of power law legend

    # We need the types of all the parameters for the config file
    BoolList = ['twoD', 'set_ylim', 'set_xlim',
                'spatial_x', 'spatial_y', 'normalize_spectra',
                'show_ions', 'show_electrons', 'rest_frame','show_legend',
                'PrtlIntegrationRelative', 'BoostedIons'
                'SetTe', 'SetTi','MeasureEpsP', 'MeasureEpsE',
                'DoPowerLawFitElectron', 'DoPowerLawFitIon']
    IntList = ['spectral_type']
    FloatList = ['y_min', 'y_max', 'x_min', 'x_max',
                 'DelGame', 'DelGami', 'GammaIonInjection', 'GammaElectronInjection',
                 'PowerLawElectronMin', 'PowerLawElectronMax',
                 'PowerLawIonMin', 'PowerLawIonMax', 'eNormalizer', 'iNormalizer',
                 'ElectronLeft', 'ElectronRight', 'IonLeft', 'IonRight',]
    StrList = ['T_legend_loc', 'PL_legend_loc']

    def __init__(self, parent, figwrapper):
        self.settings_window = None
        self.FigWrap = figwrapper
        self.parent = parent
        self.ChartTypes = self.FigWrap.PlotTypeDict.keys()
        self.chartType = self.FigWrap.chartType
        self.figure = self.FigWrap.figure
        self.PowerlawPworked = False
        self.PowerlawEworked = False
        if self.GetPlotParam('spectral_type') == 1:
            self.SetPlotParam('x_min', 0.0005, update_plot = False)
            self.SetPlotParam('x_max', 100, update_plot = False)

        # When the plot is initialized it will look to see if there are other spectral plots.
        self.spect_num = -1 #it will count itself so we have to correct for an off by one error.
        for i in range(self.parent.MainParamDict['NumOfRows']):
            for j in range(self.parent.MainParamDict['NumOfCols']):
                if self.parent.SubPlotList[i][j].chartType == 'SpectraPlot':
                    self.spect_num+=1

        self.dashes = self.parent.dashes_options[min(self.spect_num,len(self.parent.dashes_options)-1)]
        self.RegionFlag = True # Set to true if the region lines need to be written to a phase plot
        # The variables that store the eps_p & eps_e values
        self.eps_pVar = Tk.StringVar()
        self.eps_pVar.set('N/A')

        self.eps_eVar = Tk.StringVar()
        self.eps_eVar.set('N/A')

    def ChangePlotType(self, str_arg):
        self.FigWrap.ChangeGraph(str_arg)

    def set_plot_keys(self):
        '''A helper function that will insure that each hdf5 file will only be
        opened once per time step'''
        # First make sure that omega_plasma & istep is loaded so we can fix the
        # integrate over the specified regions

        self.arrs_needed = ['c_omp', 'istep', 'gamma', 'xsl', 'mi', 'me', 'gamma0']
        if self.GetPlotParam('rest_frame'):
            # Set the loading of the rest frame spectra
            self.arrs_needed.append('specerest')
            self.arrs_needed.append('specprest')
        else:
            # Load the normal spectra
            self.arrs_needed.append('spece')
            self.arrs_needed.append('specp')


        return self.arrs_needed

    def LoadData(self):
        self.c_omp = self.FigWrap.LoadKey('c_omp')[0]
        self.istep = self.FigWrap.LoadKey('istep')[0]
#        self.gamma0 = self.FigWrap.LoadKey('gamma0')[0]
        self.xsl = self.FigWrap.LoadKey('xsl')/self.c_omp
        # Load gamma-1 of the spectra
        self.gamma = self.FigWrap.LoadKey('gamma')
        self.massRatio = self.FigWrap.LoadKey('mi')[0]/self.FigWrap.LoadKey('me')[0]
        self.keyname = 'spectra_data'
        if self.GetPlotParam('normalize_spectra'):
            self.keyname +='_normalized_'

        if self.GetPlotParam('rest_frame'):
            self.keyname +='in_rest_frame'


        # Select the x-range from which to take the spectra
        self.e_left_loc = self.GetPlotParam('ElectronLeft')
        self.e_right_loc = self.GetPlotParam('ElectronRight')

        if self.GetPlotParam('PrtlIntegrationRelative'):
            self.e_left_loc += self.parent.shock_loc
            self.e_right_loc += self.parent.shock_loc


        self.i_left_loc = self.GetPlotParam('IonLeft')
        self.i_right_loc = self.GetPlotParam('IonRight')

        if self.GetPlotParam('PrtlIntegrationRelative'):
            self.i_left_loc += self.parent.shock_loc
            self.i_right_loc += self.parent.shock_loc

        # A list that will make sure that the data has the same int region
        self.region_args = list([self.e_left_loc, self.e_right_loc, self.i_left_loc, self.i_right_loc])

        is_loaded = False
        if self.keyname in self.parent.DataDict.keys():
            # make sure it is integrating over the correct region.
            if np.all(self.region_args == self.parent.DataDict[self.keyname][-1]):
                is_loaded = True

                #unpack the values
                self.gamma, self.edist, self.pdist, self.momentum, self.momedist, self.mompdist, region_args = self.parent.DataDict[self.keyname]

        if not is_loaded:

            # Calculate the momentum := gamma*beta
            #  NOTE: self.gamma is actually the real lorentz factor, gamma, minus 1

            self.momentum=np.sqrt((self.gamma+1)**2-1.)

            if self.GetPlotParam('rest_frame'):
                self.spece = np.copy(self.FigWrap.LoadKey('specerest'))
                self.specp = np.copy(self.FigWrap.LoadKey('specprest'))
            else:
                self.spece = np.copy(self.FigWrap.LoadKey('spece'))
                self.specp = np.copy(self.FigWrap.LoadKey('specp'))


            # In output.F90, spece (specp) is defined by the number of electons (ions)
            # divided by gamma in each logrithmic energy bin. So we multiply by gamma.

            for i in range(len(self.xsl)):
                self.spece[:,i] *= self.gamma
                self.specp[:,i] *= self.gamma

            ###############################
            ###### energy spectra, f=(dN/dE)/N
            ###############################
            self.dgamma = np.empty(len(self.gamma))
            delta=np.log10(self.gamma[-1]/self.gamma[0])/len(self.gamma)
            for i in range(len(self.dgamma)):
                self.dgamma[i]=self.gamma[i]*(10**delta-1.)

            eL = self.xsl.searchsorted(self.e_left_loc)
            eR = self.xsl.searchsorted(self.e_right_loc, side='right')

            iL = self.xsl.searchsorted(self.i_left_loc)
            iR = self.xsl.searchsorted(self.i_right_loc, side='right')

            if iL >= iR:
                iL = iR
                iR += 1
            if eL >= eR:
                eL = eR
                eR += 1

            # energy distribution, f(E)=(dn/dE)/N
            self.fe=np.empty(len(self.gamma))
            self.fp=np.empty(len(self.gamma))


            norme = np.ones(len(self.xsl))
            normp = np.ones(len(self.xsl))

            # total particles in each linear x bin
            for i in range(len(norme)):
                norme[i]=sum(self.spece[:,i])
                normp[i]=sum(self.specp[:,i])

            for k in range(len(self.fe)):
                self.fe[k]=sum(self.spece[k][eL:eR])
                self.fp[k]=sum(self.specp[k][iL:iR])

            if sum(norme[eL:eR]) > 0:
                if self.GetPlotParam('normalize_spectra'):
                    self.fe *= 1.0/sum(norme[eL:eR])

                self.edist=np.copy(self.fe)
                self.fe *= self.dgamma**(-1)
                self.femom=self.fe/(4*np.pi*self.momentum)/(self.gamma+1)
                self.momedist=self.femom*self.momentum**4

            else:
                print('RUNTIME WARNING: spectra.py can\'t find the electrons in the integration region, spectra will be incorrect')
                self.edist = 1E-200*np.ones(len(self.fe))
                self.momedist = 1E-200*np.ones(len(self.fe))

            if sum(normp[iL:iR]) > 0:
                if self.GetPlotParam('normalize_spectra'):
                    self.fp *= 1.0/sum(normp[iL:iR])

                self.pdist=np.copy(self.fp)
                self.fp *= self.dgamma**(-1)
                self.fpmom=self.fp/(4*np.pi*self.momentum)/(self.gamma+1)
                self.mompdist=self.fpmom*self.momentum**4

            else:
                print('RUNTIME WARNING: spectra.py can\'t find ions in the integration region, spectra will be incorrect')
                self.pdist = 1E-200*np.ones(len(self.fp))
                self.mompdist = 1E-200*np.ones(len(self.fp))

            self.parent.DataDict[self.keyname] = self.gamma, self.edist, self.pdist, self.momentum, self.momedist, self.mompdist, self.region_args

    def draw(self):
        self.RegionFlag = True
        # Set the tick color
        tick_color = 'black'
        # Create a gridspec to handle spacing better
        self.gs = gridspec.GridSpecFromSubplotSpec(100,100, subplot_spec = self.parent.gs0[self.FigWrap.pos])#, bottom=0.2,left=0.1,right=0.95, top = 0.95)

        self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]])

        self.ion_spect = self.axes.plot(self.momentum, self.mompdist, color = self.parent.ion_color)
        if not self.GetPlotParam('show_ions'):
            self.ion_spect[0].set_visible(False)
        self.boostedIonSpect = self.axes.plot(self.momentum*self.massRatio, self.mompdist*self.massRatio, color = self.parent.ion_color, ls='-.')
        if not self.GetPlotParam('BoostedIons'):
            self.boostedIonSpect[0].set_visible(False)
        self.ion_temp = self.axes.plot(self.momentum, self.mompdist,
                                       color = self.parent.ion_fit_color,
                                       linestyle = '--', linewidth = 1.5) # a placeholder

        # See if we want to plot the electron temperature
        if not self.GetPlotParam('SetTi'):
            self.ion_temp[0].set_visible(False)

         # a placeholder
        self.PLP = self.axes.plot(self.momentum, self.mompdist,
                                        color = self.parent.ion_fit_color,
                                        linestyle = ':', linewidth = 1.5)
        self.PLP[0].set_visible(False)


        self.electron_spect = self.axes.plot(self.momentum, self.momedist, color = self.parent.electron_color)
        if not self.GetPlotParam('show_electrons'):
            self.electron_spect[0].set_visible(False)

        # a placeholder
        self.electron_temp =  self.axes.plot(self.momentum, self.momedist,
                                            color = self.parent.electron_fit_color,
                                            linestyle = '--', linewidth = 1.5)
        if not self.GetPlotParam('SetTe'):
            self.electron_temp[0].set_visible(False)

        # a placeholder
        self.electron_Einj =  self.axes.axvline(0, color = self.parent.electron_fit_color,
                                            linestyle = '-.', linewidth = 1.5)
        if not self.GetPlotParam('MeasureEpsE'):
            self.electron_Einj.set_visible(False)


        # a placeholder
        self.ion_Einj =  self.axes.axvline(0, color = self.parent.ion_fit_color,
                                            linestyle = '-.', linewidth = 1.5)
        if not self.GetPlotParam('MeasureEpsP'):
            self.ion_Einj.set_visible(False)

        self.PLE = self.axes.plot(self.momentum, self.momedist,
                                  color = self.parent.electron_fit_color,
                                  linestyle = ':', linewidth = 1.5)
        self.PLE[0].set_visible(False)

        self.axes.set_xscale("log")
        self.axes.set_yscale("log")

        if int(matplotlib.__version__[0]) < 2:
            self.axes.set_axis_bgcolor('lightgrey')
        else:
            self.axes.set_facecolor('lightgrey')

        if self.GetPlotParam('set_xlim'):
            self.axes.set_xlim(self.GetPlotParam('x_min'), self.GetPlotParam('x_max'))
        if self.GetPlotParam('set_ylim'):
            self.axes.set_ylim(10**self.GetPlotParam('y_min'), 10**self.GetPlotParam('y_max'))

        self.axes.tick_params(labelsize = self.parent.MainParamDict['NumFontSize'], color=tick_color)

        if self.GetPlotParam('spectral_type') == 0:
            self.axes.set_xlabel(r'$\gamma\beta$', labelpad = self.parent.MainParamDict['xLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
            self.axes.set_ylabel(r'$p^4f(p)$', labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
        else:
            self.axes.set_xlabel(r'$\gamma-1$', labelpad = -2, color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
            self.axes.set_ylabel(r'$E(dn/dE)/n$', labelpad = 0, color = 'black', size = self.parent.MainParamDict['AxLabelSize'])

        self.refresh()

    def refresh(self):
        self.RegionFlag = True
        # First let set the dashes
        self.ion_spect[0].set_dashes(self.dashes)
        self.electron_spect[0].set_dashes(self.dashes)


        ############
        #
        # First we'll calculate eps_e or eps_e
        #
        ############

        if self.GetPlotParam('MeasureEpsP'):
            eps_p = self.measure_eps(self.pdist, self.GetPlotParam('GammaIonInjection'), 'proton')
            self.eps_pVar.set('%.6f' % eps_p)
        else:
            self.eps_pVar.set('N/A')

        if self.GetPlotParam('MeasureEpsE'):
            eps_e = self.measure_eps(self.edist, self.GetPlotParam('GammaElectronInjection'), 'electron')
            self.eps_eVar.set('%.6f' % eps_e)
        else:
            self.eps_eVar.set('N/A')
        if self.GetPlotParam('spectral_type') == 0: #Show the momentum dist
            if self.GetPlotParam('show_ions'):
                self.ion_spect[0].set_data(self.momentum, self.mompdist)
                if self.GetPlotParam('BoostedIons'):
                    self.boostedIonSpect[0].set_data(self.momentum*self.massRatio, self.mompdist*self.massRatio)
                if self.GetPlotParam('MeasureEpsP'):
                    self.ion_Einj.set_visible(True)
                    self.ion_Einj.set_xdata([np.sqrt((self.GetPlotParam('GammaIonInjection')+1)**2-1.),np.sqrt((self.GetPlotParam('GammaIonInjection')+1)**2-1.)])
                else:
                    self.ion_Einj.set_visible(False)
                if self.GetPlotParam('SetTi'):
                    self.ion_temp[0].set_visible(True)
                    if self.GetPlotParam('DelGami') >= 0.013:
                        aconst = 1/(self.GetPlotParam('DelGami')*np.exp(1.0/self.GetPlotParam('DelGami')*kn(2, 1.0/self.GetPlotParam('DelGami'))))
                    else:
                        aconst = np.sqrt(2/np.pi)/self.GetPlotParam('DelGami')**1.5
                        aconst -= 15.0/(4.*np.sqrt(self.GetPlotParam('DelGami'))*np.sqrt(2*np.pi))
                        aconst += (345*np.sqrt(self.GetPlotParam('DelGami')))/(64.*np.sqrt(2*np.pi))
                        aconst -= (3285*self.GetPlotParam('DelGami')**1.5)/(512.*np.sqrt(2*np.pi))
                        aconst += (95355*self.GetPlotParam('DelGami')**2.5)/(16384.*np.sqrt(2*np.pi))
                        aconst -= (232065*self.GetPlotParam('DelGami')**3.5)/(131072.*np.sqrt(2*np.pi))

                    fpmommax = self.momentum**4*aconst*(self.gamma+1.0)*np.sqrt((self.gamma+1.0)**2-1)
                    fpmommax *= np.exp(-self.gamma/self.GetPlotParam('DelGami'))/(4*np.pi*self.momentum)/(self.gamma+1.0)
                    fpmommax *= 10**self.GetPlotParam('iNormalizer')
                    self.ion_temp[0].set_data(self.momentum, fpmommax)

                else:
                    self.ion_temp[0].set_visible(False)

                self.PowerlawPworked = False
                self.PLP[0].set_visible(False)
                if self.GetPlotParam('DoPowerLawFitIon'):

                    # Find the part of the spectrum chosen by the fitting
                    #first convert to
                    mompleft = np.sqrt((self.GetPlotParam('PowerLawIonMin')+1)**2-1)
                    mompright = np.sqrt((self.GetPlotParam('PowerLawIonMax')+1.)**2-1.)
                    impLeft = self.momentum.searchsorted(mompleft)
                    impRight = self.momentum.searchsorted(mompright, side='Right')

                    # a while loop to make sure that the program won't mess up if mompdist == 0
                    while np.abs(self.mompdist[impLeft]) <= 1E-10 and impLeft < len(self.mompdist)-1:
                        impLeft += 1
                    while np.abs(self.mompdist[impRight]) <= 1E-10 and impRight > 0:
                        impRight -=1

                    if impRight>impLeft:
                        self.pslope, self.pintercept, pr_value, pp_value, pstderr = linregress(np.log(self.momentum[impLeft:impRight]), np.log(self.mompdist[impLeft:impRight]))

                        if np.isnan(self.pslope) == False and np.isnan(self.pintercept) == False:
                            self.PowerlawPworked = True
                            self.PLP[0].set_data(self.momentum, np.exp(self.pintercept)*self.momentum**self.pslope)
                            self.PLP[0].set_visible(True)

            if self.GetPlotParam('show_electrons'):
                self.electron_spect[0].set_data(self.momentum, self.momedist)
                if self.GetPlotParam('MeasureEpsE'):
                    self.electron_Einj.set_visible(True)
                    self.electron_Einj.set_xdata([np.sqrt((self.GetPlotParam('GammaElectronInjection')+1)**2-1.),np.sqrt((self.GetPlotParam('GammaElectronInjection')+1)**2-1.)])
                else:
                    self.electron_Einj.set_visible(False)

                if self.GetPlotParam('SetTe'):
                    self.electron_temp[0].set_visible(True)

                    self.delgame0=self.GetPlotParam('DelGame')*self.FigWrap.LoadKey('mi')[0]/self.FigWrap.LoadKey('me')[0]
                    if self.delgame0 >= 0.013:
                        aconst = 1/(self.delgame0*np.exp(1.0/self.delgame0)*kn(2, 1.0/self.delgame0))
                    else:
                        aconst = np.sqrt(2/np.pi)/self.delgame0**1.5
                        aconst -= 15.0/(4.*np.sqrt(self.delgame0)*np.sqrt(2*np.pi))
                        aconst += (345*np.sqrt(self.delgame0))/(64.*np.sqrt(2*np.pi))
                        aconst -= (3285*self.delgame0**1.5)/(512.*np.sqrt(2*np.pi))
                        aconst += (95355*self.delgame0**2.5)/(16384.*np.sqrt(2*np.pi))
                        aconst -= (232065*self.delgame0**3.5)/(131072.*np.sqrt(2*np.pi))

                    femommax = self.momentum**4*aconst*(self.gamma+1.0)*np.sqrt((self.gamma+1.0)**2-1)
                    femommax *= np.exp(-self.gamma/self.delgame0)/(4*np.pi*self.momentum)/(self.gamma+1.0)
                    femommax *= 10**self.GetPlotParam('eNormalizer')
                    self.electron_temp[0].set_data(self.momentum, femommax)


                else:
                    self.electron_temp[0].set_visible(False)

                # Now the power-law
                self.PowerlawEworked = False
                self.PLE[0].set_visible(False)
                if self.GetPlotParam('DoPowerLawFitElectron'):


                    # Find the part of the spectrum chosen by the fitting
                    #first convert to
                    momeleft = np.sqrt((self.GetPlotParam('PowerLawElectronMin')+1)**2-1)
                    momeright = np.sqrt((self.GetPlotParam('PowerLawElectronMax')+1.)**2-1.)
                    imeLeft = self.momentum.searchsorted(momeleft)
                    imeRight = self.momentum.searchsorted(momeright, side='Right')

                    # a while loop to make sure that the program won't mess up if momedist == 0
                    while np.abs(self.momedist[imeLeft]) <= 1E-14 and imeLeft < len(self.momedist)-1:
                        imeLeft += 1

                    while np.abs(self.momedist[imeRight]) <= 1E-14 and imeRight > 0:
                        imeRight -=1


                    if imeRight>imeLeft:
                        self.eslope, self.eintercept, r_value, p_value, stderr = linregress(np.log(self.momentum[imeLeft:imeRight]), np.log(self.momedist[imeLeft:imeRight]))

                        if np.isnan(self.eslope) == False and np.isnan(self.eintercept) == False:
                            self.PowerlawEworked = True
                            self.PLE[0].set_data(self.momentum, np.exp(self.eintercept)*self.momentum**self.eslope)
                            self.PLE[0].set_visible(True)

            self.axes.set_xlim(self.GetPlotParam('x_min'), self.GetPlotParam('x_max'))
            self.axes.set_ylim(10**self.GetPlotParam('y_min'), 10**self.GetPlotParam('y_max'))


        #####
        #
        # The momentum section is done, now we see what happens if we use energy spectra instead
        #
        #####


        if self.GetPlotParam('spectral_type') == 1: #Show the energy dist
            if self.GetPlotParam('show_electrons'):
                self.electron_spect[0].set_data(self.gamma, self.edist)
                if self.GetPlotParam('MeasureEpsE'):
                    self.electron_Einj.set_visible(True)
                    self.electron_Einj.set_xdata([self.GetPlotParam('GammaElectronInjection'),self.GetPlotParam('GammaElectronInjection')])
                else:
                    self.electron_Einj.set_visible(False)


                # The temperature
                if self.GetPlotParam('SetTe'):
                    self.electron_temp[0].set_visible(True)
                    self.delgame0=self.GetPlotParam('DelGame')*self.FigWrap.LoadKey('mi')[0]/self.FigWrap.LoadKey('me')[0]
                    if self.delgame0 >= 0.013:
                        aconst = 1/(self.delgame0*np.exp(1.0/self.delgame0)*kn(2, 1.0/self.delgame0))
                    else:
                        aconst = np.sqrt(2/np.pi)/self.delgame0**1.5
                        aconst -= 15.0/(4.*np.sqrt(self.delgame0)*np.sqrt(2*np.pi))
                        aconst += (345*np.sqrt(self.delgame0))/(64.*np.sqrt(2*np.pi))
                        aconst -= (3285*self.delgame0**1.5)/(512.*np.sqrt(2*np.pi))
                        aconst += (95355*self.delgame0**2.5)/(16384.*np.sqrt(2*np.pi))
                        aconst -= (232065*self.delgame0**3.5)/(131072.*np.sqrt(2*np.pi))

                    femax = aconst*self.gamma*(self.gamma+1.0)*np.sqrt((self.gamma+1.0)**2-1)
                    femax *= np.exp(-self.gamma/self.delgame0)
                    femax *= 10**self.GetPlotParam('eNormalizer')
                    self.electron_temp[0].set_data(self.gamma, femax)

                else:
                    self.electron_temp[0].set_visible(False)

                # the power-law
                self.PowerlawEworked = False
                self.PLE[0].set_visible(False)

                if self.GetPlotParam('DoPowerLawFitElectron'):
                    # Find the part of the spectrum chosen by the fitting
                    #first convert to
                    ieLeft = self.gamma.searchsorted(self.GetPlotParam('PowerLawElectronMin'))
                    ieRight = self.gamma.searchsorted(self.GetPlotParam('PowerLawElectronMax'), side='Right')
                    if ieLeft == len(self.gamma):
                        ieLeft -= 1

                    if ieRight == len(self.gamma):
                        ieRight -= 1
                    # a while loop to make sure that the program won't mess up if momedist == 0
                    while np.abs(self.edist[ieLeft]) <= 1E-14 and ieLeft < len(self.edist)-1:
                        ieLeft += 1

                    while np.abs(self.edist[ieRight]) <= 1E-14 and ieRight > 0:
                        ieRight -=1

                    if ieRight>ieLeft:
                        self.eslope, self.eintercept, er_value, ep_value, estderr = linregress(np.log(self.gamma[ieLeft:ieRight]), np.log(self.edist[ieLeft:ieRight]))

                        if np.isnan(self.eslope) == False and np.isnan(self.eintercept) == False:
                            self.PowerlawEworked = True

                            self.PLE[0].set_data(self.gamma, np.exp(self.eintercept)*self.gamma**self.eslope)
                            self.PLE[0].set_visible(True)

            if self.GetPlotParam('show_ions'):
                self.ion_spect[0].set_data(self.gamma, self.pdist)
                if self.GetPlotParam('BoostedIons'):
                    self.boostedIonSpect[0].set_data(self.gamma*self.massRatio, self.pdist)
                if self.GetPlotParam('MeasureEpsP'):
                    self.ion_Einj.set_visible(True)
                    self.ion_Einj.set_xdata([self.GetPlotParam('GammaIonInjection'),self.GetPlotParam('GammaIonInjection') ])
                else:
                    self.ion_Einj.set_visible(False)

                if self.GetPlotParam('SetTi'):
                    self.ion_temp[0].set_visible(True)
                    if self.GetPlotParam('DelGami') >= 0.013:
                        aconst = 1/(self.GetPlotParam('DelGami')*np.exp(1.0/self.GetPlotParam('DelGami'))*kn(2, 1.0/self.GetPlotParam('DelGami')))
                    else:
                        aconst = np.sqrt(2/np.pi)/self.GetPlotParam('DelGami')**1.5
                        aconst -= 15.0/(4.*np.sqrt(self.GetPlotParam('DelGami'))*np.sqrt(2*np.pi))
                        aconst += (345*np.sqrt(self.GetPlotParam('DelGami')))/(64.*np.sqrt(2*np.pi))
                        aconst -= (3285*self.GetPlotParam('DelGami')**1.5)/(512.*np.sqrt(2*np.pi))
                        aconst += (95355*self.GetPlotParam('DelGami')**2.5)/(16384.*np.sqrt(2*np.pi))
                        aconst -= (232065*self.GetPlotParam('DelGami')**3.5)/(131072.*np.sqrt(2*np.pi))

                    fpmax = aconst*self.gamma*(self.gamma+1.0)*np.sqrt((self.gamma+1.0)**2-1)
                    fpmax *= np.exp(-self.gamma/self.GetPlotParam('DelGami'))
                    fpmax *= 10**self.GetPlotParam('iNormalizer')
                    self.ion_temp[0].set_data(self.gamma, fpmax)
                else:
                    self.ion_temp[0].set_visible(False)
                self.PowerlawPworked = False
                if not self.GetPlotParam('DoPowerLawFitIon'):
                    self.PLP[0].set_visible(False)
                else:
                    # Find the part of the spectrum chosen by the fitting
                    iepLeft = self.gamma.searchsorted(self.GetPlotParam('PowerLawIonMin'))
                    iepRight = self.gamma.searchsorted(self.GetPlotParam('PowerLawIonMax'), side='Right')

                    if iepLeft == len(self.gamma):
                        iepLeft -= 1
                    if iepRight == len(self.gamma):
                        iepRight -=1
                    # a while loop to make sure that the program won't mess up if mompdist == 0
                    while np.abs(self.pdist[iepLeft]) <= 1E-14 and iepLeft < len(self.pdist)-1:
                        iepLeft += 1
                    while np.abs(self.pdist[iepRight]) <= 1E-14 and iepRight > 0:
                        iepRight -=1

                    if iepRight>iepLeft:
                        self.pslope, self.pintercept, per_value, pep_value, pestderr = linregress(np.log(self.gamma[iepLeft:iepRight]), np.log(self.pdist[iepLeft:iepRight]))
                        if np.isnan(self.pslope) == False and np.isnan(self.pintercept) == False:
                            self.PowerlawPworked = True
                            self.PLP[0].set_data(self.gamma, np.exp(self.pintercept)*self.gamma**self.pslope)
                            self.PLP[0].set_visible(True)

                self.axes.set_xlim(self.GetPlotParam('x_min'), self.GetPlotParam('x_max'))
                self.axes.set_ylim(10**self.GetPlotParam('y_min'), 10**self.GetPlotParam('y_max'))

        if self.GetPlotParam('show_legend'):
            self.MakeLegend()

    def measure_eps(self, energy_dist, e_inj, prtl_type):
        injloc=self.gamma.searchsorted(e_inj)

        if injloc >= len(energy_dist)-1:
            return 0.0

        else:
            '''
            HERE is a trapezoidal integration in logspace not using anymore
            because it goes bad when energy_dist is zero. I'm leaving it in, commented out in case it is
            valuable later.
            # do a cummulative
            LogY = np.log(energy_dist[e_injloc-1:])
            # Make all values where LogY ==-inf == -1E4 effectively set them to zero
            LogY[np.isinf(LogY)] = -1E4
            LogX = np.log(self.gamma[e_injloc-1:])
            dLogY = np.diff(LogY)
            dLogX = np.diff(LogX)
            dF = np.exp(LogY[:-1]+LogX[:-1])
            dF *= dLogX*(np.exp(dLogX+dLogY)-1)
            dF *= (dLogX+dLogY)**(-1)
            endSumF = np.cumsum(dF[::-1])[::-1]
            '''
            # Do a reverse cummulative trapz integration over the energy_dist
            endSum = cumtrapz(energy_dist[::-1],self.gamma[::-1])[::-1]
            endSum *= -1

            logF = np.nan

            if injloc ==0:
                logF = np.log(endSum[0])
            else:
                # The answer we want will be between the [injloc-1] and [injloc] place in
                # endSumF We do linear interpolation in logspace to find the answer
                logY_left = np.log(endSum[injloc-1])
                logY_right = np.log(endSum[injloc])
                logX_left = np.log(self.gamma[injloc-1])
                logX_right = np.log(self.gamma[injloc])
                # do the linear interpolation in logspace
                logF = np.log(e_inj)-logX_left
                logF *= (logX_right-logX_left)**(-1)
                logF *= logY_right-logY_left
                logF += logY_left

            # Now calculate the eps
            if prtl_type == 'electron':
                psum = np.trapz(self.pdist,self.gamma)
                return np.exp(logF)/psum * self.FigWrap.LoadKey('me')[0]/self.FigWrap.LoadKey('mi')[0]
            else:
                return np.exp(logF)/endSum[0]

    def GetPlotParam(self, keyname):
        return self.FigWrap.GetPlotParam(keyname)

    def SetPlotParam(self, keyname, value, update_plot = True, NeedsRedraw = False):#, refresh = False):
        self.FigWrap.SetPlotParam(keyname, value, update_plot = update_plot, NeedsRedraw = NeedsRedraw)
        #if refresh == True:
        #    self.refresh()
        #    self.parent.canvas.draw()
        #    self.parent.canvas.get_tk_widget().update_idletasks()

    def MakeLegend(self):
        ''' A helper function to make the legend'''
        # First make the temperature legend
        Tlegend_handles = []
        Tlegend_labels = []
        # check if the legend already exists. It it does, we need to remove it.

        try:
            self.legT.remove()
        except:
            pass


        if self.GetPlotParam('show_electrons') and self.GetPlotParam('SetTe'):
            Tlegend_handles.append(self.electron_temp[0])
            tmpstr = '%.3f' % self.delgame0
            Tlegend_labels.append(r'$T_e\ = $' +  ' ' + tmpstr + ' ' + r'$m_e c^2$')
        if self.GetPlotParam('show_ions') and self.GetPlotParam('SetTi'):
            Tlegend_handles.append(self.ion_temp[0])
            tmpcon =self.GetPlotParam('DelGami')*self.FigWrap.LoadKey('mi')[0]/self.FigWrap.LoadKey('me')[0]
            tmpstr = '%.3f' % tmpcon
            Tlegend_labels.append(r'$T_p\ = $' +  ' ' + tmpstr + ' ' + r'$m_e c^2$')


        # now make the power-fit legend
        legend_handles = []
        legend_labels = []
        if self.GetPlotParam('show_electrons') and self.PowerlawEworked:
            legend_handles.append(self.PLE[0])
            if self.GetPlotParam('spectral_type') == 0:
                tmpnum = 4.0-self.eslope
                tmpstr = '%.1f' % tmpnum
            else:
                tmpnum = 1-self.eslope
                tmpstr = '%.1f' % tmpnum
            legend_labels.append(r'$\delta_e\ = $' +  ' ' + tmpstr)
        if self.GetPlotParam('show_ions') and self.PowerlawPworked:
            legend_handles.append(self.PLP[0])
            if self.GetPlotParam('spectral_type') == 0:
                tmpnum = 4.0 - self.pslope
                tmpstr = '%.1f' % tmpnum
            else:
                tmpnum = 1-self.pslope
                tmpstr = '%.1f' % tmpnum

            legend_labels.append(r'$\delta_p\ = $' +  ' ' + tmpstr)


        # Draw the legend, there is a complication here because each plot can only have one legend.
        #if len(Tlegend_handles)> 0 and len(legend_handles) == 0:
        #    self.legT = self.axes.legend(Tlegend_handles, Tlegend_labels, framealpha = .05, fontsize = 11, loc = 'upper left')
        #    self.legT.get_frame().set_facecolor('k')
        #    self.legT.get_frame().set_linewidth(0.0)

        if len(Tlegend_handles) > 0:
            self.legT = self.axes.legend(Tlegend_handles, Tlegend_labels, framealpha = .05, fontsize = self.parent.MainParamDict['legendLabelSize'], loc = 2)
            self.legT.get_frame().set_facecolor('k')
            self.legT.get_frame().set_linewidth(0.0)
            self.legT.draggable(update = 'loc')
            if self.GetPlotParam('T_legend_loc') != 'N/A':
                tmp_tup = float(self.GetPlotParam('T_legend_loc').split()[0]),float(self.GetPlotParam('T_legend_loc').split()[1])
                self.legT._set_loc(tmp_tup)

            if len(legend_handles) != 0:
                self.axes.add_artist(self.legT)

        if len(legend_handles)> 0:
            self.legDelta = self.axes.legend(legend_handles, legend_labels,
            framealpha = .05, fontsize = self.parent.MainParamDict['legendLabelSize'], loc = 'upper right')
            self.legDelta.get_frame().set_facecolor('k')
            self.legDelta.get_frame().set_linewidth(0.0)
            self.legDelta.draggable()
            if self.GetPlotParam('PL_legend_loc') != 'N/A':
                tmp_tup = float(self.GetPlotParam('PL_legend_loc').split()[0]),float(self.GetPlotParam('PL_legend_loc').split()[1])
                self.legDelta._set_loc(tmp_tup)


    def OpenSettings(self):
        if self.settings_window is None:
            self.settings_window = SpectraSettings(self)
        else:
            self.settings_window.destroy()
            self.settings_window = SpectraSettings(self)


class SpectraSettings(Tk.Toplevel):
    def __init__(self, parent):
        self.parent = parent
        Tk.Toplevel.__init__(self)

        self.wm_title('Spectrum (%d,%d) Settings' % self.parent.FigWrap.pos)
        nb = ttk.Notebook(self) # This allows to make our subplot settings tabbed.

        frm = ttk.Frame(nb) # frm is the Spectral Chart settings
        frm2 = ttk.Frame(nb) # frm2 holds the measurement info
        nb.add(frm, text = 'Settings')
        nb.add(frm2, text = 'Measurements')
        nb.pack(fill=Tk.BOTH, expand = True)

        self.protocol('WM_DELETE_WINDOW', self.OnClosing)
        #Create some sizers

        self.bind('<Return>', self.TxtEnter)

        #####
        #
        # FIRST INITIALIZE THE STUFF FOR THE CHART SETTINGS!!!
        #
        #####

        # Create the OptionMenu to chooses the Chart Type:
        self.ctypevar = Tk.StringVar(self)
        self.ctypevar.set(self.parent.chartType) # default value
        self.ctypevar.trace('w', self.ctypeChanged)

        ttk.Label(frm, text="Choose Chart Type:").grid(row=0, column = 0)
        cmapChooser = ttk.OptionMenu(frm, self.ctypevar, self.parent.chartType, *tuple(self.parent.ChartTypes))
        cmapChooser.grid(row =0, column = 1, sticky = Tk.W + Tk.E)


        # the Radiobox Control to choose the Field Type
        self.SpectList = ['Momentum', 'Energy']
        self.SpectTypeVar  = Tk.IntVar()
        self.SpectTypeVar.set(self.parent.GetPlotParam('spectral_type'))

        ttk.Label(frm, text='Choose Spectrum Type:').grid(row = 2, sticky = Tk.W)

        for i in range(len(self.SpectList)):
            ttk.Radiobutton(frm,
                text=self.SpectList[i],
                variable=self.SpectTypeVar,
                command = self.RadioSpect,
                value=i).grid(row = 3+i, sticky =Tk.W)


        # show ions
        self.IonVar = Tk.IntVar()
        self.IonVar.set(self.parent.GetPlotParam('show_ions'))
        cb = ttk.Checkbutton(frm, text = "Show ions",
                        variable = self.IonVar,
                        command = self.IonVarHandler)
        cb.grid(row = 6, column = 0, sticky = Tk.W)

        # show electrons
        self.eVar = Tk.IntVar()
        self.eVar.set(self.parent.GetPlotParam('show_electrons'))
        cb = ttk.Checkbutton(frm, text = "Show electrons",
                        variable = self.eVar,
                        command = self.eVarHandler)
        cb.grid(row = 6, column = 1, sticky = Tk.W)

        # Normalize spectra
        self.NormalizeVar = Tk.IntVar()
        self.NormalizeVar.set(self.parent.GetPlotParam('normalize_spectra'))
        cb = ttk.Checkbutton(frm, text = "Normalize Spectra",
                        variable = self.NormalizeVar,
                        command = self.NormalizeVarHandler)
        cb.grid(row = 7, column = 1, sticky = Tk.W)

        # show in rest frame
        self.RestVar = Tk.IntVar()
        self.RestVar.set(self.parent.GetPlotParam('rest_frame'))
        cb = ttk.Checkbutton(frm, text = "Show in rest frame",
                        variable = self.RestVar,
                        command = lambda:
                        self.parent.SetPlotParam('rest_frame', self.RestVar.get()))
        cb.grid(row = 7, column = 0, sticky = Tk.W)

        # show in rest frame
        self.BoostVar = Tk.IntVar()
        self.BoostVar.set(self.parent.GetPlotParam('BoostedIons'))
        cb = ttk.Checkbutton(frm, text = "Show ions spect shifted by mi/me",
                        variable = self.BoostVar,
                        command =self.BoostHandler)
        cb.grid(row = 8, column = 0, sticky = Tk.W)

        self.xLimVar = Tk.IntVar()
        self.xLimVar.set(self.parent.GetPlotParam('set_xlim'))
        self.xLimVar.trace('w', self.xLimChanged)



        self.xmin = Tk.StringVar()
        self.xmin.set(str(self.parent.GetPlotParam('x_min')))
        self.xmax = Tk.StringVar()
        self.xmax.set(str(self.parent.GetPlotParam('x_max')))


        cb = ttk.Checkbutton(frm, text ='Set xlim',
                        variable = self.xLimVar)
        cb.grid(row = 4, column =3,sticky = Tk.W)
#        ttk.Label(frm, text = 'Set xlim').grid(row = 4, column = 3, sticky = Tk.W)
        self.eLEnter = ttk.Entry(frm, textvariable=self.xmin, width=7)
        self.eLEnter.grid(row = 4, column =4)
        self.eREnter = ttk.Entry(frm, textvariable=self.xmax, width=7)
        self.eREnter.grid(row = 4, column =5)


        self.yLimVar = Tk.IntVar()
        self.yLimVar.set(self.parent.GetPlotParam('set_ylim'))
        self.yLimVar.trace('w', self.yLimChanged)



        self.ymin = Tk.StringVar()
        self.ymin.set(str(self.parent.GetPlotParam('y_min')))
        self.ymax = Tk.StringVar()
        self.ymax.set(str(self.parent.GetPlotParam('y_max')))


        cb = ttk.Checkbutton(frm, text ='Set log(y) lim',
                        variable = self.yLimVar)
        cb.grid(row = 5,  column =3, sticky = Tk.W)
#        ttk.Label(frm, text = 'Set log(y) lim').grid(row = 5, column =3, sticky = Tk.W)
        self.eLEnter = ttk.Entry(frm, textvariable=self.ymin, width=7)
        self.eLEnter.grid(row = 5, column =4)
        self.eREnter = ttk.Entry(frm, textvariable=self.ymax, width=7)
        self.eREnter.grid(row = 5, column =5)

        #####
        #
        # NOW WE INITALIZE STUFF FOR THE MEASUREMENTS
        #
        #####

        # Make an entry to change the integration region
        # A StringVar for a box to type in a value for the left ion region
        self.ileft = Tk.StringVar()
        # set it to the left value
        self.ileft.set(str(self.parent.GetPlotParam('IonLeft')))

        # A StringVar for a box to type in a value for the right ion region
        self.iright = Tk.StringVar()
        # set it to the right value
        self.iright.set(str(self.parent.GetPlotParam('IonRight')))

        # Now the electrons
        self.eleft = Tk.StringVar()
        self.eleft.set(str(self.parent.GetPlotParam('ElectronLeft')))
        self.eright = Tk.StringVar()
        self.eright.set(str(self.parent.GetPlotParam('ElectronRight')))

        ttk.Label(frm2, text='Energy Int region:').grid(row = 0, sticky = Tk.W)
        ttk.Label(frm2, text='left').grid(row = 0, column = 1, sticky = Tk.N)
        ttk.Label(frm2, text='right').grid(row = 0, column = 2, sticky = Tk.N)

        # the ion row
        ttk.Label(frm2, text='ions').grid(row= 1, sticky = Tk.N)
        # Make an button to change the wait time

        self.iLEnter = ttk.Entry(frm2, textvariable=self.ileft, width=7)
        self.iLEnter.grid(row =1, column = 1, sticky = Tk.W + Tk.E)

        self.iREnter = ttk.Entry(frm2, textvariable=self.iright, width=7)
        self.iREnter.grid(row = 1, column =2, sticky = Tk.W + Tk.E)

        ttk.Label(frm2, text='electrons').grid(row= 2, sticky = Tk.N)
        self.eLEnter = ttk.Entry(frm2, textvariable=self.eleft, width=7)
        self.eLEnter.grid(row = 2, column =1, sticky = Tk.W + Tk.E)
        self.eREnter = ttk.Entry(frm2, textvariable=self.eright, width=7)
        self.eREnter.grid(row = 2, column =2, sticky = Tk.W + Tk.E)

        self.RelVar = Tk.IntVar()
        self.RelVar.set(self.parent.GetPlotParam('PrtlIntegrationRelative'))
        self.RelVar.trace('w', self.RelChanged)
        cb = ttk.Checkbutton(frm2, text = "Energy Region relative to shock?",
                        variable = self.RelVar)
        cb.grid(row = 3, columnspan = 3, sticky = Tk.W)

        self.SetTeVar = Tk.IntVar()
        self.SetTeVar.set(self.parent.GetPlotParam('SetTe'))
        self.SetTeVar.trace('w', self.SetTeChanged)
        cb = ttk.Checkbutton(frm2, text='Show T_e', variable =  self.SetTeVar)
        cb.grid(row = 5, sticky = Tk.W)

        ttk.Label(frm2, text=u'\u0394'+u'\u0263' + ' =').grid(row= 5, column =1, sticky = Tk.N)

        self.SetTpVar = Tk.IntVar()
        self.SetTpVar.set(self.parent.GetPlotParam('SetTi'))
        self.SetTpVar.trace('w', self.SetTpChanged)

        cb = ttk.Checkbutton(frm2, text='Show T_i', variable =  self.SetTpVar)
        cb.grid(row = 6, sticky = Tk.W)
        ttk.Label(frm2, text=u'\u0394'+u'\u0263' + ' =').grid(row= 6, column =1, sticky = Tk.N)

        self.delgameVar = Tk.StringVar()
        self.delgameVar.set(str(self.parent.GetPlotParam('DelGame')))
        self.delgampVar = Tk.StringVar()
        self.delgampVar.set(str(self.parent.GetPlotParam('DelGami')))


        ttk.Entry(frm2, textvariable=self.delgameVar, width = 7).grid(row = 5, column = 2, sticky = Tk.N)
        ttk.Entry(frm2, textvariable=self.delgampVar, width = 7).grid(row = 6, column =2, sticky = Tk.N)

        #####
        #
        # ADD NORMALIZERS
        #
        #####

        self.eTempNormVar = Tk.StringVar()
        self.eTempNormVar.set(self.parent.GetPlotParam('eNormalizer'))
        labele = ttk.Label(frm2, text='log10(Norm)')#
        labele.grid(row = 5, column = 3)#, expand=0)


        # A slider that will select the 2D slice in the simulation


        self.txtEnterTe = ttk.Entry(frm2, textvariable=self.eTempNormVar, width=6)
        self.txtEnterTe.grid(row = 5, column = 4)
        self.sliderTe = ttk.Scale(frm2, from_=-6, to=6, command = self.TeScaleHandler)
        self.sliderTe.set(self.eTempNormVar.get())
        self.sliderTe.grid(row = 5, column = 5)#, expand=1)

        self.iTempNormVar = Tk.StringVar()
        self.iTempNormVar.set(self.parent.GetPlotParam('iNormalizer'))
        labeli = ttk.Label(frm2, text='log10(Norm)')#
        labeli.grid(row = 6, column = 3)#, expand=0)


        # A slider that will select the 2D slice in the simulation


        self.txtEnterTi = ttk.Entry(frm2, textvariable=self.iTempNormVar, width=6)
        self.txtEnterTi.grid(row = 6, column = 4)
        self.sliderTi = ttk.Scale(frm2, from_=-6, to=6, command = self.TiScaleHandler)
        self.sliderTi.set(self.iTempNormVar.get())
        self.sliderTi.grid(row = 6, column = 5)#, expand=1)

        # bind releasing the moust button to updating the plots.
        #self.slidery.bind("<ButtonRelease-1>", self.yUpdateValue)

        ttk.Label(frm2, text='Powerlaw fits:').grid(row = 8, sticky = Tk.W)
        ttk.Label(frm2, text='E_min [mc^2]').grid(row = 8, column = 1, sticky = Tk.N)
        ttk.Label(frm2, text='E_max [mc^2]').grid(row = 8, column = 2, sticky = Tk.N)

        self.PLFitEVar = Tk.IntVar()
        self.PLFitEVar.set(self.parent.GetPlotParam('DoPowerLawFitElectron'))
        self.PLFitEVar.trace('w', self.PLFitEChanged)
        ttk.Checkbutton(frm2, text='Electrons', variable =  self.PLFitEVar).grid(row = 9, sticky = Tk.W)

        self.E1Var = Tk.StringVar()
        self.E1Var.set(str(self.parent.GetPlotParam('PowerLawElectronMin')))
        self.E2Var = Tk.StringVar()
        self.E2Var.set(str(self.parent.GetPlotParam('PowerLawElectronMax')))


        ttk.Entry(frm2, textvariable=self.E1Var, width = 7).grid(row = 9, column = 1, sticky = Tk.N)
        ttk.Entry(frm2, textvariable=self.E2Var, width = 7).grid(row = 9, column =2, sticky = Tk.N)


        self.PLFitPVar = Tk.IntVar()
        self.PLFitPVar.set(self.parent.GetPlotParam('DoPowerLawFitIon'))
        self.PLFitPVar.trace('w', self.PLFitPChanged)
        ttk.Checkbutton(frm2, text='Ions', variable =  self.PLFitPVar).grid(row = 10, sticky = Tk.W)

        self.P1Var = Tk.StringVar()
        self.P1Var.set(str(self.parent.GetPlotParam('PowerLawIonMin')))
        self.P2Var = Tk.StringVar()
        self.P2Var.set(str(self.parent.GetPlotParam('PowerLawIonMax')))

        ttk.Entry(frm2, textvariable=self.P1Var, width = 7).grid(row = 10, column = 1, sticky = Tk.N)
        ttk.Entry(frm2, textvariable=self.P2Var, width = 7).grid(row = 10, column =2, sticky = Tk.N)


        ttk.Label(frm2, text='Measure eps:').grid(row = 11, column = 0, sticky = Tk.W)
        ttk.Label(frm2, text='E_inj [mc^2]').grid(row = 11, column = 1, sticky = Tk.N)
        ttk.Label(frm2, text='eps').grid(row = 11, column = 2, sticky = Tk.N)

        self.eps_p_fitVar = Tk.IntVar()
        self.eps_p_fitVar.set(self.parent.GetPlotParam('MeasureEpsP'))
        self.eps_p_fitVar.trace('w', self.eps_pFitChanged)
        ttk.Checkbutton(frm2, text='protons', variable =  self.eps_p_fitVar).grid(row = 12, sticky = Tk.W)

        self.EinjPVar = Tk.StringVar()
        self.EinjPVar.set(str(self.parent.GetPlotParam('GammaIonInjection')))
        ttk.Entry(frm2, textvariable=self.EinjPVar, width = 7).grid(row = 12, column = 1, sticky = Tk.N)
        ttk.Entry(frm2, textvariable=self.parent.eps_pVar, width = 7, state = 'readonly').grid(row = 12, column =2, sticky = Tk.N)

        self.eps_e_fitVar = Tk.IntVar()
        self.eps_e_fitVar.set(self.parent.GetPlotParam('MeasureEpsE'))
        self.eps_e_fitVar.trace('w', self.eps_eFitChanged)
        ttk.Checkbutton(frm2, text='electrons', variable =  self.eps_e_fitVar).grid(row = 13, sticky = Tk.W)

        self.EinjEVar = Tk.StringVar()
        self.EinjEVar.set(str(self.parent.GetPlotParam('GammaElectronInjection')))
        ttk.Entry(frm2, textvariable=self.EinjEVar, width = 7).grid(row = 13, column = 1, sticky = Tk.N)
        ttk.Entry(frm2, textvariable=self.parent.eps_eVar, width = 7, state = 'readonly').grid(row = 13, column =2, sticky = Tk.N)

    def TeScaleHandler(self, e):
        # if changing the scale will change the value of the parameter, do so
        if np.abs(float(self.eTempNormVar.get()) - self.sliderTe.get()) > 1E-4:
            self.eTempNormVar.set(self.sliderTe.get())
            self.parent.SetPlotParam('eNormalizer', self.sliderTe.get(), update_plot= False)
            self.parent.refresh()
            self.parent.parent.canvas.draw()
            self.parent.parent.canvas.get_tk_widget().update_idletasks()
    def BoostHandler(self):
        # if changing the scale will change the value of the parameter, do so
        if self.BoostVar.get() != self.parent.GetPlotParam('BoostedIons'):
            self.parent.boostedIonSpect[0].set_visible(self.BoostVar.get())
            self.parent.SetPlotParam('BoostedIons', self.BoostVar.get(), update_plot= False)
            self.parent.refresh()
            self.parent.parent.canvas.draw()
            self.parent.parent.canvas.get_tk_widget().update_idletasks()

    def TiScaleHandler(self, e):
        # if changing the scale will change the value of the parameter, do so
        if np.abs(float(self.iTempNormVar.get()) - self.sliderTi.get()) > 1E-4:
            self.iTempNormVar.set(self.sliderTi.get())
            self.parent.SetPlotParam('iNormalizer', self.sliderTi.get(), update_plot= False)
            self.parent.refresh()
            self.parent.parent.canvas.draw()
            self.parent.parent.canvas.get_tk_widget().update_idletasks()

    def ctypeChanged(self, *args):
        if self.ctypevar.get() == self.parent.chartType:
            pass
        else:
            self.parent.ChangePlotType(self.ctypevar.get())
            self.destroy()
    def IonVarHandler(self, *args):
        if self.IonVar.get() == self.parent.GetPlotParam('show_ions'):
            pass
        else:
            self.parent.ion_spect[0].set_visible(self.IonVar.get())
            # I could do a better job here...
            self.parent.SetPlotParam('show_ions', self.IonVar.get(), NeedsRedraw = True)

    def eVarHandler(self, *args):
        if self.eVar.get() == self.parent.GetPlotParam('show_electrons'):
            pass
        else:
            self.parent.electron_spect[0].set_visible(self.eVar.get())
            self.parent.SetPlotParam('show_electrons', self.eVar.get(), NeedsRedraw = True)


    def NormalizeVarHandler(self, *args):
        if self.NormalizeVar.get() == self.parent.GetPlotParam('normalize_spectra'):
            pass
        else:
            self.parent.SetPlotParam('normalize_spectra', self.NormalizeVar.get())

    def RadioSpect(self):
        if self.SpectTypeVar.get() == self.parent.GetPlotParam('spectral_type'):
            pass
        else:
            if self.SpectTypeVar.get() == 1:
                self.parent.axes.set_xlabel(r'$\gamma-1$', size = self.parent.parent.MainParamDict['AxLabelSize'])
                self.parent.axes.set_ylabel(r'$E(dn/dE)/n$', size = self.parent.parent.MainParamDict['AxLabelSize'])
                self.parent.SetPlotParam('x_min', 0.0005, update_plot = False)
                self.parent.SetPlotParam('x_max', 100, update_plot = False)

            else:
                self.parent.axes.set_xlabel(r'$\gamma\beta$', size = self.parent.parent.MainParamDict['AxLabelSize'])
                self.parent.axes.set_ylabel(r'$p^4f(p)$', size = self.parent.parent.MainParamDict['AxLabelSize'])
                self.parent.SetPlotParam('x_min', 0.05, update_plot = False)
                self.parent.SetPlotParam('x_max', 200, update_plot = False)

            self.xmin.set(str(self.parent.GetPlotParam('x_min')))
            self.xmax.set(str(self.parent.GetPlotParam('x_max')))
            self.parent.SetPlotParam('spectral_type', self.SpectTypeVar.get())


    def TxtEnter(self, e):
        to_reload = False
        to_reload += self.FieldsCallback()
        to_reload += self.MeasuresCallback()
        if to_reload:
            self.parent.SetPlotParam('x_min', self.parent.GetPlotParam('x_min'))
    def eps_pFitChanged(self, *args):
        if self.eps_p_fitVar.get() == self.parent.GetPlotParam('MeasureEpsP'):
            pass
        else:
            self.parent.SetPlotParam('MeasureEpsP',self.eps_p_fitVar.get())

    def eps_eFitChanged(self, *args):
        if self.eps_e_fitVar.get() == self.parent.GetPlotParam('MeasureEpsE'):
            pass
        else:
            self.parent.SetPlotParam('MeasureEpsE', self.eps_e_fitVar.get())


    def PLFitEChanged(self, *args):
        if self.PLFitEVar.get() == self.parent.GetPlotParam('DoPowerLawFitElectron'):
            pass
        else:
            self.parent.SetPlotParam('DoPowerLawFitElectron', self.PLFitEVar.get())

    def PLFitPChanged(self, *args):
        if self.PLFitPVar.get() == self.parent.GetPlotParam('DoPowerLawFitIon'):
            pass
        else:
            self.parent.SetPlotParam('DoPowerLawFitIon', self.PLFitPVar.get())

    def CheckIfTeChanged(self):
        to_reload = False
        try:
            # make sure the user types in a float
            if np.abs(float(self.delgameVar.get()) - self.parent.GetPlotParam('DelGame')) > 1E-12:
                self.parent.SetPlotParam('DelGame', float(self.delgameVar.get()), update_plot = False)
                self.parent.refresh()
                self.parent.parent.canvas.draw()
                self.parent.parent.canvas.get_tk_widget().update_idletasks()

        except ValueError:
            #if they type in random stuff, just set it ot the param value
            self.delgameVar.set(str(self.parent.GetPlotParam('DelGame')))
        return to_reload

    def CheckIfTpChanged(self):
        to_reload = False
        try:
            # make sure the user types in a float
            if np.abs(float(self.delgampVar.get()) - self.parent.GetPlotParam('DelGami')) > 1E-12:
                self.parent.SetPlotParam('DelGami', float(self.delgampVar.get()), update_plot=False)
                self.parent.refresh()
                self.parent.parent.canvas.draw()
                self.parent.parent.canvas.get_tk_widget().update_idletasks()



        except ValueError:
            #if they type in random stuff, just set it ot the param value
            self.delgampVar.set(str(self.parent.GetPlotParam('DelGami')))
        return to_reload

    def CheckIfEpsChanged(self):
        to_reload = False

        # The protons first
        try:
            # First check if the injection energy changed
            if np.abs(float(self.EinjPVar.get()) -self.parent.GetPlotParam('GammaIonInjection'))>1E-8:
                # Set the parent value to the var value
                self.parent.SetPlotParam('GammaIonInjection', float(self.EinjPVar.get()), update_plot = False )
                to_reload += self.parent.GetPlotParam('MeasureEpsP')
        except ValueError:
            #if they type in random stuff, just set it to the value
            self.EinjPVar.set(str(self.parent.GetPlotParam('GammaIonInjection')))

        # Now the electrons
        try:
            # First check if the injection energy changed
            if np.abs(float(self.EinjEVar.get()) -self.parent.GetPlotParam('GammaElectronInjection'))>1E-8:
                # Set the parent value to the var value
                self.parent.SetPlotParam('GammaElectronInjection', float(self.EinjEVar.get()), update_plot = False)
                to_reload += self.parent.GetPlotParam('MeasureEpsE')
        except ValueError:
            #if they type in random stuff, just set it to the value
            self.EinjEVar.set(str(self.parent.GetPlotParam('GammaElectronInjection')))


        return to_reload

    def CheckIfPLChanged(self):
        to_reload = False
        VarList = [[self.E1Var, self.E2Var], [self.P1Var, self.P2Var]]
        KeyList = [['PowerLawElectronMin', 'PowerLawElectronMax'], ['PowerLawIonMin', 'PowerLawIonMax']]
        PLList = ['DoPowerLawFitElectron', 'DoPowerLawFitIon']

        for j in range(2):
            try:
                # First check if the left index changed
                if np.abs(float(VarList[j][0].get())- self.parent.GetPlotParam(KeyList[j][0]))>1E-8:
                    # See if the left index is larger than the right index
                    if float(VarList[j][0].get()) > float(VarList[j][1].get()):
                        # it is, so make it larger:
                        VarList[j][1].set(str(float(VarList[j][0].get())*2))
                        #set the parent value to the right var value
                        self.parent.SetPlotParam(KeyList[j][1], float(VarList[j][1].get()), update_plot = False)

                    # Set the parent value to the left var value
                    self.parent.SetPlotParam(KeyList[j][0], float(VarList[j][0].get()), update_plot = False)
                    to_reload += self.parent.GetPlotParam(PLList[j])
            except ValueError:
                #if they type in random stuff, just set it to the value
                VarList[j][0].set(str(self.parent.GetPlotParam(KeyList[j][0])))

            try:
                # Now see if the right index changed
                if np.abs(float(VarList[j][1].get())- self.parent.GetPlotParam(KeyList[j][1]))>1E-8:
                    # See if the left index is smaller than the right index
                    if float(VarList[j][1].get()) < float(VarList[j][0].get()):
                        # it is, so make it smaller:
                        VarList[j][0].set(str(float(VarList[j][1].get())*.5))
                        #set the parent value to the left var value
                        self.parent.SetPlotParam(KeyList[j][0], float(VarList[j][0].get()), update_plot = False)

                    # Set the parent value to the right var value
                    self.parent.SetPlotParam(KeyList[j][1], float(VarList[j][1].get()), update_plot = False)
                    to_reload += self.parent.GetPlotParam(PLList[j])

            except ValueError:
                #if they type in random stuff, just set it to the value
                VarList[j][1].set(str(self.parent.GetPlotParam(KeyList[j][1])))
        return to_reload

    def CheckIfFloatChanged(self, tkVar, paramKey):
        to_reload = False
        try:
            #make sure the user types in a float
            if np.abs(float(tkVar.get() ) - self.parent.GetPlotParam(paramKey)) >1E-8:
                self.parent.SetPlotParam(paramKey, float(tkVar.get()), update_plot = False)
                to_reload = True
            return to_reload

        except ValueError:
            #if they type in random stuff, just set it to the param value
            tkVar.set(str(self.parent.GetPlotParam(paramKey)))
            return to_reload

    def SetTeChanged(self, *args):
        if self.SetTeVar.get()==self.parent.GetPlotParam('SetTe'):
            pass
        else:
            self.parent.SetPlotParam('SetTe', self.SetTeVar.get())

    def SetTpChanged(self, *args):
        if self.SetTpVar.get()==self.parent.GetPlotParam('SetTi'):
            pass
        else:
            self.parent.SetPlotParam('SetTi', self.SetTpVar.get())


    def RelChanged(self, *args):
        if self.RelVar.get()==self.parent.GetPlotParam('PrtlIntegrationRelative'):
            pass
        else:
            self.parent.SetPlotParam('PrtlIntegrationRelative', self.RelVar.get())


    def CheckIfNormsChanged(self):
        to_reload = False
        if np.abs(float(self.eTempNormVar.get()) - self.parent.GetPlotParam('eNormalizer'))>1E-8:
            self.parent.SetPlotParam('eNormalizer', float(self.eTempNormVar.get()), update_plot = False)
            self.sliderTe.set(float(self.eTempNormVar.get()))
            to_reload += self.parent.GetPlotParam('SetTe')
        if np.abs(float(self.iTempNormVar.get()) - self.parent.GetPlotParam('iNormalizer'))>1E-8:
            self.parent.SetPlotParam('iNormalizer', float(self.iTempNormVar.get()), update_plot = False)
            self.sliderTi.set(float(self.iTempNormVar.get()))
            to_reload += self.parent.GetPlotParam('SetTi')
        return to_reload

    def MeasuresCallback(self):
        tkvarIntList = [self.ileft, self.iright, self.eleft, self.eright]
        IntValList = ['IonLeft', 'IonRight', 'ElectronLeft', 'ElectronRight']

        to_reload = False



        for j in range(len(tkvarIntList)):
            to_reload += self.CheckIfFloatChanged(tkvarIntList[j], IntValList[j])

        to_reload += self.CheckIfTeChanged()
        to_reload += self.CheckIfTpChanged()
        to_reload += self.CheckIfNormsChanged()
        to_reload += self.CheckIfPLChanged()
        to_reload += self.CheckIfEpsChanged()

        return to_reload

    def FieldsCallback(self):
        tkvarLimList = [self.xmin, self.xmax, self.ymin, self.ymax]
        plot_param_List = ['x_min', 'x_max', 'y_min', 'y_max']
        tkvarSetList = [self.xLimVar, self.xLimVar, self.yLimVar, self.yLimVar]
        to_reload = False
        for j in range(len(tkvarLimList)):
            try:
            #make sure the user types in a int
                if np.abs(float(tkvarLimList[j].get()) - self.parent.GetPlotParam(plot_param_List[j])) > 1E-8:
                    self.parent.SetPlotParam(plot_param_List[j], float(tkvarLimList[j].get()), update_plot = False)
                    to_reload += True*tkvarSetList[j].get()

            except ValueError:
                #if they type in random stuff, just set it ot the param value
                tkvarLimList[j].set(str(self.parent.GetPlotParam(plot_param_List[j])))
        return to_reload


    def xLimChanged(self, *args):
        if self.xLimVar.get() == self.parent.GetPlotParam('set_xlim'):
            pass
        else:
            self.parent.SetPlotParam('set_xlim', self.xLimVar.get())

    def yLimChanged(self, *args):
        if self.yLimVar.get() == self.parent.GetPlotParam('set_ylim'):
            pass
        else:
            self.parent.SetPlotParam('set_ylim', self.yLimVar.get())

    def OnClosing(self):
        self.parent.settings_window = None
        self.destroy()
