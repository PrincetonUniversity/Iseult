#!/usr/bin/env pythonw
import matplotlib, sys
sys.path.append('../')
import numpy as np
import numpy.ma as ma
import new_cmaps
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as PathEffects
from scipy.special import kn # Modified Bessel function
from scipy.stats import linregress
from scipy.integrate import cumulative_trapezoid
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
                       'PL_legend_loc': 'N/A',
                       'face_color': 'gainsboro'} # lcation of power law legend

    def __init__(self, parent, pos, param_dict):
        self.param_dict = {}
        for key, val in self.plot_param_dict.items():
            self.param_dict[key] = val
        for key, val in param_dict.items():
            self.param_dict[key] = val
        self.pos = pos
        self.parent = parent
        self.chartType = 'SpectraPlot'
        self.figure = self.parent.figure
        self.PowerlawPworked = False
        self.PowerlawEworked = False

        # When the plot is initialized it will look to see if there are other spectral plots.
        #self.spect_num = -1 #it will count itself so we have to correct for an off by one error.
        self.spect_num = self.parent.spect_plot_counter #it will count itself so we have to correct for an off by one error.
        self.parent.spect_plot_counter += 1
        #for i in range(self.parent.MainParamDict['NumOfRows']):
        #    for j in range(self.parent.MainParamDict['NumOfCols']):
        #        if self.parent.SubPlotList[i][j].chartType == 'SpectraPlot':
        #            self.spect_num+=1

        self.dashes = self.parent.dashes_options[min(self.spect_num,len(self.parent.dashes_options)-1)]
        self.RegionFlag = True # Set to true if the region lines need to be written to a phase plot
        # The variables that store the eps_p & eps_e values

    def update_data(self, output):
        self.c_omp = getattr(output, 'c_omp')
        self.istep = getattr(output, 'istep')
#        self.gamma0 = getattr(output, 'gamma0')[0]
        self.xsl = getattr(output, 'xsl')/self.c_omp
        # Load gamma-1 of the spectra
        self.gamma = getattr(output, 'gamma')
        self.massRatio = getattr(output, 'mi')/getattr(output, 'me')
        self.mi = output.mi
        self.me = output.me

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

        # Calculate the momentum := gamma*beta
        #  NOTE: self.gamma is actually the real lorentz factor, gamma, minus 1

        self.momentum=np.sqrt((self.gamma+1)**2-1.)

        if self.GetPlotParam('rest_frame'):
            self.spece = np.copy(getattr(output, 'specerest'))
            self.specp = np.copy(getattr(output, 'specprest'))
        else:
            self.spece = np.copy(getattr(output, 'spece'))
            self.specp = np.copy(getattr(output, 'specp'))


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


    def draw(self):
        self.RegionFlag = True
        # Set the tick color
        tick_color = 'black'
        # Create a gridspec to handle spacing better
        self.gs = gridspec.GridSpecFromSubplotSpec(100,100, subplot_spec = self.parent.gs0[self.pos])#, bottom=0.2,left=0.1,right=0.95, top = 0.95)

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
        self.electron_temp = self.axes.plot(self.momentum, self.momedist,
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
        self.ion_Einj = self.axes.axvline(0, color = self.parent.ion_fit_color,
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
            self.axes.set_axis_bgcolor(self.GetPlotParam('face_color'))
        else:
            self.axes.set_facecolor(self.GetPlotParam('face_color'))

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

                    self.delgame0=self.GetPlotParam('DelGame')*self.mi/self.me
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
        #        # The momentum section is done, now we see what happens if we use energy spectra instead
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
                    self.delgame0=self.GetPlotParam('DelGame')*self.mi/self.me
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
            endSum =  cumulative_trapezoid(energy_dist[::-1],self.gamma[::-1])[::-1]
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
                return np.exp(logF)/psum * self.me/self.mi
            else:
                return np.exp(logF)/endSum[0]

    def GetPlotParam(self, keyname):
        return self.param_dict[keyname]

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
            tmpcon =self.GetPlotParam('DelGami')*self.mi/self.me
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
            self.legT.set_draggable(True, update = 'loc')
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
            self.legDelta.set_draggable(True, update='loc')
            if self.GetPlotParam('PL_legend_loc') != 'N/A':
                tmp_tup = float(self.GetPlotParam('PL_legend_loc').split()[0]),float(self.GetPlotParam('PL_legend_loc').split()[1])
                self.legDelta._set_loc(tmp_tup)
