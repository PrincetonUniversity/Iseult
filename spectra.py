#!/usr/bin/env pythonw
import Tkinter as Tk
import ttk as ttk
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
                       'twoD': 0,
                       'show_ions': 1,
                       'show_electrons': 1,
                       'rest_frame': False,
                       'set_ylim': True,
                       'set_xlim': True,
                       'x_min': 0.05,
                       'x_max': 200,
                       'spatial_x': False,
                       'spatial_y': False,
                       'show_legend': True,
                       'y_min': -6,
                       'y_max': 0}

    # We need the types of all the parameters for the config file
    BoolList = ['twoD', 'set_ylim', 'set_xlim',
                'spatial_x', 'spatial_y',
                'show_ions', 'show_electrons', 'rest_frame','show_legend']
    IntList = ['spectral_type']
    FloatList = ['y_min', 'y_max', 'x_min', 'x_max']
    StrList = []

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
        self.gamma0 = self.FigWrap.LoadKey('gamma0')[0]
        self.xsl = self.FigWrap.LoadKey('xsl')/self.c_omp
        self.gamma = self.FigWrap.LoadKey('gamma')

        self.keyname = 'spectra_data'
        if self.GetPlotParam('rest_frame'):
            self.keyname +='in_rest_frame'

        # A list that will make sure that the data has the same int region
        self.region_args = [self.parent.MainParamDict['PrtlIntegrationRelative'], self.parent.MainParamDict['ElectronLeft'], self.parent.MainParamDict['ElectronRight'], self.parent.MainParamDict['IonLeft'], self.parent.MainParamDict['IonRight']]

        is_loaded = False
        if self.keyname in self.parent.DataDict.keys():

            # make sure it is integrating over the correct region.
            if np.all(self.region_args == self.parent.DataDict[self.keyname][-1]):
                is_loaded = True

                #unpack the values
                self.gamma, self.edist, self.pdist, self.momentum, self.momedist, self.mompdist, region_args = self.parent.DataDict[self.keyname]

        if not is_loaded:
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

            # Select the x-range from which to take the spectra
            e_left_loc = self.parent.MainParamDict['ElectronLeft']
            e_right_loc = self.parent.MainParamDict['ElectronRight']

            if self.parent.MainParamDict['PrtlIntegrationRelative']:
                e_left_loc = self.parent.shock_loc+self.parent.MainParamDict['ElectronLeft']
                e_right_loc = self.parent.shock_loc+self.parent.MainParamDict['ElectronRight']

            eL = self.xsl.searchsorted(e_left_loc)
            eR = self.xsl.searchsorted(e_right_loc, side='right')

            i_left_loc = self.parent.MainParamDict['IonLeft']
            i_right_loc = self.parent.MainParamDict['IonRight']

            if self.parent.MainParamDict['PrtlIntegrationRelative']:
                i_left_loc = self.parent.shock_loc+self.parent.MainParamDict['IonLeft']
                i_right_loc = self.parent.shock_loc+self.parent.MainParamDict['IonRight']

            iL = self.xsl.searchsorted(i_left_loc)
            iR = self.xsl.searchsorted(i_right_loc, side='right')

            if iL >= iR:
                iL = iR
                iR += 1
            if eL >= eR:
                eL = eR
                eR += 1

                # total particles in each linear x bin
            norme = np.copy(self.xsl)
            normp = np.copy(self.xsl)
            for i in range(len(norme)):
                norme[i]=sum(self.spece[:,i])
                normp[i]=sum(self.specp[:,i])

            # energy distribution, f(E)=(dn/dE)/N
            self.fe=np.empty(len(self.gamma))
            self.fp=np.empty(len(self.gamma))


            for k in range(len(self.fe)):
                if sum(norme[eL:eR])*self.dgamma[k] > 1E-100:
                    self.fe[k]=sum(self.spece[k][eL:eR])/(sum(norme[eL:eR])*self.dgamma[k])
                else:
                    print 'RUNTIME WARNING: spectra.py can\'t find the electrons in the integration region, plot may be wrong'
                    self.fe[k] = 1E-100
                if sum(norme[iL:iR])*self.dgamma[k] > 1E-100:
                    self.fp[k]=sum(self.specp[k][iL:iR])/(sum(normp[iL:iR])*self.dgamma[k])
                else:
                    print 'RUNTIME WARNING: spectra.py can\'t find ions in the integration region, plot may be wrong'
                    self.fp[k] = 1E-100

            self.fe[self.fe <= 0] = 1E-100
            self.fp[self.fe <= 0] = 1E-100
            #  NOTE: self.gamma is actually the real lorentz factor, gamma, minus 1 ***
            self.edist=self.gamma*self.fe
            self.pdist=self.gamma*self.fp

            self.momentum=np.sqrt((self.gamma+1)**2-1.)
            self.femom=self.fe/(4*np.pi*self.momentum)/(self.gamma+1)
            self.momedist=self.femom*self.momentum**4
            self.fpmom=self.fp/(4*np.pi*self.momentum)/(self.gamma+1)
            self.mompdist=self.fpmom*self.momentum**4
            self.parent.DataDict[self.keyname] = self.gamma, self.edist, self.pdist, self.momentum, self.momedist, self.mompdist, self.region_args

    def draw(self):

        # Set the tick color
        tick_color = 'black'
        # Create a gridspec to handle spacing better
        self.gs = gridspec.GridSpecFromSubplotSpec(100,100, subplot_spec = self.parent.gs0[self.FigWrap.pos])#, bottom=0.2,left=0.1,right=0.95, top = 0.95)

        self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]])

        self.ion_spect = self.axes.plot(self.momentum, self.mompdist, color = self.parent.ion_color)
        if not self.GetPlotParam('show_ions'):
            self.ion_spect[0].set_visible(False)
        self.ion_temp = self.axes.plot(self.momentum, self.mompdist,
                                       color = self.parent.ion_fit_color,
                                       linestyle = '--', linewidth = 1.5) # a placeholder
        # See if we want to plot the electron temperature
        if not self.parent.MainParamDict['SetTi']:
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
        if not self.parent.MainParamDict['SetTe']:
            self.electron_temp[0].set_visible(False)

        # a placeholder
        self.electron_Einj =  self.axes.axvline(0, color = self.parent.electron_fit_color,
                                            linestyle = '-.', linewidth = 1.5)
        if not self.parent.MainParamDict['MeasureEpsE']:
            self.electron_Einj.set_visible(False)


        # a placeholder
        self.ion_Einj =  self.axes.axvline(0, color = self.parent.ion_fit_color,
                                            linestyle = '-.', linewidth = 1.5)
        if not self.parent.MainParamDict['MeasureEpsP']:
            self.ion_Einj.set_visible(False)

        self.PLE = self.axes.plot(self.momentum, self.momedist,
                                  color = self.parent.electron_fit_color,
                                  linestyle = ':', linewidth = 1.5)
        self.PLE[0].set_visible(False)

        self.axes.set_xscale("log")
        self.axes.set_yscale("log")
        self.axes.set_axis_bgcolor('lightgrey')
        if self.GetPlotParam('set_xlim'):
            self.axes.set_xlim(self.GetPlotParam('x_min'), self.GetPlotParam('x_max'))
        if self.GetPlotParam('set_ylim'):
            self.axes.set_ylim(10**self.GetPlotParam('y_min'), 10**self.GetPlotParam('y_max'))

        self.axes.tick_params(labelsize = self.parent.MainParamDict['NumFontSize'], color=tick_color)

        if self.GetPlotParam('spectral_type') == 0:
            self.axes.set_xlabel(r'$\gamma\beta$', labelpad = self.parent.MainParamDict['xLabelPad'], color = 'black')
            self.axes.set_ylabel(r'$p^4f(p)$', labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black')
        else:
            self.axes.set_xlabel(r'$\gamma-1$', labelpad = -2, color = 'black')
            self.axes.set_ylabel(r'$E(dn/dE)/n$', labelpad = 0, color = 'black')

        self.refresh()

    def refresh(self):
        ############
        #
        # First we'll calculate eps_e or eps_e
        #
        ############

        if self.parent.MainParamDict['MeasureEpsP']:
            eps_p = self.measure_eps(self.pdist, self.parent.MainParamDict['GammaIonInjection'], 'protron')
            self.parent.eps_pVar.set('%.6f' % eps_p)
        else:
            self.parent.eps_pVar.set('N/A')

        if self.parent.MainParamDict['MeasureEpsE']:
            eps_e = self.measure_eps(self.edist, self.parent.MainParamDict['GammaElectronInjection'], 'electron')
            self.parent.eps_eVar.set('%.6f' % eps_e)
        else:
            self.parent.eps_eVar.set('N/A')
        if self.GetPlotParam('spectral_type') == 0: #Show the momentum dist
            if self.GetPlotParam('show_ions'):
                self.ion_spect[0].set_data(self.momentum, self.mompdist)

                if self.parent.MainParamDict['MeasureEpsP']:
                    self.ion_Einj.set_visible(True)
                    self.ion_Einj.set_xdata([np.sqrt((self.parent.MainParamDict['GammaIonInjection']+1)**2-1.),np.sqrt((self.parent.MainParamDict['GammaIonInjection']+1)**2-1.)])
                else:
                    self.ion_Einj.set_visible(False)
                if self.parent.MainParamDict['SetTi']:
                    self.ion_temp[0].set_visible(True)
                    if self.parent.MainParamDict['DelGami'] >= 0.013:
                        aconst = 1/(self.parent.MainParamDict['DelGami']*np.exp(1.0/self.parent.MainParamDict['DelGami'])*kn(2, 1.0/self.parent.MainParamDict['DelGami']))
                    else:
                        aconst = np.sqrt(2/np.pi)/self.parent.MainParamDict['DelGami']**1.5
                        aconst -= 15.0/(4.*np.sqrt(self.parent.MainParamDict['DelGami'])*np.sqrt(2*np.pi))
                        aconst += (345*np.sqrt(self.parent.MainParamDict['DelGami']))/(64.*np.sqrt(2*np.pi))
                        aconst -= (3285*self.parent.MainParamDict['DelGami']**1.5)/(512.*np.sqrt(2*np.pi))
                        aconst += (95355*self.parent.MainParamDict['DelGami']**2.5)/(16384.*np.sqrt(2*np.pi))
                        aconst -= (232065*self.parent.MainParamDict['DelGami']**3.5)/(131072.*np.sqrt(2*np.pi))

                    fpmommax = self.momentum**4*aconst*(self.gamma+1.0)*np.sqrt((self.gamma+1.0)**2-1)
                    fpmommax *= np.exp(-self.gamma/self.parent.MainParamDict['DelGami'])/(4*np.pi*self.momentum)/(self.gamma+1.0)
                    self.ion_temp[0].set_data(self.momentum, fpmommax)

                else:
                    self.ion_temp[0].set_visible(False)

                self.PowerlawPworked = False
                self.PLP[0].set_visible(False)
                if self.parent.MainParamDict['DoPowerLawFitIon']:

                    # Find the part of the spectrum chosen by the fitting
                    #first convert to
                    mompleft = np.sqrt((self.parent.MainParamDict['PowerLawIonMin']+1)**2-1)
                    mompright = np.sqrt((self.parent.MainParamDict['PowerLawIonMax']+1.)**2-1.)
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
                if self.parent.MainParamDict['MeasureEpsE']:
                    self.electron_Einj.set_visible(True)
                    self.electron_Einj.set_xdata([np.sqrt((self.parent.MainParamDict['GammaElectronInjection']+1)**2-1.),np.sqrt((self.parent.MainParamDict['GammaElectronInjection']+1)**2-1.)])
                else:
                    self.electron_Einj.set_visible(False)

                if self.parent.MainParamDict['SetTe']:
                    self.electron_temp[0].set_visible(True)

                    self.delgame0=self.parent.MainParamDict['DelGame']*self.FigWrap.LoadKey('mi')[0]/self.FigWrap.LoadKey('me')[0]
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
                    self.electron_temp[0].set_data(self.momentum, femommax)


                else:
                    self.electron_temp[0].set_visible(False)

                # Now the power-law
                self.PowerlawEworked = False
                self.PLE[0].set_visible(False)
                if self.parent.MainParamDict['DoPowerLawFitElectron']:


                    # Find the part of the spectrum chosen by the fitting
                    #first convert to
                    momeleft = np.sqrt((self.parent.MainParamDict['PowerLawElectronMin']+1)**2-1)
                    momeright = np.sqrt((self.parent.MainParamDict['PowerLawElectronMax']+1.)**2-1.)
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


        if self.GetPlotParam('spectral_type') == 1: #Show the energy dist
            if self.GetPlotParam('show_electrons'):
                self.electron_spect[0].set_data(self.gamma, self.edist)
                if self.parent.MainParamDict['MeasureEpsE']:
                    self.electron_Einj.set_visible(True)
                    self.electron_Einj.set_xdata([self.parent.MainParamDict['GammaElectronInjection'],self.parent.MainParamDict['GammaElectronInjection']])
                else:
                    self.electron_Einj.set_visible(False)


                # The temperature
                if self.parent.MainParamDict['SetTe']:
                    self.electron_temp[0].set_visible(True)
                    self.delgame0=self.parent.MainParamDict['DelGame']*self.FigWrap.LoadKey('mi')[0]/self.FigWrap.LoadKey('me')[0]
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
                    self.electron_temp[0].set_data(self.gamma, femax)


                # the power-law
                self.PowerlawEworked = False
                self.PLE[0].set_visible(False)

                if self.parent.MainParamDict['DoPowerLawFitElectron']:
                    # Find the part of the spectrum chosen by the fitting
                    #first convert to
                    ieLeft = self.gamma.searchsorted(self.parent.MainParamDict['PowerLawElectronMin'])
                    ieRight = self.gamma.searchsorted(self.parent.MainParamDict['PowerLawElectronMax'], side='Right')
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
                if self.parent.MainParamDict['MeasureEpsP']:
                    self.ion_Einj.set_visible(True)
                    self.ion_Einj.set_xdata([self.parent.MainParamDict['GammaIonInjection'],self.parent.MainParamDict['GammaIonInjection']])
                else:
                    self.ion_Einj.set_visible(False)

                if self.parent.MainParamDict['SetTi']:
                    self.ion_temp[0].set_visible(True)
                    if self.parent.MainParamDict['DelGami'] >= 0.013:
                        aconst = 1/(self.parent.MainParamDict['DelGami']*np.exp(1.0/self.parent.MainParamDict['DelGami'])*kn(2, 1.0/self.parent.MainParamDict['DelGami']))
                    else:
                        aconst = np.sqrt(2/np.pi)/self.parent.MainParamDict['DelGami']**1.5
                        aconst -= 15.0/(4.*np.sqrt(self.parent.MainParamDict['DelGami'])*np.sqrt(2*np.pi))
                        aconst += (345*np.sqrt(self.parent.MainParamDict['DelGami']))/(64.*np.sqrt(2*np.pi))
                        aconst -= (3285*self.parent.MainParamDict['DelGami']**1.5)/(512.*np.sqrt(2*np.pi))
                        aconst += (95355*self.parent.MainParamDict['DelGami']**2.5)/(16384.*np.sqrt(2*np.pi))
                        aconst -= (232065*self.parent.MainParamDict['DelGami']**3.5)/(131072.*np.sqrt(2*np.pi))

                    fpmax = aconst*self.gamma*(self.gamma+1.0)*np.sqrt((self.gamma+1.0)**2-1)
                    fpmax *= np.exp(-self.gamma/self.parent.MainParamDict['DelGami'])
                    self.ion_temp[0].set_data(self.gamma, fpmax)

                self.PowerlawPworked = False
                if not self.parent.MainParamDict['DoPowerLawFitIon']:
                    self.PLP[0].set_visible(False)
                else:
                    # Find the part of the spectrum chosen by the fitting
                    iepLeft = self.gamma.searchsorted(self.parent.MainParamDict['PowerLawIonMin'])
                    iepRight = self.gamma.searchsorted(self.parent.MainParamDict['PowerLawIonMax'], side='Right')

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
            because it goes when energy_dist is zero. I'm leaving it in, commented out in case it is
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

    def SetPlotParam(self, keyname, value, update_plot = True):
        self.FigWrap.SetPlotParam(keyname, value, update_plot = update_plot)

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


        if self.GetPlotParam('show_electrons') and self.parent.MainParamDict['SetTe']:
            Tlegend_handles.append(self.electron_temp[0])
            tmpstr = '%.3f' % self.delgame0
            Tlegend_labels.append(r'$T_e\ = $' +  ' ' + tmpstr + ' ' + r'$m_e c^2$')
        if self.GetPlotParam('show_ions') and self.parent.MainParamDict['SetTi']:
            Tlegend_handles.append(self.ion_temp[0])
            tmpcon =self.parent.MainParamDict['DelGami']*self.FigWrap.LoadKey('mi')[0]/self.FigWrap.LoadKey('me')[0]
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
        if len(Tlegend_handles)> 0 and len(legend_handles) == 0:
            self.legT = self.axes.legend(Tlegend_handles, Tlegend_labels, framealpha = .05, fontsize = 11, loc = 'upper left')
            self.legT.get_frame().set_facecolor('k')
            self.legT.get_frame().set_linewidth(0.0)

        elif len(Tlegend_handles) > 0:
            self.legT = self.axes.legend(Tlegend_handles, Tlegend_labels, framealpha = .05, fontsize = 11, loc = 'upper left')
            self.legT.get_frame().set_facecolor('k')
            self.legT.get_frame().set_linewidth(0.0)
            self.axes.add_artist(self.legT)

        if len(legend_handles)> 0:
            self.legDelta = self.axes.legend(legend_handles, legend_labels,
            framealpha = .05, fontsize = 11, loc = 'upper right')
            self.legDelta.get_frame().set_facecolor('k')
            self.legDelta.get_frame().set_linewidth(0.0)

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
        self.parent = parent
        frm = ttk.Frame(self)
        frm.pack(fill=Tk.BOTH, expand=True)
        self.protocol('WM_DELETE_WINDOW', self.OnClosing)
        #Create some sizers

        self.bind('<Return>', self.TxtEnter)

        # Create the OptionMenu to chooses the Chart Type:
        self.ctypevar = Tk.StringVar(self)
        self.ctypevar.set(self.parent.chartType) # default value
        self.ctypevar.trace('w', self.ctypeChanged)

        ttk.Label(frm, text="Choose Chart Type:").grid(row=0, column = 0)
        cmapChooser = apply(ttk.OptionMenu, (frm, self.ctypevar, self.parent.chartType) + tuple(self.parent.ChartTypes))
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
        # show in rest frame
        self.RestVar = Tk.IntVar()
        self.RestVar.set(self.parent.GetPlotParam('rest_frame'))
        cb = ttk.Checkbutton(frm, text = "Show in rest frame",
                        variable = self.RestVar,
                        command = lambda:
                        self.parent.SetPlotParam('rest_frame', self.RestVar.get()))
        cb.grid(row = 7, column = 0, sticky = Tk.W)

        self.xLimVar = Tk.IntVar()
        self.xLimVar.set(self.parent.GetPlotParam('set_xlim'))
        self.xLimVar.trace('w', self.xLimChanged)



        self.xmin = Tk.StringVar()
        self.xmin.set(str(self.parent.GetPlotParam('x_min')))
        self.xmax = Tk.StringVar()
        self.xmax.set(str(self.parent.GetPlotParam('x_max')))


#        cb = ttk.Checkbutton(frm, text ='Set xlim',
#                        variable = self.xLimVar)
#        cb.grid(row = 4, column =3,sticky = Tk.W)
        ttk.Label(frm, text = 'Set xlim').grid(row = 4, column = 3, sticky = Tk.W)
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


#        cb = ttk.Checkbutton(frm, text ='Set log(y) lim',
#                        variable = self.yLimVar)
#        cb.grid(row = 5,  column =3, sticky = Tk.W)
        ttk.Label(frm, text = 'Set log(y) lim').grid(row = 5, column =3, sticky = Tk.W)
        self.eLEnter = ttk.Entry(frm, textvariable=self.ymin, width=7)
        self.eLEnter.grid(row = 5, column =4)
        self.eREnter = ttk.Entry(frm, textvariable=self.ymax, width=7)
        self.eREnter.grid(row = 5, column =5)

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
            self.parent.SetPlotParam('show_ions', self.IonVar.get())

    def eVarHandler(self, *args):
        if self.eVar.get() == self.parent.GetPlotParam('show_electrons'):
            pass
        else:
            self.parent.electron_spect[0].set_visible(self.eVar.get())
            self.parent.SetPlotParam('show_electrons', self.eVar.get())

    def RadioSpect(self):
        if self.SpectTypeVar.get() == self.parent.GetPlotParam('spectral_type'):
            pass
        else:
            if self.SpectTypeVar.get() == 1:
                self.parent.axes.set_xlabel(r'$\gamma-1$')
                self.parent.axes.set_ylabel(r'$E(dn/dE)/n$')
                self.parent.SetPlotParam('x_min', 0.0005, update_plot = False)
                self.parent.SetPlotParam('x_max', 100, update_plot = False)

            else:
                self.parent.axes.set_xlabel(r'$\gamma\beta$')
                self.parent.axes.set_ylabel(r'$p^4f(p)$')
                self.parent.SetPlotParam('x_min', 0.05, update_plot = False)
                self.parent.SetPlotParam('x_max', 200, update_plot = False)

            self.xmin.set(str(self.parent.GetPlotParam('x_min')))
            self.xmax.set(str(self.parent.GetPlotParam('x_max')))
            self.parent.SetPlotParam('spectral_type', self.SpectTypeVar.get())


    def TxtEnter(self, e):
        self.FieldsCallback()

    def FieldsCallback(self):
        tkvarLimList = [self.xmin, self.xmax, self.ymin, self.ymax]
        plot_param_List = ['x_min', 'x_max', 'y_min', 'y_max']
        tkvarSetList = [self.xLimVar, self.xLimVar, self.yLimVar, self.yLimVar]
        to_reload = False
        for j in range(len(tkvarLimList)):
            try:
            #make sure the user types in a int
                if np.abs(float(tkvarLimList[j].get()) - self.parent.GetPlotParam(plot_param_List[j])) > 1E-4:
                    self.parent.SetPlotParam(plot_param_List[j], float(tkvarLimList[j].get()), update_plot = False)
                    to_reload += True*tkvarSetList[j].get()

            except ValueError:
                #if they type in random stuff, just set it ot the param value
                tkvarLimList[j].set(str(self.parent.GetPlotParam(plot_param_List[j])))
        if to_reload:
            self.parent.SetPlotParam('x_min', self.parent.GetPlotParam('x_min'))
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
