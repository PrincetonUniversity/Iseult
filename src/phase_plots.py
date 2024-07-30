#!/usr/bin/env pythonw
import tkinter as Tk
from tkinter import ttk
import matplotlib
import numpy as np
import numpy.ma as ma
import new_cmaps
from new_cnorms import PowerNormWithNeg
from Numba2DHist import Fast2DHist, Fast2DWeightedHist, vecLog10Norm
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as PathEffects

class PhasePanel:
    # A dictionary of all of the parameters for this plot with the default parameters

    plot_param_dict = {'twoD' : 1,
                       'mom_dim': 0,
                       'masked': 1,
                       'cnorm_type': 'Log', #Colormap normalization. Opts are Log or Linear
                       'prtl_type': 0,
                       'cpow_num': 0.6,
                       'show_cbar': True,
                       'weighted': False,
                       'show_shock': False,
                       'show_int_region': True,
                       'xbins' : 200,
                       'pbins' : 200,
                       'v_min': -2.0,
                       'v_max' : 0,
                       'set_v_min': False,
                       'set_v_max': False,
                       'p_min': -2.0,
                       'p_max' : 2,
                       'set_E_min' : False,
                       'E_min': 1.0,
                       'set_E_max': False,
                       'E_max': 200.0,
                       'set_p_min': False,
                       'set_p_max': False,
                       'spatial_x': True,
                       'spatial_y': False,
                       'symmetric': False,
                       'interpolation': 'nearest',
                       'face_color': 'gainsboro'}


    prtl_opts = ['proton_p', 'electron_p']
    direction_opts = ['x-x', 'y-x', 'z-x']
    # Old labels:
    #ylabel_list =[
    #              [[r'$P_{px}\ [m_i c]$', r'$P_{py}\ [m_i c]$',r'$P_{pz}\ [m_i c]$'],
    #              [r'$P_{ex}\ [m_e c]$', r'$P_{ey}\ [m_e c]$',r'$P_{ez}\ [m_e c]$']],
    #              [[r'$P\prime_{px}\ [m_i c]$', r'$P\prime_{py}\ [m_i c]$',r'$P\prime_{pz}\ [m_i c]$'],
    #              [r'$P\prime_{ex}\ [m_e c]$', r'$P\prime_{ey}\ [m_e c]$',r'$P\prime_{ez}\ [m_e c]$']]
    #             ]
    ylabel_list =[
                 [[r'$\gamma_i\beta_{x,i}$',r'$\gamma_i\beta_{y,i}$',r'$\gamma_i\beta_{z,i}$'],
                  [r'$\gamma_e\beta_{x,e}$',r'$\gamma_e\beta_{y,e}$',r'$\gamma_e\beta_{z,e}$']],
                 [[r'$\gamma\prime_i\beta\prime_{x,i}$',r'$\gamma\prime_i\beta\prime_{y,i}$',r'$\gamma\prime_i\beta\prime_{z,i}$'],
                  [r'$\gamma\prime_e\beta\prime_{x,e}$',r'$\gamma\prime_e\beta\prime_{y,e}$',r'$\gamma\prime_e\beta\prime_{z,e}$']]
                 ]

    gradient =  np.linspace(0, 1, 256)# A way to make the colorbar display better
    gradient = np.vstack((gradient, gradient))

    def __init__(self, parent, figwrapper):

        self.settings_window = None
        self.FigWrap = figwrapper
        self.parent = parent
        self.ChartTypes = self.FigWrap.PlotTypeDict.keys()
        self.chartType = self.FigWrap.chartType
        self.figure = self.FigWrap.figure
        self.InterpolationMethods = ['none','nearest', 'bilinear', 'bicubic', 'spline16',
            'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric',
            'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']

        # A variable that controls whether the energy integration region
        # is shown
        self.IntRegVar = Tk.IntVar()
        self.IntRegVar.set(self.GetPlotParam('show_int_region'))
        self.IntRegVar.trace('w', self.IntVarHandler)
        # Figure out the energy color the intergration region
        if self.GetPlotParam('prtl_type') == 1: #electons
            self.energy_color = self.parent.electron_color
        else:
            self.energy_color = self.parent.ion_color
        # A list that will hold any lines for the integration region


    def IntVarHandler(self, *args):
        # This should only be called by the user-interaction when all the plots already exist...
        # so we can take some shortcut  s and assume a lot of things are already created.
        self.SetPlotParam('show_int_region', self.IntRegVar.get(), update_plot = False)
        if self.IntRegVar.get() == True:
            # We need to show the integration region.

            # Look for all the spectra plots and plot the lines.
            for i in range(self.parent.MainParamDict['NumOfRows']):
                for j in range(self.parent.MainParamDict['NumOfCols']):
                    if self.parent.SubPlotList[i][j].chartType == 'SpectraPlot':
                        k = min(self.parent.SubPlotList[i][j].graph.spect_num, len(self.parent.dashes_options)-1)
                        # figure out if we are as ion phase diagram or an electron one
                        if self.GetPlotParam('prtl_type') == 0:
                            # Append the left line to the list
                            self.IntRegionLines.append(self.axes.axvline(
                            max(self.parent.SubPlotList[i][j].graph.i_left_loc, self.xmin+1),
                            linewidth = 1.5, linestyle = '-', color = self.energy_color))
                            # Choose the right dashes pattern
                            self.IntRegionLines[-1].set_dashes(self.parent.dashes_options[k])
                            # Append the left line to the list
                            self.IntRegionLines.append(self.axes.axvline(
                            min(self.parent.SubPlotList[i][j].graph.i_right_loc, self.xmax-1),
                            linewidth = 1.5, linestyle = '-', color = self.energy_color))
                            # Choose the right dashes pattern
                            self.IntRegionLines[-1].set_dashes(self.parent.dashes_options[k])
                        else:
                            # Append the left line to the list
                            self.IntRegionLines.append(self.axes.axvline(
                            max(self.parent.SubPlotList[i][j].graph.e_left_loc, self.xmin+1),
                            linewidth = 1.5, linestyle = '-', color = self.energy_color))
                            # Choose the right dashes pattern
                            self.IntRegionLines[-1].set_dashes(self.parent.dashes_options[k])
                            # Append the left line to the list
                            self.IntRegionLines.append(self.axes.axvline(
                            min(self.parent.SubPlotList[i][j].graph.e_right_loc, self.xmax-1),
                            linewidth = 1.5, linestyle = '-', color = self.energy_color))
                            # Choose the right dashes pattern
                            self.IntRegionLines[-1].set_dashes(self.parent.dashes_options[k])

        # CLOSES IF. NOW IF WE TURN OFF THE INTEGRATION REGIONS, we have to delete all the lines.
        else:
            for i in xrange(len(self.IntRegionLines)):
                self.IntRegionLines.pop(0).remove()
        # Update the canvas
        self.parent.canvas.draw()
        self.parent.canvas.get_tk_widget().update_idletasks()

    def ChangePlotType(self, str_arg):
        self.FigWrap.ChangeGraph(str_arg)

    def norm(self, vmin=None, vmax=None):
        if self.GetPlotParam('cnorm_type') == "Log":
            return  mcolors.LogNorm(vmin, vmax)

        else:
            return mcolors.Normalize(vmin, vmax)


    def set_plot_keys(self):
        '''A helper function that will insure that each hdf5 file will only be
        opened once per time step'''
        self.arrs_needed = ['c_omp', 'bx', 'istep', 'me', 'mi']
        # First see if we will need to know the energy of the particle
        # (requied for lorentz boosts and setting e_min and e_max)
        Need_Energy = self.parent.MainParamDict['DoLorentzBoost'] and np.abs(self.parent.MainParamDict['GammaBoost'])>1E-8
        Need_Energy = Need_Energy or self.GetPlotParam('set_E_min')
        Need_Energy = Need_Energy or self.GetPlotParam('set_E_max')

        if self.GetPlotParam('prtl_type') == 0:
            self.arrs_needed.append('xi')
            if self.GetPlotParam('weighted'):
                self.arrs_needed.append('chi')
            if Need_Energy:
                self.arrs_needed.append('ui')
                self.arrs_needed.append('vi')
                self.arrs_needed.append('wi')
            elif self.GetPlotParam('mom_dim') == 0:
                self.arrs_needed.append('ui')
            elif self.GetPlotParam('mom_dim') == 1:
                self.arrs_needed.append('vi')
            elif self.GetPlotParam('mom_dim') == 2:
                self.arrs_needed.append('wi')

        if self.GetPlotParam('prtl_type') == 1:
            self.arrs_needed.append('xe')
            if self.GetPlotParam('weighted'):
                self.arrs_needed.append('che')
            if Need_Energy:
                self.arrs_needed.append('ue')
                self.arrs_needed.append('ve')
                self.arrs_needed.append('we')
            elif self.GetPlotParam('mom_dim') == 0:
                self.arrs_needed.append('ue')
            elif self.GetPlotParam('mom_dim') == 1:
                self.arrs_needed.append('ve')
            elif self.GetPlotParam('mom_dim') == 2:
                self.arrs_needed.append('we')
        return self.arrs_needed

    def LoadData(self):
        ''' A helper function that checks if the histogram has
        already been calculated and if it hasn't, it calculates
        it then stores it.'''
        self.key_name = str(self.GetPlotParam('pbins')) + 'x' + str(self.GetPlotParam('xbins'))

        if self.GetPlotParam('masked'):
            self.key_name += 'masked_'

        if self.GetPlotParam('weighted'):
            self.key_name += 'weighted_'

        if self.GetPlotParam('set_E_min'):
            self.key_name += 'Emin_'+str(self.GetPlotParam('E_min')) + '_'

        if self.GetPlotParam('set_E_max'):
            self.key_name += 'Emax_'+str(self.GetPlotParam('E_max')) + '_'

        if self.parent.MainParamDict['DoLorentzBoost'] and np.abs(self.parent.MainParamDict['GammaBoost'])>1E-8:
            self.key_name += 'boosted_'+ str(self.parent.MainParamDict['GammaBoost'])+'_'

        self.key_name += self.prtl_opts[self.GetPlotParam('prtl_type')]
        self.key_name += self.direction_opts[self.GetPlotParam('mom_dim')]
        self.key_name += str(int(self.parent.MainParamDict['PrtlStride']))

        if self.key_name in self.parent.DataDict.keys():
            self.hist2d = self.parent.DataDict[self.key_name]

        elif self.parent.MainParamDict['DoLorentzBoost'] and np.abs(self.parent.MainParamDict['GammaBoost'])>1E-8:
            # Gotta boost it
            self.c_omp = self.FigWrap.LoadKey('c_omp')
            self.istep = self.FigWrap.LoadKey('istep')
            self.weights = None
            self.x_values = None
            self.y_values = None

            # x_min & x_max before boostin'
            self.xmin = 0
            self.xmax = self.FigWrap.LoadKey('bx').shape[2]/self.c_omp*self.istep
            self.xmax = self.xmax if (self.xmax != self.xmin) else self.xmin + 1

            # First calculate beta and gamma
            if self.parent.MainParamDict['GammaBoost'] >=1:
                self.GammaBoost = self.parent.MainParamDict['GammaBoost']
                self.betaBoost = np.sqrt(1-1/self.parent.MainParamDict['GammaBoost']**2)
            elif self.parent.MainParamDict['GammaBoost'] >-1:
                self.betaBoost = self.parent.MainParamDict['GammaBoost']
                self.GammaBoost = np.sqrt(1-self.betaBoost**2)**(-1)

            else:
                self.GammaBoost = -self.parent.MainParamDict['GammaBoost']
                self.betaBoost = -np.sqrt(1-1/self.parent.MainParamDict['GammaBoost']**2)



            # Now load the data. We require all 3 dimensions to determine
            # the velociy and LF in the boosted frame.
            if self.GetPlotParam('prtl_type') == 0:
                # first load everything downstream frame
                self.x_values = self.FigWrap.LoadKey('xi')/self.c_omp

                u = self.FigWrap.LoadKey('ui')
                v = self.FigWrap.LoadKey('vi')
                w = self.FigWrap.LoadKey('wi')
                if self.GetPlotParam('weighted'):
                    self.weights = self.FigWrap.LoadKey('chi')

            if self.GetPlotParam('prtl_type') == 1: #electons
                self.x_values = self.FigWrap.LoadKey('xe')/self.c_omp
                u = self.FigWrap.LoadKey('ue')
                v = self.FigWrap.LoadKey('ve')
                w = self.FigWrap.LoadKey('we')

                if self.GetPlotParam('weighted'):
                    self.weights = self.FigWrap.LoadKey('che')


            # Now calculate gamma of the particles in downstream restframe
            gamma_ds = np.sqrt(u**2+v**2+w**2+1)



            # calculate the velocities from the momenta
            vx = u/gamma_ds
            vy = v/gamma_ds
            vz = w/gamma_ds

            # Now calulate the velocities in the boosted frames
            tmp_helper = 1-vx*self.betaBoost
            vx_prime = (vx-self.betaBoost)/tmp_helper
            vy_prime = vy/self.GammaBoost/tmp_helper
            vz_prime = vz/self.GammaBoost/tmp_helper

            # Now calculate the LF in the boosted frames using rapidity
            # Initial rapidity
            rap_prtl = np.arccosh(gamma_ds)
            rap_boost = np.arccosh(self.GammaBoost)

            #v_tot_sq = vx_prime**2 + vy_prime**2 + vz_prime**2
            #gamma_old_way = 1/np.sqrt(1-v_tot_sq)

            gamma_prime = gamma_ds*self.GammaBoost-np.sign(u)*np.sign(self.betaBoost)*np.sinh(rap_prtl)*np.sinh(rap_boost)/np.sqrt(1+(v/u)**2+(w/u)**2)

            if self.GetPlotParam('mom_dim') == 0:
                self.y_values  = vx_prime*gamma_prime
            if self.GetPlotParam('mom_dim') == 1:
                self.y_values  = vy_prime*gamma_prime
            if self.GetPlotParam('mom_dim') == 2:
                self.y_values  = vz_prime*gamma_prime

            # Some of the values are becoming NaN.
            # ignore them, but I don't think this should be happening anymore....
            nan_ind = np.isnan(self.y_values)


            self.pmin = 0.0 if len(self.y_values) == 0 else min(self.y_values)
            self.pmax = 0.0 if len(self.y_values) == 0 else max(self.y_values)
            self.pmax = self.pmax if (self.pmax != self.pmin) else self.pmin + 1


            if self.GetPlotParam('set_E_min') or self.GetPlotParam('set_E_max'):
                # We need to calculate the total energy in units m_e c^2
                if self.GetPlotParam('prtl_type')==0:
                    energy = gamma_ds*self.FigWrap.LoadKey('mi')/self.FigWrap.LoadKey('me')
                else:
                    energy = np.copy(gamma_ds)

                # Now find the particles that fall in our range
                if self.GetPlotParam('set_E_min'):
                    inRange = energy >= self.FigWrap.GetPlotParam('E_min')
                    if self.GetPlotParam('set_E_max'):
                        inRange *= energy <= self.GetPlotParam('E_max')
                elif self.GetPlotParam('set_E_max'):
                    inRange = energy <= self.GetPlotParam('E_max')
                inRange *= ~nan_ind
                if self.GetPlotParam('weighted'):
                    self.hist2d = Fast2DWeightedHist(self.y_values[inRange], self.x_values[inRange], self.weights[inRange], self.pmin,self.pmax, self.GetPlotParam('pbins'), self.xmin,self.xmax, self.GetPlotParam('xbins')), [self.pmin, self.pmax], [self.xmin, self.xmax]

                else:
                    self.hist2d = Fast2DHist(self.y_values[inRange], self.x_values[inRange], self.pmin,self.pmax, self.GetPlotParam('pbins'), self.xmin,self.xmax, self.GetPlotParam('xbins')), [self.pmin, self.pmax], [self.xmin, self.xmax]

            else:
                if self.GetPlotParam('weighted'):
                    self.hist2d = Fast2DWeightedHist(self.y_values, self.x_values, self.weights, self.pmin,self.pmax, self.GetPlotParam('pbins'), self.xmin,self.xmax, self.GetPlotParam('xbins')), [self.pmin, self.pmax], [self.xmin, self.xmax]
                else:
                    self.hist2d = Fast2DHist(self.y_values, self.x_values, self.pmin,self.pmax, self.GetPlotParam('pbins'), self.xmin,self.xmax, self.GetPlotParam('xbins')), [self.pmin, self.pmax], [self.xmin, self.xmax]
            try:
                if self.GetPlotParam('masked'):
                    zval = ma.masked_array(self.hist2d[0])
                    zval[zval == 0] = ma.masked
                    zval *= float(zval.max())**(-1)
                    tmplist = [zval[~zval.mask].min(), zval.max()]
                else:
                    zval = np.copy(self.hist2d[0])
                    zval[zval==0] = 0.5
                    zval *= float(zval.max())**(-1)
                    tmplist = [zval.min(), zval.max()]
            except ValueError:
                tmplist=[0.1,1]
            self.hist2d = zval, self.hist2d[1], self.hist2d[2], tmplist

            self.parent.DataDict[self.key_name] = self.hist2d


        else:
            # Generate the X-axis values
            self.c_omp = self.FigWrap.LoadKey('c_omp')
            self.istep = self.FigWrap.LoadKey('istep')
            self.weights = None
            self.x_values = None
            self.y_values = None

            # Choose the particle type and px, py, or pz
            if self.GetPlotParam('prtl_type') == 0: #protons
                self.x_values = self.FigWrap.LoadKey('xi')/self.c_omp
                if self.GetPlotParam('weighted'):
                    self.weights = self.FigWrap.LoadKey('chi')
                if self.GetPlotParam('mom_dim') == 0:
                    self.y_values = self.FigWrap.LoadKey('ui')
                if self.GetPlotParam('mom_dim') == 1:
                    self.y_values = self.FigWrap.LoadKey('vi')
                if self.GetPlotParam('mom_dim') == 2:
                    self.y_values = self.FigWrap.LoadKey('wi')

            if self.GetPlotParam('prtl_type') == 1: #electons
                self.energy_color = self.parent.electron_color
                self.x_values = self.FigWrap.LoadKey('xe')/self.c_omp
                if self.GetPlotParam('weighted'):
                    self.weights = self.FigWrap.LoadKey('che')
                if self.GetPlotParam('mom_dim') == 0:
                    self.y_values = self.FigWrap.LoadKey('ue')
                if self.GetPlotParam('mom_dim') == 1:
                    self.y_values = self.FigWrap.LoadKey('ve')
                if self.GetPlotParam('mom_dim') == 2:
                    self.y_values = self.FigWrap.LoadKey('we')

            self.pmin = 0.0 if len(self.y_values) == 0 else min(self.y_values)
            self.pmax = 0.0 if len(self.y_values) == 0 else max(self.y_values)
            self.pmax = self.pmax if (self.pmax != self.pmin) else self.pmin + 1

            self.xmin = 0
            self.xmax = self.FigWrap.LoadKey('bx').shape[2]/self.c_omp*self.istep
            self.xmax = self.xmax if (self.xmax != self.xmin) else self.xmin + 1

            if self.GetPlotParam('set_E_min') or self.GetPlotParam('set_E_max'):
                # We need to calculate the total energy of each particle in
                # units m_e c^2

                # First load the data. We require all 3 dimensions of momentum
                # to determine the energy in the downstream frame
                if self.GetPlotParam('prtl_type') == 0:
                    u = self.FigWrap.LoadKey('ui')
                    v = self.FigWrap.LoadKey('vi')
                    w = self.FigWrap.LoadKey('wi')

                if self.GetPlotParam('prtl_type') == 1: #electons
                    self.x_values = self.FigWrap.LoadKey('xe')/self.c_omp
                    u = self.FigWrap.LoadKey('ue')
                    v = self.FigWrap.LoadKey('ve')
                    w = self.FigWrap.LoadKey('we')

                # Now calculate LF of the particles in downstream restframe
                energy = np.sqrt(u**2+v**2+w**2+1)
                # If they are electrons this already the energy in units m_e c^2.
                # Otherwise...
                if self.GetPlotParam('prtl_type')==0:
                    energy *= self.FigWrap.LoadKey('mi')/self.FigWrap.LoadKey('me')

                # Now find the particles that fall in our range
                if self.GetPlotParam('set_E_min'):
                    inRange = energy >= self.FigWrap.GetPlotParam('E_min')
                    if self.GetPlotParam('set_E_max'):
                        inRange *= energy <= self.GetPlotParam('E_max')
                elif self.GetPlotParam('set_E_max'):
                    inRange = energy <= self.GetPlotParam('E_max')
                if self.GetPlotParam('weighted'):
                    self.hist2d = Fast2DWeightedHist(self.y_values[inRange], self.x_values[inRange], self.weights[inRange], self.pmin,self.pmax, self.GetPlotParam('pbins'), self.xmin,self.xmax, self.GetPlotParam('xbins')), [self.pmin, self.pmax], [self.xmin, self.xmax]
                else:
                    self.hist2d = Fast2DHist(self.y_values[inRange], self.x_values[inRange], self.pmin,self.pmax, self.GetPlotParam('pbins'), self.xmin,self.xmax, self.GetPlotParam('xbins')), [self.pmin, self.pmax], [self.xmin, self.xmax]
            else:
                if self.GetPlotParam('weighted'):
                    self.hist2d = Fast2DWeightedHist(self.y_values, self.x_values, self.weights, self.pmin,self.pmax, self.GetPlotParam('pbins'), self.xmin,self.xmax, self.GetPlotParam('xbins')), [self.pmin, self.pmax], [self.xmin, self.xmax]
                else:
                    self.hist2d = Fast2DHist(self.y_values, self.x_values, self.pmin,self.pmax, self.GetPlotParam('pbins'), self.xmin,self.xmax, self.GetPlotParam('xbins')), [self.pmin, self.pmax], [self.xmin, self.xmax]

            try:
                if self.GetPlotParam('masked'):
                    zval = ma.masked_array(self.hist2d[0])
                    zval[zval == 0] = ma.masked
                    zval *= float(zval.max())**(-1)
                    tmplist = [zval[~zval.mask].min(), zval.max()]
                else:
                    zval = np.copy(self.hist2d[0])
                    zval[zval==0] = 0.5
                    zval *= float(zval.max())**(-1)
                    tmplist = [zval.min(), zval.max()]
            except ValueError:
                tmplist = [0.1,1]
            self.hist2d = zval, self.hist2d[1], self.hist2d[2], tmplist
            self.parent.DataDict[self.key_name] = self.hist2d

    def UpdateLabelsandColors(self):
        # set the colors
        if self.GetPlotParam('prtl_type') == 0: #protons
            self.energy_color = self.parent.ion_color
        else: #electons
            self.energy_color = self.parent.electron_color

        for line in self.IntRegionLines:
            line.set_color(self.energy_color)
        #set the xlabels
        if self.parent.MainParamDict['DoLorentzBoost'] and np.abs(self.parent.MainParamDict['GammaBoost'])>1E-8:
            self.x_label = r'$x\prime\ [c/\omega_{\rm pe}]$'
        else:
            self.x_label = r'$x\ [c/\omega_{\rm pe}]$'
        #set the ylabel
        self.y_label  = self.ylabel_list[self.parent.MainParamDict['DoLorentzBoost']][self.GetPlotParam('prtl_type')][self.GetPlotParam('mom_dim')]

    def draw(self):
        # In order to speed up the plotting, we only recalculate everything
        # if necessary.
        self.IntRegionLines = []
        # Figure out the color and ylabel
        # Choose the particle type and px, py, or pz
        self.UpdateLabelsandColors()

        self.xmin = self.hist2d[2][0]
        self.xmax = self.hist2d[2][-1]

        self.ymin = self.hist2d[1][0]
        self.ymax = self.hist2d[1][-1]


        if self.GetPlotParam('masked'):
            self.tick_color = 'k'
        else:
            self.tick_color = 'white'


        self.clim = list(self.hist2d[3])

        if self.GetPlotParam('set_v_min'):
            self.clim[0] = 10**self.GetPlotParam('v_min')
        if self.GetPlotParam('set_v_max'):
            self.clim[1] = 10**self.GetPlotParam('v_max')


        self.gs = gridspec.GridSpecFromSubplotSpec(100,100, subplot_spec = self.parent.gs0[self.FigWrap.pos])#, bottom=0.2,left=0.1,right=0.95, top = 0.95)

        if self.parent.MainParamDict['LinkSpatial'] == 1:
            if self.FigWrap.pos == self.parent.first_x:
                self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]])
            else:
                self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]], sharex = self.parent.SubPlotList[self.parent.first_x[0]][self.parent.first_x[1]].graph.axes)
        else:
            self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]])

        self.cax = self.axes.imshow(self.hist2d[0],
                                    cmap = new_cmaps.cmaps[self.parent.MainParamDict['ColorMap']],
                                    norm = self.norm(), origin = 'lower',
                                    aspect = 'auto',
                                    interpolation=self.GetPlotParam('interpolation'))

        self.cax.set_extent([self.xmin, self.xmax, self.ymin, self.ymax])

        self.cax.set_clim(self.clim)

        self.shock_line = self.axes.axvline(self.parent.shock_loc, linewidth = 1.5, linestyle = '--', color = self.parent.shock_color, path_effects=[PathEffects.Stroke(linewidth=2, foreground='k'),
                   PathEffects.Normal()])
        if not self.GetPlotParam('show_shock'):
            self.shock_line.set_visible(False)




        self.axC = self.figure.add_subplot(self.gs[self.parent.cbar_extent[0]:self.parent.cbar_extent[1], self.parent.cbar_extent[2]:self.parent.cbar_extent[3]])
        self.parent.cbarList.append(self.axC)

        # Technically I should use the colorbar class here,
        # but I found it annoying in some of it's limitations.
        if self.parent.MainParamDict['HorizontalCbars']:
            self.cbar = self.axC.imshow(self.gradient, aspect='auto',
                                    cmap=new_cmaps.cmaps[self.parent.MainParamDict['ColorMap']])
            # Make the colobar axis more like the real colorbar
            self.axC.tick_params(axis='x',
                                which = 'both', # bothe major and minor ticks
                                top = False, # turn off top ticks
                                labelsize=self.parent.MainParamDict['NumFontSize'])

            self.axC.tick_params(axis='y',          # changes apply to the y-axis
                                which='both',      # both major and minor ticks are affected
                                left=False,      # ticks along the bottom edge are off
                                right=False,         # ticks along the top edge are off
                                labelleft=False)

        else:
            self.cbar = self.axC.imshow(np.transpose(self.gradient)[::-1], aspect='auto',
                                    cmap=new_cmaps.cmaps[self.parent.MainParamDict['ColorMap']])
            # Make the colobar axis more like the real colorbar
            self.axC.tick_params(axis='x',
                                which = 'both', # bothe major and minor ticks
                                top = False, # turn off top ticks
                                bottom = False,
                                labelbottom = False,
                                labelsize=self.parent.MainParamDict['NumFontSize'])

            self.axC.tick_params(axis='y',          # changes apply to the y-axis
                                which='both',      # both major and minor ticks are affected
                                left=False,      # ticks along the bottom edge are off
                                right=True,         # ticks along the top edge are off
                                labelleft=False,
                                labelright=True,
                                labelsize=self.parent.MainParamDict['NumFontSize'])

        self.cbar.set_extent([0, 1.0, 0, 1.0])

        if not self.GetPlotParam('show_cbar'):
            self.axC.set_visible(False)

        if int(matplotlib.__version__[0]) < 2:
            self.axes.set_axis_bgcolor(self.GetPlotParam('face_color'))
        else:
            self.axes.set_facecolor(self.GetPlotParam('face_color'))
        self.axes.tick_params(labelsize = self.parent.MainParamDict['NumFontSize'], color=self.tick_color)
        self.axes.set_xlabel(self.x_label, labelpad = self.parent.MainParamDict['xLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
        self.axes.set_ylabel(self.y_label, labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])

        self.refresh()

    def refresh(self):
        '''This is a function that will be called only if self.axes already
        holds a density type plot. We only update things that have shown. If
        hasn't changed, or isn't viewed, don't touch it. The difference between this and last
        time, is that we won't actually do any drawing in the plot. The plot
        will be redrawn after all subplots data is changed. '''


        # Main goal, only change what is showing..
        self.xmin = self.hist2d[2][0]
        self.xmax = self.hist2d[2][-1]
        self.ymin = self.hist2d[1][0]
        self.ymax = self.hist2d[1][-1]
        self.clim = list(self.hist2d[3])

        self.cax.set_data(self.hist2d[0])

        self.cax.set_extent([self.xmin,self.xmax, self.ymin, self.ymax])


        if self.GetPlotParam('set_v_min'):
            self.clim[0] =  10**self.GetPlotParam('v_min')
        if self.GetPlotParam('set_v_max'):
            self.clim[1] =  10**self.GetPlotParam('v_max')

        self.cax.set_clim(self.clim)
        if self.GetPlotParam('show_cbar'):
            self.CbarTickFormatter()


        if self.GetPlotParam('show_shock'):
            self.shock_line.set_xdata([self.parent.shock_loc,self.parent.shock_loc])

        self.UpdateLabelsandColors()
        self.axes.set_xlabel(self.x_label, labelpad = self.parent.MainParamDict['xLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
        self.axes.set_ylabel(self.y_label, labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])

        if self.GetPlotParam('set_p_min'):
            self.ymin = self.GetPlotParam('p_min')
        if self.GetPlotParam('set_p_max'):
            self.ymax = self.GetPlotParam('p_max')
        if self.GetPlotParam('symmetric'):
            self.ymin = -max(abs(self.ymin), abs(self.ymax))
            self.ymax = abs(self.ymin)
        self.axes.set_ylim(self.ymin, self.ymax)

        if self.parent.MainParamDict['SetxLim'] and self.parent.MainParamDict['LinkSpatial'] == 1:
            if self.parent.MainParamDict['xLimsRelative']:
                self.axes.set_xlim(self.parent.MainParamDict['xLeft'] + self.parent.shock_loc,
                                   self.parent.MainParamDict['xRight'] + self.parent.shock_loc)
            else:
                self.axes.set_xlim(self.parent.MainParamDict['xLeft'], self.parent.MainParamDict['xRight'])

        else:
            self.axes.set_xlim(self.xmin,self.xmax)

    def CbarTickFormatter(self):
        ''' A helper function that sets the cbar ticks & labels. This used to be
        easier, but because I am no longer using the colorbar class i have to do
        stuff manually.'''
        clim = np.copy(self.cax.get_clim())
        if self.GetPlotParam('show_cbar'):
            if self.GetPlotParam('cnorm_type') == "Log":
                if self.parent.MainParamDict['HorizontalCbars']:
                    self.cbar.set_extent([np.log10(clim[0]),np.log10(clim[1]),0,1])
                    self.axC.set_xlim(np.log10(clim[0]),np.log10(clim[1]))
                    self.axC.xaxis.set_label_position("top")
                    if self.GetPlotParam('prtl_type') ==0:
                        self.axC.set_xlabel(r'$\log{\ \ f_i(p)}$', size = self.parent.MainParamDict['AxLabelSize'])#, labelpad =15, rotation = -90)
                    else:
                        self.axC.set_xlabel(r'$\log{\ \ f_e(p)}$', size = self.parent.MainParamDict['AxLabelSize'])#, size = 12,labelpad =15, rotation = -90)

                else:
                    self.cbar.set_extent([0,1,np.log10(clim[0]),np.log10(clim[1])])
                    self.axC.set_ylim(np.log10(clim[0]),np.log10(clim[1]))
                    self.axC.locator_params(axis='y',nbins=6)
                    self.axC.yaxis.set_label_position("right")
                    if self.GetPlotParam('prtl_type') ==0:
                        self.axC.set_ylabel(r'$\log{\ \ f_i(p)}$', labelpad = self.parent.MainParamDict['cbarLabelPad'], rotation = -90, size = self.parent.MainParamDict['AxLabelSize'])
                    else:
                        self.axC.set_ylabel(r'$\log{\ \ f_e(p)}$', labelpad = self.parent.MainParamDict['cbarLabelPad'], rotation = -90, size = self.parent.MainParamDict['AxLabelSize'])

            else:# self.GetPlotParam('cnorm_type') == "Linear":
                if self.parent.MainParamDict['HorizontalCbars']:
                    self.cbar.set_extent([clim[0], clim[1], 0, 1])
                    self.axC.set_xlim(clim[0], clim[1])
                    if self.GetPlotParam('prtl_type') ==0:
                        self.axC.set_xlabel(r'$f_i(p)$', size = self.parent.MainParamDict['AxLabelSize'])
                    else:
                        self.axC.set_xlabel(r'$f_e(p)$', size = self.parent.MainParamDict['AxLabelSize'])

                else:
                    self.cbar.set_extent([0, 1, clim[0], clim[1]])
                    self.axC.set_ylim(clim[0], clim[1])
                    self.axC.locator_params(axis='y', nbins=6)
                    self.axC.yaxis.set_label_position("right")
                    if self.GetPlotParam('prtl_type') ==0:
                        self.axC.set_ylabel(r'$f_i(p)$', labelpad = self.parent.MainParamDict['cbarLabelPad'], rotation = -90, size = self.parent.MainParamDict['AxLabelSize'])
                    else:
                        self.axC.set_ylabel(r'$f_e(p)$', labelpad = self.parent.MainParamDict['cbarLabelPad'], rotation = -90, size = self.parent.MainParamDict['AxLabelSize'])



    def GetPlotParam(self, keyname):
        return self.FigWrap.GetPlotParam(keyname)

    def SetPlotParam(self, keyname, value,  update_plot = True):
        self.FigWrap.SetPlotParam(keyname, value,  update_plot = update_plot)

    def OpenSettings(self):
        if self.settings_window is None:
            self.settings_window = PhaseSettings(self)
        else:
            self.settings_window.destroy()
            self.settings_window = PhaseSettings(self)


class PhaseSettings(Tk.Toplevel):
    def __init__(self, parent):
        self.parent = parent
        Tk.Toplevel.__init__(self)

        self.wm_title('Phase Plot (%d,%d) Settings' % self.parent.FigWrap.pos)
        self.parent = parent
        frm = ttk.Frame(self)
        frm.pack(fill=Tk.BOTH, expand=True)
        self.protocol('WM_DELETE_WINDOW', self.OnClosing)
        self.bind('<Return>', self.TxtEnter)

        # Create the OptionMenu to chooses the Interpolation Type:
        self.InterpolVar = Tk.StringVar(self)
        self.InterpolVar.set(self.parent.GetPlotParam('interpolation')) # default value
        self.InterpolVar.trace('w', self.InterpolChanged)

        ttk.Label(frm, text="Interpolation Method:").grid(row=0, column = 2)
        InterplChooser = ttk.OptionMenu(frm, self.InterpolVar, self.parent.GetPlotParam('interpolation'), *tuple(self.parent.InterpolationMethods))
        InterplChooser.grid(row =0, column = 3, sticky = Tk.W + Tk.E)

        # Create the OptionMenu to chooses the Chart Type:
        self.ctypevar = Tk.StringVar(self)
        self.ctypevar.set(self.parent.chartType) # default value
        self.ctypevar.trace('w', self.ctypeChanged)

        ttk.Label(frm, text="Choose Chart Type:").grid(row=0, column = 0)
        cmapChooser = ttk.OptionMenu(frm, self.ctypevar, self.parent.chartType, *tuple(self.parent.ChartTypes))
        cmapChooser.grid(row =0, column = 1, sticky = Tk.W + Tk.E)


        # the Radiobox Control to choose the particle
        self.prtlList = ['ion', 'electron']
        self.pvar = Tk.IntVar()
        self.pvar.set(self.parent.GetPlotParam('prtl_type'))

        ttk.Label(frm, text='Particle:').grid(row = 1, sticky = Tk.W)

        for i in range(len(self.prtlList)):
            ttk.Radiobutton(frm,
                text=self.prtlList[i],
                variable=self.pvar,
                command = self.RadioPrtl,
                value=i).grid(row = 2+i, sticky =Tk.W)

        # the Radiobox Control to choose the momentum dim
        self.dimList = ['x-px', 'x-py', 'x-pz']
        self.dimvar = Tk.IntVar()
        self.dimvar.set(self.parent.GetPlotParam('mom_dim'))

        ttk.Label(frm, text='Dimenison:').grid(row = 1, column = 1, sticky = Tk.W)

        for i in range(len(self.dimList)):
            ttk.Radiobutton(frm,
                text=self.dimList[i],
                variable=self.dimvar,
                command = self.RadioDim,
                value=i).grid(row = 2+i, column = 1, sticky = Tk.W)


        # Control whether or not Cbar is shown
        self.CbarVar = Tk.IntVar()
        self.CbarVar.set(self.parent.GetPlotParam('show_cbar'))
        cb = ttk.Checkbutton(frm, text = "Show Color bar",
                        variable = self.CbarVar,
                        command = self.CbarHandler)
        cb.grid(row = 6, sticky = Tk.W)

        # show shock
        self.ShockVar = Tk.IntVar()
        self.ShockVar.set(self.parent.GetPlotParam('show_shock'))
        cb = ttk.Checkbutton(frm, text = "Show Shock",
                        variable = self.ShockVar,
                        command = self.ShockVarHandler)
        cb.grid(row = 6, column = 1, sticky = Tk.W)
        # Use full div cmap
        self.SymVar = Tk.IntVar()
        self.SymVar.set(self.parent.GetPlotParam('symmetric'))
        cb = ttk.Checkbutton(frm, text = "Symmetric about zero",
                        variable = self.SymVar,
                        command = self.SymmetricHandler)
        cb.grid(row = 8, column = 1, sticky = Tk.W)


        # Control if the plot is weightedd
        self.WeightVar = Tk.IntVar()
        self.WeightVar.set(self.parent.GetPlotParam('weighted'))
        cb = ttk.Checkbutton(frm, text = "Weight by charge",
                        variable = self.WeightVar,
                        command = lambda:
                        self.parent.SetPlotParam('weighted', self.WeightVar.get()))
        cb.grid(row = 7, sticky = Tk.W)

        # Show energy integration region
        #self.IntRegVar = Tk.IntVar()
        #self.IntRegVar.set(self.parent.GetPlotParam('show_int_region'))
        cb = ttk.Checkbutton(frm, text = "Show Energy Region",
                        variable = self.parent.IntRegVar)#,
        #               command = self.ShowIntRegionHandler)
        cb.grid(row = 7, column = 1, sticky = Tk.W)

        # control mask
        self.MaskVar = Tk.IntVar()
        self.MaskVar.set(self.parent.GetPlotParam('masked'))
        cb = ttk.Checkbutton(frm, text = "Mask Zeros",
                        variable = self.MaskVar,
                        command = lambda:
                        self.parent.SetPlotParam('masked', self.MaskVar.get()))
        cb.grid(row = 8, sticky = Tk.W)

        self.TrueVar = Tk.IntVar()
        self.TrueVar.set(1)
        self.pBins = Tk.StringVar()
        self.pBins.set(str(self.parent.GetPlotParam('pbins')))
        ttk.Label(frm, text ='# of pbins').grid(row = 9, column = 0, sticky = Tk.W)
        ttk.Entry(frm, textvariable=self.pBins, width=7).grid(row = 9, column = 1)

        self.xBins = Tk.StringVar()
        self.xBins.set(str(self.parent.GetPlotParam('xbins')))
        ttk.Label(frm, text ='# of xbins').grid(row = 10, column = 0, sticky = Tk.W)
        ttk.Entry(frm, textvariable=self.xBins, width=7).grid(row = 10, column = 1)


#        ttk.Label(frm, text = 'If the zero values are not masked they are set to z_min/2').grid(row =9, columnspan =2)
    # Define functions for the events
        # Now the field lim
        self.setVminVar = Tk.IntVar()
        self.setVminVar.set(self.parent.GetPlotParam('set_v_min'))
        self.setVminVar.trace('w', self.setVminChanged)

        self.setVmaxVar = Tk.IntVar()
        self.setVmaxVar.set(self.parent.GetPlotParam('set_v_max'))
        self.setVmaxVar.trace('w', self.setVmaxChanged)



        self.Vmin = Tk.StringVar()
        self.Vmin.set(str(self.parent.GetPlotParam('v_min')))

        self.Vmax = Tk.StringVar()
        self.Vmax.set(str(self.parent.GetPlotParam('v_max')))


        cb = ttk.Checkbutton(frm, text ='Set log(f) min',
                        variable = self.setVminVar)
        cb.grid(row = 3, column = 2, sticky = Tk.W)
        self.VminEnter = ttk.Entry(frm, textvariable=self.Vmin, width=7)
        self.VminEnter.grid(row = 3, column = 3)

        cb = ttk.Checkbutton(frm, text ='Set log(f) max',
                        variable = self.setVmaxVar)
        cb.grid(row = 4, column = 2, sticky = Tk.W)

        self.VmaxEnter = ttk.Entry(frm, textvariable=self.Vmax, width=7)
        self.VmaxEnter.grid(row = 4, column = 3)

        # Now the y lim
        self.setPminVar = Tk.IntVar()
        self.setPminVar.set(self.parent.GetPlotParam('set_p_min'))
        self.setPminVar.trace('w', self.setPminChanged)

        self.setPmaxVar = Tk.IntVar()
        self.setPmaxVar.set(self.parent.GetPlotParam('set_p_max'))
        self.setPmaxVar.trace('w', self.setPmaxChanged)



        self.Pmin = Tk.StringVar()
        self.Pmin.set(str(self.parent.GetPlotParam('p_min')))

        self.Pmax = Tk.StringVar()
        self.Pmax.set(str(self.parent.GetPlotParam('p_max')))


        cb = ttk.Checkbutton(frm, text ='Set y_axis min',
                        variable = self.setPminVar)
        cb.grid(row = 5, column = 2, sticky = Tk.W)
        self.PminEnter = ttk.Entry(frm, textvariable=self.Pmin, width=7)
        self.PminEnter.grid(row = 5, column = 3)

        cb = ttk.Checkbutton(frm, text ='Set y_axis max',
                        variable = self.setPmaxVar)
        cb.grid(row = 6, column = 2, sticky = Tk.W)

        self.PmaxEnter = ttk.Entry(frm, textvariable=self.Pmax, width=7)
        self.PmaxEnter.grid(row = 6, column = 3)

        # Now the E lim
        self.setEminVar = Tk.IntVar()
        self.setEminVar.set(self.parent.GetPlotParam('set_E_min'))
        self.setEminVar.trace('w', self.setEminChanged)

        self.setEmaxVar = Tk.IntVar()
        self.setEmaxVar.set(self.parent.GetPlotParam('set_E_max'))
        self.setEmaxVar.trace('w', self.setEmaxChanged)


        self.Emin = Tk.StringVar()
        self.Emin.set(str(self.parent.GetPlotParam('E_min')))

        self.Emax = Tk.StringVar()
        self.Emax.set(str(self.parent.GetPlotParam('E_max')))


        cb = ttk.Checkbutton(frm, text ='Set E_min (m_e c^2)',
                        variable = self.setEminVar)
        cb.grid(row = 7, column = 2, sticky = Tk.W)
        self.EminEnter = ttk.Entry(frm, textvariable=self.Emin, width=7)
        self.EminEnter.grid(row = 7, column = 3)

        cb = ttk.Checkbutton(frm, text ='Set E_max (m_e c^2)',
                        variable = self.setEmaxVar)
        cb.grid(row = 8, column = 2, sticky = Tk.W)

        self.EmaxEnter = ttk.Entry(frm, textvariable=self.Emax, width=7)
        self.EmaxEnter.grid(row = 8, column = 3)


    def ShockVarHandler(self, *args):
        if self.parent.GetPlotParam('show_shock')== self.ShockVar.get():
            pass
        else:
            self.parent.shock_line.set_visible(self.ShockVar.get())
            self.parent.SetPlotParam('show_shock', self.ShockVar.get())


    def CbarHandler(self, *args):
        if self.parent.GetPlotParam('show_cbar')== self.CbarVar.get():
            pass
        else:
            self.parent.axC.set_visible(self.CbarVar.get())
            self.parent.SetPlotParam('show_cbar', self.CbarVar.get(), update_plot =self.parent.GetPlotParam('twoD'))


    def ctypeChanged(self, *args):
        if self.ctypevar.get() == self.parent.chartType:
            pass
        else:
            self.parent.ChangePlotType(self.ctypevar.get())
            self.destroy()

    def InterpolChanged(self, *args):
        if self.InterpolVar.get() == self.parent.GetPlotParam('interpolation'):
            pass
        else:
            self.parent.cax.set_interpolation(self.InterpolVar.get())
            self.parent.SetPlotParam('interpolation', self.InterpolVar.get())

    def RadioPrtl(self):
        if self.pvar.get() == self.parent.GetPlotParam('prtl_type'):
            pass
        else:

            self.parent.SetPlotParam('prtl_type', self.pvar.get(), update_plot =  False)
            self.parent.UpdateLabelsandColors()
            self.parent.axes.set_ylabel(self.parent.y_label, labelpad = self.parent.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.parent.MainParamDict['AxLabelSize'])
            #self.parent.lineleft.set_color(self.parent.energy_color)
            #self.parent.lineright.set_color(self.parent.energy_color)
            self.parent.SetPlotParam('prtl_type', self.pvar.get())

    def RadioDim(self):
        if self.dimvar.get() == self.parent.GetPlotParam('mom_dim'):
            pass
        else:
            self.parent.SetPlotParam('mom_dim', self.dimvar.get(), update_plot = False)
            self.parent.UpdateLabelsandColors()
            self.parent.axes.set_ylabel(self.parent.y_label, labelpad = self.parent.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.parent.MainParamDict['AxLabelSize'])
            self.parent.SetPlotParam('mom_dim', self.dimvar.get())

    def SymmetricHandler(self, *args):
        if self.parent.GetPlotParam('symmetric') == self.SymVar.get():
            pass
        else:
            self.parent.SetPlotParam('symmetric', self.SymVar.get(), update_plot = True)

    def setVminChanged(self, *args):
        if self.setVminVar.get() == self.parent.GetPlotParam('set_v_min'):
            pass
        else:
            self.parent.SetPlotParam('set_v_min', self.setVminVar.get())

    def setVmaxChanged(self, *args):
        if self.setVmaxVar.get() == self.parent.GetPlotParam('set_v_max'):
            pass
        else:
            self.parent.SetPlotParam('set_v_max', self.setVmaxVar.get())

    def setPminChanged(self, *args):
        if self.setPminVar.get() == self.parent.GetPlotParam('set_p_min'):
            pass
        else:
            self.parent.SetPlotParam('set_p_min', self.setPminVar.get())

    def setPmaxChanged(self, *args):
        if self.setPmaxVar.get() == self.parent.GetPlotParam('set_p_max'):
            pass
        else:
            self.parent.SetPlotParam('set_p_max', self.setPmaxVar.get())

    def setEminChanged(self, *args):
        if self.setEminVar.get() == self.parent.GetPlotParam('set_E_min'):
            pass
        else:
            self.parent.SetPlotParam('set_E_min', self.setEminVar.get())

    def setEmaxChanged(self, *args):
        if self.setEmaxVar.get() == self.parent.GetPlotParam('set_E_max'):
            pass
        else:
            self.parent.SetPlotParam('set_E_max', self.setEmaxVar.get())


    def TxtEnter(self, e):
        self.FieldsCallback()

    def FieldsCallback(self):
        #### First set the Float Values
        tkvarLimList = [self.Vmin, self.Vmax, self.Pmin, self.Pmax, self.Emin, self.Emax]
        plot_param_List = ['v_min', 'v_max', 'p_min', 'p_max', 'E_min', 'E_max']
        tkvarSetList = [self.setVminVar, self.setVmaxVar, self.setPminVar, self.setPmaxVar, self.setEminVar, self.setEmaxVar]
        to_reload = False
        for j in range(len(tkvarLimList)):
            try:
            #make sure the user types in a float
                if np.abs(float(tkvarLimList[j].get()) - self.parent.GetPlotParam(plot_param_List[j])) > 1E-4:
                    self.parent.SetPlotParam(plot_param_List[j], float(tkvarLimList[j].get()), update_plot = False)
                    to_reload += True*tkvarSetList[j].get()

            except ValueError:
                #if they type in random stuff, just set it ot the param value
                tkvarLimList[j].set(str(self.parent.GetPlotParam(plot_param_List[j])))

        intVarList = [self.pBins, self.xBins]
        intParamList = ['pbins', 'xbins']
        for j in range(len(intVarList)):
            try:
            #make sure the user types in a float
                intVarList[j].set(str(int(float(intVarList[j].get()))))
                if int(float(intVarList[j].get())) - int(self.parent.GetPlotParam(intParamList[j])) != 0:

                    self.parent.SetPlotParam(intParamList[j], int(float(intVarList[j].get())), update_plot = False)
                    to_reload += True

            except ValueError:
            #    print hi
                #if they type in random stuff, just set it ot the param value
                intVarList[j].set(str(self.parent.GetPlotParam(intVarList[j])))


        if to_reload:
            self.parent.SetPlotParam('v_min', self.parent.GetPlotParam('v_min'))

    def OnClosing(self):
        self.parent.settings_window = None
        self.destroy()
