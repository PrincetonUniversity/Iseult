#!/usr/bin/env python
import Tkinter as Tk
import ttk as ttk
import matplotlib
import numpy as np
import numpy.ma as ma
import new_cmaps
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as PathEffects
import matplotlib.transforms as mtransforms

class  MomentsPanel:
    # A dictionary of all of the parameters for this plot with the default parameters

    plot_param_dict = {'twoD': 0,
                       'm_type': 0, # 0 = average_velocity, 1 = average_momentum
                       'v_min': 0,
                       'v_max' : 10,
                       'set_v_min': False,
                       'set_v_max': False,
                       'show_x': True,
                       'show_y': False,
                       'show_z': False,
                       'show_ions': True,
                       'show_electrons': True,
                       'UpstreamFrame': True,
                       'weighted': False,
                       'xbins': 100,
                       'spatial_x': True,
                       'show_legend': True,
                       'spatial_y': False}

    # We need the types of all the parameters for the config file
    BoolList = ['twoD', 'set_v_min', 'set_v_max',
                'UpstreamFrame', 'show_ions', 'show_electrons', 'show_x',
                'show_y', 'show_z', 'spatial_x', 'spatial_y', 'weighted', 'show_legend']
    IntList = ['m_type','xbins']
    FloatList = ['v_min', 'v_max']
    StrList = []

    def __init__(self, parent, figwrapper):
        self.settings_window = None
        self.FigWrap = figwrapper
        self.parent = parent

        self.ylabel_list = [[r'$\langle \beta \rangle$',r'$\langle \beta \rangle$' + ' upstream frame'],[r'$\langle m\gamma\beta\rangle/m_i$', r'$\langle m\gamma\beta \rangle/m_i$' + ' upstream frame']]
        self.ChartTypes = self.FigWrap.PlotTypeDict.keys()
        self.chartType = self.FigWrap.chartType
        self.figure = self.FigWrap.figure
        self.SetPlotParam('spatial_y', self.GetPlotParam('twoD'), update_plot = False)
        self.InterpolationMethods = ['nearest', 'bilinear', 'bicubic', 'spline16',
            'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric',
            'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']

    def ChangePlotType(self, str_arg):
        self.FigWrap.ChangeGraph(str_arg)
    def stepify(self, bins, hist):
        # make it a step, there probably is a much better way to do this
        tmp_hist = np.ones(2*len(hist))
        tmp_bin = np.ones(2*len(hist))
        for j in range(len(hist)):
            tmp_hist[2*j] = hist[j]
            tmp_hist[2*j+1] = hist[j]

            tmp_bin[2*j] = bins[j]
            if j != 0:
                tmp_bin[2*j-1] = bins[j]
            if j == len(hist)-1:
                tmp_bin[2*j+1] = bins[j+1]
        return tmp_bin, tmp_hist

    def set_plot_keys(self):

        '''A helper function that will insure that each hdf5 file will only be
        opened once per time step'''
        self.arrs_needed = ['c_omp', 'bx', 'istep', 'me', 'mi', 'gamma0']

        # First see if we will need to know the energy of the particle
        # (required to Lorentz boost the momentum)

#        self.ShowingSomething = self.GetPlotParam('show_electrons') or self.GetPlotParam('show_ions')
#        tmp_bool = self.GetPlotParam('show_x') or self.GetPlotParam('show_y') or self.GetPlotParam('show_z')

#        self.ShowingSomething = self.ShowingSomething and tmp_bool

#        if self.ShowingSomething:
        if self.GetPlotParam('weighted'):
            self.arrs_needed.append('chi')
            self.arrs_needed.append('che')


        self.arrs_needed.append('xi')
        self.arrs_needed.append('xe')

        self.arrs_needed.append('ui')
        self.arrs_needed.append('vi')
        self.arrs_needed.append('wi')

        self.arrs_needed.append('ue')
        self.arrs_needed.append('ve')
        self.arrs_needed.append('we')


        '''
        Needs_All = self.GetPlotParam('m_type') == 1 * self.GetPlotParam('UpstreamFrame')
        Needs_All = Needs_All or self.GetPlotParam('m_type') == 0

        if self.GetPlotParam('show_ions'):
            self.arrs_needed.append('xi')
            if self.GetPlotParam('weighted'):
                self.arrs_needed.append('chi')
            if Needs_All:
                self.arrs_needed.append('ui')
                self.arrs_needed.append('vi')
                self.arrs_needed.append('wi')

            else:
                if self.GetPlotParam('show_x'):
                    self.arrs_needed.append('ui')
                if self.GetPlotParam('show_y'):
                    self.arrs_needed.append('vi')
                if self.GetPlotParam('show_z'):
                    self.arrs_needed.append('wi')

        if self.GetPlotParam('show_electrons'):
            self.arrs_needed.append('xe')
            if self.GetPlotParam('weighted'):
                self.arrs_needed.append('che')
            if Needs_All:
                self.arrs_needed.append('ue')
                self.arrs_needed.append('ve')
                self.arrs_needed.append('we')

            else:
                if self.GetPlotParam('show_x'):
                    self.arrs_needed.append('ue')
                if self.GetPlotParam('show_y'):
                    self.arrs_needed.append('ve')
                if self.GetPlotParam('show_z'):
                    self.arrs_needed.append('we')
            '''
        return self.arrs_needed

    def LoadData(self):
        ''' A helper function that checks if the histogram has
        already been calculated and if it hasn't, it calculates
        it then stores it.'''

        self.c_omp = self.FigWrap.LoadKey('c_omp')[0]
        self.istep = self.FigWrap.LoadKey('istep')[0]
        self.gamma0 = self.FigWrap.LoadKey('gamma0')[0]
        self.mass_ratio = self.FigWrap.LoadKey('mi')[0]/self.FigWrap.LoadKey('me')[0]
        self.iweights = None
        self.eweights = None
        self.x_bins = None
        self.ix = None
        self.iy = None
        self.iz = None
        self.ex = None
        self.ey = None
        self.ez = None

        if 'xaxis_values' in self.parent.DataDict.keys():
            self.xaxis_values = self.parent.DataDict['xaxis_values']
        else:
            # x-values haven't been calculated yet, generate them then save them to the dictionary for later.
            self.xaxis_values = np.arange(self.FigWrap.LoadKey('bx')[0,:,:].shape[1])/self.c_omp*self.istep
        #            print self.xaxis_values
            self.parent.DataDict['xaxis_values'] = np.copy(self.xaxis_values)
        self.xmin = 0
        self.xmax = self.xaxis_values[-1]
        self.key_name = str(self.GetPlotParam('xbins'))
        if self.GetPlotParam('UpstreamFrame'):
            self.key_name += 'Up'
        if self.GetPlotParam('weighted'):
            self.key_name += 'Weight'
        if self.GetPlotParam('m_type') == 0:
            self.key_name += 'vel'
        elif self.GetPlotParam('m_type') == 1:
            self.key_name += 'mom'
        # See if anything at all needs loading.

#        if self.ShowingSomething:
        if self.key_name in self.parent.DataDict.keys():
            self.x_bins = self.parent.DataDict[self.key_name+'x_bins']

            self.ex = self.parent.DataDict[self.key_name+'ex']
            self.ey = self.parent.DataDict[self.key_name+'ey']
            self.ez = self.parent.DataDict[self.key_name+'ez']


            self.ix = self.parent.DataDict[self.key_name+'ix']
            self.iy = self.parent.DataDict[self.key_name+'iy']
            self.iz = self.parent.DataDict[self.key_name+'iz']


        else:
            xi = self.FigWrap.LoadKey('xi')/self.c_omp
            ui = self.FigWrap.LoadKey('ui')
            vi = self.FigWrap.LoadKey('vi')
            wi = self.FigWrap.LoadKey('wi')

            xe = self.FigWrap.LoadKey('xe')/self.c_omp
            ue = self.FigWrap.LoadKey('ue')
            ve = self.FigWrap.LoadKey('ve')
            we = self.FigWrap.LoadKey('we')

            if self.GetPlotParam('weighted'):
                i_weights = self.FigWrap.LoadKey('chi')
                e_weights = self.FigWrap.LoadKey('che')

            else:
                i_weights = np.ones(len(xi))
                e_weights = np.ones(len(xe))

            tmp_hist, self.x_bins = np.histogram(xi, bins = self.GetPlotParam('xbins'))

            self.ex = np.zeros(len(tmp_hist))
            self.ey = np.zeros(len(tmp_hist))
            self.ez = np.zeros(len(tmp_hist))
            self.ix = np.zeros(len(tmp_hist))
            self.iy = np.zeros(len(tmp_hist))
            self.iz = np.zeros(len(tmp_hist))

            if self.GetPlotParam('m_type') == 1 and not self.GetPlotParam('UpstreamFrame'):
                # Momentum in the downstream frame. Easiest.
                for j in range(len(self.x_bins)-1):
                    i_find_loc = xi>=self.x_bins[j]
                    i_find_loc *= xi<=self.x_bins[j+1]
                    e_find_loc = xe>=self.x_bins[j]
                    e_find_loc *= xe<=self.x_bins[j+1]
                    i_ind = np.where(i_find_loc)[0]
                    e_ind = np.where(e_find_loc)[0]


                    self.ix[j]=np.sum(ui[i_ind]*i_weights[i_ind])/(len(i_ind))
                    self.iy[j]=np.sum(vi[i_ind]*i_weights[i_ind])/(len(i_ind))
                    self.iz[j]=np.sum(wi[i_ind]*i_weights[i_ind])/(len(i_ind))

                    self.ex[j]=np.sum(ue[e_ind]*e_weights[e_ind])/(len(e_ind))/self.mass_ratio
                    self.ey[j]=np.sum(ve[e_ind]*e_weights[e_ind])/(len(e_ind))/self.mass_ratio
                    self.ez[j]=np.sum(we[e_ind]*e_weights[e_ind])/(len(e_ind))/self.mass_ratio


            else:

                gi= np.sqrt(ui**2+vi**2+wi**2+1)
                ge= np.sqrt(ue**2+ve**2+we**2+1)


                if self.GetPlotParam('m_type') == 0 and not self.GetPlotParam('UpstreamFrame'):
                    # calculate the velocities in downstream frame from the momenta
                    vex = ue/ge
                    vey = ve/ge
                    vez = we/ge

                    vix = ui/gi
                    viy = vi/gi
                    viz = wi/gi

                    for j in range(len(self.x_bins)-1):
                        i_find_loc = xi>=self.x_bins[j]
                        i_find_loc *= xi<=self.x_bins[j+1]
                        e_find_loc = xe>=self.x_bins[j]
                        e_find_loc *= xe<=self.x_bins[j+1]
                        i_ind = np.where(i_find_loc)[0]
                        e_ind = np.where(e_find_loc)[0]


                        self.ix[j]=np.sum(vix[i_ind]*i_weights[i_ind])/(len(i_ind))
                        self.iy[j]=np.sum(viy[i_ind]*i_weights[i_ind])/(len(i_ind))
                        self.iz[j]=np.sum(viz[i_ind]*i_weights[i_ind])/(len(i_ind))

                        self.ex[j]=np.sum(vex[e_ind]*e_weights[e_ind])/(len(e_ind))
                        self.ey[j]=np.sum(vey[e_ind]*e_weights[e_ind])/(len(e_ind))
                        self.ez[j]=np.sum(vez[e_ind]*e_weights[e_ind])/(len(e_ind))

                if self.GetPlotParam('m_type') == 0 and self.GetPlotParam('UpstreamFrame'):
                    # Velocities in the upstream frame

                    # calculate the velocities in downstream frame from the momenta
                    vex = ue/ge
                    vey = ve/ge
                    vez = we/ge

                    vix = ui/gi
                    viy = vi/gi
                    viz = wi/gi

                    # Boost into the upstream frame:
                    betaBoost = -np.sqrt(1-1./self.gamma0**2)

                    # Now calulate the velocities in the boosted frames
                    tmp_ihelper = 1-vix*betaBoost
                    vix = (vix-betaBoost)/tmp_ihelper
                    viy = viy/self.gamma0/tmp_ihelper
                    viz = viz/self.gamma0/tmp_ihelper

                    tmp_ehelper = 1-vex*betaBoost
                    vex = (vex-betaBoost)/tmp_ehelper
                    vey = vey/self.gamma0/tmp_ehelper
                    vez = vez/self.gamma0/tmp_ehelper

                    for j in range(len(self.x_bins)-1):
                        i_find_loc = xi>=self.x_bins[j]
                        i_find_loc *= xi<=self.x_bins[j+1]
                        e_find_loc = xe>=self.x_bins[j]
                        e_find_loc *= xe<=self.x_bins[j+1]
                        i_ind = np.where(i_find_loc)[0]
                        e_ind = np.where(e_find_loc)[0]


                        self.ix[j]=np.sum(vix[i_ind]*i_weights[i_ind])/(len(i_ind))
                        self.iy[j]=np.sum(viy[i_ind]*i_weights[i_ind])/(len(i_ind))
                        self.iz[j]=np.sum(viz[i_ind]*i_weights[i_ind])/(len(i_ind))

                        self.ex[j]=np.sum(vex[e_ind]*e_weights[e_ind])/(len(e_ind))
                        self.ey[j]=np.sum(vey[e_ind]*e_weights[e_ind])/(len(e_ind))
                        self.ez[j]=np.sum(vez[e_ind]*e_weights[e_ind])/(len(e_ind))


                if self.GetPlotParam('m_type') == 1 and self.GetPlotParam('UpstreamFrame'):
                    # Momentum in the upstream frame

                    # calculate the velocities in downstream frame from the momenta
                    vex = ue/ge
                    vey = ve/ge
                    vez = we/ge

                    vix = ui/gi
                    viy = vi/gi
                    viz = wi/gi

                    # Boost into the upstream frame:
                    betaBoost = -np.sqrt(1-1./self.gamma0**2)
                    rap_boost = np.arccosh(self.gamma0) #Boost rapidity

                    # Now calulate the velocities in the boosted frames
                    tmp_ihelper = 1-vix*betaBoost
                    vix = (vix-betaBoost)/tmp_ihelper
                    viy = viy/self.gamma0/tmp_ihelper
                    viz = viz/self.gamma0/tmp_ihelper

                    tmp_ehelper = 1-vex*betaBoost
                    vex = (vex-betaBoost)/tmp_ehelper
                    vey = vey/self.gamma0/tmp_ehelper
                    vez = vez/self.gamma0/tmp_ehelper

                    # Initial rapidity
                    rap_i = np.arccosh(gi)
                    rap_e = np.arccosh(ge)

                    gi_up = gi*self.gamma0-np.sign(ui)*np.sign(betaBoost)*np.sinh(rap_i)*np.sinh(rap_boost)/np.sqrt(1+(vi/ui)**2+(wi/ui)**2)
                    ge_up = ge*self.gamma0-np.sign(ue)*np.sign(betaBoost)*np.sinh(rap_e)*np.sinh(rap_boost)/np.sqrt(1+(ve/ue)**2+(we/ue)**2)

                    vix *= gi_up
                    viy *= gi_up
                    viz *= gi_up

                    vex *= ge_up
                    vey *= ge_up
                    vez *= ge_up

                    for j in range(len(self.x_bins)-1):
                        i_find_loc = xi>=self.x_bins[j]
                        i_find_loc *= xi<=self.x_bins[j+1]
                        e_find_loc = xe>=self.x_bins[j]
                        e_find_loc *= xe<=self.x_bins[j+1]
                        i_ind = np.where(i_find_loc)[0]
                        e_ind = np.where(e_find_loc)[0]


                        self.ix[j]=np.sum(vix[i_ind]*i_weights[i_ind])/(len(i_ind))
                        self.iy[j]=np.sum(viy[i_ind]*i_weights[i_ind])/(len(i_ind))
                        self.iz[j]=np.sum(viz[i_ind]*i_weights[i_ind])/(len(i_ind))

                        self.ex[j]=np.sum(vex[e_ind]*e_weights[e_ind])/(len(e_ind))/self.mass_ratio
                        self.ey[j]=np.sum(vey[e_ind]*e_weights[e_ind])/(len(e_ind))/self.mass_ratio
                        self.ez[j]=np.sum(vez[e_ind]*e_weights[e_ind])/(len(e_ind))/self.mass_ratio



                self.parent.DataDict[self.key_name+'x_bins'] = self.x_bins

                self.parent.DataDict[self.key_name+'ex'] = self.ex
                self.parent.DataDict[self.key_name+'ey'] = self.ey
                self.parent.DataDict[self.key_name+'ez'] = self.ez

                self.parent.DataDict[self.key_name+'ix'] = self.ix
                self.parent.DataDict[self.key_name+'iy'] = self.iy
                self.parent.DataDict[self.key_name+'iz'] = self.iz

    def draw(self):

        ''' A function that draws the data. In the interest in speeding up the
        code, draw should only be called when you want to recreate the whole
        figure, i.e. it  will be slow. Most times you will only want to update
        what has changed in the figure. This will be done in a function called
        refresh, that should be much much faster.'''

        # Set the tick color
        tick_color = 'black'

        # Create a gridspec to handle spacing better
        self.gs = gridspec.GridSpecFromSubplotSpec(100,100, subplot_spec = self.parent.gs0[self.FigWrap.pos])

        if self.parent.MainParamDict['LinkSpatial'] != 0 and self.parent.MainParamDict['LinkSpatial'] != 3:
            if self.FigWrap.pos == self.parent.first_x:
                self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]])
            else:
                self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]],
                sharex = self.parent.SubPlotList[self.parent.first_x[0]][self.parent.first_x[1]].graph.axes)
        else:
            self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]])


        self.ex_plot = self.axes.plot(0,0,
                                            ls= '-', # marker = '^', markeredgecolor = self.parent.electron_color,
                                            color = self.parent.electron_color,
                                            visible = self.GetPlotParam('show_electrons')*self.GetPlotParam('show_x'))
        self.ey_plot = self.axes.plot(0,0,
                                            ls='-', #marker = 'v', markeredgecolor= self.parent.electron_color,
                                            color = self.parent.electron_color,
                                            visible = self.GetPlotParam('show_electrons')*self.GetPlotParam('show_y'))
        self.ey_plot[0].set_dashes([1, 1])
        self.ez_plot = self.axes.plot(0,0,
                                            ls= '-', #marker = 'd', markeredgecolor = self.parent.electron_color,
                                            color = self.parent.electron_color,
                                            visible = self.GetPlotParam('show_electrons')*self.GetPlotParam('show_z'))
        self.ez_plot[0].set_dashes([5, 1])

        self.ix_plot = self.axes.plot(0,0,
                                       ls= '-', #marker = '<', markeredgecolor =  self.parent.ion_color,
                                       color = self.parent.ion_color,
                                       visible =  self.GetPlotParam('show_ions')*self.GetPlotParam('show_x'))
        self.iy_plot = self.axes.plot(0,0,
                                       ls= '-',# marker = '*',  markersize = 10, markeredgecolor =  self.parent.ion_color,
                                       color = self.parent.ion_color,
                                       visible =  self.GetPlotParam('show_ions')*self.GetPlotParam('show_y'))
        self.iy_plot[0].set_dashes([1, 1])

        self.iz_plot = self.axes.plot(0,0,
                                     ls= '--', #marker = 's', markeredgecolor  = self.parent.ion_color,
                                     color = self.parent.ion_color,
                                     visible =  self.GetPlotParam('show_ions')*self.GetPlotParam('show_z'))

        self.iz_plot[0].set_dashes([5, 1])

        self.ex_plot[0].set_data(*self.stepify(self.x_bins, self.ex))
        self.ey_plot[0].set_data(*self.stepify(self.x_bins, self.ey))
        self.ez_plot[0].set_data(*self.stepify(self.x_bins, self.ez))
        self.ix_plot[0].set_data(*self.stepify(self.x_bins, self.ix))
        self.iy_plot[0].set_data(*self.stepify(self.x_bins, self.iy))
        self.iz_plot[0].set_data(*self.stepify(self.x_bins, self.iz))

        self.axes.set_axis_bgcolor('lightgrey')
        self.axes.tick_params(labelsize = self.parent.MainParamDict['NumFontSize'], color=tick_color)



        # fancy code to make sure that matplotlib sets its limits
        # only based on visible lines
        self.key_list = ['show_x', 'show_y', 'show_z']
        self.ion_plot_list = [self.ix_plot[0], self.iy_plot[0], self.iz_plot[0]]
        self.e_plot_list = [self.ex_plot[0], self.ey_plot[0], self.ez_plot[0]]
        self.label_prefix = [r'$\langle\beta', r'$\langle p']
        self.label_mid = [r'_{i',r'_{e']
        self.label_suffix = [r'x}\rangle$',r'y\rangle}$',r'z}\rangle$']



        legend_handles = []
        legend_labels = []

        self.axes.dataLim = mtransforms.Bbox.unit()
        #self.axes.dataLim.update_from_data_xy(xy = np.vstack(self.ion_plot_list[0].get_data()).T, ignore=True)
        #self.axes.dataLim.update_from_data_xy(xy = np.vstack(self.e_plot_list[0].get_data()).T, ignore=True)
        if self.GetPlotParam('show_ions'):
            for i in range(len(self.ion_plot_list)):
                if self.GetPlotParam(self.key_list[i]):
                    xy = np.vstack(self.ion_plot_list[i].get_data()).T
                    self.axes.dataLim.update_from_data_xy(xy, ignore=False)
                    legend_handles.append(self.ion_plot_list[i])
                    legend_labels.append(self.label_prefix[self.GetPlotParam('m_type')]+self.label_mid[0]+self.label_suffix[i])


        if self.GetPlotParam('show_electrons'):
            for i in range(len(self.e_plot_list)):
                if self.GetPlotParam(self.key_list[i]):
                    xy = np.vstack(self.e_plot_list[i].get_data()).T
                    self.axes.dataLim.update_from_data_xy(xy, ignore=False)
                    legend_handles.append(self.e_plot_list[i])
                    legend_labels.append(self.label_prefix[self.GetPlotParam('m_type')]+self.label_mid[1]+self.label_suffix[i])


        self.axes.autoscale()

        # now make the legend
        if self.GetPlotParam('show_ions') and self.GetPlotParam('show_electrons'):
            self.legend = self.axes.legend(legend_handles, legend_labels,
            framealpha = .05, fontsize = 11, loc = 'best', ncol = 2)

        else:
            self.legend = self.axes.legend(legend_handles, legend_labels,
            framealpha = .05, fontsize = 11, loc = 'best')
        self.legend.get_frame().set_facecolor('k')
        self.legend.get_frame().set_linewidth(0.0)
        if not self.GetPlotParam('show_legend'):
            self.legend.set_visible(False)


        if self.GetPlotParam('set_v_min'):
            self.axes.set_ylim(ymin = self.GetPlotParam('v_min'))
        if self.GetPlotParam('set_v_max'):
            self.axes.set_ylim(ymax = self.GetPlotParam('v_max'))

        if self.parent.MainParamDict['SetxLim']:
            if self.parent.MainParamDict['xLimsRelative']:
                self.axes.set_xlim(self.parent.MainParamDict['xLeft'] + self.parent.shock_loc,
                                   self.parent.MainParamDict['xRight'] + self.parent.shock_loc)
            else:
                self.axes.set_xlim(self.parent.MainParamDict['xLeft'], self.parent.MainParamDict['xRight'])
        else:
            self.axes.set_xlim(self.xmin,self.xmax)



        self.axes.set_xlabel(r'$x\  [c/\omega_{pe}]$', labelpad = self.parent.MainParamDict['xLabelPad'], color = 'black')
        self.axes.set_ylabel(self.ylabel_list[self.GetPlotParam('m_type')][self.GetPlotParam('UpstreamFrame')], labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black')

    def refresh(self):

        '''This is a function that will be called only if self.axes already
        holds a total energy type plot. We only update things that have changed & are
        shown.  If hasn't changed or isn't shown, don't touch it. The difference
        between this and last time, is that we won't actually do any drawing in
        the plot. The plot will be redrawn after all subplots are refreshed. '''

        self.ex_plot[0].set_data(*self.stepify(self.x_bins, self.ex))
        self.ey_plot[0].set_data(*self.stepify(self.x_bins, self.ey))
        self.ez_plot[0].set_data(*self.stepify(self.x_bins, self.ez))
        self.ix_plot[0].set_data(*self.stepify(self.x_bins, self.ix))
        self.iy_plot[0].set_data(*self.stepify(self.x_bins, self.iy))
        self.iz_plot[0].set_data(*self.stepify(self.x_bins, self.iz))

        # fancy code to make sure that matplotlib sets its limits
        # based only on the visible lines.

        self.axes.dataLim = mtransforms.Bbox.unit()
        #self.axes.dataLim.update_from_data_xy(xy = np.vstack(self.ion_plot_list[0].get_data()).T, ignore=True)
        if self.GetPlotParam('show_ions'):
            for i in range(len(self.ion_plot_list)):
                if self.GetPlotParam(self.key_list[i]):
                    xy = np.vstack(self.ion_plot_list[i].get_data()).T
                    self.axes.dataLim.update_from_data_xy(xy, ignore=False)
        if self.GetPlotParam('show_electrons'):
            for i in range(len(self.e_plot_list)):
                if self.GetPlotParam(self.key_list[i]):
                    xy = np.vstack(self.e_plot_list[i].get_data()).T
                    self.axes.dataLim.update_from_data_xy(xy, ignore=False)


        self.axes.autoscale()

        # Set new lims if the user chooses to do so.
        if self.GetPlotParam('set_v_min'):
            self.axes.set_ylim(ymin = self.GetPlotParam('v_min'))
        if self.GetPlotParam('set_v_max'):
            self.axes.set_ylim(ymax = self.GetPlotParam('v_max'))


        if self.parent.MainParamDict['SetxLim']:
            if self.parent.MainParamDict['xLimsRelative']:
                self.axes.set_xlim(self.parent.MainParamDict['xLeft'] + self.parent.shock_loc,
                                   self.parent.MainParamDict['xRight'] + self.parent.shock_loc)
            else:
                self.axes.set_xlim(self.parent.MainParamDict['xLeft'], self.parent.MainParamDict['xRight'])
        else:
            self.axes.set_xlim(self.xmin,self.xmax)


    def GetPlotParam(self, keyname):
        return self.FigWrap.GetPlotParam(keyname)

    def SetPlotParam(self, keyname, value, update_plot = True):
        self.FigWrap.SetPlotParam(keyname, value, update_plot = update_plot)

    def OpenSettings(self):
        if self.settings_window is None:
            self.settings_window = MomentsSettings(self)
        else:
            self.settings_window.destroy()
            self.settings_window = MomentsSettings(self)


class MomentsSettings(Tk.Toplevel):
    def __init__(self, parent):
        self.parent = parent
        Tk.Toplevel.__init__(self)

        self.wm_title('Moments (%d,%d) Settings' % self.parent.FigWrap.pos)
        self.parent = parent
        frm = ttk.Frame(self)
        frm.pack(fill=Tk.BOTH, expand=True)
        self.protocol('WM_DELETE_WINDOW', self.OnClosing)
        self.bind('<Return>', self.TxtEnter)

        # Create the OptionMenu to chooses the Chart Type:
        self.ctypevar = Tk.StringVar(self)
        self.ctypevar.set(self.parent.chartType) # default value
        self.ctypevar.trace('w', self.ctypeChanged)

        ttk.Label(frm, text="Choose Chart Type:").grid(row=0, column = 0)
        ctypeChooser = apply(ttk.OptionMenu, (frm, self.ctypevar, self.parent.chartType) + tuple(self.parent.ChartTypes))
        ctypeChooser.grid(row =0, column = 1, sticky = Tk.W + Tk.E)


        # the Radiobox Control to choose the moment
        self.momList = ['avg velocity', 'avg momentum']
        self.pvar = Tk.IntVar()
        self.pvar.set(self.parent.GetPlotParam('m_type'))

        ttk.Label(frm, text='Moment Type:').grid(row = 1, sticky = Tk.W)
        ttk.Label(frm, text='Prtls:').grid(row = 1, column = 1, sticky = Tk.W)
        ttk.Label(frm, text='Dims:').grid(row = 1, column = 2, sticky = Tk.W)

        for i in range(len(self.momList)):
            ttk.Radiobutton(frm,
                text=self.momList[i],
                variable=self.pvar,
                command = self.RadioMom,
                value=i).grid(row = 2+i, sticky =Tk.W)
        # the Radiobox Control to choose the moment
        self.refList = ['downstream', 'upstream']
        self.framevar = Tk.IntVar()
        self.framevar.set(self.parent.GetPlotParam('UpstreamFrame'))

        ttk.Label(frm, text='Reference Frame:').grid(row = 5, sticky = Tk.W)

        for i in range(len(self.refList)):
            ttk.Radiobutton(frm,
                text=self.refList[i],
                variable=self.framevar,
                command = self.RadioRefFrame,
                value=i).grid(row = 6+i, sticky =Tk.W)


        self.ShowIonsVar = Tk.IntVar(self) # Create a var to track whether or not to show electrons
        self.ShowIonsVar.set(self.parent.GetPlotParam('show_ions'))
        ttk.Checkbutton(frm, text = "Ions",
            variable = self.ShowIonsVar,
            command = self.Selector).grid(row = 2, column = 1, sticky = Tk.W)

        self.ShowElectronsVar = Tk.IntVar(self) # Create a var to track whether or not to show electrons
        self.ShowElectronsVar.set(self.parent.GetPlotParam('show_electrons'))
        ttk.Checkbutton(frm, text = "Electrons",
            variable = self.ShowElectronsVar,
            command = self.Selector).grid(row = 3, column = 1, sticky = Tk.W)


        self.ShowXVar = Tk.IntVar(self) # Create a var to track whether or not to show electrons
        self.ShowXVar.set(self.parent.GetPlotParam('show_x'))
        ttk.Checkbutton(frm, text = "Show x",
            variable = self.ShowXVar,
            command = self.Selector).grid(row = 2, column = 2, sticky = Tk.W)
        self.ShowYVar = Tk.IntVar(self) # Create a var to track whether or not to show electrons
        self.ShowYVar.set(self.parent.GetPlotParam('show_y'))
        ttk.Checkbutton(frm, text = "Show y",
            variable = self.ShowYVar,
            command = self.Selector).grid(row = 3, column = 2, sticky = Tk.W)
        self.ShowZVar = Tk.IntVar(self) # Create a var to track whether or not to show electrons
        self.ShowZVar.set(self.parent.GetPlotParam('show_z'))
        ttk.Checkbutton(frm, text = "Show z",
            variable = self.ShowZVar,
            command = self.Selector).grid(row = 4, column = 2, sticky = Tk.W)

        # Control if the plot is weightedd
        self.WeightVar = Tk.IntVar()
        self.WeightVar.set(self.parent.GetPlotParam('weighted'))
        cb = ttk.Checkbutton(frm, text = "Weight by charge",
                        variable = self.WeightVar,
                        command = lambda:
                        self.parent.SetPlotParam('weighted', self.WeightVar.get()))
        cb.grid(row = 2, column = 3, columnspan=2, sticky = Tk.W + Tk.E)

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


        cb = ttk.Checkbutton(frm, text ='Set y min',
                        variable = self.setVminVar)
        cb.grid(row = 5, column = 3, sticky = Tk.W)
        self.VminEnter = ttk.Entry(frm, textvariable=self.Vmin, width=7)
        self.VminEnter.grid(row = 5, column = 4)

        cb = ttk.Checkbutton(frm, text ='Set y max',
                        variable = self.setVmaxVar)
        cb.grid(row = 6, column = 3, sticky = Tk.W)

        self.VmaxEnter = ttk.Entry(frm, textvariable=self.Vmax, width=7)
        self.VmaxEnter.grid(row = 6, column = 4)

        self.TrueVar = Tk.IntVar()
        self.TrueVar.set(1)
        self.xBins = Tk.StringVar()
        self.xBins.set(str(self.parent.GetPlotParam('xbins')))
        ttk.Label(frm, text ='# of xbins').grid(row = 8, column = 0, sticky = Tk.W)
        ttk.Entry(frm, textvariable=self.xBins, width=8).grid(row = 8, column = 1)





    def ctypeChanged(self, *args):
        if self.ctypevar.get() == self.parent.chartType:
            pass
        else:
            self.parent.ChangePlotType(self.ctypevar.get())
            self.destroy()

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



    def RadioMom(self):
        if self.pvar.get() == self.parent.GetPlotParam('m_type'):
            pass
        else:
            self.parent.axes.set_ylabel(self.parent.ylabel_list[self.pvar.get()][self.parent.GetPlotParam('UpstreamFrame')], labelpad = self.parent.parent.MainParamDict['yLabelPad'], color = 'black')
            self.parent.SetPlotParam('m_type', self.pvar.get())
    def RadioRefFrame(self):
        if self.framevar.get() == self.parent.GetPlotParam('UpstreamFrame'):
            pass
        else:
            self.parent.axes.set_ylabel(self.parent.ylabel_list[self.parent.GetPlotParam('m_type')][self.framevar.get()], labelpad = self.parent.parent.MainParamDict['yLabelPad'], color = 'black')
            self.parent.SetPlotParam('UpstreamFrame', self.framevar.get())

    def Selector(self):
        # Repeat the lists to help me remember
        # self.parent.key_list = ['show_x', 'show_y', 'show_z']
        # self.parent.ion_plot_list = [self.ix_plot[0], self.iy_plot[0], self.iz_plot[0]]
        # self.parent.e_plot_list = [self.ex_plot[0], self.ey_plot[0], self.ez_plot[0]]
        # self.label_prefix = [r'$\langle\beta', r'$\langle p']
        # self.label_mid = [r'_{i',r'_{e']
        # self.label_suffix = [r'x}\rangle$',r'y\rangle}$',r'z}\rangle$']


        VarList = [self.ShowXVar, self.ShowYVar,  self.ShowZVar]

        legend_handles = []
        legend_labels = []

        for i in range(len(self.parent.ion_plot_list)):
            # Set the visibility to the new value
            self.parent.ion_plot_list[i].set_visible(VarList[i].get() and self.ShowIonsVar.get())
            if VarList[i].get() and self.ShowIonsVar.get():
                legend_handles.append(self.parent.ion_plot_list[i])
                legend_labels.append(self.parent.label_prefix[self.parent.GetPlotParam('m_type')]+self.parent.label_mid[0]+self.parent.label_suffix[i])

        for i in range(len(self.parent.e_plot_list)):
            # Set the visibility to the new value
            self.parent.e_plot_list[i].set_visible(VarList[i].get() and self.ShowElectronsVar.get())
            if VarList[i].get() and self.ShowElectronsVar.get():
                legend_handles.append(self.parent.e_plot_list[i])
                legend_labels.append(self.parent.label_prefix[self.parent.GetPlotParam('m_type')]+self.parent.label_mid[1]+self.parent.label_suffix[i])


        if self.ShowElectronsVar.get() and self.ShowIonsVar.get():
            self.parent.legend = self.parent.axes.legend(legend_handles, legend_labels,
            framealpha = .05, fontsize = 11, loc = 'best', ncol = 2)
        else:
            self.parent.legend = self.parent.axes.legend(legend_handles, legend_labels,
            framealpha = .05, fontsize = 11, loc = 'best')

        self.parent.legend.set_visible(self.parent.GetPlotParam('show_legend'))
        self.parent.legend.get_frame().set_facecolor('k')
        self.parent.legend.get_frame().set_linewidth(0.0)

        # Force a plot refresh
        for i in range(len(VarList)):
            self.parent.SetPlotParam(self.parent.key_list[i], VarList[i].get(), update_plot = False)
        self.parent.SetPlotParam('show_ions', self.ShowIonsVar.get(), update_plot = False)
        self.parent.SetPlotParam('show_electrons', self.ShowElectronsVar.get(), update_plot = False)
        self.parent.SetPlotParam(self.parent.key_list[0], VarList[0].get(), update_plot = True)


    def TxtEnter(self, e):
        self.FieldsCallback()

    def FieldsCallback(self):
        tkvarLimList = [self.Vmin, self.Vmax, self.xBins]
        plot_param_List = ['v_min', 'v_max', 'xbins']
        tkvarSetList = [self.setVminVar, self.setVmaxVar, self.TrueVar]
        to_reload = False
        for j in range(len(tkvarLimList)):
            try:
            #make sure the user types in a int
                if np.abs(float(tkvarLimList[j].get()) - self.parent.GetPlotParam(plot_param_List[j])) > 1E-4:
                    if j != len(tkvarLimList)-1:
                        self.parent.SetPlotParam(plot_param_List[j], float(tkvarLimList[j].get()), update_plot = False)
                        to_reload += True*tkvarSetList[j].get()
                    else:
                        self.parent.SetPlotParam(plot_param_List[j], int(tkvarLimList[j].get()), update_plot = False)
                        to_reload += True*tkvarSetList[j].get()

            except ValueError:
                #if they type in random stuff, just set it ot the param value
                tkvarLimList[j].set(str(self.parent.GetPlotParam(plot_param_List[j])))
        if to_reload:
            self.parent.SetPlotParam('v_min', self.parent.GetPlotParam('v_min'))


    def OnClosing(self):
        self.parent.settings_window = None
        self.destroy()
