
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
from NumbaMoments import stepify, CalcVxEHists, CalcVHists, CalcVxEWeightedHists, CalcVWeightedHists, CalcPHists, CalcPWeightedHists, RestFrameBoost, LorentzFactor, Total, CalcDelGamHists, CalcDelGamWeightedHists

class  MomentsPanel:
    # A dictionary of all of the parameters for this plot with the default parameters

    plot_param_dict = {'twoD': 0,
                       'm_type': 0, # 0 = average_velocity, 1 = average_momentum, 2 = Energy
                       'v_min': 0,
                       'v_max' : 10,
                       'set_v_min': False,
                       'set_v_max': False,
                       'show_x': True,
                       'show_y': False,
                       'show_z': False,
                       'show_ions': True,
                       'show_electrons': True,
                       'show_total': False,
                       'UpstreamFrame': False,
                       'weighted': False,
                       'xbins': 100,
                       'spatial_x': True,
                       'show_legend': True,
                       'spatial_y': False,
                       'logy': False,
                       'legend_loc': 'N/A'} # legend_loc is a string that stores the
                                         # location of the legend in figure pixels.
                                         # Unfortunately it is not always up to date.
                                         # N/A plots it at the 'best' location.

    # We need the types of all the parameters for the config file
    BoolList = ['twoD', 'set_v_min', 'set_v_max',
                'UpstreamFrame', 'show_ions', 'show_electrons', 'show_x',
                'show_y', 'show_z', 'spatial_x', 'spatial_y', 'weighted', 'show_legend']
    IntList = ['m_type','xbins']
    FloatList = ['v_min', 'v_max']
    StrList = ['legend_loc']

    def __init__(self, parent, figwrapper):
        self.settings_window = None
        self.FigWrap = figwrapper
        self.parent = parent

        #self.ylabel_list = [[r'$\langle \beta \rangle$',r'$\langle \beta \rangle$' + ' upstream frame'],[r'$\langle m\gamma\beta\rangle/m_i$', r'$\langle m\gamma\beta \rangle/m_i$' + ' upstream frame']]
        self.ylabel_list = [r'$\langle \beta \rangle$',r'$\langle m\gamma\beta\rangle/m_i$', r'$\langle KE \rangle/m_ic^2$']
        self.ChartTypes = self.FigWrap.PlotTypeDict.keys()
        self.chartType = self.FigWrap.chartType
        self.figure = self.FigWrap.figure
        self.legend_labels = [[[r'$\langle\beta_{ix}\rangle$',r'$\langle\beta_{iy}\rangle$',r'$\langle\beta_{iz}\rangle$'],
                               [r'$\langle\beta_{ex}\rangle$',r'$\langle\beta_{ey}\rangle$',r'$\langle\beta_{ez}\rangle$'],
                               [r'$\langle\beta_x\rangle$',r'$\langle\beta_y\rangle$',r'$\langle\beta_z\rangle$']],
                              [[r'$\langle p_{ix}\rangle$',r'$\langle p_{iy}\rangle$',r'$\langle p_{iz}\rangle$'],
                               [r'$\langle p_{ex}\rangle$',r'$\langle p_{ey}\rangle$',r'$\langle p_{ez}\rangle$'],
                               [r'$\langle p_x \rangle$',r'$\langle p_y \rangle$',r'$\langle p_z \rangle$']],
                              [[r'$\Delta \gamma_i$',r'$\langle KE_{i}\rangle$', None],
                               [r'$\Delta \gamma_e$',r'$\langle KE_{e}\rangle$', None],
                               [r'$\Delta \gamma$',r'$\langle KE \rangle$',None]]]

    def ChangePlotType(self, str_arg):
        self.FigWrap.ChangeGraph(str_arg)

    def set_plot_keys(self):

        '''A helper function that will insure that each hdf5 file will only be
        opened once per time step'''
        # To make things easier, we are just going to load and calculate all the moments
        self.key_name = str(self.GetPlotParam('xbins'))
        if self.GetPlotParam('weighted'):
            self.key_name += 'W'
        if self.GetPlotParam('m_type') == 0:
            self.key_name += 'vel'
        elif self.GetPlotParam('m_type') == 1:
            self.key_name += 'mom'
        elif self.GetPlotParam('m_type') == 2:
            self.key_name += 'Eng'
        self.key_name += str(self.parent.MainParamDict['PrtlStride'])

        #if self.key_name+'x_bins' in self.parent.DataDict.keys():
        #    print 'hi'
        #    self.preloaded = True
        #    self.arrs_needed = ['c_omp', 'bx', 'istep', 'me', 'mi']
        if self.GetPlotParam('weighted'):
            self.preloaded = False
            self.arrs_needed = ['c_omp', 'bx', 'istep', 'me', 'mi', 'gamma0','xi', 'ui', 'vi', 'wi', 'xe','ue', 've', 'we','chi', 'che']
        else:
            self.preloaded = False
            self.arrs_needed = ['c_omp', 'bx', 'istep', 'me', 'mi', 'gamma0','xi', 'ui', 'vi', 'wi', 'xe','ue', 've', 'we']
        return self.arrs_needed

    def LoadData(self):
        ''' A helper function that checks if the histogram has
        already been calculated and if it hasn't, it calculates
        it then stores it.'''

        self.c_omp = self.FigWrap.LoadKey('c_omp')[0]
        self.istep = self.FigWrap.LoadKey('istep')[0]
        self.memi = self.FigWrap.LoadKey('me')[0]/self.FigWrap.LoadKey('mi')[0]
        self.totalcolor = new_cmaps.cmaps[self.parent.MainParamDict['ColorMap']](0.0)

        if 'xaxis_values' in self.parent.DataDict.keys():
            self.xaxis_values = self.parent.DataDict['xaxis_values']
        else:
            # x-values haven't been calculated yet, generate them then save them to the dictionary for later.
            self.xaxis_values = np.arange(self.FigWrap.LoadKey('bx').shape[-1])/self.c_omp*self.istep
        #            print self.xaxis_values
            self.parent.DataDict['xaxis_values'] = np.copy(self.xaxis_values)
        self.fxmin = 0
        self.fxmax = self.xaxis_values[-1]

        if self.key_name+'x_bins' in self.parent.DataDict.keys():
            self.x_bins = self.parent.DataDict[self.key_name+'x_bins']
            self.ex = self.parent.DataDict[self.key_name+'ex']
            self.ey = self.parent.DataDict[self.key_name+'ey']
            self.ez = self.parent.DataDict[self.key_name+'ez']
            self.ecounts = self.parent.DataDict[self.key_name+'ecounts']

            self.ix = self.parent.DataDict[self.key_name+'ix']
            self.iy = self.parent.DataDict[self.key_name+'iy']
            self.iz = self.parent.DataDict[self.key_name+'iz']
            self.icounts = self.parent.DataDict[self.key_name+'icounts']

        else: #We have to calculate everything
            xbn = self.GetPlotParam('xbins')
            self.xe = np.copy(self.FigWrap.LoadKey('xe'))
            self.ue = self.FigWrap.LoadKey('ue')
            self.ve = self.FigWrap.LoadKey('ve')
            self.we = self.FigWrap.LoadKey('we')
            self.ex = np.zeros(xbn)
            self.ey = np.zeros(xbn)
            self.ez = np.zeros(xbn)
            self.ecounts = np.zeros(xbn)

            self.xi = np.copy(self.FigWrap.LoadKey('xi'))
            self.ui = self.FigWrap.LoadKey('ui')
            self.vi = self.FigWrap.LoadKey('vi')
            self.wi = self.FigWrap.LoadKey('wi')
            self.ix = np.zeros(xbn)
            self.iy = np.zeros(xbn)
            self.iz = np.zeros(xbn)
            self.icounts = np.zeros(xbn)

            if self.memi==0:
                self.xmin = np.min(self.xe)/self.c_omp
                self.xmax = np.max(self.xe)/self.c_omp
                self.memi = 1.0
                self.SetPlotParam('show_ions', False, update_plot = False)
                self.ylabel_list = [r'$\langle \beta \rangle$',r'$\langle \gamma\beta\rangle$', r'$\langle KE \rangle/m_ec^2$']
            else:
                self.xmin = min(np.min(self.xi), np.min(self.xe))/self.c_omp
                self.xmax = max(np.max(self.xi), np.max(self.xe))/self.c_omp
            bin_width = (self.xmax-self.xmin)/float(xbn)*self.c_omp
            self.x_bins = np.linspace(self.xmin, self.xmax, num = xbn+1)
            self.parent.DataDict[self.key_name+'x_bins'] = self.x_bins

            if self.GetPlotParam('m_type') == 0:
                # Calculate vx, vy, vz
                # First calculate gamma
                gi = np.empty(len(self.ui))
                LorentzFactor(self.ui, self.vi, self.wi, gi)
                ge = np.empty(len(self.ue))
                LorentzFactor(self.ue, self.ve, self.we, ge)
                if not self.GetPlotParam('weighted'):
                    CalcVHists(self.xe,self.ue, self.ve, self.we, ge, bin_width, self.xmin, self.ex, self.ey, self.ez, self.ecounts)
                    CalcVHists(self.xi,self.ui, self.vi, self.wi, gi, bin_width, self.xmin, self.ix, self.iy, self.iz, self.icounts)
                else:
                    eweights = self.FigWrap.LoadKey('che')
                    CalcVWeightedHists(self.xe,self.ue, self.ve, self.we, ge,eweights, bin_width, self.xmin, self.ex, self.ey, self.ez, self.ecounts)

                    iweights = self.FigWrap.LoadKey('chi')
                    CalcVWeightedHists(self.xi,self.ui, self.vi, self.wi, gi, iweights, bin_width, self.xmin, self.ix, self.iy, self.iz, self.icounts)

                self.parent.DataDict[self.key_name+'ex'] = self.ex
                self.parent.DataDict[self.key_name+'ey'] = self.ey
                self.parent.DataDict[self.key_name+'ez'] = self.ez
                self.parent.DataDict[self.key_name+'ecounts'] = self.ecounts

                self.parent.DataDict[self.key_name+'ix'] = self.ix
                self.parent.DataDict[self.key_name+'iy'] = self.iy
                self.parent.DataDict[self.key_name+'iz'] = self.iz
                self.parent.DataDict[self.key_name+'icounts'] = self.icounts

            if self.GetPlotParam('m_type') == 1:
                # Calculate px, py, pz
                if not self.GetPlotParam('weighted'):
                    CalcPHists(self.xe,self.ue, self.ve, self.we, bin_width, self.xmin, self.ex, self.ey, self.ez, self.ecounts)
                    CalcPHists(self.xi,self.ui, self.vi, self.wi, bin_width, self.xmin, self.ix, self.iy, self.iz, self.icounts)
                else:
                    eweights = self.FigWrap.LoadKey('che')
                    CalcPWeightedHists(self.xe,self.ue, self.ve, self.we, eweights, bin_width, self.xmin, self.ex, self.ey, self.ez, self.ecounts)
                    iweights = self.FigWrap.LoadKey('chi')
                    CalcPWeightedHists(self.xi,self.ui, self.vi, self.wi, iweights, bin_width, self.xmin, self.ix, self.iy, self.iz, self.icounts)

                self.ex *= self.memi
                self.ey *= self.memi
                self.ez *= self.memi
                self.parent.DataDict[self.key_name+'ex'] = self.ex
                self.parent.DataDict[self.key_name+'ey'] = self.ey
                self.parent.DataDict[self.key_name+'ez'] = self.ez
                self.parent.DataDict[self.key_name+'ecounts'] = self.ecounts

                self.parent.DataDict[self.key_name+'ix'] = self.ix
                self.parent.DataDict[self.key_name+'iy'] = self.iy
                self.parent.DataDict[self.key_name+'iz'] = self.iz
                self.parent.DataDict[self.key_name+'icounts'] = self.icounts

            if self.GetPlotParam('m_type') == 2:
                boost_gam = np.zeros(xbn)
                vx_avg = np.zeros(xbn)
                vix = np.zeros(xbn)
                vex = np.zeros(xbn)
                gi = np.empty(len(self.ui))
                LorentzFactor(self.ui, self.vi, self.wi, gi)
                ge = np.empty(len(self.ue))
                LorentzFactor(self.ue, self.ve, self.we, ge)

                if not self.GetPlotParam('weighted'):
                    # We'll put the Energy histograms into ex and ix, and delE intos ey and iy
                    CalcVxEHists(self.xe,self.ue, ge, bin_width, self.xmin,  vex, self.ey, self.ecounts)
                    CalcVxEHists(self.xi,self.ui, gi, bin_width, self.xmin, vix, self.iy, self.icounts)

                    CalcDelGamHists(self.xe, self.ue, self.ve, self.we, ge, vix, 1/np.sqrt(1-vix**2), bin_width, self.xmin, self.ecounts, self.ex)
                    CalcDelGamHists(self.xi, self.ui, self.vi, self.wi, gi, vex, 1/np.sqrt(1-vex**2), bin_width, self.xmin, self.icounts, self.ix)


                else:
                    eweights = self.FigWrap.LoadKey('che')
                    iweights = self.FigWrap.LoadKey('chi')
                    # We'll put the Temp histograms into ex and ix, and Energy into ey and iy
                    CalcVxEWeightedHists(self.xe,self.ue, ge, eweights, bin_width, self.xmin, vex, self.ey, self.ecounts)
                    CalcVxEWeightedHists(self.xi,self.ui, gi, iweights, bin_width, self.xmin, vix, self.iy, self.icounts)

                    RestFrameBoost(vex, self.ecounts, vix,  self.icounts, vx_avg, boost_gam)
                    CalcDelGamWeightedHists(self.xe, self.ue, self.ve, self.we, ge, eweights, vix,  1/np.sqrt(1-vix**2), bin_width, self.xmin, self.ecounts, self.ex)
                    CalcDelGamWeightedHists(self.xi, self.ui, self.vi, self.wi, gi, iweights, vex,  1/np.sqrt(1-vex**2), bin_width, self.xmin, self.icounts, self.ix)

                self.ex*=self.memi
                self.ey*=self.memi
                self.parent.DataDict[self.key_name+'ex'] = self.ex
                self.parent.DataDict[self.key_name+'ecounts'] = self.ecounts
                self.parent.DataDict[self.key_name+'ix'] = self.ix
                self.parent.DataDict[self.key_name+'icounts'] = self.icounts

                self.parent.DataDict[self.key_name+'ey'] = self.ey
                self.parent.DataDict[self.key_name+'iy'] = self.iy
                self.parent.DataDict[self.key_name+'ez'] = None
                self.parent.DataDict[self.key_name+'iz'] = None





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
                                     ls= '-', #marker = 's', markeredgecolor  = self.parent.ion_color,
                                     color = self.parent.ion_color,
                                     visible =  self.GetPlotParam('show_ions')*self.GetPlotParam('show_z'))

        self.iz_plot[0].set_dashes([5, 1])

        self.tx_plot = self.axes.plot(0,0,
                                       ls= '-', #marker = '<', markeredgecolor =  self.parent.ion_color,
                                       color = self.totalcolor,
                                       visible =  self.GetPlotParam('show_total')*self.GetPlotParam('show_x'))
        self.ty_plot = self.axes.plot(0,0,
                                       ls= '-',# marker = '*',  markersize = 10, markeredgecolor =  self.parent.ion_color,
                                       color = self.totalcolor,
                                       visible =  self.GetPlotParam('show_total')*self.GetPlotParam('show_y'))
        self.ty_plot[0].set_dashes([1, 1])

        self.tz_plot = self.axes.plot(0,0,
                                     ls= '-', #marker = 's', markeredgecolor  = self.parent.ion_color,
                                     color = self.totalcolor,
                                     visible =  self.GetPlotParam('show_total')*self.GetPlotParam('show_z'))

        self.tz_plot[0].set_dashes([5, 1])
        if self.GetPlotParam('show_electrons'):
            if self.GetPlotParam('show_x'):
                self.ex_plot[0].set_data(*stepify(self.x_bins, self.ex))
            if self.GetPlotParam('show_y'):
                self.ey_plot[0].set_data(*stepify(self.x_bins, self.ey))
            if self.GetPlotParam('show_z'):
                self.ez_plot[0].set_data(*stepify(self.x_bins, self.ez))
        if self.GetPlotParam('show_ions'):
            if self.GetPlotParam('show_x'):
                self.ix_plot[0].set_data(*stepify(self.x_bins, self.ix))
            if self.GetPlotParam('show_y'):
                self.iy_plot[0].set_data(*stepify(self.x_bins, self.iy))
            if self.GetPlotParam('show_z'):
                self.iz_plot[0].set_data(*stepify(self.x_bins, self.iz))
        if self.GetPlotParam('show_total'):
            if self.GetPlotParam('show_x'):
                self.tx_plot[0].set_data(*stepify(self.x_bins, Total(self.ix, self.icounts, self.ex, self.ecounts)))
            if self.GetPlotParam('show_y'):
                self.ty_plot[0].set_data(*stepify(self.x_bins, Total(self.iy, self.icounts, self.ey, self.ecounts)))
            if self.GetPlotParam('show_z') and self.GetPlotParam('m_type')!=2:
                self.tz_plot[0].set_data(*stepify(self.x_bins, Total(self.iz, self.icounts, self.ez, self.ecounts)))

        if int(matplotlib.__version__[0]) < 2:
            self.axes.set_axis_bgcolor('lightgrey')
        else:
            self.axes.set_facecolor('lightgrey')

        self.axes.tick_params(labelsize = self.parent.MainParamDict['NumFontSize'], color=tick_color)



        # fancy code to make sure that matplotlib sets its limits
        # only based on visible lines
        self.key_list = ['show_x', 'show_y', 'show_z']
        self.ion_plot_list = [self.ix_plot[0], self.iy_plot[0], self.iz_plot[0]]
        self.e_plot_list = [self.ex_plot[0], self.ey_plot[0], self.ez_plot[0]]
        self.t_plot_list = [self.tx_plot[0], self.ty_plot[0], self.tz_plot[0]]



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
                    legend_labels.append(self.legend_labels[self.GetPlotParam('m_type')][0][i])


        if self.GetPlotParam('show_electrons'):
            for i in range(len(self.e_plot_list)):
                if self.GetPlotParam(self.key_list[i]):
                    xy = np.vstack(self.e_plot_list[i].get_data()).T
                    self.axes.dataLim.update_from_data_xy(xy, ignore=False)
                    legend_handles.append(self.e_plot_list[i])
                    legend_labels.append(self.legend_labels[self.GetPlotParam('m_type')][1][i])

        if self.GetPlotParam('show_total'):
            for i in range(len(self.t_plot_list)):
                if self.GetPlotParam(self.key_list[i]):
                    xy = np.vstack(self.t_plot_list[i].get_data()).T
                    self.axes.dataLim.update_from_data_xy(xy, ignore=False)
                    legend_handles.append(self.t_plot_list[i])
                    legend_labels.append(self.legend_labels[self.GetPlotParam('m_type')][2][i])

        if self.GetPlotParam('logy'):
            self.axes.set_yscale('log')
        self.axes.autoscale()

        # now make the legend


        self.legend = self.axes.legend(legend_handles, legend_labels,
        framealpha = .05, fontsize = 11, loc = 1, ncol = max(1, self.GetPlotParam('show_ions')+ self.GetPlotParam('show_electrons')+self.GetPlotParam('show_total')))

        self.legend.get_frame().set_facecolor('k')
        self.legend.get_frame().set_linewidth(0.0)
        if not self.GetPlotParam('show_legend'):
            self.legend.set_visible(False)

        self.legend.draggable(update = 'loc')
        if self.GetPlotParam('legend_loc') != 'N/A':
            tmp_tup = float(self.GetPlotParam('legend_loc').split()[0]),float(self.GetPlotParam('legend_loc').split()[1])
            self.legend._set_loc(tmp_tup)

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
            self.axes.set_xlim(self.fxmin,self.fxmax)



        self.axes.set_xlabel(r'$x\  [c/\omega_{pe}]$', labelpad = self.parent.MainParamDict['xLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
        self.axes.set_ylabel(self.ylabel_list[self.GetPlotParam('m_type')], labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])

    def refresh(self):

        '''This is a function that will be called only if self.axes already
        holds a total energy type plot. We only update things that have changed & are
        shown.  If hasn't changed or isn't shown, don't touch it. The difference
        between this and last time, is that we won't actually do any drawing in
        the plot. The plot will be redrawn after all subplots are refreshed. '''

        if self.GetPlotParam('show_electrons'):
            if self.GetPlotParam('show_x'):
                self.ex_plot[0].set_data(*stepify(self.x_bins, self.ex))
            if self.GetPlotParam('show_y'):
                self.ey_plot[0].set_data(*stepify(self.x_bins, self.ey))
            if self.GetPlotParam('show_z'):
                self.ez_plot[0].set_data(*stepify(self.x_bins, self.ez))
        if self.GetPlotParam('show_ions'):
            if self.GetPlotParam('show_x'):
                self.ix_plot[0].set_data(*stepify(self.x_bins, self.ix))
            if self.GetPlotParam('show_y'):
                self.iy_plot[0].set_data(*stepify(self.x_bins, self.iy))
            if self.GetPlotParam('show_z'):
                self.iz_plot[0].set_data(*stepify(self.x_bins, self.iz))
        if self.GetPlotParam('show_total'):
            if self.GetPlotParam('show_x'):
                self.tx_plot[0].set_data(*stepify(self.x_bins, Total(self.ix, self.icounts, self.ex, self.ecounts)))
            if self.GetPlotParam('show_y'):
                self.ty_plot[0].set_data(*stepify(self.x_bins, Total(self.iy, self.icounts, self.ey, self.ecounts)))
            if self.GetPlotParam('show_z'):
                self.tz_plot[0].set_data(*stepify(self.x_bins, Total(self.iz, self.icounts, self.ez, self.ecounts)))

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
        if self.GetPlotParam('show_total'):
            for i in range(len(self.t_plot_list)):
                if self.GetPlotParam(self.key_list[i]):
                    xy = np.vstack(self.t_plot_list[i].get_data()).T
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
            self.axes.set_xlim(self.fxmin,self.fxmax)


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
        self.momList = ['avg velocity', 'avg momentum', 'Kinetic Energy']
        self.pvar = Tk.IntVar()
        self.pvar.set(self.parent.GetPlotParam('m_type'))

        ttk.Label(frm, text='Moment Type:').grid(row = 1, sticky = Tk.W)
        ttk.Label(frm, text='Prtls:').grid(row = 1, column = 1, sticky = Tk.W)
        self.Dlabel = ttk.Label(frm, text='Options:')
        self.Dlabel.grid(row = 1, column = 2, sticky = Tk.W)

        for i in range(len(self.momList)):
            ttk.Radiobutton(frm,
                text=self.momList[i],
                variable=self.pvar,
                command = self.RadioMom,
                value=i).grid(row = 2+i, sticky =Tk.W)
        # the Radiobox Control to choose the moment
        #self.refList = ['downstream', 'upstream']
        #self.framevar = Tk.IntVar()
        #self.framevar.set(self.parent.GetPlotParam('UpstreamFrame'))

        #ttk.Label(frm, text='Reference Frame:').grid(row = 5, sticky = Tk.W)

        #for i in range(len(self.refList)):
        #    ttk.Radiobutton(frm,
        #        text=self.refList[i],
        #        variable=self.framevar,
        #        command = self.RadioRefFrame,
        #        value=i).grid(row = 6+i, sticky =Tk.W)


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

        self.ShowTotalVar = Tk.IntVar(self) # Create a var to track whether or not to show electrons
        self.ShowTotalVar.set(self.parent.GetPlotParam('show_total'))
        ttk.Checkbutton(frm, text = "Total",
            variable = self.ShowTotalVar,
            command = self.Selector).grid(row = 4, column = 1, sticky = Tk.W)


        self.ShowXVar = Tk.IntVar(self) # Create a var to track whether or not to show electrons
        self.ShowXVar.set(self.parent.GetPlotParam('show_x'))
        self.Xcb = ttk.Checkbutton(frm, text = "Show x",
            variable = self.ShowXVar,
            command = self.Selector)
        self.Xcb.grid(row = 2, column = 2, sticky = Tk.W)
        self.ShowYVar = Tk.IntVar(self) # Create a var to track whether or not to show electrons
        self.ShowYVar.set(self.parent.GetPlotParam('show_y'))
        self.Ycb = ttk.Checkbutton(frm, text = "Show y",
            variable = self.ShowYVar,
            command = self.Selector)
        self.Ycb.grid(row = 3, column = 2, sticky = Tk.W)
        self.ShowZVar = Tk.IntVar(self) # Create a var to track whether or not to show electrons
        self.ShowZVar.set(self.parent.GetPlotParam('show_z'))
        self.Zcb = ttk.Checkbutton(frm, text = "Show z",
            variable = self.ShowZVar,
            command = self.Selector)
        self.Zcb.grid(row = 4, column = 2, sticky = Tk.W)
        if self.parent.GetPlotParam('m_type') == 2:
            self.Zcb.state(['disabled'])
            self.Zcb.config(text= '')
            self.Ycb.config(text= 'KE')
            self.Xcb.config(text= 'Co-moving temp (<vx>=0)')
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


        self.LogY = Tk.IntVar()
        self.LogY.set(self.parent.GetPlotParam('logy'))
        self.LogY.trace('w', self.LogYChanged)
        cb = ttk.Checkbutton(frm, text ='Y-axis logscale',
                        variable = self.LogY)
        cb.grid(row = 6, column = 0, sticky = Tk.W)

        self.VminEnter = ttk.Entry(frm, textvariable=self.Vmin, width=7)
        self.VminEnter.grid(row = 5, column = 4)

        self.TrueVar = Tk.IntVar()
        self.TrueVar.set(1)
        self.xBins = Tk.StringVar()
        self.xBins.set(str(self.parent.GetPlotParam('xbins')))
        ttk.Label(frm, text ='# of xbins').grid(row = 7, column = 0, sticky = Tk.W)
        ttk.Entry(frm, textvariable=self.xBins, width=8).grid(row = 7, column = 1)

        self.ShowLegVar = Tk.IntVar(self)
        self.ShowLegVar.set(self.parent.GetPlotParam('show_legend'))
        self.ShowLegVar.trace('w', self.ShowLegChanged)
        cb = ttk.Checkbutton(frm, text ='Show Legend',
                        variable = self.ShowLegVar)
        cb.grid(row = 7, column = 3, sticky = Tk.W)
    def ShowLegChanged(self, *args):
        if self.ShowLegVar.get() == self.parent.GetPlotParam('show_legend'):
            pass
        else:
            self.parent.SetPlotParam('show_legend', self.ShowLegVar.get(), update_plot  = 'False')
            self.parent.legend.set_visible(self.parent.GetPlotParam('show_legend'))
            self.parent.parent.canvas.draw()
            self.parent.parent.canvas.get_tk_widget().update_idletasks()
    def LogYChanged(self, *args):
        if self.LogY.get() == self.parent.GetPlotParam('logy'):
            pass
        else:
            self.parent.SetPlotParam('logy', self.LogY.get(), update_plot  = 'False')
            if self.parent.GetPlotParam('logy'):
                self.parent.axes.set_yscale('log')
            else:
                self.parent.axes.set_yscale('linear')
            self.parent.parent.canvas.draw()
            self.parent.parent.canvas.get_tk_widget().update_idletasks()

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
            self.parent.axes.set_ylabel(self.parent.ylabel_list[self.pvar.get()], labelpad = self.parent.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.parent.MainParamDict['AxLabelSize'])
            self.parent.SetPlotParam('m_type', self.pvar.get(), update_plot=False)
            if self.parent.GetPlotParam('m_type') == 2:
                self.Zcb.state(['disabled'])
                self.Zcb.config(text= '')
                self.ShowZVar.set(0)
                self.Dlabel.config(text= 'options')
                self.Ycb.config(text= 'KE')
                self.Xcb.config(text= 'Co-moving temp (<vx>=0)')
            else:
                self.Zcb.state(['!disabled'])
                self.Zcb.config(text= 'Show Z')
                self.Dlabel.config(text= 'options')
                self.Ycb.config(text= 'Show Y')
                self.Xcb.config(text= 'Show X')

            self.Selector()

    def RadioRefFrame(self):
        if self.framevar.get() == self.parent.GetPlotParam('UpstreamFrame'):
            pass
        else:

            self.parent.axes.set_ylabel(self.parent.ylabel_list[self.parent.GetPlotParam('m_type')][self.framevar.get()], labelpad = self.parent.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.parent.MainParamDict['AxLabelSize'])
            self.parent.SetPlotParam('UpstreamFrame', self.framevar.get())

    def Selector(self):
        # Repeat the lists to help me remember
        # self.parent.key_list = ['show_x', 'show_y', 'show_z']
        # self.parent.ion_plot_list = [self.ix_plot[0], self.iy_plot[0], self.iz_plot[0]]
        # self.parent.e_plot_list = [self.ex_plot[0], self.ey_plot[0], self.ez_plot[0]]
        # self.label_prefix = [r'$\langle\beta', r'$\langle p']
        # self.label_mid = [r'_{i',r'_{e']
        # self.label_suffix = [r'x}\rangle$',r'y\rangle}$',r'z}\rangle$']

        # Update current legend position
        if self.parent.legend._get_loc() != 1:
            self.parent.SetPlotParam('legend_loc', ' '.join(str(x) for x in self.parent.legend._get_loc()), update_plot = False)

        VarList = [self.ShowXVar, self.ShowYVar,  self.ShowZVar]

        legend_handles = []
        legend_labels = []

        for i in range(len(self.parent.ion_plot_list)):
            # Set the visibility to the new value
            self.parent.ion_plot_list[i].set_visible(VarList[i].get() and self.ShowIonsVar.get())
            if VarList[i].get() and self.ShowIonsVar.get():
                legend_handles.append(self.parent.ion_plot_list[i])
                legend_labels.append(self.parent.legend_labels[self.parent.GetPlotParam('m_type')][0][i])

        for i in range(len(self.parent.e_plot_list)):
            # Set the visibility to the new value
            self.parent.e_plot_list[i].set_visible(VarList[i].get() and self.ShowElectronsVar.get())
            if VarList[i].get() and self.ShowElectronsVar.get():
                legend_handles.append(self.parent.e_plot_list[i])
                legend_labels.append(self.parent.legend_labels[self.parent.GetPlotParam('m_type')][1][i])
        for i in range(len(self.parent.t_plot_list)):
            # Set the visibility to the new value
            self.parent.t_plot_list[i].set_visible(VarList[i].get() and self.ShowTotalVar.get())
            if VarList[i].get() and self.ShowTotalVar.get():
                legend_handles.append(self.parent.t_plot_list[i])
                legend_labels.append(self.parent.legend_labels[self.parent.GetPlotParam('m_type')][2][i])



        for i in range(len(VarList)):
            self.parent.SetPlotParam(self.parent.key_list[i], VarList[i].get(), update_plot = False)
        self.parent.SetPlotParam('show_ions', self.ShowIonsVar.get(), update_plot = False)
        self.parent.SetPlotParam('show_electrons', self.ShowElectronsVar.get(), update_plot = False)
        self.parent.SetPlotParam('show_total', self.ShowTotalVar.get(), update_plot = False)



        self.parent.legend = self.parent.axes.legend(legend_handles, legend_labels,
        framealpha = .05, fontsize = 11, loc = 1, ncol = max(1, self.ShowElectronsVar.get()+self.ShowTotalVar.get()+self.ShowIonsVar.get()))

        self.parent.legend.set_visible(self.parent.GetPlotParam('show_legend'))
        self.parent.legend.get_frame().set_facecolor('k')
        self.parent.legend.get_frame().set_linewidth(0.0)
        self.parent.legend.draggable(update = 'loc')
        if self.parent.GetPlotParam('legend_loc') != 'N/A':
            tmp_tup = float(self.parent.GetPlotParam('legend_loc').split()[0]),float(self.parent.GetPlotParam('legend_loc').split()[1])
            self.parent.legend._set_loc(tmp_tup)
        # Force a plot refresh
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
