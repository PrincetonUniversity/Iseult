#!/usr/bin/env python
import matplotlib, sys
sys.path.append('../')

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
                       'legend_loc': 'N/A',
                       'face_color': 'gainsboro'
                       } # legend_loc is a string that stores the
                         # location of the legend in figure pixels.
                         # Unfortunately it is not always up to date.
                         # N/A plots it at the 'best' location.

    # We need the types of all the parameters for the config file

    def __init__(self, parent, pos, param_dict):
        self.param_dict = {}
        for key, val in self.plot_param_dict.items():
            self.param_dict[key] = val
        for key, val in param_dict.items():
            self.param_dict[key] = val
        self.pos = pos
        self.parent = parent
        self.chartType = 'Moments'
        self.figure = self.parent.figure
        #self.ylabel_list = [[r'$\langle \beta \rangle$',r'$\langle \beta \rangle$' + ' upstream frame'],[r'$\langle m\gamma\beta\rangle/m_i$', r'$\langle m\gamma\beta \rangle/m_i$' + ' upstream frame']]
        self.ylabel_list = [r'$\langle \beta \rangle$',r'$\langle m\gamma\beta\rangle/m_i$', r'$\langle KE \rangle/m_ic^2$']
        self.legend_labels = [[[r'$\langle\beta_{ix}\rangle$',r'$\langle\beta_{iy}\rangle$',r'$\langle\beta_{iz}\rangle$'],
                               [r'$\langle\beta_{ex}\rangle$',r'$\langle\beta_{ey}\rangle$',r'$\langle\beta_{ez}\rangle$'],
                               [r'$\langle\beta_x\rangle$',r'$\langle\beta_y\rangle$',r'$\langle\beta_z\rangle$']],
                              [[r'$\langle p_{ix}\rangle$',r'$\langle p_{iy}\rangle$',r'$\langle p_{iz}\rangle$'],
                               [r'$\langle p_{ex}\rangle$',r'$\langle p_{ey}\rangle$',r'$\langle p_{ez}\rangle$'],
                               [r'$\langle p_x \rangle$',r'$\langle p_y \rangle$',r'$\langle p_z \rangle$']],
                              [[r'$\Delta \gamma_i$',r'$\langle KE_{i}\rangle$', None],
                               [r'$\Delta \gamma_e$',r'$\langle KE_{e}\rangle$', None],
                               [r'$\Delta \gamma$',r'$\langle KE \rangle$',None]]]








    def draw(self, output):

        ''' A function that draws the data. In the interest in speeding up the
        code, draw should only be called when you want to recreate the whole
        figure, i.e. it  will be slow. Most times you will only want to update
        what has changed in the figure. This will be done in a function called
        refresh, that should be much much faster.'''

        self.c_omp = getattr(output, 'c_omp')
        self.istep = getattr(output, 'istep')
        self.memi = getattr(output, 'me')/getattr(output, 'mi')
        self.totalcolor = new_cmaps.cmaps[self.parent.MainParamDict['ColorMap']](0.0)

        # x-values haven't been calculated yet, generate them then save them to the dictionary for later.
        self.xaxis_values = np.arange(getattr(output, 'bx').shape[-1])/self.c_omp*self.istep
        self.fxmin = 0
        self.fxmax = self.xaxis_values[-1]

        xbn = self.GetPlotParam('xbins')
        self.xe = np.copy(getattr(output, 'xe'))
        self.ue = getattr(output, 'ue')
        self.ve = getattr(output, 've')
        self.we = getattr(output, 'we')
        self.ex = np.zeros(xbn)
        self.ey = np.zeros(xbn)
        self.ez = np.zeros(xbn)
        self.ecounts = np.zeros(xbn)

        self.xi = np.copy(getattr(output, 'xi'))
        self.ui = getattr(output, 'ui')
        self.vi = getattr(output, 'vi')
        self.wi = getattr(output, 'wi')
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
                eweights = getattr(output, 'che')
                CalcVWeightedHists(self.xe,self.ue, self.ve, self.we, ge,eweights, bin_width, self.xmin, self.ex, self.ey, self.ez, self.ecounts)

                iweights = getattr(output, 'chi')
                CalcVWeightedHists(self.xi,self.ui, self.vi, self.wi, gi, iweights, bin_width, self.xmin, self.ix, self.iy, self.iz, self.icounts)

        if self.GetPlotParam('m_type') == 1:
            # Calculate px, py, pz
            if not self.GetPlotParam('weighted'):
                CalcPHists(self.xe,self.ue, self.ve, self.we, bin_width, self.xmin, self.ex, self.ey, self.ez, self.ecounts)
                CalcPHists(self.xi,self.ui, self.vi, self.wi, bin_width, self.xmin, self.ix, self.iy, self.iz, self.icounts)
            else:
                eweights = getattr(output, 'che')
                CalcPWeightedHists(self.xe,self.ue, self.ve, self.we, eweights, bin_width, self.xmin, self.ex, self.ey, self.ez, self.ecounts)
                iweights = getattr(output, 'chi')
                CalcPWeightedHists(self.xi,self.ui, self.vi, self.wi, iweights, bin_width, self.xmin, self.ix, self.iy, self.iz, self.icounts)

            self.ex *= self.memi
            self.ey *= self.memi
            self.ez *= self.memi

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
                eweights = getattr(output, 'che')
                iweights = getattr(output, 'chi')
                # We'll put the Temp histograms into ex and ix, and Energy into ey and iy
                CalcVxEWeightedHists(self.xe,self.ue, ge, eweights, bin_width, self.xmin, vex, self.ey, self.ecounts)
                CalcVxEWeightedHists(self.xi,self.ui, gi, iweights, bin_width, self.xmin, vix, self.iy, self.icounts)

                CalcDelGamWeightedHists(self.xe, self.ue, self.ve, self.we, ge, eweights, vex,  1/np.sqrt(1-vex**2), bin_width, self.xmin, self.ecounts, self.ex)
                CalcDelGamWeightedHists(self.xi, self.ui, self.vi, self.wi, gi, iweights, vix,  1/np.sqrt(1-vix**2), bin_width, self.xmin, self.icounts, self.ix)

            self.ex*=self.memi
            self.ey*=self.memi

        # Set the tick color
        tick_color = 'black'

        # Create a gridspec to handle spacing better
        self.gs = gridspec.GridSpecFromSubplotSpec(100,100, subplot_spec = self.parent.gs0[self.pos])

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
            self.axes.set_axis_bgcolor(self.GetPlotParam('face_color'))
        else:
            self.axes.set_facecolor(self.GetPlotParam('face_color'))

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
        framealpha = .05, fontsize = self.parent.MainParamDict['legendLabelSize'], loc = 1, ncol = max(1, self.GetPlotParam('show_ions')+ self.GetPlotParam('show_electrons')+self.GetPlotParam('show_total')))

        self.legend.get_frame().set_facecolor('k')
        self.legend.get_frame().set_linewidth(0.0)
        if not self.GetPlotParam('show_legend'):
            self.legend.set_visible(False)

        self.legend.set_draggable(True,update = 'loc')
        if self.GetPlotParam('legend_loc') != 'N/A':
            tmp_tup = float(self.GetPlotParam('legend_loc').split()[0]),float(self.GetPlotParam('legend_loc').split()[1])
            self.legend._set_loc(tmp_tup)

        if self.GetPlotParam('set_v_min'):
            self.axes.set_ylim(bottom = self.GetPlotParam('v_min'))
        if self.GetPlotParam('set_v_max'):
            self.axes.set_ylim(top = self.GetPlotParam('v_max'))

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

    def GetPlotParam(self, keyname):
        return self.param_dict[keyname]
