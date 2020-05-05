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

class TotEnergyPanel:
    # A dictionary of all of the parameters for this plot with the default parameters,
    # For the config file loading to work, the dictionary can only contain strings,
    # bools, ints, or floats.
    plot_param_dict = {'twoD': 0,
                       'show_prtl_KE': True,
                       'show_field_E': True,
                       'show_ion_E': False,
                       'show_electron_E': False,
                       'show_total_E': True,
                       'show_Bx_energy': False,
                       'show_By_energy': False,
                       'show_Bz_energy': False,
                       'show_Ex_energy': False,
                       'show_Ey_energy': False,
                       'show_Ez_energy': False,

                       'show_B_E': False,
                       'show_E_E': False,
                       'y_min': 0,
                       'y_max' : 10,
                       'set_y_min': False,
                       'set_y_max': False,
                       'show_legend': True,
                       'show_current_time': True,
                       'x_min': 0,
                       'x_max' : 10,
                       'set_x_min': False,
                       'set_x_max': False,
                       'yLog': True,
                       'spatial_x': False,
                       'spatial_y': False,
                       'legend_loc': 'N/A',
                       'face_color': 'gainsboro'}

    # We need the types of all the parameters for the config file

    def __init__(self, parent, pos, param_dict):
        self.param_dict = {}
        for key, val in self.plot_param_dict.items():
            self.param_dict[key] = val
        for key, val in param_dict.items():
            self.param_dict[key] = val
        self.pos = pos
        self.parent = parent
        self.chartType = 'TotalEnergyPlot'
        self.figure = self.parent.figure
    def update_data(self, output):
        self.time = output.time
    def draw(self):
        ''' A function that draws the data. In the interest in speeding up the
        code, draw should only be called when you want to recreate the whole
        figure, i.e. it  will be slow. Most times you will only want to update
        what has changed in the figure. This will be done in a function called
        refresh, that should be much much faster.'''

        # Set the tick color
        tick_color = 'black'

        # Create a gridspec to handle spacing better
        self.gs = gridspec.GridSpecFromSubplotSpec(100,100, subplot_spec = self.parent.gs0[self.pos])
        self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]])

        self.prtlcolor = new_cmaps.cmaps[self.parent.MainParamDict['ColorMap']](0.2)
        self.totalcolor = new_cmaps.cmaps[self.parent.MainParamDict['ColorMap']](0.0)
        self.fieldcolor = new_cmaps.cmaps[self.parent.MainParamDict['ColorMap']](0.8)

        self.electron_plot = self.axes.plot(self.parent.TotalEnergyTimes, self.parent.TotalElectronEnergy,
                                            ls= ':', marker = '^', markeredgecolor = self.prtlcolor,
                                            color = self.prtlcolor, visible = self.GetPlotParam('show_electron_E'))
        self.ion_plot = self.axes.plot(self.parent.TotalEnergyTimes,  self.parent.TotalIonEnergy,
                                       ls= ':', marker = 'v', markeredgecolor = self.prtlcolor,
                                       color = self.prtlcolor, visible = self.GetPlotParam('show_ion_E'))
        self.prtl_plot = self.axes.plot(self.parent.TotalEnergyTimes, self.parent.TotalElectronEnergy + self.parent.TotalIonEnergy,
                                        ls= ':', marker = 'd', markeredgecolor = self.prtlcolor,
                                        color = self.prtlcolor, visible = self.GetPlotParam('show_prtl_KE'))

        self.Bx_plot = self.axes.plot(self.parent.TotalEnergyTimes, self.parent.TotalBxEnergy,
                                       ls= ':', marker = '1', markeredgecolor = self.fieldcolor,
                                       color = self.fieldcolor, visible = self.GetPlotParam('show_Bx_energy'))
        self.By_plot = self.axes.plot(self.parent.TotalEnergyTimes, self.parent.TotalByEnergy,
                                       ls= ':', marker = '2', markeredgecolor = self.fieldcolor,
                                       color = self.fieldcolor, visible = self.GetPlotParam('show_By_energy'))
        self.Bz_plot = self.axes.plot(self.parent.TotalEnergyTimes, self.parent.TotalBzEnergy,
                                       ls= ':', marker = '3', markeredgecolor = self.fieldcolor,
                                       color = self.fieldcolor, visible = self.GetPlotParam('show_Bz_energy'))

        self.Ex_plot = self.axes.plot(self.parent.TotalEnergyTimes, self.parent.TotalExEnergy,
                                       ls= ':', marker = '4', markeredgecolor = 'green',
                                       color = 'green', visible = self.GetPlotParam('show_Ex_energy'))
        self.Ey_plot = self.axes.plot(self.parent.TotalEnergyTimes, self.parent.TotalEyEnergy,
                                       ls= ':', marker = 'CARETLEFT', markeredgecolor = 'green',
                                       color = 'green', visible = self.GetPlotParam('show_Ey_energy'))
        self.Ez_plot = self.axes.plot(self.parent.TotalEnergyTimes, self.parent.TotalEzEnergy,
                                       ls= ':', marker = 'CARETRIGHT', markeredgecolor = 'green',
                                       color = 'green', visible = self.GetPlotParam('show_Ez_energy'))

        self.mag_plot = self.axes.plot(self.parent.TotalEnergyTimes, self.parent.TotalMagEnergy,
                                       ls= ':', marker = '*',  markersize = 10, markeredgecolor = self.fieldcolor,
                                       color = self.fieldcolor, visible = self.GetPlotParam('show_B_E'))
        self.e_plot = self.axes.plot(self.parent.TotalEnergyTimes, self.parent.TotalElectricEnergy,
                                     ls= ':', marker = 's', markeredgecolor = self.fieldcolor,
                                     color = self.fieldcolor, visible = self.GetPlotParam('show_E_E'))
        self.field_plot = self.axes.plot(self.parent.TotalEnergyTimes, self.parent.TotalMagEnergy + self.parent.TotalElectricEnergy,
                                         ls= ':', marker = 'o', markeredgecolor = self.fieldcolor,
                                         color = self.fieldcolor, visible = self.GetPlotParam('show_field_E'))

        self.total_plot = self.axes.plot(self.parent.TotalEnergyTimes, self.parent.TotalMagEnergy + self.parent.TotalElectricEnergy + self.parent.TotalElectronEnergy + self.parent.TotalIonEnergy,
                                         ls= ':', marker = 'x', markeredgecolor = self.totalcolor,
                                         color = self.totalcolor, visible = self.GetPlotParam('show_total_E'))


        self.cur_time = self.axes.axvline(self.time, linewidth = 1.5, linestyle = '--',
                                          color = 'k', alpha = .4,
                                          visible = self.GetPlotParam('show_current_time'))


        if int(matplotlib.__version__[0]) < 2:
            self.axes.set_axis_bgcolor(self.GetPlotParam('face_color'))
        else:
            self.axes.set_facecolor(self.GetPlotParam('face_color'))

        self.axes.tick_params(labelsize = self.parent.MainParamDict['NumFontSize'], color=tick_color)


        if self.GetPlotParam('yLog'):
            self.axes.set_yscale('log')

        # fancy code to make sure that matplotlib sets its limits
        # only based on visible lines
        self.key_list = ['show_total_E', 'show_prtl_KE', 'show_ion_E', 'show_electron_E', 'show_field_E', 'show_E_E',
                                'show_B_E', 'show_Bx_energy', 'show_By_energy', 'show_Bz_energy',
                                'show_Ex_energy', 'show_Ey_energy', 'show_Ez_energy']
        self.plot_list = [self.total_plot[0], self.prtl_plot[0], self.ion_plot[0], self.electron_plot[0], self.field_plot[0],
                                self.e_plot[0], self.mag_plot[0],
                                self.Bx_plot[0], self.By_plot[0], self.Bz_plot[0],
                                self.Ex_plot[0], self.Ey_plot[0], self.Ez_plot[0]]
        self.label_names = ['Prtl+Field', 'Particles', 'Ions', 'Electrons', 'EM Field', 'Electric Field', 'Magnetic Field',
                            r'$B_x^2$', r'$B_y^2$', r'$B_z^2$', r'$E_x^2$', r'$E_y^2$', r'$E_z^2$']

        for i in range(len(self.plot_list)):
            line = self.plot_list[i]
            if self.GetPlotParam(self.key_list[i]):
                xmin_max[0] = min(xmin_max[0], np.min(line.get_xdata()))
                xmin_max[1] = max(xmin_max[1], np.max(line.get_xdata()))
                ymin_max[0] = min(ymin_max[0], np.min(line.get_ydata()))
                ymin_max[1] = max(ymin_max[1], np.max(line.get_ydata()))

        if np.isinf(xmin_max[0]):
            xmin_max=[None, None]
            ymin_max=[None, None]

        if self.GetPlotParam('yLog'):
            for i in range(2):
                if ymin_max[i] <= 0:
                    ymin_max[i] = None
                if ymin_max[i] is not None:
                    ymin_max[i] = np.log10(ymin_max[i])


        if ymin_max[0] is not None and ymin_max[1] is not None:
            dist = ymin_max[1]-ymin_max[0]
            ymin_max[0] -= 0.04*dist
            ymin_max[1] += 0.04*dist
            if self.GetPlotParam('yLog'):
                ymin_max = [10**elm for elm in ymin_max]

        self.axes.set_xlim(xmin_max)
        self.axes.set_ylim(ymin_max)

        if self.GetPlotParam('set_y_min'):
            self.axes.set_ylim(bottom = self.GetPlotParam('y_min'))
        if self.GetPlotParam('set_y_max'):
            self.axes.set_ylim(top = self.GetPlotParam('y_max'))

        if self.GetPlotParam('set_x_min'):
            self.axes.set_xlim(left = self.GetPlotParam('x_min'))
        if self.GetPlotParam('set_x_max'):
            self.axes.set_xlim(right = self.GetPlotParam('x_max'))

        # now make the total energy legend
        legend_handles = []
        legend_labels = []
        for i in range(len(self.key_list)):
            if self.GetPlotParam(self.key_list[i]):
                legend_handles.append(self.plot_list[i])
                legend_labels.append(self.label_names[i])
        self.legend = self.axes.legend(legend_handles, legend_labels,
        framealpha = .05, fontsize = self.parent.MainParamDict['legendLabelSize'], loc = 1)
        self.legend.get_frame().set_facecolor('k')
        self.legend.get_frame().set_linewidth(0.0)
        if not self.GetPlotParam('show_legend'):
            self.legend.set_visible(False)

        self.legend.set_draggable(True, update = 'loc')
        if self.GetPlotParam('legend_loc') != 'N/A':
            tmp_tup = float(self.GetPlotParam('legend_loc').split()[0]),float(self.GetPlotParam('legend_loc').split()[1])
            self.legend._set_loc(tmp_tup)

        self.axes.set_xlabel(r'$t\ \  [\omega^{-1}_{pe}]$', labelpad = self.parent.MainParamDict['xLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
        self.axes.set_ylabel('Energy [arb. unit]', labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
    def refresh(self):

        '''This is a function that will be called only if self.axes already
        holds a total energy type plot. We only update things that have changed & are
        shown.  If hasn't changed or isn't shown, don't touch it. The difference
        between this and last time, is that we won't actually do any drawing in
        the plot. The plot will be redrawn after all subplots are refreshed. '''

        self.cur_time.set_xdata([self.time,self.time])

    def GetPlotParam(self, keyname):
        return self.param_dict[keyname]
