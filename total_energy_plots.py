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
                       'show_Bz_energy': False,
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
                       'legend_loc': 'N/A'}

    # We need the types of all the parameters for the config file
    BoolList = ['twoD', 'set_y_min', 'set_y_max','show_prtl_KE', 'show_field_E','show_B_E', 'show_E_E',
                'show_ion_E', 'show_electron_E', 'show_total_E', 'yLog', 'set_x_min',
                'set_x_max', 'show_legend','show_current_time', 'show_Bz_energy'] # spatial_x and spatial_y should never be true, remove from boollist
    IntList = ['E_type']
    FloatList = ['y_min', 'y_max','x_min', 'x_max']
    StrList = ['legend_loc']

    def __init__(self, parent, figwrapper):
        self.settings_window = None
        self.FigWrap = figwrapper
        self.parent = parent
        self.ChartTypes = self.FigWrap.PlotTypeDict.keys()
        self.chartType = self.FigWrap.chartType
        self.figure = self.FigWrap.figure
    def ChangePlotType(self, str_arg):
        self.FigWrap.ChangeGraph(str_arg)

    def set_plot_keys(self):
        '''A helper function that will insure that each hdf5 file will only be
        opened once per time step'''
        # First make sure that time is loaded so we can show the current time on the plot.
        self.arrs_needed = ['time']
        return self.arrs_needed

    def LoadData(self):
        ''' A Helper function that loads the data for the plot'''
        self.time = self.FigWrap.LoadKey('time')[0]


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

        self.Bz_plot = self.axes.plot(self.parent.TotalEnergyTimes, self.parent.TotalBzEnergy,
                                       ls= ':', marker = '<', markeredgecolor = self.fieldcolor,
                                       color = self.fieldcolor, visible = self.GetPlotParam('show_Bz_energy'))
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
            self.axes.set_axis_bgcolor('lightgrey')
        else:
            self.axes.set_facecolor('lightgrey')

        self.axes.tick_params(labelsize = self.parent.MainParamDict['NumFontSize'], color=tick_color)


        if self.GetPlotParam('yLog'):
            self.axes.set_yscale('log')

        # fancy code to make sure that matplotlib sets its limits
        # only based on visible lines
        self.key_list = ['show_total_E', 'show_prtl_KE', 'show_ion_E', 'show_electron_E', 'show_field_E', 'show_E_E', 'show_B_E', 'show_Bz_energy']
        self.plot_list = [self.total_plot[0], self.prtl_plot[0], self.ion_plot[0], self.electron_plot[0], self.field_plot[0], self.e_plot[0], self.mag_plot[0], self.Bz_plot[0]]
        self.label_names = ['Prtl+Field', 'Particles', 'Ions', 'Electrons', 'EM Field', 'Electric Field', 'Magnetic Field', r'$B_z^2$']

        self.axes.dataLim = mtransforms.Bbox.unit()
        self.axes.dataLim.update_from_data_xy(xy = np.vstack(self.field_plot[0].get_data()).T, ignore=True)
        for i in range(len(self.plot_list)):
            line = self.plot_list[i]
            if self.GetPlotParam(self.key_list[i]):
                xy = np.vstack(line.get_data()).T
                self.axes.dataLim.update_from_data_xy(xy, ignore=False)
        self.axes.autoscale()

        if self.GetPlotParam('set_y_min'):
            self.axes.set_ylim(ymin = self.GetPlotParam('y_min'))
        if self.GetPlotParam('set_y_max'):
            self.axes.set_ylim(ymax = self.GetPlotParam('y_max'))

        if self.GetPlotParam('set_x_min'):
            self.axes.set_xlim(xmin = self.GetPlotParam('x_min'))
        if self.GetPlotParam('set_x_max'):
            self.axes.set_xlim(xmax = self.GetPlotParam('x_max'))

        # now make the total energy legend
        legend_handles = []
        legend_labels = []
        for i in range(len(self.key_list)):
            if self.GetPlotParam(self.key_list[i]):
                legend_handles.append(self.plot_list[i])
                legend_labels.append(self.label_names[i])
        self.legend = self.axes.legend(legend_handles, legend_labels,
        framealpha = .05, fontsize = 11, loc = 1)
        self.legend.get_frame().set_facecolor('k')
        self.legend.get_frame().set_linewidth(0.0)
        if not self.GetPlotParam('show_legend'):
            self.legend.set_visible(False)

        self.legend.draggable(update = 'loc')
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

        self.electron_plot[0].set_data(self.parent.TotalEnergyTimes, self.parent.TotalElectronEnergy)
        self.ion_plot[0].set_data(self.parent.TotalEnergyTimes, self.parent.TotalIonEnergy)
        self.prtl_plot[0].set_data(self.parent.TotalEnergyTimes, self.parent.TotalElectronEnergy + self.parent.TotalIonEnergy)
        self.mag_plot[0].set_data(self.parent.TotalEnergyTimes, self.parent.TotalMagEnergy)
        self.e_plot[0].set_data(self.parent.TotalEnergyTimes, self.parent.TotalElectricEnergy)
        self.field_plot[0].set_data(self.parent.TotalEnergyTimes, self.parent.TotalMagEnergy + self.parent.TotalElectricEnergy)
        self.total_plot[0].set_data(self.parent.TotalEnergyTimes, self.parent.TotalMagEnergy + self.parent.TotalElectricEnergy + self.parent.TotalElectronEnergy + self.parent.TotalIonEnergy)
        self.Bz_plot[0].set_data(self.parent.TotalEnergyTimes, self.parent.TotalBzEnergy)

        self.cur_time.set_xdata([self.time,self.time])
        # fancy code to make sure that matplotlib sets its limits
        # based only on the visible lines.
        self.axes.dataLim = mtransforms.Bbox.unit()
        self.axes.dataLim.update_from_data_xy(xy = np.vstack(self.field_plot[0].get_data()).T, ignore=True)
        for i in range(len(self.plot_list)):
            line = self.plot_list[i]
            if self.GetPlotParam(self.key_list[i]):
                xy = np.vstack(line.get_data()).T
                self.axes.dataLim.update_from_data_xy(xy, ignore=False)
        self.axes.autoscale()

        # Set a new lims if the user chooses to do so.
        if self.GetPlotParam('set_y_min'):
            self.axes.set_ylim(ymin = self.GetPlotParam('y_min'))
        if self.GetPlotParam('set_y_max'):
            self.axes.set_ylim(ymax = self.GetPlotParam('y_max'))

#        self.axes.set_xlim(self.LimFinder(self.parent.TotalEnergyTimes))
        if self.GetPlotParam('set_x_min'):
            self.axes.set_xlim(xmin = self.GetPlotParam('x_min'))
        if self.GetPlotParam('set_x_max'):
            self.axes.set_xlim(xmax = self.GetPlotParam('x_max'))

    def GetPlotParam(self, keyname):
        return self.FigWrap.GetPlotParam(keyname)

    def SetPlotParam(self, keyname, value, update_plot = True):
        self.FigWrap.SetPlotParam(keyname, value, update_plot = update_plot)

    def OpenSettings(self):
        if self.settings_window is None:
            self.settings_window = TotEnergySettings(self)
        else:
            self.settings_window.destroy()
            self.settings_window = TotEnergySettings(self)


class TotEnergySettings(Tk.Toplevel):
    def __init__(self, parent):
        self.parent = parent
        Tk.Toplevel.__init__(self)

        self.wm_title('TotEnergy (%d,%d) Settings' % self.parent.FigWrap.pos)
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

        # the Check boxes for the dimension
        ttk.Label(frm, text='Show Prtl Energy:').grid(row = 1, column = 0, sticky = Tk.W)
        ttk.Label(frm, text='Show Field Energy:').grid(row = 1, column = 1, sticky = Tk.W)

        self.ShowElectronVar = Tk.IntVar(self) # Create a var to track whether or not to show electrons
        self.ShowElectronVar.set(self.parent.GetPlotParam('show_electron_E'))
        ttk.Checkbutton(frm, text = "Electrons",
            variable = self.ShowElectronVar,
            command = self.Selector).grid(row = 2, column = 0, sticky = Tk.W)

        self.ShowIonVar = Tk.IntVar(self) # Create a var to track whether or not to show ions
        self.ShowIonVar.set(self.parent.GetPlotParam('show_ion_E'))
        ttk.Checkbutton(frm, text = "Ions",
            variable = self.ShowIonVar,
            command = self.Selector).grid(row = 3, column = 0, sticky = Tk.W)

        self.ShowPrtlVar = Tk.IntVar(self) # Create a var to track whether or not to show tot Prtl
        self.ShowPrtlVar.set(self.parent.GetPlotParam('show_prtl_KE'))
        ttk.Checkbutton(frm, text = "Total Prtls",
            variable = self.ShowPrtlVar,
            command = self.Selector).grid(row = 4, column = 0, sticky = Tk.W)

        self.ShowMagVar = Tk.IntVar(self) # Create a var to track whether or not to plot Mag Field
        self.ShowMagVar.set(self.parent.GetPlotParam('show_B_E'))
        cb = ttk.Checkbutton(frm, text = "Magnetic Field",
            variable = self.ShowMagVar,
            command = self.Selector)
        cb.grid(row = 2, column = 1, sticky = Tk.W)

        self.ShowEVar = Tk.IntVar(self) # Create a var to track whether or not to plot E field
        self.ShowEVar.set(self.parent.GetPlotParam('show_E_E'))
        cb = ttk.Checkbutton(frm, text = "Electric Field",
            variable = self.ShowEVar,
            command = self.Selector)
        cb.grid(row = 3, column = 1, sticky = Tk.W)


        self.ShowFieldVar = Tk.IntVar(self) # Create a var to track whether or not to plot poynting energy
        self.ShowFieldVar.set(self.parent.GetPlotParam('show_field_E'))
        cb = ttk.Checkbutton(frm, text = "Total E&M Fields",
            variable = self.ShowFieldVar,
            command = self.Selector)
        cb.grid(row = 4, column = 1, sticky = Tk.W)

        self.ShowTotalVar = Tk.IntVar(self) # Create a var to track whether or not to plot poynting energy
        self.ShowTotalVar.set(self.parent.GetPlotParam('show_total_E'))
        cb = ttk.Checkbutton(frm, text = "E&M + Prtls",
            variable = self.ShowTotalVar,
            command = self.Selector)
        cb.grid(row = 5, column = 0, columnspan =1, sticky = Tk.W)

        self.ShowBzVar = Tk.IntVar(self) # Create a var to track whether or not to plot poynting energy
        self.ShowBzVar.set(self.parent.GetPlotParam('show_Bz_energy'))
        cb = ttk.Checkbutton(frm, text = "B_z*B_z",
            variable = self.ShowBzVar,
            command = self.Selector)
        cb.grid(row = 5, column = 1, columnspan =1, sticky = Tk.W)



        # Now the x & y lim
        self.setZminVar = Tk.IntVar()
        self.setZminVar.set(self.parent.GetPlotParam('set_y_min'))
        self.setZminVar.trace('w', self.setZminChanged)

        self.setZmaxVar = Tk.IntVar()
        self.setZmaxVar.set(self.parent.GetPlotParam('set_y_max'))
        self.setZmaxVar.trace('w', self.setZmaxChanged)



        self.Zmin = Tk.StringVar()
        self.Zmin.set(str(self.parent.GetPlotParam('y_min')))

        self.Zmax = Tk.StringVar()
        self.Zmax.set(str(self.parent.GetPlotParam('y_max')))


        cb = ttk.Checkbutton(frm, text ='Set y min',
                        variable = self.setZminVar)
        cb.grid(row = 2, column = 2, sticky = Tk.W)
        self.ZminEnter = ttk.Entry(frm, textvariable=self.Zmin, width=7)
        self.ZminEnter.grid(row = 2, column = 3)

        cb = ttk.Checkbutton(frm, text ='Set y max',
                        variable = self.setZmaxVar)
        cb.grid(row = 3, column = 2, sticky = Tk.W)

        self.ZmaxEnter = ttk.Entry(frm, textvariable=self.Zmax, width=7)
        self.ZmaxEnter.grid(row = 3, column = 3)

        self.setXminVar = Tk.IntVar()
        self.setXminVar.set(self.parent.GetPlotParam('set_x_min'))
        self.setXminVar.trace('w', self.setXminChanged)

        self.setXmaxVar = Tk.IntVar()
        self.setXmaxVar.set(self.parent.GetPlotParam('set_x_max'))
        self.setXmaxVar.trace('w', self.setXmaxChanged)



        self.Xmin = Tk.StringVar()
        self.Xmin.set(str(self.parent.GetPlotParam('x_min')))

        self.Xmax = Tk.StringVar()
        self.Xmax.set(str(self.parent.GetPlotParam('x_max')))


        cb = ttk.Checkbutton(frm, text ='Set x min',
                        variable = self.setXminVar)
        cb.grid(row = 4, column = 2, sticky = Tk.W)
        self.XminEnter = ttk.Entry(frm, textvariable=self.Xmin, width=7)
        self.XminEnter.grid(row = 4, column = 3)

        cb = ttk.Checkbutton(frm, text ='Set x max',
                        variable = self.setXmaxVar)
        cb.grid(row = 5, column = 2, sticky = Tk.W)

        self.XmaxEnter = ttk.Entry(frm, textvariable=self.Xmax, width=7)
        self.XmaxEnter.grid(row = 5, column = 3)


        # Now whether or not the y axes should be in logspace

        self.yLogVar = Tk.IntVar()
        self.yLogVar.set(self.parent.GetPlotParam('yLog'))
        self.yLogVar.trace('w', self.yLogChanged)



        cb = ttk.Checkbutton(frm, text ='y-axis logscale',
                        variable = self.yLogVar)
        cb.grid(row = 6, column = 0, sticky = Tk.W)

        self.LegendVar = Tk.IntVar()
        self.LegendVar.set(self.parent.GetPlotParam('show_legend'))
        self.LegendVar.trace('w', self.showLegendChanged)



        cb = ttk.Checkbutton(frm, text ='Show legend',
                        variable = self.LegendVar)
        cb.grid(row = 6, column = 1, sticky = Tk.W)


        self.CurTimeVar = Tk.IntVar()
        self.CurTimeVar.set(self.parent.GetPlotParam('show_current_time'))
        self.LegendVar.trace('w', self.showtimeChanged)

        cb = ttk.Checkbutton(frm, text ='Show current time',
                        variable = self.CurTimeVar)
        cb.grid(row = 7, column = 0, sticky = Tk.W)



    def ctypeChanged(self, *args):
        if self.ctypevar.get() == self.parent.chartType:
            pass
        else:
            self.parent.ChangePlotType(self.ctypevar.get())
            self.destroy()

    def setZminChanged(self, *args):
        if self.setZminVar.get() == self.parent.GetPlotParam('set_y_min'):
            pass
        else:
            self.parent.SetPlotParam('set_y_min', self.setZminVar.get())

    def setZmaxChanged(self, *args):
        if self.setZmaxVar.get() == self.parent.GetPlotParam('set_y_max'):
            pass
        else:
            self.parent.SetPlotParam('set_y_max', self.setZmaxVar.get())

    def setXminChanged(self, *args):
        if self.setXminVar.get() == self.parent.GetPlotParam('set_x_min'):
            pass
        else:
            self.parent.SetPlotParam('set_x_min', self.setXminVar.get())

    def setXmaxChanged(self, *args):
        if self.setXmaxVar.get() == self.parent.GetPlotParam('set_x_max'):
            pass
        else:
            self.parent.SetPlotParam('set_x_max', self.setXmaxVar.get())


    def yLogChanged(self, *args):
        if self.yLogVar.get() == self.parent.GetPlotParam('yLog'):
            pass
        else:
            if self.yLogVar.get():
                self.parent.axes.set_yscale('log')
            else:
                self.parent.axes.set_yscale('linear')

            self.parent.SetPlotParam('yLog', self.yLogVar.get())

    def showLegendChanged(self, *args):
        if self.LegendVar.get() == self.parent.GetPlotParam('show_legend'):
            pass
        else:
            self.parent.legend.set_visible(self.LegendVar.get())
            self.parent.SetPlotParam('show_legend', self.LegendVar.get())

    def showtimeChanged(self, *args):
        if self.CurTimeVar.get() == self.parent.GetPlotParam('show_current_time'):
            pass
        else:
            self.parent.cur_time.set_visible(self.CurTimeVar.get())
            self.parent.SetPlotParam('show_current_time', self.CurTimeVar.get())


    def Selector(self):
        # Repeat the lists to remember the order
        # self.parent.key_list = ['show_prtl_KE', 'show_ion_E', 'show_electron_E', 'show_field_E', 'show_E_E', 'show_B_E']
        # self.parent.plot_list = [self.prtl_plot[0], self.ion_plot[0], self.electron_plot[0], self.field_plot[0], self.e_plot[0], self.mag_plot[0]]
        # self.parent.label_names = ['Particles', 'Ions', 'Electrons', 'EM Field', 'Electric Field', 'Magnetic Field']
        VarList = [self.ShowTotalVar, self.ShowPrtlVar,  self.ShowIonVar, self.ShowElectronVar, self.ShowFieldVar, self.ShowEVar, self.ShowMagVar, self.ShowBzVar]

        # Update current legend position
        if self.parent.legend._get_loc() != 1:
            self.parent.SetPlotParam('legend_loc', ' '.join(str(x) for x in self.parent.legend._get_loc()), update_plot = False)

        # First set the visibility of the plots to their new value
        for i in range(len(self.parent.key_list)):
            if self.parent.GetPlotParam(self.parent.key_list[i]) != VarList[i].get():
                self.parent.plot_list[i].set_visible(VarList[i].get())
                self.parent.SetPlotParam(self.parent.key_list[i], VarList[i].get(), update_plot = False)

        # Now create a new legend with only the visible lines:
        legend_handles = []
        legend_labels = []
        for i in range(len(self.parent.key_list)):
            if VarList[i].get():
                legend_handles.append(self.parent.plot_list[i])
                legend_labels.append(self.parent.label_names[i])
        self.parent.legend = self.parent.axes.legend(legend_handles, legend_labels,
        framealpha = .05, fontsize = 11, loc = 1)
        self.parent.legend.set_visible(self.parent.GetPlotParam('show_legend'))
        self.parent.legend.get_frame().set_facecolor('k')
        self.parent.legend.get_frame().set_linewidth(0.0)
        self.parent.legend.draggable()
        tmp_tup = 1
        if self.parent.GetPlotParam('legend_loc') != 'N/A':
            tmp_tup = float(self.parent.GetPlotParam('legend_loc').split()[0]),float(self.parent.GetPlotParam('legend_loc').split()[1])
            self.parent.legend._set_loc(tmp_tup)
        self.parent.legend._set_loc(tmp_tup)

        # Force a plot refresh
        self.parent.SetPlotParam(self.parent.key_list[0], VarList[0].get())

    def TxtEnter(self, e):
        self.FieldsCallback()

    def FieldsCallback(self):
        tkvarLimList = [self.Zmin, self.Zmax, self.Xmin, self.Xmax]
        plot_param_List = ['y_min', 'y_max', 'x_min', 'x_max']
        tkvarSetList = [self.setZminVar, self.setZmaxVar, self.setXminVar, self.setXmaxVar]
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
            self.parent.SetPlotParam('y_min', self.parent.GetPlotParam('y_min'))


    def OnClosing(self):
        self.parent.settings_window = None
        self.destroy()
