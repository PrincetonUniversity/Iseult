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

class DensPanel:
    # A dictionary of all of the parameters for this plot with the default parameters

    plot_param_dict = {'twoD': 0,
                       'dens_type': 0, #0 = n, 1 = rho
                       'show_cbar': True,
                       'set_color_limits': False,
                       'z_min': 0,
                       'z_max' : 10,
                       'set_z_min': False,
                       'set_z_max': False,
                       'show_labels' : True,
                       'show_shock' : False,
                       'OutlineText': True,
                       'spatial_x': True,
                       'spatial_y': None,
                       'interpolation': 'nearest'}

    def __init__(self, parent, figwrapper):
        self.settings_window = None
        self.FigWrap = figwrapper
        self.parent = parent
        self.ChartTypes = self.FigWrap.PlotTypeDict.keys()
        self.chartType = self.FigWrap.chartType
        self.figure = self.FigWrap.figure
        self.SetPlotParam('spatial_y', self.GetPlotParam('twoD'), update_plot = False)
        self.InterpolationMethods = ['nearest', 'bilinear', 'bicubic', 'spline16',
            'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric',
            'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']


    def ChangePlotType(self, str_arg):
        self.FigWrap.ChangeGraph(str_arg)

    def set_plot_keys(self):
        '''A helper function that will insure that each hdf5 file will only be
        opened once per time step'''
        # First make sure that omega_plasma & xi is loaded so we can fix the
        # x & y distances.

        self.arrs_needed = ['c_omp', 'istep', 'dens', 'densi']
        # To plot rho we need both dens and densi
#        if self.GetPlotParam('dens_type') == 1: # Load the ion density
#            self.arrs_needed.append('densi')

        return self.arrs_needed
    def LoadData(self):
        ''' A Helper function that loads the data for the plot'''

        self.dens_color = new_cmaps.cmaps[self.parent.cmap](0.5)
        # get c_omp and istep to convert cells to physical units
        self.c_omp = self.FigWrap.LoadKey('c_omp')[0]
        self.istep = self.FigWrap.LoadKey('istep')[0]

        self.dens = self.FigWrap.LoadKey('dens')[0,:,:]
        # see if this time has already been checked
        if 'my_rho' in self.parent.DataDict.keys():
            self.rho = self.parent.DataDict['my_rho']
        else:
            self.rho = 2*self.FigWrap.LoadKey('densi')[0,:,:] -self.FigWrap.LoadKey('dens')[0,:,:]
            self.parent.DataDict['my_rho'] = self.rho

        # see if the min/max of all the arrays has aready been calculated.
        if 'dens_min_max' in self.parent.DataDict.keys():
            self.dens_min_max = self.parent.DataDict['dens_min_max']
        else:
            self.dens_min_max = [self.dens.min(), self.dens.max()]
            self.parent.DataDict['dens_min_max'] = self.dens_min_max

        if 'rho_min_max' in self.parent.DataDict.keys():
            self.rho_min_max = self.parent.DataDict['rho_min_max']
        else:
            self.rho_min_max = [self.rho.min(), self.rho.max()]
            self.parent.DataDict['rho_min_max'] = self.rho_min_max

        if 'xaxis_values' in self.parent.DataDict.keys():
            # Generate the x and y axes
            self.xaxis_values = self.parent.DataDict['xaxis_values']
        else:
            self.xaxis_values = np.arange(self.dens.shape[1])/self.c_omp*self.istep
            self.parent.DataDict['xaxis_values'] = self.xaxis_values
        # y values not needed so commenting out
        # self.y_values =  np.arange(self.zval.shape[0])/self.c_omp*self.istep


    def draw(self):
        if self.GetPlotParam('OutlineText'):
            self.annotate_kwargs = {'horizontalalignment': 'right',
            'verticalalignment': 'top',
            'size' : 18,
            'path_effects' : [PathEffects.withStroke(linewidth=1.5,foreground="k")]
            }
        else:
            self.annotate_kwargs = {'horizontalalignment' : 'right',
            'verticalalignment' : 'top',
            'size' : 18}

        # Set the tick color
        tick_color = 'black'

        # Create a gridspec to handle spacing better
        self.gs = gridspec.GridSpecFromSubplotSpec(100,100, subplot_spec = self.parent.gs0[self.FigWrap.pos])#, bottom=0.2,left=0.1,right=0.95, top = 0.95)

        self.LoadData()

        # Now that the data is loaded, start making the plots
        if self.GetPlotParam('twoD'):

            # Link up the spatial axes if desired
            if self.parent.LinkSpatial != 0:
                if self.FigWrap.pos == self.parent.first_x and self.FigWrap.pos == self.parent.first_y:
                    self.axes = self.figure.add_subplot(self.gs[18:92,:])
                elif self.FigWrap.pos == self.parent.first_x:
                    self.axes = self.figure.add_subplot(self.gs[18:92,:],
                    sharey = self.parent.SubPlotList[self.parent.first_y[0]][self.parent.first_y[1]].graph.axes)
                elif self.FigWrap.pos == self.parent.first_y:
                    self.axes = self.figure.add_subplot(self.gs[18:92,:],
                    sharex = self.parent.SubPlotList[self.parent.first_x[0]][self.parent.first_x[1]].graph.axes)
                else:
                    self.axes = self.figure.add_subplot(self.gs[18:92,:],
                    sharex = self.parent.SubPlotList[self.parent.first_x[0]][self.parent.first_x[1]].graph.axes,
                    sharey = self.parent.SubPlotList[self.parent.first_y[0]][self.parent.first_y[1]].graph.axes)


            else:
                self.axes = self.figure.add_subplot(self.gs[18:92,:])

            # First choose the 'zval' to plot, we can only do one because it is 2-d.
            if self.FigWrap.GetPlotParam('dens_type') == 0:
                self.zval = self.dens
                self.two_d_label = r'$n_e$'
            if self.FigWrap.GetPlotParam('dens_type') == 1:
                self.zval = self.rho
                self.two_d_label = r'$\rho$'

            self.ymin = 0
            self.ymax =  self.zval.shape[0]/self.c_omp*self.istep
            self.xmin = 0
            self.xmax =  self.zval.shape[1]/self.c_omp*self.istep

            self.vmin = None
            if self.GetPlotParam('set_z_min'):
                self.vmin = self.GetPlotParam('z_min')
            self.vmax = None
            if self.GetPlotParam('set_z_max'):
                self.vmax = self.GetPlotParam('z_max')

            if self.parent.plot_aspect:
                self.cax = self.axes.imshow(self.zval,
                    cmap = new_cmaps.cmaps[self.parent.cmap],
                    origin = 'lower',
                    extent = (self.xmin,self.xmax, self.ymin, self.ymax),
                    vmin = self.vmin, vmax = self.vmax,
                    interpolation=self.GetPlotParam('interpolation'))

            else:
                self.cax = self.axes.imshow(self.zval,
                    cmap = new_cmaps.cmaps[self.parent.cmap],
                    origin = 'lower', aspect = 'auto',
                    extent = (self.xmin,self.xmax, self.ymin, self.ymax),
                    vmin = self.vmin, vmax = self.vmax,
                    interpolation=self.GetPlotParam('interpolation'))

            self.shockline_2d = self.axes.axvline(self.parent.shock_loc,
                                                    linewidth = 1.5,
                                                    linestyle = '--',
                                                    color = self.parent.shock_color,
                                                    path_effects=[PathEffects.Stroke(linewidth=2, foreground='k'),
                                                    PathEffects.Normal()])
            self.shockline_2d.set_visible(self.GetPlotParam('show_shock'))

            self.an_2d = self.axes.annotate(self.two_d_label,
                                            xy = (0.9,.9),
                                            xycoords= 'axes fraction',
                                            color = 'white',
                                            **self.annotate_kwargs)
            self.an_2d.set_visible(self.GetPlotParam('show_labels'))

            self.axes.set_axis_bgcolor('lightgrey')

            self.axC = self.figure.add_subplot(self.gs[:4,:])
            self.cbar = self.figure.colorbar(self.cax, ax = self.axes, cax = self.axC, orientation = 'horizontal')
            if self.GetPlotParam('show_cbar'):
                cmin = self.zval.min()
                if self.vmin:
                    cmin = self.vmin
                cmax = self.zval.max()
                if self.vmax:
                    cmax = self.vmax
                self.cbar.set_ticks(np.linspace(cmin, cmax, 5))
                self.cbar.ax.tick_params(labelsize=self.parent.num_font_size)
            else:
                self.axC.set_visible(False)

            self.axes.set_axis_bgcolor('lightgrey')
            self.axes.tick_params(labelsize = self.parent.num_font_size, color=tick_color)

            if self.parent.xlim[0] and self.parent.LinkSpatial != 0:
                self.axes.set_xlim(self.parent.xlim[1],self.parent.xlim[2])
            else:
                self.axes.set_xlim(self.xmin,self.xmax)

            if self.parent.ylim[0]:
                self.axes.set_ylim(self.parent.ylim[1], self.parent.ylim[2])

            self.axes.set_xlabel(r'$x\ [c/\omega_{\rm pe}]$', labelpad = self.parent.xlabel_pad, color = 'black')
            self.axes.set_ylabel(r'$y\ [c/\omega_{\rm pe}]$', labelpad = self.parent.ylabel_pad, color = 'black')

        else:
            # Do the 1D Plots
            if self.parent.LinkSpatial != 0 and self.parent.LinkSpatial != 3:
                if self.FigWrap.pos == self.parent.first_x:
                    self.axes = self.figure.add_subplot(self.gs[18:92,:])
                else:
                    self.axes = self.figure.add_subplot(self.gs[18:92,:],
                    sharex = self.parent.SubPlotList[self.parent.first_x[0]][self.parent.first_x[1]].graph.axes)
            else:
                self.axes = self.figure.add_subplot(self.gs[18:92,:])

            # Make the 1-D plots
            self.linedens = self.axes.plot(self.xaxis_values, self.dens[self.dens.shape[0]/2,:], color = self.dens_color)
            self.linedens[0].set_visible(not self.GetPlotParam('dens_type')) #visible if dens_type == 0

            self.linerho = self.axes.plot(self.xaxis_values, self.dens[self.dens.shape[0]/2,:], color = self.dens_color)
            self.linerho[0].set_visible(self.GetPlotParam('dens_type'))


            self.shock_line =self.axes.axvline(self.parent.shock_loc, linewidth = 1.5, linestyle = '--', color = self.parent.shock_color, path_effects=[PathEffects.Stroke(linewidth=2, foreground='k'),
                    PathEffects.Normal()])
            self.shock_line.set_visible(self.GetPlotParam('show_shock'))

            self.axes.set_axis_bgcolor('lightgrey')
            self.axes.tick_params(labelsize = self.parent.num_font_size, color=tick_color)

            if self.parent.xlim[0]:
                self.axes.set_xlim(self.parent.xlim[1],self.parent.xlim[2])
            else:
                self.axes.set_xlim(self.xaxis_values[0],self.xaxis_values[-1])

            if self.GetPlotParam('set_z_min'):
                self.axes.set_ylim(ymin = self.GetPlotParam('z_min'))
            if self.GetPlotParam('set_z_max'):
                self.axes.set_ylim(ymax = self.GetPlotParam('z_max'))

            # Handle the axes labeling
            tmp_str = r'$\rm density$'
            if self.GetPlotParam('dens_type') == 1:
                tmp_str = r'$\rho$'
            self.axes.set_xlabel(r'$x\ [c/\omega_{\rm pe}]$', labelpad = self.parent.xlabel_pad, color = 'black')
            self.axes.set_ylabel(tmp_str, labelpad = self.parent.ylabel_pad, color = 'black')


    def refresh(self):
        '''This is a function that will be called only if self.axes already
        holds a fields type plot. We only update things that have shown.  If
        hasn't changed, or isn't viewed, don't touch it. The difference between this and last
        time, is that we won't actually do any drawing in the plot. The plot
        will be redrawn after all subplots data is changed. '''
        self.LoadData()

        # Main goal, only change what is showing..
        # First do the 1D plots, because it is simpler
        if self.GetPlotParam('twoD') == 0:
            if self.GetPlotParam('dens_type') == 0:
                self.linedens[0].set_data(self.xaxis_values, self.dens[self.dens.shape[0]/2,:])
                self.axes.set_ylim(self.dens_min_max)
            else:
                self.linerho[0].set_data(self.xaxis_values, self.rho[self.rho.shape[0]/2,:])
                self.axes.set_ylim(self.rho_min_max)
            if self.GetPlotParam('show_shock'):
                self.shock_line.set_xdata([self.parent.shock_loc,self.parent.shock_loc])



            if self.parent.xlim[0]:
                self.axes.set_xlim(self.parent.xlim[1],self.parent.xlim[2])
            else:
                self.axes.set_xlim(self.xaxis_values[0], self.xaxis_values[-1])


            if self.GetPlotParam('set_z_min'):
                self.axes.set_ylim(ymin = self.GetPlotParam('z_min'))
            if self.GetPlotParam('set_z_max'):
                self.axes.set_ylim(ymax = self.GetPlotParam('z_max'))

        else: # Now refresh the plot if it is 2D
            if self.GetPlotParam('dens_type')==0:
                self.cax.set_data(self.dens)
                self.ymin = 0
                self.ymax =  self.dens.shape[0]/self.c_omp*self.istep
                self.xmin = 0
                self.xmax =  self.xaxis_values[-1]
                self.clims = self.dens_min_max

            else:
                self.cax.set_data(self.rho)
                self.ymin = 0
                self.ymax =  self.rho.shape[0]/self.c_omp*self.istep
                self.xmin = 0
                self.xmax =  self.xaxis_values[-1]
                self.clims = self.rho_min_max

            self.cax.set_extent([self.xmin,self.xmax, self.ymin, self.ymax])
            self.axes.set_xlim(self.xmin,self.xmax)
            self.axes.set_ylim(self.ymin,self.ymax)

            self.cax.set_clim(self.clims)

            self.climArgs = {}
            if self.GetPlotParam('set_z_min'):
                self.climArgs['vmin'] =  self.GetPlotParam('z_min')
            if self.GetPlotParam('set_z_max'):
                self.climArgs['vmax'] =  self.GetPlotParam('z_max')
            if len(self.climArgs)>0:
                self.cax.set_clim(**self.climArgs)
            if self.GetPlotParam('show_cbar'):
                self.cbar.set_ticks(np.linspace(self.cax.get_clim()[0],self.cax.get_clim()[1], 5))

            if self.GetPlotParam('show_shock'):
                self.shockline_2d.set_xdata([self.parent.shock_loc,self.parent.shock_loc])

    def GetPlotParam(self, keyname):
        return self.FigWrap.GetPlotParam(keyname)

    def SetPlotParam(self, keyname, value, update_plot = True):
        self.FigWrap.SetPlotParam(keyname, value, update_plot = update_plot)

    def OpenSettings(self):
        if self.settings_window is None:
            self.settings_window = DensSettings(self)
        else:
            self.settings_window.destroy()
            self.settings_window = DensSettings(self)


class DensSettings(Tk.Toplevel):
    def __init__(self, parent):
        self.parent = parent
        Tk.Toplevel.__init__(self)

        self.wm_title('Dens Plot (%d,%d) Settings' % self.parent.FigWrap.pos)
        self.parent = parent
        frm = ttk.Frame(self)
        frm.pack(fill=Tk.BOTH, expand=True)
        self.protocol('WM_DELETE_WINDOW', self.OnClosing)
        self.bind('<Return>', self.TxtEnter)


        # Create the OptionMenu to chooses the Chart Type:
        self.InterpolVar = Tk.StringVar(self)
        self.InterpolVar.set(self.parent.GetPlotParam('interpolation')) # default value
        self.InterpolVar.trace('w', self.InterpolChanged)

        ttk.Label(frm, text="Interpolation Method:").grid(row=0, column = 2)
        InterplChooser = apply(ttk.OptionMenu, (frm, self.InterpolVar, self.parent.GetPlotParam('interpolation')) + tuple(self.parent.InterpolationMethods))
        InterplChooser.grid(row =0, column = 3, sticky = Tk.W + Tk.E)

        # Create the OptionMenu to chooses the Chart Type:
        self.ctypevar = Tk.StringVar(self)
        self.ctypevar.set(self.parent.chartType) # default value
        self.ctypevar.trace('w', self.ctypeChanged)

        ttk.Label(frm, text="Choose Chart Type:").grid(row=0, column = 0)
        cmapChooser = apply(ttk.OptionMenu, (frm, self.ctypevar, self.parent.chartType) + tuple(self.parent.ChartTypes))
        cmapChooser.grid(row =0, column = 1, sticky = Tk.W + Tk.E)


        self.TwoDVar = Tk.IntVar(self) # Create a var to track whether or not to plot in 2-D
        self.TwoDVar.set(self.parent.GetPlotParam('twoD'))
        cb = ttk.Checkbutton(frm, text = "Show in 2-D",
                variable = self.TwoDVar,
                command = self.Change2d)
        cb.grid(row = 1, sticky = Tk.W)

        # the Radiobox Control to choose the Field Type
        self.DensList = ['dens_e', 'rho']
        self.DensTypeVar  = Tk.IntVar()
        self.DensTypeVar.set(self.parent.GetPlotParam('dens_type'))

        ttk.Label(frm, text='Choose Density:').grid(row = 2, sticky = Tk.W)

        for i in range(len(self.DensList)):
            ttk.Radiobutton(frm,
                text=self.DensList[i],
                variable=self.DensTypeVar,
                command = self.RadioField,
                value=i).grid(row = 3+i, sticky =Tk.W)


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

        # Now the field lim
        self.setZminVar = Tk.IntVar()
        self.setZminVar.set(self.parent.GetPlotParam('set_z_min'))
        self.setZminVar.trace('w', self.setZminChanged)

        self.setZmaxVar = Tk.IntVar()
        self.setZmaxVar.set(self.parent.GetPlotParam('set_z_max'))
        self.setZmaxVar.trace('w', self.setZmaxChanged)



        self.Zmin = Tk.StringVar()
        self.Zmin.set(str(self.parent.GetPlotParam('z_min')))

        self.Zmax = Tk.StringVar()
        self.Zmax.set(str(self.parent.GetPlotParam('z_max')))


        cb = ttk.Checkbutton(frm, text ='Set dens min',
                        variable = self.setZminVar)
        cb.grid(row = 3, column = 2, sticky = Tk.W)
        self.ZminEnter = ttk.Entry(frm, textvariable=self.Zmin, width=7)
        self.ZminEnter.grid(row = 3, column = 3)

        cb = ttk.Checkbutton(frm, text ='Set dens max',
                        variable = self.setZmaxVar)
        cb.grid(row = 4, column = 2, sticky = Tk.W)

        self.ZmaxEnter = ttk.Entry(frm, textvariable=self.Zmax, width=7)
        self.ZmaxEnter.grid(row = 4, column = 3)

    def ShockVarHandler(self, *args):
        if self.parent.GetPlotParam('show_shock')== self.ShockVar.get():
            pass
        else:
            if self.parent.GetPlotParam('twoD'):
                self.parent.shockline_2d.set_visible(self.ShockVar.get())
            else:
                self.parent.shock_line.set_visible(self.ShockVar.get())

            self.parent.SetPlotParam('show_shock', self.ShockVar.get())

    def CbarHandler(self, *args):
        if self.parent.GetPlotParam('show_cbar')== self.CbarVar.get():
            pass
        else:
            if self.parent.GetPlotParam('twoD'):
                self.parent.axC.set_visible(self.CbarVar.get())

            self.parent.SetPlotParam('show_cbar', self.CbarVar.get(), update_plot =self.parent.GetPlotParam('twoD'))


    def Change2d(self):
        if self.TwoDVar.get() == self.parent.GetPlotParam('twoD'):
            pass
        else:
            self.parent.SetPlotParam('spatial_y', self.TwoDVar.get(), update_plot = False)
            self.parent.SetPlotParam('twoD', self.TwoDVar.get())


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
            self.parent.SetPlotParam('interpolation', self.InterpolVar.get())

    def setZminChanged(self, *args):
        if self.setZminVar.get() == self.parent.GetPlotParam('set_z_min'):
            pass
        else:
            self.parent.SetPlotParam('set_z_min', self.setZminVar.get())

    def setZmaxChanged(self, *args):
        if self.setZmaxVar.get() == self.parent.GetPlotParam('set_z_max'):
            pass
        else:
            self.parent.SetPlotParam('set_z_max', self.setZmaxVar.get())

    def TxtEnter(self, e):
        self.FieldsCallback()

    def FieldsCallback(self):
        tkvarLimList = [self.Zmin, self.Zmax]
        plot_param_List = ['z_min', 'z_max']
        tkvarSetList = [self.setZminVar, self.setZmaxVar]
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
            self.parent.SetPlotParam('z_min', self.parent.GetPlotParam('z_min'))


    def RadioField(self):
        if self.DensTypeVar.get() == self.parent.GetPlotParam('dens_type'):
            pass
        else:
            if self.parent.GetPlotParam('twoD'):
                if self.DensTypeVar.get() == 0:
                    self.parent.an_2d.set_text(r'$n_e$')
                else:
                    self.parent.an_2d.set_text(r'$\rho$')
            else:
                if self.DensTypeVar.get() == 0:
                    self.parent.linerho[0].set_visible(False)
                    self.parent.linedens[0].set_visible(True)
                    self.parent.axes.set_ylabel('density')
                else:
                    self.parent.linerho[0].set_visible(True)
                    self.parent.linedens[0].set_visible(False)
                    self.parent.axes.set_ylabel(r'$\rho$')

            self.parent.SetPlotParam('dens_type', self.DensTypeVar.get())


    def OnClosing(self):
        self.parent.settings_window = None
        self.destroy()
