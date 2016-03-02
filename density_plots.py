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
    # A diction of all of the parameters for this plot with the default parameters

    plot_param_dict = {'twoD': 0,
                       'dens_type': 0, #0 = n, 1 = rho
                       'show_cbar': True,
                       'set_color_limits': False,
                       'z_min': None,
                       'z_max' : None,
                       'show_labels' : True,
                       'show_shock' : False,
                       'OutlineText': True,
                       'interpolation': 'hermite'}

    def __init__(self, parent, figwrapper):
        self.settings_window = None
        self.FigWrap = figwrapper
        self.parent = parent
        self.ChartTypes = self.FigWrap.PlotTypeDict.keys()
        self.chartType = self.FigWrap.chartType
        self.figure = self.FigWrap.figure
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

        self.arrs_needed = ['c_omp', 'istep', 'sizex', 'dens']
        # To plot rho we need both dens and densi
        if self.GetPlotParam('dens_type') == 1: # Load the ion density
            self.arrs_needed.append('densi')

        return self.arrs_needed

    def draw(self):
        # Get the color from the colormap
        self.dens_color = new_cmaps.cmaps[self.parent.cmap](0.5)
        # get c_omp and istep to convert cells to physical units
        self.c_omp = self.FigWrap.LoadKey('c_omp')[0]
        self.istep = self.FigWrap.LoadKey('istep')[0]


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

        # load the density values
        self.zval= self.FigWrap.LoadKey('dens')[0,:,:]


        if self.GetPlotParam('dens_type') == 1: # Load calculate rho
            self.zval = 2*self.FigWrap.LoadKey('densi')[0,:,:]-self.zval

        # Generate the x and y axes
        self.y_values =  np.arange(self.zval.shape[0])/self.c_omp*self.istep
        self.x_values =  np.arange(self.zval.shape[1])/self.c_omp*self.istep

        self.axes = self.figure.add_subplot(self.gs[18:92,:])

        # Now that the data is loaded, start making the plots
        if self.GetPlotParam('twoD'):
            # First choose the 'zval' to plot, we can only do one because it is 2-d.
            self.two_d_label = r'$n_e$'
            if self.FigWrap.GetPlotParam('dens_type') == 1:
                self.two_d_label = r'$\rho$'

            self.ymin = 0
            self.ymax =  self.zval.shape[0]/self.c_omp*self.istep
            self.xmin = 0
            self.xmax =  self.zval.shape[1]/self.c_omp*self.istep

            self.cax = self.axes.imshow(self.zval,
                cmap = new_cmaps.cmaps[self.parent.cmap],
                origin = 'lower', aspect = 'auto',
                extent = (self.xmin,self.xmax, self.ymin, self.ymax),
                interpolation=self.GetPlotParam('interpolation'))


            if self.GetPlotParam('show_shock'):
                self.axes.axvline(self.parent.shock_loc, linewidth = 1.5, linestyle = '--', color = self.parent.shock_color, path_effects=[PathEffects.Stroke(linewidth=2, foreground='k'),
                                    PathEffects.Normal()])

            if self.FigWrap.GetPlotParam('show_labels'):
                self.axes.annotate(self.two_d_label,
                                xy = (0.9,.9),
                                xycoords= 'axes fraction',
                                color = 'white',
                                **self.annotate_kwargs)
            self.axes.set_axis_bgcolor('lightgrey')

            if self.GetPlotParam('show_cbar'):
                self.axC = self.figure.add_subplot(self.gs[:4,:])
                self.cbar = self.figure.colorbar(self.cax, ax = self.axes, cax = self.axC, orientation = 'horizontal')

                self.cbar.set_ticks(np.linspace(self.zval.min(),self.zval.max(), 5))
                self.cbar.ax.tick_params(labelsize=10)

            self.axes.set_axis_bgcolor('lightgrey')
            self.axes.tick_params(labelsize = 10, color=tick_color)
#        self.axes.set_xlim(self.xmin,self.xmax)
            self.axes.set_xlabel(r'$x\ [c/\omega_{\rm pe}]$', labelpad = self.parent.xlabel_pad, color = 'black')
            self.axes.set_ylabel(r'$y\ [c/\omega_{\rm pe}]$', labelpad = self.parent.ylabel_pad, color = 'black')

        else:
            # Make the 1-D plots
            self.axes.plot(self.x_values, self.zval[self.zval.shape[0]/2,:], color = self.dens_color)
            tmp_str = r'$\rm density$'
            if self.GetPlotParam('dens_type') == 1:
                tmp_str = r'$\rho$'


            if self.GetPlotParam('show_shock'):
                self.axes.axvline(self.parent.shock_loc, linewidth = 1.5, linestyle = '--', color = self.parent.shock_color, path_effects=[PathEffects.Stroke(linewidth=2, foreground='k'),
                        PathEffects.Normal()])


            self.axes.set_axis_bgcolor('lightgrey')
            self.axes.tick_params(labelsize = 10, color=tick_color)
            self.axes.set_xlim(self.x_values[0],self.x_values[-1])

            self.axes.set_xlabel(r'$x\ [c/\omega_{\rm pe}]$', labelpad = self.parent.xlabel_pad, color = 'black')
            self.axes.set_ylabel(tmp_str, labelpad = self.parent.ylabel_pad, color = 'black')

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
        #Create some sizers




#                   'field_type': 0, #0 = B-Field, 1 = E-field
#                   'show_x' : 1,
#                   'show_y' : 1,
#                   'show_z' : 1,
#                   'show_cbar': True,
#                   'set_color_limits': False,
#                   'v_min': None,
#                   'OutlineText': True,
#                   'interpolation': 'hermite',
#                   'v_max': None
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
                        command = lambda:
                        self.parent.SetPlotParam('show_cbar', self.CbarVar.get()))
        cb.grid(row = 6, sticky = Tk.W)

        # show shock
        self.ShockVar = Tk.IntVar()
        self.ShockVar.set(self.parent.GetPlotParam('show_shock'))
        cb = ttk.Checkbutton(frm, text = "Show Shock",
                        variable = self.ShockVar,
                        command = lambda:
                        self.parent.SetPlotParam('show_shock', self.ShockVar.get()))
        cb.grid(row = 6, column = 1, sticky = Tk.W)




    def Change2d(self):
        if self.TwoDVar.get() == self.parent.GetPlotParam('twoD'):
            pass
        else:
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

    def RadioField(self):
        if self.DensTypeVar.get() == self.parent.GetPlotParam('dens_type'):
            pass
        else:
            self.parent.SetPlotParam('dens_type', self.DensTypeVar.get())


    def OnClosing(self):
        self.parent.settings_window = None
        self.destroy()