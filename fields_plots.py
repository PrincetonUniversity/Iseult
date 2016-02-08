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

class FieldsPanel:
    # A diction of all of the parameters for this plot with the default parameters

    plot_param_dict = {'twoD': 0,
                       'field_type': 0, #0 = B-Field, 1 = E-field
                       'show_x' : 1,
                       'show_y' : 1,
                       'show_z' : 1,
                       'show_cbar': True,
                       'set_color_limits': False,
                       'z_min': None,
                       'z_max' : None,
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

        self.arrs_needed = ['c_omp']
        # Then see if we are plotting E-field or B-Field
        if self.GetPlotParam('field_type') == 0: # Load the B-Field
            if self.GetPlotParam('twoD'):
                #If it is 2-D we can only show one dim of the field.
                if self.GetPlotParam('show_x'):
                    self.arrs_needed.append('bx')
                elif self.GetPlotParam('show_y'):
                    self.arrs_needed.append('by')
                elif self.GetPlotParam('show_z'):
                    self.arrs_needed.append('bz')
            else:
                if self.GetPlotParam('show_x'):
                    self.arrs_needed.append('bx')
                if self.GetPlotParam('show_y'):
                    self.arrs_needed.append('by')
                if self.GetPlotParam('show_z'):
                    self.arrs_needed.append('bz')

        if self.GetPlotParam('field_type') == 1: # Load the E-Field
            if self.GetPlotParam('twoD'):
                if self.GetPlotParam('show_x'):
                    self.arrs_needed.append('ex')
                elif self.GetPlotParam('show_y'):
                    self.arrs_needed.append('ey')
                elif self.GetPlotParam('show_z'):
                    self.arrs_needed.append('ez')

            else:
                if self.GetPlotParam('show_x'):
                    self.arrs_needed.append('ex')
                if self.GetPlotParam('show_y'):
                    self.arrs_needed.append('ey')
                if self.GetPlotParam('show_z'):
                    self.arrs_needed.append('ez')

        return self.arrs_needed

    def draw(self):
        # Get the x, y, and z colors from the colormap
        self.xcolor = new_cmaps.cmaps[self.parent.cmap](0.2)
        self.ycolor = new_cmaps.cmaps[self.parent.cmap](0.5)
        self.zcolor = new_cmaps.cmaps[self.parent.cmap](0.8)


        self.c_omp = self.FigWrap.LoadKey('c_omp')[0]

        # Set the tick color
        tick_color = 'white'

        # Create a gridspec to handle spacing better
        self.gs = gridspec.GridSpecFromSubplotSpec(100,100, subplot_spec = self.parent.gs0[self.FigWrap.pos])#, bottom=0.2,left=0.1,right=0.95, top = 0.95)

        if self.GetPlotParam('field_type') == 0: # Load the B-Field
            if self.GetPlotParam('twoD'):
                if self.GetPlotParam('show_x'):
                    self.zval = self.FigWrap.LoadKey('bx')[0,:,:]
                elif self.GetPlotParam('show_y'):
                    self.zval = self.FigWrap.LoadKey('by')[0,:,:]
                elif self.GetPlotParam('show_z'):
                    self.zval = self.FigWrap.LoadKey('bz')[0,:,:]

            else:
                if self.GetPlotParam('show_x'):
                    self.fx = self.FigWrap.LoadKey('bx')[0,:,:]
                if self.GetPlotParam('show_y'):
                    self.fy = self.FigWrap.LoadKey('by')[0,:,:]
                if self.GetPlotParam('show_z'):
                    self.fz = self.FigWrap.LoadKey('bz')[0,:,:]

        if self.GetPlotParam('field_type') == 1: # Load the e-Field
            if self.GetPlotParam('twoD'):
                if self.GetPlotParam('show_x'):
                    self.zval = self.FigWrap.LoadKey('ex')[0,:,:]
                elif self.GetPlotParam('show_y'):
                    self.zval = self.FigWrap.LoadKey('ey')[0,:,:]
                elif self.GetPlotParam('show_z'):
                    self.zval = self.FigWrap.LoadKey('ez')[0,:,:]

            else:
                if self.GetPlotParam('show_x'):
                    self.fx = self.FigWrap.LoadKey('ex')[0,:,:]
                if self.GetPlotParam('show_y'):
                    self.fy = self.FigWrap.LoadKey('ey')[0,:,:]
                if self.GetPlotParam('show_z'):
                    self.fz = self.FigWrap.LoadKey('ez')[0,:,:]

        self.axes = self.figure.add_subplot(self.gs[18:92,:])
        if self.GetPlotParam('twoD'):
            self.cax = self.axes.imshow(self.zval, cmap = new_cmaps.cmaps[self.parent.cmap],  origin = 'lower', aspect = 'auto', interpolation=self.GetPlotParam('interpolation'))
            self.axes.set_axis_bgcolor('lightgrey')

            if self.GetPlotParam('show_cbar'):
                self.axC = self.figure.add_subplot(self.gs[:4,:])
                self.cbar = self.figure.colorbar(self.cax, ax = self.axes, cax = self.axC, orientation = 'horizontal')

                self.cbar.set_ticks(np.linspace(self.zval.min(),self.zval.max(), 5))

            self.axes.set_axis_bgcolor('lightgrey')
            self.axes.tick_params(labelsize = 10, color=tick_color)
#        self.axes.set_xlim(self.xmin,self.xmax)
            self.axes.set_xlabel(r'$x\ [c/\omega_{\rm pe}]$', labelpad = -2, color = 'black')
#        self.axes.set_ylabel(self.y_label, labelpad = -2, color = 'black')
        else:
            self.annotate_pos = [0.8,0.9]
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
            if self.GetPlotParam('show_x'):
                self.axes.plot(self.fx[self.fx.shape[0]/2,:], color = self.xcolor)
                if self.GetPlotParam('field_type') == 0:
                    tmp_str = r'$B_x$'
                else:
                    tmp_str = r'$E_x$'
                self.axes.annotate(tmp_str, xy = self.annotate_pos,
                                xycoords = 'axes fraction',
                                color = self.xcolor,
                                **self.annotate_kwargs)
                self.annotate_pos[0] += .08
            if self.GetPlotParam('show_y'):
                self.axes.plot(self.fy[self.fy.shape[0]/2,:], color = self.ycolor)
                if self.GetPlotParam('field_type') == 0:
                    tmp_str = r'$B_y$'
                else:
                    tmp_str = r'$E_y$'
                self.axes.annotate(tmp_str, xy = self.annotate_pos,
                                xycoords= 'axes fraction',
                                color = self.ycolor,
                                **self.annotate_kwargs)

                self.annotate_pos[0] += .08

            if self.GetPlotParam('show_z'):
                self.axes.plot(self.fz[self.fz.shape[0]/2,:], color = self.zcolor)
                if self.GetPlotParam('field_type') == 0:
                    tmp_str = r'$B_z$'
                else:
                    tmp_str = r'$E_z$'
                self.axes.annotate(tmp_str, xy = self.annotate_pos,
                                xycoords= 'axes fraction',
                                color = self.zcolor,
                                **self.annotate_kwargs
                                )


            self.axes.set_axis_bgcolor('lightgrey')


    def GetPlotParam(self, keyname):
        return self.FigWrap.GetPlotParam(keyname)

    def SetPlotParam(self, keyname, value):
        self.FigWrap.SetPlotParam(keyname, value)

    def OpenSettings(self):
        if self.settings_window is None:
            self.settings_window = FieldSettings(self)
        else:
            self.settings_window.destroy()
            self.settings_window = FieldSettings(self)


class FieldSettings(Tk.Toplevel):
    def __init__(self, parent):
        self.parent = parent
        Tk.Toplevel.__init__(self)

        self.wm_title('Phase Plot (%d,%d) Settings' % self.parent.FigWrap.pos)
        self.parent = parent
        frm = ttk.Frame(self)
        frm.pack(fill=Tk.BOTH, expand=True)
        self.protocol('WM_DELETE_WINDOW', self.OnClosing)
        #Create some sizers

        self.TwoDVar = Tk.IntVar(self) # Create a var to track whether or not to plot in 2-D
        self.TwoDVar.set(self.parent.GetPlotParam('twoD'))
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


        # the Radiobox Control to choose the Field Type
        self.FieldList = ['B Field', 'E field']
        self.FieldTypeVar  = Tk.IntVar()
        self.FieldTypeVar.set(self.parent.GetPlotParam('field_type'))

        ttk.Label(frm, text='Choose Field:').grid(row = 1, sticky = Tk.W)

        for i in range(len(self.FieldList)):
            ttk.Radiobutton(frm,
                text=self.FieldList[i],
                variable=self.FieldTypeVar,
                command = self.RadioField,
                value=i).grid(row = 2+i, sticky =Tk.W)

        # the Radiobox Control to choose the momentum dim
        self.dimList = ['x', 'y', 'z']
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


        # Control if the plot is weightedd
        self.WeightVar = Tk.IntVar()
        self.WeightVar.set(self.parent.GetPlotParam('weighted'))
        cb = ttk.Checkbutton(frm, text = "Weight by charge",
                        variable = self.WeightVar,
                        command = lambda:
                        self.parent.SetPlotParam('weighted', self.WeightVar.get()))
        cb.grid(row = 7, sticky = Tk.W)

        # Show energy integration region
        self.IntRegVar = Tk.IntVar()
        self.IntRegVar.set(self.parent.GetPlotParam('show_int_region'))
        cb = ttk.Checkbutton(frm, text = "Show Energy Region",
                        variable = self.IntRegVar,
                        command = lambda:
                        self.parent.SetPlotParam('show_int_region', self.IntRegVar.get()))
        cb.grid(row = 7, column = 1, sticky = Tk.W)

        # control mask
        self.MaskVar = Tk.IntVar()
        self.MaskVar.set(self.parent.GetPlotParam('masked'))
        cb = ttk.Checkbutton(frm, text = "Mask Zeros",
                        variable = self.MaskVar,
                        command = lambda:
                        self.parent.SetPlotParam('masked', self.MaskVar.get()))
        cb.grid(row = 8, sticky = Tk.W)


#        ttk.Label(frm, text = 'If the zero values are not masked they are set to z_min/2').grid(row =9, columnspan =2)
    # Define functions for the events


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
        if self.FieldTypeVar.get() == self.parent.GetPlotParam('field_type'):
            pass
        else:
            self.parent.SetPlotParam('field_type', self.FieldTypeVar.get())


    def RadioPrtl(self):
        if self.pvar.get() == self.parent.GetPlotParam('prtl_type'):
            pass
        else:
            self.parent.SetPlotParam('prtl_type', self.pvar.get())

    def RadioDim(self):
        if self.dimvar.get() == self.parent.GetPlotParam('mom_dim'):
            pass
        else:
            self.parent.SetPlotParam('mom_dim', self.dimvar.get())

    def RadioNorm(self):
        if self.cnormList[self.normvar.get()] == self.parent.GetPlotParam('norm_type'):
            pass
        else:
            self.parent.SetPlotParam('norm_type', self.cnormList[self.normvar.get()])

    def OnClosing(self):
        self.parent.settings_window = None
        self.destroy()
