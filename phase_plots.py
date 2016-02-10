#!/usr/bin/env pythonw
import Tkinter as Tk
import ttk as ttk
import matplotlib
import numpy as np
import numpy.ma as ma
import new_cmaps
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec

class PhasePanel:
    # A diction of all of the parameters for this plot with the default parameters

    plot_param_dict = {'mom_dim': 0,
                       'masked': 1,
                       'norm_type': 'LogNorm',
                       'prtl_type': 0,
                       'pow_num': 0.4,
                       'show_cbar': True,
                       'weighted': False,
                       'show_shock': True,
                       'show_int_region': True,
                       'set_color_limits': False,
                       'xbins' : 200,
                       'pbins' : 200,
                       'v_min': None,
                       'interpolation': 'hermite',
                       'v_max': None}

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


    def UpdatePowNum(self, value):
        self.SetPlotParam('pow_num', value)
        if self.GetPlotParam('norm_type') == "PowerNorm":
            self.draw()

    def ChangePlotType(self, str_arg):
        self.FigWrap.ChangeGraph(str_arg)

    def norm(self, vmin=None,vmax=None):
        if self.GetPlotParam('norm_type') =="Linear":
            return mcolors.Normalize(vmin, vmax)
        elif self.GetPlotParam('norm_type') == "LogNorm":
            return  mcolors.LogNorm(vmin, vmax)
        else:
            return  mcolors.PowerNorm(self.FigWrap.GetPlotParam('pow_num'), vmin, vmax)

    def set_plot_keys(self):
        '''A helper function that will insure that each hdf5 file will only be
        opened once per time step'''
        self.arrs_needed = ['c_omp', 'bx', 'istep']
        if self.GetPlotParam('prtl_type') == 0:
            self.arrs_needed.append('xi')
            if self.GetPlotParam('weighted'):
                self.arrs_needed.append('chi')
            if self.GetPlotParam('mom_dim') == 0:
                self.arrs_needed.append('ui')
            elif self.GetPlotParam('mom_dim') == 1:
                self.arrs_needed.append('vi')
            elif self.GetPlotParam('mom_dim') == 2:
                self.arrs_needed.append('wi')

        if self.GetPlotParam('prtl_type') == 1:
            self.arrs_needed.append('xe')
            if self.GetPlotParam('weighted'):
                self.arrs_needed.append('che')

            elif self.GetPlotParam('mom_dim') == 0:
                self.arrs_needed.append('ue')
            elif self.GetPlotParam('mom_dim') == 1:
                self.arrs_needed.append('ve')
            elif self.GetPlotParam('mom_dim') == 2:
                self.arrs_needed.append('we')
        return self.arrs_needed

    def draw(self):
        # Choose the normalization

        # Generate the X-axis values
        self.c_omp = self.FigWrap.LoadKey('c_omp')[0]
        self.weights = None
        # Choose the particle type and px, py, or pz
        if self.GetPlotParam('prtl_type') == 0:
            self.x_values = self.FigWrap.LoadKey('xi')/self.c_omp
            if self.GetPlotParam('weighted'):
                self.weights = self.FigWrap.LoadKey('chi')
            if self.GetPlotParam('mom_dim') == 0:
                self.y_values = self.FigWrap.LoadKey('ui')
                self.y_label  = r'$P_{px}\ [c]$'
            if self.GetPlotParam('mom_dim') == 1:
                self.y_values = self.FigWrap.LoadKey('vi')
                self.y_label  = r'$P_{py}\ [c]$'
            if self.GetPlotParam('mom_dim') == 2:
                self.y_values = self.FigWrap.LoadKey('wi')
                self.y_label  = r'$P_{pz}\ [c]$'
        if self.GetPlotParam('prtl_type') == 1:
            self.x_values = self.FigWrap.LoadKey('xe')/self.c_omp
            if self.GetPlotParam('weighted'):
                self.weights = self.FigWrap.LoadKey('che')

            if self.GetPlotParam('mom_dim') == 0:
                self.y_values = self.FigWrap.LoadKey('ue')
                self.y_label  = r'$P_{ex}\ [c]$'
            if self.GetPlotParam('mom_dim') == 1:
                self.y_values = self.FigWrap.LoadKey('ve')
                self.y_label  = r'$P_{ey}\ [c]$'
            if self.GetPlotParam('mom_dim') == 2:
                self.y_values = self.FigWrap.LoadKey('we')
                self.y_label  = r'$P_{ez}\ [c]$'

        self.pmin = min(self.y_values)
        self.pmax = max(self.y_values)
        self.xmin = 0
        self.xmax = self.FigWrap.LoadKey('bx').shape[2]/self.c_omp*self.FigWrap.LoadKey('istep')[0]
        self.hist2d = np.histogram2d(self.y_values, self.x_values, bins = [self.GetPlotParam('pbins'), self.GetPlotParam('xbins')], range = [[self.pmin,self.pmax],[0,192.2]], weights = self.weights)

        self.zval = ma.masked_array(self.hist2d[0])
        tick_color = 'white'
        if self.GetPlotParam('masked'):
            self.zval[self.zval == 0] = ma.masked
            tick_color = 'k'
        else:
            self.zval[self.zval==0] = 1
        self.zval *= self.zval.max()**(-1)
        self.gs = gridspec.GridSpecFromSubplotSpec(100,100, subplot_spec = self.parent.gs0[self.FigWrap.pos])#, bottom=0.2,left=0.1,right=0.95, top = 0.95)

        self.axes = self.figure.add_subplot(self.gs[18:92,:])

        self.cax = self.axes.imshow(self.zval, cmap = new_cmaps.cmaps[self.parent.cmap], norm = self.norm(), origin = 'lower', aspect = 'auto', extent=[self.xmin,self.xmax,self.hist2d[1][-1],self.hist2d[1][0]], interpolation=self.GetPlotParam('interpolation'))

        if self.GetPlotParam('show_cbar'):
            self.axC = self.figure.add_subplot(self.gs[:4,:])
            self.cbar = self.figure.colorbar(self.cax, ax = self.axes, cax = self.axC, orientation = 'horizontal')
            if self.GetPlotParam('norm_type')== 'PowerNorm':
                self.cbar.set_ticks(np.linspace(self.zval.min(),self.zval.max(), 5)**(1./self.FigWrap.GetPlotParam('pow_num')))

            if self.GetPlotParam('norm_type') == 'LogNorm':

                ctick_range = np.logspace(np.log10(self.zval.min()),np.log10(self.zval.max()), 5)
                self.cbar.set_ticks(ctick_range)
                ctick_labels = []
                for elm in ctick_range:
                    if np.abs(np.log10(elm))<1E-2:
                        tmp_s = '0'
                    else:
                        tmp_s = '%.2f' % np.log10(elm)
                    ctick_labels.append(tmp_s)

                self.cbar.set_ticklabels(ctick_labels)
                self.cbar.ax.tick_params(labelsize=10)
            if self.GetPlotParam('norm_type')== 'Linear':
                self.cbar.set_ticks(np.linspace(self.zval.min(),self.zval.max(), 5))


        self.axes.set_axis_bgcolor('lightgrey')
        self.axes.tick_params(labelsize = 10, color=tick_color)
        self.axes.set_xlim(self.xmin,self.xmax)
        self.axes.set_xlabel(r'$x\ [c/\omega_{\rm pe}]$', labelpad = -2, color = 'black')
        self.axes.set_ylabel(self.y_label, labelpad = -2, color = 'black')


    def GetPlotParam(self, keyname):
        return self.FigWrap.GetPlotParam(keyname)

    def SetPlotParam(self, keyname, value):
        self.FigWrap.SetPlotParam(keyname, value)

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
        #Create some sizers

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

        ''' Commenting out some lines that allows you to change the color norm,
        No longer needed
        self.cnormList = ['Linear', 'LogNorm', 'PowerNorm']
        self.normvar = Tk.IntVar()
        self.normvar.set(self.cnormList.index(self.parent.GetPlotParam('norm_type')))

        ttk.Label(frm, text = 'Cmap Norm:').grid(row = 6, sticky = Tk.W)

        for i in range(3):
            ttk.Radiobutton(frm,
                            text = self.cnormList[i],
                            variable = self.normvar,
                            command = self.RadioNorm,
                            value = i).grid(row = 7+i, sticky = Tk.W)

        '''
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
