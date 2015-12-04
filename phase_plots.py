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
    plot_param_dict = {'mom_dim': 0, 'masked': 0, 'norm_type': 'LogNorm', 'prtl_type': 0, 'pow_num': 0.4, 'show_cbar': True, 'weighted': False}
    def __init__(self, parent, figwrapper):
        self.FigWrap = figwrapper
        self.parent = parent
        self.ChartTypes = self.FigWrap.PlotTypeDict.keys()
        self.chartType = self.FigWrap.chartType
        self.subplotlist = []
        self.figure = self.FigWrap.figure


    def sizeHandler(self, *args, **kwargs):
        '''Make it so the plot scales with resizing of the window'''
        self.canvas.SetSize(self.GetSize())

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
        self.arrs_needed = ['c_omp']
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
            if self.GetPlotParam('mom_dim') == 1:
                self.y_values = self.FigWrap.LoadKey('vi')
            if self.GetPlotParam('mom_dim') == 2:
                self.y_values = self.FigWrap.LoadKey('wi')

        if self.GetPlotParam('prtl_type') == 1:
            self.x_values = self.FigWrap.LoadKey('xe')/self.c_omp
            if self.GetPlotParam('weighted'):
                self.weights = self.FigWrap.LoadKey('che')

            if self.GetPlotParam('mom_dim') == 0:
                self.y_values = self.FigWrap.LoadKey('ue')
            if self.GetPlotParam('mom_dim') == 1:
                self.y_values = self.FigWrap.LoadKey('ve')
            if self.GetPlotParam('mom_dim') == 2:
                self.y_values = self.FigWrap.LoadKey('we')

        self.hist2d = np.histogram2d(self.y_values, self.x_values, bins = [200,200], weights = self.weights)
        self.zval = ma.masked_array(self.hist2d[0])

        if self.GetPlotParam('masked'):
            self.zval[self.zval == 0] += ma.masked
        else:
            self.zval[self.zval==0] = 0.5
        self.zval *= self.zval.max()**(-1)
        self.gs = gridspec.GridSpecFromSubplotSpec(100,100, subplot_spec = self.parent.gs0[self.FigWrap.pos])#, bottom=0.2,left=0.1,right=0.95, top = 0.95)

        if self.GetPlotParam('show_cbar'):
            self.axes = self.figure.add_subplot(self.gs[20:,:])
            self.axC = self.figure.add_subplot(self.gs[:5,:])
            self.cax = self.axes.pcolormesh(self.hist2d[2], self.hist2d[1], self.zval, cmap = new_cmaps.cmaps[self.parent.cmap], norm = self.norm())
            self.axes.set_xlim(self.hist2d[2].min(), self.hist2d[2].max())
            self.axes.set_axis_bgcolor('lightgrey')

            self.axes.set_ylim(self.hist2d[1].min(), self.hist2d[1].max())
            self.cbar = self.figure.colorbar(self.cax, ax = self.axes, cax = self.axC, orientation = 'horizontal')
            if self.GetPlotParam('norm_type')== 'PowerNorm':
                self.cbar.set_ticks(np.linspace(self.zval.min(),self.zval.max(), 5)**(1./self.FigWrap.GetPlotParam('pow_num')))

            if self.GetPlotParam('norm_type') == 'LogNorm':

                self.cbar.set_ticks(np.logspace(np.log10(self.zval.min()), 0, 5))
#                print tuple(np.logspace(np.log10(self.zval.min()), 0, 5))

#                self.cbar.set_ticklabels(str.split('%f, %f, %f, %f, %f', ',') % tuple(np.logspace(np.log10(self.zval.min()), 0, 5)))
            if self.GetPlotParam('norm_type')== 'Linear':
                self.cbar.set_ticks(np.linspace(self.zval.min(),self.zval.max(), 5))

            self.axes.set_xlabel(r'$x/\omega_{\rm pe}$')
        else:
            self.axes = self.figure.add_subplot(self.gs[5:,:])
            self.cax = self.axes.hist2d(self.x_values,self.y_values, bins = [200,200],cmap = new_cmaps.cmaps[self.Parent.cmap], norm = self.norm(1))
            self.axes.set_xlabel(r'$x/\omega_{\rm pe}$')


    def GetPlotParam(self, keyname):
        return self.FigWrap.GetPlotParam(keyname)

    def SetPlotParam(self, keyname, value):
        self.FigWrap.SetPlotParam(keyname, value)

    def OpenSettings(self):
        PhaseSettings(self)


class PhaseSettings(Tk.Toplevel):
    def __init__(self, parent):
        self.parent = parent
        Tk.Toplevel.__init__(self)

        self.wm_title('Phase Plot (%d,%d) Settings' % self.parent.FigWrap.pos)
        self.parent = parent
        frm = ttk.Frame(self)
        frm.pack(fill=Tk.BOTH, expand=True)

        #Create some sizers

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
        self.cbColorbar = wx.CheckBox(self, -1, "Show Cbar")
        self.cbColorbar.SetValue(self.parent.GetPlotParam('show_cbar'))
        self.Bind(wx.EVT_CHECKBOX, self.EvtCheckCbar, self.cbColorbar)
        grid.Add(self.cbColorbar, pos=(3,0))

        #Adding WEIGHT *TO DO*
        self.cbWeight = wx.CheckBox(self, -1, "Weight")
        self.cbWeight.SetValue(self.parent.GetPlotParam('weighted'))
        self.Bind(wx.EVT_CHECKBOX, self.EvtCheckWeight, self.cbWeight)
        grid.Add(self.cbWeight, pos = (3,1))
        '''

    # Define functions for the events

    def ctypeChanged(self, *args):
        if self.ctypevar.get() ==self.parent.chartType:
            pass
        else:
            self.parent.ChangePlotType(self.ctypevar.get())
            self.Destroy()


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


    def EvtCheckCbar(self, evt):
        self.parent.SetPlotParam('show_cbar', evt.IsChecked())
        self.parent.draw()

    def EvtCheckWeight(self, evt):
        self.parent.SetPlotParam('weighted', evt.IsChecked())
        self.parent.draw()

    def EvtRadioDim(self, evt):
        self.parent.SetPlotParam('mom_dim', evt.GetInt())
        self.parent.draw()

    def EvtRadioNorm(self, evt):
        self.parent.SetPlotParam('norm_type', evt.GetEventObject().GetLabel())
        self.parent.draw()

    def EvtNormNumSet(self, evt):
        tmp_num = evt.GetString()
        if not tmp_num:
            tmp_num = 1E-2
        elif float(tmp_num)<=1E-2:
            tmp_num = 1E-2
        else:
            tmp_num = float(tmp_num)
        self.parent.SetPlotParam('pow_num', tmp_num)
        if self.parent.GetPlotParam('norm_type') == "PowerNorm":
            self.parent.draw()


    def OnCloseMe(self, event):
        self.Close(True)

    def OnCloseWindow(self, event):
        self.Destroy()
