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
from mpl_toolkits.mplot3d import Axes3D

class ThreeDBPanel:
    # A dictionary of all of the parameters for this plot with the default parameters

    plot_param_dict = {'twoD': 0,
                       'show_shock' : False,
                       'spatial_x': None,
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

        ### Commenting out this because loading the HDF5 file is more expensive
        ### than just storing it in RAM. Therefore I should just load everything
        '''
        self.arrs_needed = ['c_omp', 'istep']
        # Then see if we are plotting E-field or B-Field
        if self.GetPlotParam('field_type') == 0: # Load the B-Field
            if self.GetPlotParam('show_x'):
                self.arrs_needed.append('bx')
            if self.GetPlotParam('show_y'):
                self.arrs_needed.append('by')
            if self.GetPlotParam('show_z'):
                self.arrs_needed.append('bz')

        if self.GetPlotParam('field_type') == 1: # Load the E-Field
            if self.GetPlotParam('show_x'):
                self.arrs_needed.append('ex')
            if self.GetPlotParam('show_y'):
                self.arrs_needed.append('ey')
            if self.GetPlotParam('show_z'):
                self.arrs_needed.append('ez')
        '''
        self.arrs_needed = ['c_omp', 'istep', 'bx', 'by', 'bz']
        return self.arrs_needed

    def LoadData(self):
        ''' A Helper function that loads the data for the plot'''
        # First see of the x_axis and y_axis values have already been calculated
        # and stored in the DataDict for this time step
        self.c_omp = self.FigWrap.LoadKey('c_omp')[0]
        self.istep = self.FigWrap.LoadKey('istep')[0]
        self.mag_color = new_cmaps.cmaps[self.parent.cmap](0.5)


        # see if the axis values are saved in the data dict
        if 'xaxis_values' in self.parent.DataDict.keys():
            self.xaxis_values = self.parent.DataDict['xaxis_values']
        else:
            # x-values haven't been calculated yet, generate them then save them to the dictionary for later.
            self.xaxis_values = np.arange(self.FigWrap.LoadKey('bx')[0,:,:].shape[1])/self.c_omp*self.istep
        #            print self.xaxis_values
            self.parent.DataDict['xaxis_values'] = np.copy(self.xaxis_values)

        if 'delta_b_perp' in self.parent.DataDict.keys():
            self.fy = self.parent.DataDict['delta_b_perp']
        elif np.isnan(self.parent.btheta):
            print 'hi'
            self.fy = 1.0
            self.fz = 1.0
        else:
            print self.parent.bx0
            if np.abs(self.parent.bz0) <= 1E-6:
                # Decompose the perpendicular components of the magnetic fields into two parts
                # one is bz, the other is the in plane components

                sign = np.sign(self.FigWrap.LoadKey('bx')[0,:,:]*np.tan(self.parent.btheta)-self.FigWrap.LoadKey('by')[0,:,:])
                delta_bx_perp = (self.FigWrap.LoadKey('by')[0,:,:]-self.parent.bx0)*np.sin(self.parent.btheta)
                delta_by_perp = (self.FigWrap.LoadKey('by')[0,:,:]-self.parent.by0)*np.cos(self.parent.btheta)

                self.fy = sign*np.sqrt(delta_bx_perp**2+delta_by_perp**2)
                self.parent.DataDict['delta_b_perp'] = self.fy
        self.fz = self.FigWrap.LoadKey('bz')[0,:,:]/self.parent.b0
        self.oneDslice = self.fz.shape[0]/2

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

        # Now that the data is loaded, start making the plots
        self.axes = self.figure.add_subplot(self.gs[18:92,:], projection = '3d')

        shock_ind = self.xaxis_values.searchsorted(self.parent.shock_loc)

        self.line = self.axes.plot(self.xaxis_values[shock_ind:], self.fy[self.oneDslice,shock_ind:], self.fz[self.oneDslice,shock_ind:], color = self.mag_color)
        self.axes.set_xlim(self.parent.shock_loc, None)
    def refresh(self):

        '''This is a function that will be called only if self.axes already
        holds a fields type plot. We only update things that have changed & are
        shown.  If hasn't changed or isn't shown, don't touch it. The difference
        between this and last time, is that we won't actually do any drawing in
        the plot. The plot will be redrawn after all subplots are refreshed. '''


        # Main goal, only change what is showing..
        # First do the 1D plots, because it is simpler
        shock_ind = self.xaxis_values.searchsorted(self.parent.shock_loc)
        self.line[0].set_data(self.xaxis_values[shock_ind:], self.fy[self.oneDslice,shock_ind:])
        self.line[0].set_3d_properties(self.fz[self.oneDslice,shock_ind:])
        self.axes.set_xlim(self.parent.shock_loc, None)

            #self.axes.draw_artist(self.axes.patch)
            #self.axes.draw_artist(self.cax)
            #self.axes.draw_artist(self.axes.xaxis)


    def GetPlotParam(self, keyname):
        return self.FigWrap.GetPlotParam(keyname)

    def SetPlotParam(self, keyname, value, update_plot = True):
        self.FigWrap.SetPlotParam(keyname, value, update_plot = update_plot)

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

        self.wm_title('ThetaB (%d,%d) Settings' % self.parent.FigWrap.pos)
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
        self.MagList = ['Theta B', 'Delta B/B_0']
        self.MagTypeVar  = Tk.IntVar()
        self.MagTypeVar.set(self.parent.GetPlotParam('mag_plot_type'))

        ttk.Label(frm, text='Choose Quantity:').grid(row = 2, sticky = Tk.W)

        for i in range(len(self.MagList)):
            ttk.Radiobutton(frm,
                text=self.MagList[i],
                variable=self.MagTypeVar,
                command = self.RadioMag,
                value=i).grid(row = 3+i, sticky =Tk.W)

        # Control whether or not Cbar is shown
        self.CbarVar = Tk.IntVar()
        self.CbarVar.set(self.parent.GetPlotParam('show_cbar'))
        cb = ttk.Checkbutton(frm, text = "Show Color bar",
                        variable = self.CbarVar,
                        command = self.CbarHandler)
        cb.grid(row = 6, sticky = Tk.W)

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


        cb = ttk.Checkbutton(frm, text ='Set Theta or delta B min',
                        variable = self.setZminVar)
        cb.grid(row = 3, column = 2, sticky = Tk.W)
        self.ZminEnter = ttk.Entry(frm, textvariable=self.Zmin, width=7)
        self.ZminEnter.grid(row = 3, column = 3)

        cb = ttk.Checkbutton(frm, text ='Set Theta or delta B max',
                        variable = self.setZmaxVar)
        cb.grid(row = 4, column = 2, sticky = Tk.W)

        self.ZmaxEnter = ttk.Entry(frm, textvariable=self.Zmax, width=7)
        self.ZmaxEnter.grid(row = 4, column = 3)

        self.ShockVar = Tk.IntVar()
        self.ShockVar.set(self.parent.GetPlotParam('show_shock'))
        cb = ttk.Checkbutton(frm, text = "Show Shock",
                        variable = self.ShockVar,
                        command = self.ShockVarHandler)
        cb.grid(row = 6, column = 1, sticky = Tk.W)

    def CbarHandler(self, *args):
        if self.parent.GetPlotParam('show_cbar')== self.CbarVar.get():
            pass
        else:
            if self.parent.GetPlotParam('twoD'):
                self.parent.axC.set_visible(self.CbarVar.get())

            self.parent.SetPlotParam('show_cbar', self.CbarVar.get(), update_plot =self.parent.GetPlotParam('twoD'))

    def ShockVarHandler(self, *args):
        if self.parent.GetPlotParam('show_shock')== self.ShockVar.get():
            pass
        else:
            if self.parent.GetPlotParam('twoD'):
                self.parent.shockline_2d.set_visible(self.ShockVar.get())
            else:
                self.parent.shock_line.set_visible(self.ShockVar.get())

            self.parent.SetPlotParam('show_shock', self.ShockVar.get())

    def Change2d(self):

        if self.TwoDVar.get() == self.parent.GetPlotParam('twoD'):
            pass
        else:
            self.parent.SetPlotParam('spatial_y', self.TwoDVar.get(), update_plot=False)
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
            self.parent.cax.set_interpolation(self.InterpolVar.get())
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

    def RadioMag(self):
        if self.MagTypeVar.get() == self.parent.GetPlotParam('mag_plot_type'):
            pass
        else:
            if self.parent.GetPlotParam('twoD') == 0:
                if self.MagTypeVar.get() == 0:
                    self.parent.axes.set_ylabel(r'$\theta_B$', labelpad = self.parent.parent.ylabel_pad, color = 'black')

                elif np.isnan(self.parent.parent.btheta):
                    self.parent.axes.set_ylabel(r'$|\delta B|$', labelpad = self.parent.parent.ylabel_pad, color = 'black')
                else:
                    self.parent.axes.set_ylabel(r'$|\delta B|/B_0$', labelpad = self.parent.parent.ylabel_pad, color = 'black')

            self.parent.SetPlotParam('mag_plot_type', self.MagTypeVar.get())


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


    def OnClosing(self):
        self.parent.settings_window = None
        self.destroy()
