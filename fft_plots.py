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

class FFTPanel:
    # A dictionary of all of the parameters for this plot with the default parameters

    plot_param_dict = {'twoD': 0,
                       'FFT_type':0, # 0 = FFT_Bz, 1 = FFT_perp_in_plane, 2 = Chi -> sin(2*Chi) = V/I where V & I are the Stokes parameters
                       'y_min': 0,
                       'y_max' : 10,
                       'set_y_min': False,
                       'set_y_max': False,
                       'xLog': True,
                       'yLog': True,
                       'spatial_x': False,
                       'spatial_y': False}

    # We need the types of all the parameters for the config file
    BoolList = ['twoD', 'set_y_min', 'set_y_max',
                'xLog', 'yLog', 'spatial_x', 'spatial_y']
    IntList = ['FFT_type']
    FloatList = ['y_min', 'y_max']
    StrList = []
    
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
        self.ylabel_list = [r'$|\delta B_z(k)|/B_0$', r'$|\delta B_\perp(k)|/B_0$', r'$|E_x(k)|/E_0$', r'$\chi(k)\ [\ \! ^\circ]$']
    def ChangePlotType(self, str_arg):
        self.FigWrap.ChangeGraph(str_arg)

    def set_plot_keys(self):
        '''A helper function that will insure that each hdf5 file will only be
        opened once per time step'''
        # First make sure that omega_plasma & xi is loaded so we can fix the
        # distances. We need all 3 magnetic field directions to calculate the FFT
        # or Chi
        self.arrs_needed = ['c_omp', 'istep', 'bx', 'by', 'bz', 'ex']
        return self.arrs_needed

    def LoadData(self):
        ''' A Helper function that loads the data for the plot'''
        # Load c_omp, istep, and the line color
        self.c_omp = self.FigWrap.LoadKey('c_omp')[0]
        self.istep = self.FigWrap.LoadKey('istep')[0]
        self.line_color = new_cmaps.cmaps[self.parent.MainParamDict['ColorMap']](0.5)

        # Load the x-axis
        if 'xaxis_values' in self.parent.DataDict.keys():
            self.xaxis_values = self.parent.DataDict['xaxis_values']
        else:
            # x-values haven't been calculated yet, generate them then save them to the dictionary for later.
            self.xaxis_values = np.arange(self.FigWrap.LoadKey('bx')[0,:,:].shape[1])/self.c_omp*self.istep
            self.parent.DataDict['xaxis_values'] = np.copy(self.xaxis_values)

        # Load the left & right locations
        self.left_loc = self.parent.MainParamDict['FFTLeft'] + self.parent.shock_loc*self.parent.MainParamDict['FFTRelative']
        self.left_loc = max(self.left_loc, self.xaxis_values[0])

        self.right_loc = self.parent.MainParamDict['FFTRight'] + self.parent.shock_loc*self.parent.MainParamDict['FFTRelative']
        self.right_loc = min(self.right_loc, self.xaxis_values[-1])


        # A list that will make sure that the data has the same int region
        self.region_args = [self.parent.MainParamDict['FFTRelative'], self.parent.MainParamDict['FFTLeft'], self.parent.MainParamDict['FFTRight']]
        # Check if the region is the same in the DataDict
        is_loaded = False
        if 'FFTs' in self.parent.DataDict.keys():
            # Check that the save value has the same region for the FFT
            if np.all(self.region_args == self.parent.DataDict['FFTs'][-1]):
                is_loaded = True
                self.k_axis, self.BzFFT, self.BperpFFT, self.ExFFT, self.StokesChi, self.all_min_max, self.klims, region_args = self.parent.DataDict['FFTs']
                self.min_max = self.all_min_max[self.GetPlotParam('FFT_type')]
        if not is_loaded:
            ####
            # First find the left & right indices
            ####
            iL = self.xaxis_values.searchsorted(self.left_loc)
            iR = self.xaxis_values.searchsorted(self.right_loc)
            self.all_min_max = []
            # Calculate K_axis
            self.k_axis = np.arange(iR-iL)*(2*np.pi/(self.xaxis_values[1]-self.xaxis_values[0]))/(iR-iL-1)-(2*np.pi/(self.xaxis_values[1]-self.xaxis_values[0]))*.5
            # Calculate all of the FFTs, just simpler to do it this way...
            bz = self.FigWrap.LoadKey('bz')[0,:,:]
            self.oneDslice = bz.shape[0]/2
            self.BzFFT = np.fft.fft(bz[self.oneDslice,iL:iR]*self.parent.b0**(-1.0))
            # Shift the fft so it is centered
            self.BzFFT = np.fft.fftshift(self.BzFFT)
            self.all_min_max.append(self.LimFinder(np.abs(self.BzFFT)))

            bx = self.FigWrap.LoadKey('bx')[0,:,:]
            by = self.FigWrap.LoadKey('by')[0,:,:]
            b_perp_in_plane = (by*np.cos(self.parent.btheta)-bx*np.sin(self.parent.btheta))*self.parent.b0**(-1.0)
            self.BperpFFT = np.fft.fft(b_perp_in_plane[self.oneDslice,iL:iR])
            # Shift the fft so it is centered
            self.BperpFFT = np.fft.fftshift(self.BperpFFT)

            self.all_min_max.append(self.LimFinder(np.abs(self.BperpFFT)))

	    ex = self.FigWrap.LoadKey('ex')[0,:,:]
	    self.ExFFT = np.fft.fft(ex[self.oneDslice,iL:iR]*self.parent.e0**(-1.0))
 	    self.ExFFT = np.fft.fftshift(self.ExFFT)
	    self.all_min_max.append(self.LimFinder(np.abs(self.ExFFT)))

            # Calculate Stokes I
            I = self.BperpFFT*np.conjugate(self.BperpFFT)
            I += self.BzFFT*np.conjugate(self.BzFFT)

            # Calculate Stokes V
            V = self.BperpFFT*np.conjugate(self.BzFFT)
            V += -self.BzFFT*np.conjugate(self.BperpFFT)
            V *= -1j

            # Calculate StokesChi
            self.StokesChi = 0.5 * np.rad2deg(np.arcsin(V.real/I.real))
            self.all_min_max.append(self.LimFinder(self.StokesChi))


            max_place, = np.where(np.abs(self.BzFFT)>=np.abs(self.BzFFT).max())
            kmaxval = np.abs(self.k_axis[max_place[0]])

#            self.klims = [[0, 5*kmaxval], [kmaxval/5, 50*kmaxval]]
            self.klims = [[0, 5*kmaxval], [kmaxval/5., 50*kmaxval]]
            self.parent.DataDict['FFTs'] = self.k_axis, self.BzFFT, self.BperpFFT, self.ExFFT, self.StokesChi, self.all_min_max, self.klims, self.region_args

        # figure out the y lims for the plot
        self.min_max = self.all_min_max[self.GetPlotParam('FFT_type')]

        # Get the labels, get the y values
        if self.GetPlotParam('FFT_type') == 0:
            self.y = np.abs(self.BzFFT)

        if self.GetPlotParam('FFT_type') == 1:
            self.y = np.abs(self.BperpFFT)

        if self.GetPlotParam('FFT_type') == 2:
            self.y = np.abs(self.ExFFT)

        if self.GetPlotParam('FFT_type') == 3:
            self.y = self.StokesChi

        self.ylabel = self.ylabel_list[self.GetPlotParam('FFT_type')]


    def LimFinder(self, arr):
        oneD_lims = [arr.min(), arr.max()]
        # now give it a bit of spacing, a 4% percent difference of the distance
        dist = oneD_lims[1]-oneD_lims[0]
        oneD_lims[0] -= 0.04*dist
        oneD_lims[1] += 0.04*dist
        return oneD_lims

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


        # See if there are other k-axes, and share with them
        if self.parent.MainParamDict['LinkK'] and self.FigWrap.pos != self.parent.first_k:
            self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]],
                    sharex = self.parent.SubPlotList[self.parent.first_k[0]][self.parent.first_k[1]].graph.axes)
        else:
            self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]])


        self.line = self.axes.plot(self.k_axis, self.y, color = self.line_color)
        self.axes.set_axis_bgcolor('lightgrey')
        self.axes.tick_params(labelsize = self.parent.MainParamDict['NumFontSize'], color=tick_color)

        if self.parent.MainParamDict['SetkLim']:
            self.axes.set_xlim(self.parent.MainParamDict['kLeft'],self.parent.MainParamDict['kRight'])
        elif self.GetPlotParam('xLog'):
            self.axes.set_xscale('log')
            self.axes.set_xlim(self.klims[1])
        else:
            self.axes.set_xlim(self.klims[0])

        if self.GetPlotParam('yLog'):
            self.axes.set_yscale('log')

        if self.GetPlotParam('set_y_min'):
            self.axes.set_ylim(ymin = self.GetPlotParam('y_min'))
        if self.GetPlotParam('set_y_max'):
            self.axes.set_ylim(ymax = self.GetPlotParam('y_max'))

        self.axes.set_xlabel(r'$k \ [\omega_{pe}/c]$', labelpad = self.parent.MainParamDict['xLabelPad'], color = 'black')
        self.axes.set_ylabel(self.ylabel, labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black')

    def refresh(self):

        '''This is a function that will be called only if self.axes already
        holds a fields type plot. We only update things that have changed & are
        shown.  If hasn't changed or isn't shown, don't touch it. The difference
        between this and last time, is that we won't actually do any drawing in
        the plot. The plot will be redrawn after all subplots are refreshed. '''


        # Main goal, only change what is showing..

        self.line[0].set_data(self.k_axis, self.y)
        tmp_lims = np.copy(self.min_max)
        if self.GetPlotParam('set_y_min'):
            tmp_lims[0]=self.GetPlotParam('y_min')
        if self.GetPlotParam('set_y_max'):
            tmp_lims[1]=self.GetPlotParam('y_max')

        self.axes.set_ylim(tmp_lims)
        if self.parent.MainParamDict['SetkLim']:
            self.axes.set_xlim(self.parent.MainParamDict['kLeft'],self.parent.MainParamDict['kRight'])
        elif self.GetPlotParam('xLog'):
            self.axes.set_xlim(self.klims[1])
        else:
            self.axes.set_xlim(self.klims[0])



    def GetPlotParam(self, keyname):
        return self.FigWrap.GetPlotParam(keyname)

    def SetPlotParam(self, keyname, value, update_plot = True):
        self.FigWrap.SetPlotParam(keyname, value, update_plot = update_plot)

    def OpenSettings(self):
        if self.settings_window is None:
            self.settings_window = FFTSettings(self)
        else:
            self.settings_window.destroy()
            self.settings_window = FFTSettings(self)


class FFTSettings(Tk.Toplevel):
    def __init__(self, parent):
        self.parent = parent
        Tk.Toplevel.__init__(self)

        self.wm_title('FFT (%d,%d) Settings' % self.parent.FigWrap.pos)
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


        # the Radiobox Control to choose the Field Type
        self.FFTList = ['FFT_Bz', 'FFT_perp_in_plane', 'FFT_Ex', 'Chi']
        self.FFTTypeVar  = Tk.IntVar()
        self.FFTTypeVar.set(self.parent.GetPlotParam('FFT_type'))

        ttk.Label(frm, text='Choose FFT Plot:').grid(row = 2, sticky = Tk.W)

        for i in range(len(self.FFTList)):
            ttk.Radiobutton(frm,
                text=self.FFTList[i],
                variable=self.FFTTypeVar,
                command = self.RadioFFT,
                value=i).grid(row = 3+i, sticky =Tk.W)

        # Now the field lim
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
        cb.grid(row = 3, column = 2, sticky = Tk.W)
        self.ZminEnter = ttk.Entry(frm, textvariable=self.Zmin, width=7)
        self.ZminEnter.grid(row = 3, column = 3)

        cb = ttk.Checkbutton(frm, text ='Set y max',
                        variable = self.setZmaxVar)
        cb.grid(row = 4, column = 2, sticky = Tk.W)

        self.ZmaxEnter = ttk.Entry(frm, textvariable=self.Zmax, width=7)
        self.ZmaxEnter.grid(row = 4, column = 3)

        # Now whether or not the x and y axes should be in logspace
        self.xLogVar = Tk.IntVar()
        self.xLogVar.set(self.parent.GetPlotParam('xLog'))
        self.xLogVar.trace('w', self.xLogChanged)

        self.yLogVar = Tk.IntVar()
        self.yLogVar.set(self.parent.GetPlotParam('yLog'))
        self.yLogVar.trace('w', self.yLogChanged)


        cb = ttk.Checkbutton(frm, text ='k-axis logscale',
                        variable = self.xLogVar)
        cb.grid(row = 7, column = 0, sticky = Tk.W)

        cb = ttk.Checkbutton(frm, text ='y-axis logscale',
                        variable = self.yLogVar)
        cb.grid(row = 8, column = 0, sticky = Tk.W)





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


    def xLogChanged(self, *args):
        if self.xLogVar.get() == self.parent.GetPlotParam('xLog'):
            pass
        else:
            if self.xLogVar.get():
                self.parent.axes.set_xscale('log')
            else:
                self.parent.axes.set_xscale('linear')

            self.parent.SetPlotParam('xLog', self.xLogVar.get())

    def yLogChanged(self, *args):
        if self.yLogVar.get() == self.parent.GetPlotParam('yLog'):
            pass
        else:
            if self.yLogVar.get():
                self.parent.axes.set_yscale('log')
            else:
                self.parent.axes.set_yscale('linear')

            self.parent.SetPlotParam('yLog', self.yLogVar.get())


    def RadioFFT(self):
        if self.FFTTypeVar.get() == self.parent.GetPlotParam('FFT_type'):
            pass
        else:
            # Change the labels
            self.parent.ylabel = self.parent.ylabel_list[self.FFTTypeVar.get()]

            if self.FFTTypeVar.get() == 2:
                self.parent.SetPlotParam('yLog', False, update_plot = False)
                self.parent.axes.set_yscale('linear')
                self.yLogVar.set(False)

            self.parent.axes.set_ylabel(self.parent.ylabel)
            self.parent.SetPlotParam('FFT_type', self.FFTTypeVar.get(), update_plot = True)


    def TxtEnter(self, e):
        self.FieldsCallback()

    def FieldsCallback(self):
        tkvarLimList = [self.Zmin, self.Zmax]
        plot_param_List = ['y_min', 'y_max']
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
            self.parent.SetPlotParam('y_min', self.parent.GetPlotParam('y_min'))


    def OnClosing(self):
        self.parent.settings_window = None
        self.destroy()
