#!/usr/bin/env python
import Tkinter as Tk
import ttk as ttk
import matplotlib
import numpy as np
import numpy.ma as ma
import new_cmaps
from new_cnorms import PowerNormWithNeg, PowerNormFunc
from modest_image import imshow
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as PathEffects

class BPanel:
    # A dictionary of all of the parameters for this plot with the default parameters

    plot_param_dict = {'twoD': 0,
                       'mag_plot_type':0, # 0 = theta_b, 2 = |deltaB|/B0, 3 = |deltaB_perp|/B0, 4 = |deltaB_para|/b0
                       'show_cbar': True,
                       'v_min': 0,
                       'v_max' : 10,
                       'set_v_min': False,
                       'set_v_max': False,
                       'show_shock' : False,
                       'OutlineText': True,
                       'spatial_x': True,
                       'spatial_y': False,
                       'show_FFT_region': False,
                       'interpolation': 'none',
                       'cnorm_type': 'Linear', # Colormap norm;  options are Log, Pow or Linear
                       'cpow_num': 0.6, # Used in the PowerNorm,
                       'div_midpoint': 0.0, # The cpow color norm normalizes data to [0,1] using np.sign(x-midpoint)*np.abs(x-midpoint)**(-cpow_num) -> [0,midpoint,1] if it is a divering cmap or [0,1] if it is not a divering cmap
                       'UseDivCmap': True,
                       'stretch_colors': False,
                       'cmap': 'None', # If cmap is none, the plot will inherit the parent's cmap,
                       'show_cpu_domains': False # plots lines showing how the CPUs are divvying up the computational region
                       }
    # We need the types of all the parameters for the config file
    BoolList = ['twoD', 'show_cbar', 'set_v_min', 'set_v_max',
             'show_shock', 'OutlineText', 'spatial_x', 'spatial_y', 'Show_FFT_region',
             'UseDivCmap', 'stretch_colors']
    IntList = ['mag_plot_type']
    FloatList = ['v_min', 'v_max', 'cpow_num', 'div_midpoint']
    #StrList = ['interpolation', 'cnorm_type', 'cmap'] # Not loading interpolation anymore
    StrList = ['interpolation','cnorm_type', 'cmap']

    gradient =  np.linspace(0, 1, 256)# A way to make the colorbar display better
    gradient = np.vstack((gradient, gradient))

    def __init__(self, parent, figwrapper):
        self.settings_window = None
        self.FigWrap = figwrapper
        self.parent = parent
        self.ChartTypes = self.FigWrap.PlotTypeDict.keys()
        self.chartType = self.FigWrap.chartType
        self.figure = self.FigWrap.figure
        self.SetPlotParam('spatial_y', self.GetPlotParam('twoD'), update_plot = False)
        self.InterpolationMethods = ['none','nearest', 'bilinear', 'bicubic', 'spline16',
            'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric',
            'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']

    def norm(self, vmin=None, vmax=None):
        if self.GetPlotParam('cnorm_type') =="Linear":
            if self.GetPlotParam('UseDivCmap'):
                return PowerNormWithNeg(1.0, vmin, vmax, midpoint = self.GetPlotParam('div_midpoint'), stretch_colors = self.GetPlotParam('stretch_colors'))
            else:
                return mcolors.Normalize(vmin, vmax)
        else:
            return PowerNormWithNeg(self.GetPlotParam('cpow_num'), vmin, vmax, div_cmap = self.GetPlotParam('UseDivCmap'),midpoint = self.GetPlotParam('div_midpoint'), stretch_colors = self.GetPlotParam('stretch_colors'))

    def ChangePlotType(self, str_arg):
        self.FigWrap.ChangeGraph(str_arg)

    def set_plot_keys(self):
        '''A helper function that will insure that each hdf5 file will only be
        opened once per time step'''
        # First make sure that omega_plasma & xi is loaded so we can fix the
        # x & y distances.
        self.arrs_needed = ['c_omp', 'istep', 'bx', 'by', 'bz']
        return self.arrs_needed

    def LoadData(self):
        ''' A Helper function that loads the data for the plot'''
        # First see of the x_axis and y_axis values have already been calculated
        # and stored in the DataDict for this time step
        self.c_omp = self.FigWrap.LoadKey('c_omp')[0]
        self.istep = self.FigWrap.LoadKey('istep')[0]
        self.mag_color = new_cmaps.cmaps[self.parent.MainParamDict['ColorMap']](0.5)

        if self.GetPlotParam('cmap') == 'None':
            if self.GetPlotParam('mag_plot_type')==1:
                self.cmap = 'phase'
            elif self.GetPlotParam('UseDivCmap'):
                self.cmap = self.parent.MainParamDict['DivColorMap']
            else:
                self.cmap = self.parent.MainParamDict['ColorMap']
        else:
            self.cmap = self.GetPlotParam('cmap')
        # see if the axis values are saved in the data dict
        if 'xaxis_values' in self.parent.DataDict.keys():
            self.xaxis_values = self.parent.DataDict['xaxis_values']
        else:
            # x-values haven't been calculated yet, generate them then save them to the dictionary for later.
            self.xaxis_values = np.arange(self.FigWrap.LoadKey('bx')[0,:,:].shape[1])/self.c_omp*self.istep
            self.parent.DataDict['xaxis_values'] = np.copy(self.xaxis_values)

        if self.GetPlotParam('mag_plot_type') == 0: # set f to thetaB
            self.ylabel = r'$\theta_B$'
            self.ann_label = r'$\theta_B$'
            if 'thetaB' in self.parent.DataDict.keys():
                self.f = self.parent.DataDict['thetaB']
            else:
                bx = self.FigWrap.LoadKey('bx')
                by = self.FigWrap.LoadKey('by')
                #bz = self.FigWrap.LoadKey('bz')
                #self.f = np.rad2deg(np.arctan2(np.sqrt(by**2+bz**2),bx))
                self.f = np.rad2deg(np.arctan2(by,bx))
                self.parent.DataDict['thetaB'] = self.f
        #removing this field
        #if self.GetPlotParam('mag_plot_type') == 1: # set f to xi_perp the phase of the perpendicular magnetic field
        #    self.ylabel = r'$\xi_\perp$'
        #    self.ann_label = r'$\xi_\perp$'
        #
        #    if 'xi_perp'+str(self.ind) in self.parent.DataDict.keys():
        #        self.f = self.parent.DataDict['xi_perp'+str(self.ind)]
        #    else:
        #        if np.abs(self.parent.bz0) <= 1E-6:

        #            sign = np.sign(self.FigWrap.LoadKey('bx')[self.ind,:,:]*np.tan(self.parent.btheta)-self.FigWrap.LoadKey('by')[self.ind,:,:])
        #            delta_bx_perp = (self.FigWrap.LoadKey('bx')[self.ind,:,:]-self.parent.bx0)*np.sin(self.parent.btheta)
        #            delta_by_perp = (self.FigWrap.LoadKey('by')[self.ind,:,:]-self.parent.by0)*np.cos(self.parent.btheta)
        #            delta_b_perp = sign*np.sqrt(delta_bx_perp**2+delta_by_perp**2)
        #            self.f = np.rad2deg(np.arctan2(self.FigWrap.LoadKey('bz')[0,:,:],delta_b_perp))
        #            self.parent.DataDict['xi_perp'] = self.f
        #        else:
        #            self.f = 0
        #    self.oneDslice = self.f.shape[0]/2
        #    # Have we already calculated min/max?
        #    if 'xi_perpmin_max'+str(self.ind) in self.parent.DataDict.keys():
        #    self.min_max = self.parent.DataDict['xi_perpmin_max'+str(self.ind)]

        #    else:
        #        self.min_max =  self.min_max_finder(self.f)
        #        self.parent.DataDict['xi_perpmin_max'+str(self.ind)] = list(self.min_max)

        if self.GetPlotParam('mag_plot_type') == 1: # Set f to deltaB/B0
            if ~np.isnan(self.parent.btheta):
                self.ylabel = r'$|\delta B|/B_0$'
                self.ann_label = r'$|\delta B|/B_0$'
            else:
                self.ylabel = r'$|\delta B|$'
                self.ann_label = r'$|\delta B|$'

            bx = self.FigWrap.LoadKey('bx')
            by = self.FigWrap.LoadKey('by')
            bz = self.FigWrap.LoadKey('bz')

            deltaB = (bx-self.parent.bx0)**2
            deltaB += (by-self.parent.by0)**2
            deltaB += (bz-self.parent.bz0)**2
            self.f = np.sqrt(deltaB)/self.parent.b0


        if self.GetPlotParam('mag_plot_type') == 2: # Set f to deltaB_perp/B0
            if ~np.isnan(self.parent.btheta):
                self.ylabel = r'$|\delta B_\perp|/B_0$'
                self.ann_label = r'$|\delta B_\perp|/B_0$'
            else:
                self.ylabel = r'$|\delta B_\perp|$'
                self.ann_label = r'$|\delta B_\perp|$'
            if 'deltaB_perp' in self.parent.DataDict.keys():
                self.f = self.parent.DataDict['deltaB_perp']
            elif np.isnan(self.parent.btheta):
                self.f = 1.0
            else:
                if np.abs(self.parent.bz0) <= 1E-6:
                    # Decompose the perpendicular components of the magnetic fields into two parts
                    # one is bz, the other is the in plane components
                    delta_by_perp = (self.FigWrap.LoadKey('by')-self.parent.by0)#*np.cos(self.parent.btheta)
                    self.f = np.sqrt(delta_by_perp**2+(self.FigWrap.LoadKey('bz')-self.parent.bz0)**2)/self.parent.b0
                    self.parent.DataDict['delta_b_perp'] = self.f

                self.parent.DataDict['deltaB_perp'] = self.f

        if self.GetPlotParam('mag_plot_type') == 3: # Set f to deltaB_para/B0
            if ~np.isnan(self.parent.btheta):
                self.ylabel = r'$\delta B_\parallel/B_0$'
                self.ann_label = r'$\delta B_\parallel/B_0$'
            else:
                self.ylabel = r'$\delta B_\parallel$'
                self.ann_label = r'$\delta B_\parallel$'

            if 'deltaB_para' in self.parent.DataDict.keys():
                self.f = self.parent.DataDict['deltaB_para']
            elif np.isnan(self.parent.btheta):
                self.f = 1.0
            else:
                if np.abs(self.parent.bz0) <= 1E-6:
                # Decompose the perpendicular components of the magnetic fields into two parts
                # one is bz, the other is the in plane components

                    # First take the dot product:
                    b_para = self.FigWrap.LoadKey('bx')#*self.parent.bx0+self.FigWrap.LoadKey('by')[0,:,:]*self.parent.by0
                    #b_para *= self.parent.b0**(-1)
                    self.f = (b_para-self.parent.b0)/self.parent.b0
                    self.parent.DataDict['deltaB_para'] = self.f

    def draw(self):

        ''' A function that draws the data. In the interest in speeding up the
        code, draw should only be called when you want to recreate the whole
        figure, i.e. it  will be slow. Most times you will only want to update
        what has changed in the figure. This will be done in a function called
        refresh, that should be much much faster.'''

        if self.GetPlotParam('OutlineText'):
            self.annotate_kwargs = {'horizontalalignment': 'right',
            'verticalalignment': 'top',
            'size' : self.parent.MainParamDict['annotateTextSize'],
            'path_effects' : [PathEffects.withStroke(linewidth=1.5,foreground="k")]
            }
        else:
            self.annotate_kwargs = {'horizontalalignment' : 'right',
            'verticalalignment' : 'top',
            'size' : self.parent.MainParamDict['annotateTextSize']}


        # Set the tick color
        tick_color = 'black'

        # Create a gridspec to handle spacing better
        self.gs = gridspec.GridSpecFromSubplotSpec(100,100, subplot_spec = self.parent.gs0[self.FigWrap.pos])

        # Now that the data is loaded, start making the plots
        if self.GetPlotParam('twoD'):
            if self.parent.MainParamDict['LinkSpatial'] != 0:
                if self.FigWrap.pos == self.parent.first_x and self.FigWrap.pos == self.parent.first_y:
                    self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]])
                elif self.FigWrap.pos == self.parent.first_x:
                    self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]],
                    sharey = self.parent.SubPlotList[self.parent.first_y[0]][self.parent.first_y[1]].graph.axes)
                elif self.FigWrap.pos == self.parent.first_y:
                    self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]],
                    sharex = self.parent.SubPlotList[self.parent.first_x[0]][self.parent.first_x[1]].graph.axes)
                else:
                    self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]],
                    sharex = self.parent.SubPlotList[self.parent.first_x[0]][self.parent.first_x[1]].graph.axes,
                    sharey = self.parent.SubPlotList[self.parent.first_y[0]][self.parent.first_y[1]].graph.axes)

            else:
                self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]])

            if self.parent.MainParamDict['2DSlicePlane'] == 0: # x-y plane
                if self.parent.MainParamDict['ImageAspect']:
                    self.cax = imshow(self.axes, self.f[self.parent.zSlice,:,:], origin = 'lower', norm = self.norm())
                else:
                    self.cax = imshow(self.axes, self.f[self.parent.zSlice,:,:], origin = 'lower', norm = self.norm(),
                                            aspect= 'auto')
            elif self.parent.MainParamDict['2DSlicePlane'] == 1: # x-y plane
                if self.parent.MainParamDict['ImageAspect']:
                    self.cax = imshow(self.axes, self.f[:, self.parent.ySlice,:], origin = 'lower', norm = self.norm())
                else:
                    self.cax = imshow(self.axes,self.f[:,self.parent.ySlice,:], origin = 'lower', norm = self.norm(),
                                            aspect= 'auto')

            self.ymin = 0
            self.ymax =  self.cax.get_array().shape[0]/self.c_omp*self.istep
            self.xmin = 0
            self.xmax =  self.cax.get_array().shape[1]/self.c_omp*self.istep

            self.vmin = self.cax.get_array().min()
            if self.GetPlotParam('set_v_min'):
                self.vmin = self.GetPlotParam('v_min')
            self.vmax = self.cax.get_array().max()
            if self.GetPlotParam('set_v_max'):
                self.vmax = self.GetPlotParam('v_max')
            if self.GetPlotParam('UseDivCmap') and not self.GetPlotParam('stretch_colors'):
                self.vmax = max(np.abs(self.vmin), self.vmax)
                self.vmin = -self.vmax
            self.cax.norm.vmin = self.vmin
            self.cax.norm.vmax = self.vmax
            self.cax.set_interpolation(self.GetPlotParam('interpolation'))
            self.cax.set_cmap(new_cmaps.cmaps[self.cmap])
            self.cax.set_extent([self.xmin,self.xmax, self.ymin, self.ymax])
            self.axes.add_artist(self.cax)


            self.TwoDan = self.axes.annotate(self.ann_label,
                            xy = (0.9,.9),
                            xycoords= 'axes fraction',
                            color = 'white',
                            **self.annotate_kwargs)

            self.axC = self.figure.add_subplot(self.gs[self.parent.cbar_extent[0]:self.parent.cbar_extent[1], self.parent.cbar_extent[2]:self.parent.cbar_extent[3]])
            if self.parent.MainParamDict['HorizontalCbars']:
                self.cbar = self.axC.imshow(self.gradient, aspect='auto',
                                            cmap=new_cmaps.cmaps[self.cmap])

                # Make the colobar axis more like the real colorbar
                self.cbar.set_extent([0, 1.0, 0, 1.0])
                self.axC.tick_params(axis='x',
                                which = 'both', # bothe major and minor ticks
                                top = 'off', # turn off top ticks
                                labelsize=self.parent.MainParamDict['NumFontSize'])

                self.axC.tick_params(axis='y',          # changes apply to the y-axis
                                which='both',      # both major and minor ticks are affected
                                left='off',      # ticks along the bottom edge are off
                                right='off',         # ticks along the top edge are off
                                labelleft='off')
            else:
                self.cbar = self.axC.imshow(np.transpose(self.gradient)[::-1], aspect='auto',
                                            cmap=new_cmaps.cmaps[self.cmap])

                # Make the colobar axis more like the real colorbar
                self.cbar.set_extent([0, 1.0, 0, 1.0])
                self.axC.tick_params(axis='x',
                                which = 'both', # bothe major and minor ticks
                                top = 'off', # turn off top ticks
                                bottom = 'off',
                                labelbottom = 'off',
                                labelsize=self.parent.MainParamDict['NumFontSize'])

                self.axC.tick_params(axis='y',          # changes apply to the y-axis
                                which='both',      # both major and minor ticks are affected
                                left='off',      # ticks along the bottom edge are off
                                right='on',         # ticks along the top edge are off
                                labelleft = 'off',
                                labelright  = 'on',
                                labelsize=self.parent.MainParamDict['NumFontSize'])

            if self.GetPlotParam('show_cbar') == 0:
                self.axC.set_visible = False
            else:
                self.CbarTickFormatter()

            self.shockline_2d = self.axes.axvline(self.parent.shock_loc, linewidth = 1.5, linestyle = '--', color = self.parent.shock_color, path_effects=[PathEffects.Stroke(linewidth=2, foreground='k'),
                                    PathEffects.Normal()])
            self.shockline_2d.set_visible(self.GetPlotParam('show_shock'))

            if int(matplotlib.__version__[0]) < 2:
                self.axes.set_axis_bgcolor('lightgrey')
            else:
                self.axes.set_facecolor('lightgrey')

            self.axes.tick_params(labelsize = self.parent.MainParamDict['NumFontSize'], color=tick_color)

            if self.parent.MainParamDict['SetxLim']:
                if self.parent.MainParamDict['xLimsRelative']:
                    self.axes.set_xlim(self.parent.MainParamDict['xLeft'] + self.parent.shock_loc,
                                       self.parent.MainParamDict['xRight'] + self.parent.shock_loc)
                else:
                    self.axes.set_xlim(self.parent.MainParamDict['xLeft'], self.parent.MainParamDict['xRight'])
            else:
                self.axes.set_xlim(self.xmin, self.xmax)
            if self.parent.MainParamDict['SetyLim']:
                self.axes.set_ylim(self.parent.MainParamDict['yBottom'],self.parent.MainParamDict['yTop'])
            else:
                self.axes.set_ylim(self.ymin, self.ymax)
            self.axes.set_xlabel(r'$x\ [c/\omega_{\rm pe}]$', labelpad = self.parent.MainParamDict['xLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
            if self.parent.MainParamDict['2DSlicePlane'] == 0:
                self.axes.set_ylabel(r'$y\ [c/\omega_{\rm pe}]$', labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
            if self.parent.MainParamDict['2DSlicePlane'] == 1:
                self.axes.set_ylabel(r'$z\ [c/\omega_{\rm pe}]$', labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])

        else:
            if self.parent.MainParamDict['LinkSpatial'] != 0 and self.parent.MainParamDict['LinkSpatial'] != 3:
                if self.FigWrap.pos == self.parent.first_x:
                    self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]])
                else:
                    self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]],
                    sharex = self.parent.SubPlotList[self.parent.first_x[0]][self.parent.first_x[1]].graph.axes)
            else:
                self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]])

            self.annotate_pos = [0.8,0.9]
            # Make the 1-D plots
            if self.parent.MainParamDict['Average1D']:
                self.line = self.axes.plot(self.xaxis_values, np.average(self.f.reshape(-1,self.f.shape[-1]), axis = 0), color = self.mag_color)
            else:
                self.line = self.axes.plot(self.xaxis_values, self.f[self.parent.zSlice,self.parent.ySlice,:], color = self.mag_color)

            # Set the Ymin/Ymax. Unnecessary.
            min_max = [self.line[0].get_data()[1].min(), self.line[0].get_data()[1].max()]
            dist = min_max[1]-min_max[0]
            min_max[0] -= 0.04*dist
            min_max[1] += 0.04*dist
            self.axes.set_ylim(min_max)

            self.shock_line = self.axes.axvline(self.parent.shock_loc, linewidth = 1.5, linestyle = '--', color = self.parent.shock_color, path_effects=[PathEffects.Stroke(linewidth=2, foreground='k'),
                        PathEffects.Normal()])

            self.shock_line.set_visible(self.GetPlotParam('show_shock'))

            if int(matplotlib.__version__[0]) < 2:
                self.axes.set_axis_bgcolor('lightgrey')
            else:
                self.axes.set_facecolor('lightgrey')

            self.axes.tick_params(labelsize = self.parent.MainParamDict['NumFontSize'], color=tick_color)#, tick1On= False, tick2On= False)

            if self.parent.MainParamDict['SetxLim']:
                if self.parent.MainParamDict['xLimsRelative']:
                    self.axes.set_xlim(self.parent.MainParamDict['xLeft'] + self.parent.shock_loc,
                                       self.parent.MainParamDict['xRight'] + self.parent.shock_loc)
                else:
                    self.axes.set_xlim(self.parent.MainParamDict['xLeft'], self.parent.MainParamDict['xRight'])
            else:
                self.axes.set_xlim(self.xaxis_values[0],self.xaxis_values[-1])


            if self.GetPlotParam('set_v_min'):
                self.axes.set_ylim(ymin = self.GetPlotParam('v_min'))
            if self.GetPlotParam('set_v_max'):
                self.axes.set_ylim(ymax = self.GetPlotParam('v_max'))

            self.axes.set_xlabel(r'$x\ [c/\omega_{\rm pe}]$', labelpad = self.parent.MainParamDict['xLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
            self.axes.set_ylabel(self.ylabel, labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])

        ####
        # FFT REGION PLOTTING CODE
        ####


        self.lineleft = self.axes.axvline(0, linewidth = 1.5, linestyle = ':', color = self.parent.FFT_color)
        self.lineright = self.axes.axvline(0, linewidth = 1.5, linestyle = ':', color = self.parent.FFT_color)
        self.lineleft.set_visible(self.GetPlotParam('show_FFT_region'))
        self.lineright.set_visible(self.GetPlotParam('show_FFT_region'))

        if self.GetPlotParam('show_FFT_region'):
            self.left_loc = self.parent.MainParamDict['FFTLeft'] + self.parent.shock_loc*self.parent.MainParamDict['FFTRelative']
            self.left_loc = max(self.left_loc, self.xaxis_values[0])
            self.lineleft.set_xdata([self.left_loc,self.left_loc])

            self.right_loc = self.parent.MainParamDict['FFTRight'] + self.parent.shock_loc*self.parent.MainParamDict['FFTRelative']
            self.right_loc = min(self.right_loc, self.xaxis_values[-1])
            self.lineright.set_xdata([self.right_loc,self.right_loc])
        if self.GetPlotParam('show_cpu_domains'):
            self.FigWrap.SetCpuDomainLines()
    def refresh(self):

        '''This is a function that will be called only if self.axes already
        holds a fields type plot. We only update things that have changed & are
        shown.  If hasn't changed or isn't shown, don't touch it. The difference
        between this and last time, is that we won't actually do any drawing in
        the plot. The plot will be redrawn after all subplots are refreshed. '''


        # Main goal, only change what is showing..

        self.lineleft.set_visible(self.GetPlotParam('show_FFT_region'))
        self.lineright.set_visible(self.GetPlotParam('show_FFT_region'))

        if self.GetPlotParam('show_FFT_region'):
            self.left_loc = self.parent.MainParamDict['FFTLeft'] + self.parent.shock_loc*self.parent.MainParamDict['FFTRelative']
            self.left_loc = max(self.left_loc, self.xaxis_values[0])
            self.lineleft.set_xdata([self.left_loc,self.left_loc])

            self.right_loc = self.parent.MainParamDict['FFTRight'] + self.parent.shock_loc*self.parent.MainParamDict['FFTRelative']
            self.right_loc = min(self.right_loc, self.xaxis_values[-1])
            self.lineright.set_xdata([self.right_loc,self.right_loc])


        # First do the 1D plots, because it is simpler
        if self.GetPlotParam('twoD') == 0:
            if self.parent.MainParamDict['Average1D']:
                self.line[0].set_data(self.xaxis_values, np.average(self.f.reshape(-1,self.f.shape[-1]), axis = 0))
            else:
                self.line[0].set_data(self.xaxis_values, self.f[self.parent.zSlice,self.parent.ySlice,:])

            min_max = [self.line[0].get_data()[1].min(), self.line[0].get_data()[1].max()]
            dist = min_max[1]-min_max[0]
            min_max[0] -= 0.04*dist
            min_max[1] += 0.04*dist
            self.axes.set_ylim(min_max)
            if self.GetPlotParam('show_shock'):
                self.shock_line.set_xdata([self.parent.shock_loc,self.parent.shock_loc])
            # xlims
            if self.parent.MainParamDict['SetxLim']:
                if self.parent.MainParamDict['xLimsRelative']:
                    self.axes.set_xlim(self.parent.MainParamDict['xLeft'] + self.parent.shock_loc,
                                       self.parent.MainParamDict['xRight'] + self.parent.shock_loc)
                else:
                    self.axes.set_xlim(self.parent.MainParamDict['xLeft'], self.parent.MainParamDict['xRight'])
            else:
                self.axes.set_xlim(self.xaxis_values[0], self.xaxis_values[-1])

            if self.GetPlotParam('set_v_min'):
                self.axes.set_ylim(ymin = self.GetPlotParam('v_min'))
            if self.GetPlotParam('set_v_max'):
                self.axes.set_ylim(ymax = self.GetPlotParam('v_max'))
            self.axes.set_ylabel(self.ylabel, size = self.parent.MainParamDict['AxLabelSize'])

        else: # Now refresh the plot if it is 2D
            if self.parent.MainParamDict['2DSlicePlane'] == 0: # x-y plane
                self.cax.set_data(self.f[self.parent.zSlice,:,:])
            elif self.parent.MainParamDict['2DSlicePlane'] == 1: # x-y plane
                self.cax.set_data(self.f[:,self.parent.ySlice,:])

            self.ymin = 0
            self.ymax =  self.cax.get_array().shape[0]/self.c_omp*self.istep
            self.xmin = 0
            self.xmax = self.xaxis_values[-1]
            self.TwoDan.set_text(self.ann_label)
            self.clims = np.copy([self.cax.get_array().min(), self.cax.get_array().max()])

            if self.parent.MainParamDict['SetxLim']:
                if self.parent.MainParamDict['xLimsRelative']:
                    self.axes.set_xlim(self.parent.MainParamDict['xLeft'] + self.parent.shock_loc,
                                       self.parent.MainParamDict['xRight'] + self.parent.shock_loc)
                else:
                    self.axes.set_xlim(self.parent.MainParamDict['xLeft'], self.parent.MainParamDict['xRight'])
            else:
                self.axes.set_xlim(self.xmin,self.xmax)
            if self.parent.MainParamDict['SetyLim']:
                self.axes.set_ylim(self.parent.MainParamDict['yBottom'],self.parent.MainParamDict['yTop'])
            else:
                self.axes.set_ylim(self.ymin,self.ymax)

            self.cax.set_extent([self.xmin, self.xmax, self.ymin, self.ymax])

            self.vmin = self.cax.get_array().min()
            if self.GetPlotParam('set_v_min'):
                self.vmin = self.GetPlotParam('v_min')
            self.vmax = self.cax.get_array().max()
            if self.GetPlotParam('set_v_max'):
                self.vmax = self.GetPlotParam('v_max')
            if self.GetPlotParam('UseDivCmap') and not self.GetPlotParam('stretch_colors'):
                self.vmax = max(np.abs(self.vmin), self.vmax)
                self.vmin = -self.vmax
            self.cax.norm.vmin = self.vmin
            self.cax.norm.vmax = self.vmax

            self.CbarTickFormatter()
            if self.parent.MainParamDict['2DSlicePlane'] == 0:
                self.axes.set_ylabel(r'$y\ [c/\omega_{\rm pe}]$', labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
            if self.parent.MainParamDict['2DSlicePlane'] == 1:
                self.axes.set_ylabel(r'$z\ [c/\omega_{\rm pe}]$', labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])

            if self.GetPlotParam('show_shock'):
                self.shockline_2d.set_xdata([self.parent.shock_loc,self.parent.shock_loc])
            #self.axes.draw_artist(self.axes.patch)
            #self.axes.draw_artist(self.cax)
            #self.axes.draw_artist(self.axes.xaxis)
        if self.GetPlotParam('show_cpu_domains'):
            self.FigWrap.UpdateCpuDomainLines()
    def CbarTickFormatter(self):
        ''' A helper function that sets the cbar ticks & labels. This used to be
        easier, but because I am no longer using the colorbar class i have to do
        stuff manually.'''
        clim = np.copy(self.cax.get_clim())
        if self.GetPlotParam('show_cbar'):
            if self.GetPlotParam('cnorm_type') == "Log":
                self.cbar.set_extent([np.log10(clim[0]),np.log10(clim[1]),0,1])
                self.axC.set_xlim(np.log10(clim[0]),np.log10(clim[1]))

            elif self.GetPlotParam('cnorm_type') == "Pow":
                # re-create the gradient with the data values
                # First make a colorbar in the negative region that is linear in the pow_space
                data_range = np.linspace(clim[0],clim[1],512)

                cbardata = PowerNormFunc(data_range, vmin = data_range[0], vmax = data_range[-1], gamma = self.GetPlotParam('cpow_num'), midpoint = self.GetPlotParam('div_midpoint'), div_cmap = self.GetPlotParam('UseDivCmap'), stretch_colors = self.GetPlotParam('stretch_colors'))
                cbardata = np.vstack((cbardata,cbardata))
                if self.parent.MainParamDict['HorizontalCbars']:
                    self.cbar.set_data(cbardata)
                    self.cbar.set_extent([clim[0],clim[1],0,1])
                    self.axC.set_xlim(clim[0],clim[1])
                else:
                    self.cbar.set_data(np.transpose(cbardata)[::-1])
                    self.cbar.set_extent([0,1,clim[0],clim[1]])
                    self.axC.set_ylim(clim[0],clim[1])
                    self.axC.locator_params(axis='y',nbins=6)

            elif self.GetPlotParam('cnorm_type') == "Linear" and self.GetPlotParam('UseDivCmap'):
                # re-create the gradient with the data values
                # First make a colorbar in the negative region that is linear in the pow_space
                data_range = np.linspace(clim[0],clim[1],512)

                cbardata = PowerNormFunc(data_range, vmin = data_range[0], vmax = data_range[-1], gamma = 1.0, div_cmap = self.GetPlotParam('UseDivCmap'), midpoint = self.GetPlotParam('div_midpoint'), stretch_colors = self.GetPlotParam('stretch_colors'))
                cbardata = np.vstack((cbardata,cbardata))
                if self.parent.MainParamDict['HorizontalCbars']:
                    self.cbar.set_data(cbardata)
                    self.cbar.set_extent([clim[0],clim[1],0,1])
                    self.axC.set_xlim(clim[0],clim[1])
                else:
                    self.cbar.set_data(np.transpose(cbardata)[::-1])
                    self.cbar.set_extent([0,1,clim[0],clim[1]])
                    self.axC.set_ylim(clim[0],clim[1])
                    self.axC.locator_params(axis='y',nbins=6)

            else:# self.GetPlotParam('cnorm_type') == "Linear":
                if self.parent.MainParamDict['HorizontalCbars']:

                    self.cbar.set_extent([clim[0],clim[1],0,1])
                    self.axC.set_xlim(clim[0],clim[1])
                else:
#                    self.cbar.set_data(np.transpose(self.gradient)[::-1])
                    self.cbar.set_extent([0,1,clim[0],clim[1]])
                    self.axC.set_ylim(clim[0],clim[1])
                    self.axC.locator_params(axis='y',nbins=6)


    def GetPlotParam(self, keyname):
        return self.FigWrap.GetPlotParam(keyname)

    def SetPlotParam(self, keyname, value, update_plot = True, NeedsRedraw = False):
        self.FigWrap.SetPlotParam(keyname, value, update_plot = update_plot, NeedsRedraw = NeedsRedraw)

    def OpenSettings(self):
        if self.settings_window is None:
            self.settings_window = BSettings(self)
        else:
            self.settings_window.destroy()
            self.settings_window = BSettings(self)


class BSettings(Tk.Toplevel):
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
        ctypeChooser = apply(ttk.OptionMenu, (frm, self.ctypevar, self.parent.chartType) + tuple(self.parent.ChartTypes))
        ctypeChooser.grid(row =0, column = 1, sticky = Tk.W + Tk.E)


        self.TwoDVar = Tk.IntVar(self) # Create a var to track whether or not to plot in 2-D
        self.TwoDVar.set(self.parent.GetPlotParam('twoD'))
        cb = ttk.Checkbutton(frm, text = "Show in 2-D",
                variable = self.TwoDVar,
                command = self.Change2d)
        cb.grid(row = 1, sticky = Tk.W)

        # the Radiobox Control to choose the Field Type
        self.MagList = ['Theta B','Delta B/B_0','Delta B_perp', 'Delta B_para']
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
        cb.grid(row = 6, column =2, sticky = Tk.W)

        # Now the field lim
        self.setZminVar = Tk.IntVar()
        self.setZminVar.set(self.parent.GetPlotParam('set_v_min'))
        self.setZminVar.trace('w', self.setZminChanged)

        self.setZmaxVar = Tk.IntVar()
        self.setZmaxVar.set(self.parent.GetPlotParam('set_v_max'))
        self.setZmaxVar.trace('w', self.setZmaxChanged)



        self.Zmin = Tk.StringVar()
        self.Zmin.set(str(self.parent.GetPlotParam('v_min')))

        self.Zmax = Tk.StringVar()
        self.Zmax.set(str(self.parent.GetPlotParam('v_max')))


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
        cb.grid(row = 6, column = 3, sticky = Tk.W)

        self.CPUVar = Tk.IntVar()
        self.CPUVar.set(self.parent.GetPlotParam('show_cpu_domains'))
        cb = ttk.Checkbutton(frm, text = "Show CPU domains",
                        variable = self.CPUVar,
                        command = self.CPUVarHandler)
        cb.grid(row = 8, column = 0, sticky = Tk.W)

        self.FFTVar = Tk.IntVar()
        self.FFTVar.set(self.parent.GetPlotParam('show_FFT_region'))
        cb = ttk.Checkbutton(frm, text = "Show FFT Region",
                        variable = self.FFTVar,
                        command = self.FFTVarHandler)
        cb.grid(row = 5, column =2, sticky = Tk.W)

        # Control whether or not diverging cmap is used
        self.DivVar = Tk.IntVar()
        self.DivVar.set(self.parent.GetPlotParam('UseDivCmap'))
        cb = ttk.Checkbutton(frm, text = "Use Diverging Cmap",
                        variable = self.DivVar,
                        command = self.DivHandler)
        cb.grid(row = 7, column = 2, sticky = Tk.W)

        # Use full div cmap
        self.StretchVar = Tk.IntVar()
        self.StretchVar.set(self.parent.GetPlotParam('stretch_colors'))
        cb = ttk.Checkbutton(frm, text = "Asymmetric Color Space",
                        variable = self.StretchVar,
                        command = self.StretchHandler)
        cb.grid(row = 7, column = 3, sticky = Tk.W)

        # Create the OptionMenu to chooses the cnorm_type:
        self.cnormvar = Tk.StringVar(self)
        self.cnormvar.set(self.parent.chartType) # default value
        self.cnormvar.trace('w', self.cnormChanged)

        ttk.Label(frm, text="Choose Color Norm:").grid(row=8, column = 2)
        cnormChooser = apply(ttk.OptionMenu, (frm, self.cnormvar, self.parent.GetPlotParam('cnorm_type')) + tuple(['Pow', 'Linear']))
        cnormChooser.grid(row =8, column = 3, sticky = Tk.W + Tk.E)



        # Now the gamma of the pow norm
        self.powGamma = Tk.StringVar()
        self.powGamma.set(str(self.parent.GetPlotParam('cpow_num')))
        ttk.Label(frm, text ='gamma =').grid(row = 9, column = 2, sticky =Tk.E)
        self.GammaEnter = ttk.Entry(frm, textvariable=self.powGamma, width=7)
        self.GammaEnter.grid(row = 9, column = 3)

        ttk.Label(frm, text ='If cnorm is Pow => sign(data)*|data|**gamma').grid(row = 10, column = 1,columnspan = 3, sticky =Tk.E)




    def CbarHandler(self, *args):
        if self.parent.GetPlotParam('show_cbar')== self.CbarVar.get():
            pass
        else:
            if self.parent.GetPlotParam('twoD'):
                self.parent.axC.set_visible(self.CbarVar.get())

            self.parent.SetPlotParam('show_cbar', self.CbarVar.get(), update_plot =self.parent.GetPlotParam('twoD'))

    def DivHandler(self, *args):
        if self.parent.GetPlotParam('UseDivCmap')== self.DivVar.get():
            pass
        elif self.parent.GetPlotParam('twoD'):
            self.parent.SetPlotParam('UseDivCmap', self.DivVar.get(),  NeedsRedraw = True)
        else:
            self.parent.SetPlotParam('UseDivCmap', self.DivVar.get(),  update_plot = False)


    def StretchHandler(self, *args):
        if self.parent.GetPlotParam('stretch_colors')== self.StretchVar.get():
            pass
        elif self.parent.GetPlotParam('twoD'):
            self.parent.SetPlotParam('stretch_colors', self.StretchVar.get(), NeedsRedraw = True)
        else:
            self.parent.SetPlotParam('stretch_colors', self.StretchVar.get(), update_plot = False)

    def cnormChanged(self, *args):
        if self.parent.GetPlotParam('cnorm_type')== self.cnormvar.get():
            pass
        elif self.parent.GetPlotParam('twoD'):
            self.parent.SetPlotParam('cnorm_type', self.cnormvar.get(), NeedsRedraw = True)
        else:
            self.parent.SetPlotParam('cnorm_type', self.cnormvar.get(), update_plot = False)


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

    def FFTVarHandler(self, *args):
        if self.parent.GetPlotParam('show_FFT_region')== self.FFTVar.get():
            pass
        else:
            self.parent.SetPlotParam('show_FFT_region', self.FFTVar.get())

    def InterpolChanged(self, *args):
        if self.InterpolVar.get() == self.parent.GetPlotParam('interpolation'):
            pass
        else:
            if self.parent.GetPlotParam('twoD'):
                self.parent.cax.set_interpolation(self.InterpolVar.get())
            self.parent.SetPlotParam('interpolation', self.InterpolVar.get())

    def setZminChanged(self, *args):
        if self.setZminVar.get() == self.parent.GetPlotParam('set_v_min'):
            pass
        else:
            self.parent.SetPlotParam('set_v_min', self.setZminVar.get())

    def setZmaxChanged(self, *args):
        if self.setZmaxVar.get() == self.parent.GetPlotParam('set_v_max'):
            pass
        else:
            self.parent.SetPlotParam('set_v_max', self.setZmaxVar.get())

    def RadioMag(self):
        if self.MagTypeVar.get() == self.parent.GetPlotParam('mag_plot_type'):
            pass
        else:
            self.parent.SetPlotParam('mag_plot_type', self.MagTypeVar.get())

    def TxtEnter(self, e):
        self.FieldsCallback()
        self.GammaCallback()

    def GammaCallback(self):
        try:
        #make sure the user types in a float
            if np.abs(float(self.powGamma.get()) - self.parent.GetPlotParam('cpow_num')) > 1E-4:
                if self.parent.GetPlotParam('twoD') and self.parent.GetPlotParam('cnorm_type')=='Pow':
                    self.parent.SetPlotParam('cpow_num', float(self.powGamma.get()), NeedsRedraw = True)

                else:
                    self.parent.SetPlotParam('cpow_num', float(self.powGamma.get()), update_plot = False)

        except ValueError:
            #if they type in random stuff, just set it ot the param value
            self.powGamma.set(str(self.parent.GetPlotParam('cpow_num')))

    def FieldsCallback(self):
        tkvarLimList = [self.Zmin, self.Zmax]
        plot_param_List = ['v_min', 'v_max']
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
            self.parent.SetPlotParam('v_min', self.parent.GetPlotParam('v_min'))

    def CPUVarHandler(self, *args):
        if self.parent.GetPlotParam('show_cpu_domains')== self.CPUVar.get():
            pass
        else:
            self.parent.SetPlotParam('show_cpu_domains', self.CPUVar.get(), update_plot = False)
            if self.parent.GetPlotParam('show_cpu_domains'):
                self.parent.FigWrap.SetCpuDomainLines()

            else: # We need to get remove of the cpu lines. Pop them out of the array and remove them from the list.
                self.parent.FigWrap.RemoveCpuDomainLines()
            self.parent.parent.canvas.draw()
            self.parent.parent.canvas.get_tk_widget().update_idletasks()

    def OnClosing(self):
        self.parent.settings_window = None
        self.destroy()
