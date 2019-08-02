#!/usr/bin/env python
#import tkinter as Tk
#from tkinter import ttk
import matplotlib
import numpy as np
import numpy.ma as ma
import sys
sys.path.append('../')

import new_cmaps
from new_cnorms import PowerNormWithNeg, PowerNormFunc
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
                       'show_cpu_domains': False, # plots lines showing how the CPUs are divvying up the computational region
                       'face_color': 'gainsboro'
                       }
    # We need the types of all the parameters for the config file

    gradient =  np.linspace(0, 1, 256)# A way to make the colorbar display better
    gradient = np.vstack((gradient, gradient))

    def __init__(self, parent, pos, param_dict):
        self.param_dict = {}
        for key, val in self.plot_param_dict.items():
            self.param_dict[key] = val
        for key, val in param_dict.items():
            self.param_dict[key] = val
        self.pos = pos
        self.parent = parent
        self.chartType = 'MagPlots'
        self.figure = self.parent.figure
        self.InterpolationMethods = ['none','nearest', 'bilinear', 'bicubic', 'spline16',
            'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric',
            'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']
        self.mag_color = new_cmaps.cmaps[self.parent.MainParamDict['ColorMap']](0.5)
    def norm(self, vmin=None, vmax=None):
        if self.GetPlotParam('cnorm_type') =="Linear":
            if self.GetPlotParam('UseDivCmap'):
                return PowerNormWithNeg(1.0, vmin, vmax, midpoint = self.GetPlotParam('div_midpoint'), stretch_colors = self.GetPlotParam('stretch_colors'))
            else:
                return mcolors.Normalize(vmin, vmax)
        else:
            return PowerNormWithNeg(self.GetPlotParam('cpow_num'), vmin, vmax, div_cmap = self.GetPlotParam('UseDivCmap'),midpoint = self.GetPlotParam('div_midpoint'), stretch_colors = self.GetPlotParam('stretch_colors'))
    def update_data(self, output):
        ''' A Helper function that loads the data for the plot'''
        # First see of the x_axis and y_axis values have already been calculated
        # and stored in the DataDict for this time step
        self.c_omp = getattr(output, 'c_omp')
        self.istep = getattr(output, 'istep')


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
        # x-values haven't been calculated yet, generate them then save them to the dictionary for later.
        self.xaxis_values = np.arange(getattr(output, 'bx')[0,:,:].shape[1])/self.c_omp*self.istep


        if self.GetPlotParam('mag_plot_type') == 0: # set f to thetaB
            self.ylabel = r'$\theta_B$'
            self.ann_label = r'$\theta_B$'
            bx = getattr(output, 'bx')
            by = getattr(output, 'by')
            bz = getattr(output, 'bz')
            self.f = np.rad2deg(np.arctan2(np.sqrt(by**2+bz**2),np.abs(bx)))


        if self.GetPlotParam('mag_plot_type') == 1: # Set f to deltaB/B0
            if ~np.isnan(self.parent.btheta):
                self.ylabel = r'$|\delta B|/B_0$'
                self.ann_label = r'$|\delta B|/B_0$'
            else:
                self.ylabel = r'$|\delta B|$'
                self.ann_label = r'$|\delta B|$'

            bx = getattr(output, 'bx')
            by = getattr(output, 'by')
            bz = getattr(output, 'bz')

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
            if np.isnan(self.parent.btheta):
                self.f = 1.0
            else:
                if np.abs(self.parent.b0) >= 1E-10:
                    # Decompose the perpendicular components of the magnetic fields into two parts
                    # one is bz, the other is the in plane components
                    delta_by_perp = (getattr(output, 'by')-self.parent.by0)#*np.cos(self.parent.btheta)
                    self.f = np.sqrt(delta_by_perp**2+(getattr(output, 'bz')-self.parent.bz0)**2)/self.parent.b0
                else:
                    # Decompose the perpendicular components of the magnetic fields into two parts
                    # one is bz, the other is the in plane components
                    delta_by_perp = (getattr(output, 'by')-self.parent.by0)#*np.cos(self.parent.btheta)
                    self.f = np.sqrt(delta_by_perp**2+(getattr(output, 'bz')-self.parent.bz0)**2)



        if self.GetPlotParam('mag_plot_type') == 3: # Set f to deltaB_para/B0
            if ~np.isnan(self.parent.btheta):
                self.ylabel = r'$\delta B_\parallel/B_0$'
                self.ann_label = r'$\delta B_\parallel/B_0$'
            else:
                self.ylabel = r'$\delta B_\parallel$'
                self.ann_label = r'$\delta B_\parallel$'

            if np.isnan(self.parent.btheta):
                self.f = 1.0
            else:
                # First take the dot product:
                b_para = getattr(output, 'bx')#*self.parent.bx0+getattr(output, 'by')[0,:,:]*self.parent.by0
                #b_para *= self.parent.b0**(-1)
                self.f = (b_para-self.parent.b0)/self.parent.b0

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
        self.gs = gridspec.GridSpecFromSubplotSpec(100,100, subplot_spec = self.parent.gs0[self.pos])

        # Now that the data is loaded, start making the plots
        if self.GetPlotParam('twoD'):
            self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]])

            if self.parent.MainParamDict['2DSlicePlane'] == 0: # x-y plane
                if self.parent.MainParamDict['ImageAspect']:
                    self.cax = self.axes.imshow( self.f[self.parent.zSlice,:,:], origin = 'lower', norm = self.norm())
                else:
                    self.cax = self.axes.imshow( self.f[self.parent.zSlice,:,:], origin = 'lower', norm = self.norm(),
                                            aspect= 'auto')
            elif self.parent.MainParamDict['2DSlicePlane'] == 1: # x-y plane
                if self.parent.MainParamDict['ImageAspect']:
                    self.cax = self.axes.imshow( self.f[:, self.parent.ySlice,:], origin = 'lower', norm = self.norm())
                else:
                    self.cax = self.axes.imshow(self.f[:,self.parent.ySlice,:], origin = 'lower', norm = self.norm(),
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
                                top = False, # turn off top ticks
                                labelsize=self.parent.MainParamDict['NumFontSize'])

                self.axC.tick_params(axis='y',          # changes apply to the y-axis
                                which='both',      # both major and minor ticks are affected
                                left=False,      # ticks along the bottom edge are off
                                right=False,         # ticks along the top edge are off
                                labelleft=False)
            else:
                self.cbar = self.axC.imshow(np.transpose(self.gradient)[::-1], aspect='auto', origin='upper',
                                            cmap=new_cmaps.cmaps[self.cmap])

                # Make the colobar axis more like the real colorbar
                self.cbar.set_extent([0, 1.0, 0, 1.0])
                self.axC.tick_params(axis='x',
                                which = 'both', # bothe major and minor ticks
                                top = False, # turn off top ticks
                                bottom = False,
                                labelbottom = False,
                                labelsize=self.parent.MainParamDict['NumFontSize'])

                self.axC.tick_params(axis='y',          # changes apply to the y-axis
                                which='both',      # both major and minor ticks are affected
                                left=False,      # ticks along the bottom edge are off
                                right= True,         # ticks along the top edge are off
                                labelleft = False,
                                labelright  = True,
                                labelsize=self.parent.MainParamDict['NumFontSize'])

            if self.GetPlotParam('show_cbar') == 0:
                self.axC.set_visible = False
            else:
                self.CbarTickFormatter()

            self.shockline_2d = self.axes.axvline(self.parent.shock_loc, linewidth = 1.5, linestyle = '--', color = self.parent.shock_color, path_effects=[PathEffects.Stroke(linewidth=2, foreground='k'),
                                    PathEffects.Normal()])
            self.shockline_2d.set_visible(self.GetPlotParam('show_shock'))

            if int(matplotlib.__version__[0]) < 2:
                self.axes.set_axis_bgcolor(self.GetPlotParam('face_color'))
            else:
                self.axes.set_facecolor(self.GetPlotParam('face_color'))

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
                self.axes.set_axis_bgcolor(self.GetPlotParam('face_color'))
            else:
                self.axes.set_facecolor(self.GetPlotParam('face_color'))

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
                self.axes.set_ylim(bottom = self.GetPlotParam('v_min'))
            if self.GetPlotParam('set_v_max'):
                self.axes.set_ylim(top = self.GetPlotParam('v_max'))

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
                self.axes.set_ylim(bottom = self.GetPlotParam('v_min'))
            if self.GetPlotParam('set_v_max'):
                self.axes.set_ylim(top = self.GetPlotParam('v_max'))
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
        return self.param_dict[keyname]
