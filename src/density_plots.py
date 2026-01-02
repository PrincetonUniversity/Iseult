#!/usr/bin/env pythonw
import tkinter as Tk
from tkinter import ttk
import matplotlib
import numpy as np
import numpy.ma as ma
import new_cmaps
from new_cnorms import PowerNormWithNeg, PowerNormFunc
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as PathEffects
from matplotlib.ticker import FuncFormatter
import streamlines

class DensPanel:
    # A dictionary of all of the parameters for this plot with the default parameters

    plot_param_dict = {'twoD': 0,
                       'dens_type': 0, #0 = n, 1 = n_i, 2 =n_e, 3=rho
                       'show_cbar': True,
                       'set_color_limits': False,
                       'v_min': 0,
                       'v_max' : 10,
                       'set_v_min': False,
                       'set_v_max': False,
                       'show_labels' : True,
                       'show_shock' : False,
                       'OutlineText': True,
                       'spatial_x': True,
                       'spatial_y': False,
                       'interpolation': 'none',
                       'normalize_density': True, # Normalize density to it's upstream values
                       'cnorm_type': 'Linear', # Colormap norm;  options are Pow or Linear
                       'cpow_num': 0.6, # Used in the PowerNorm
                       'UseDivCmap': False,
                       'div_midpoint': 0.0,
                       'stretch_colors': False,
                       'cmap': 'None', # If cmap is none, the plot will inherit the parent's cmap
                       'show_cpu_domains': False, # plots lines showing how the CPUs are divvying up the computational region
                       'face_color': 'gainsboro'}

    streamlines.add_streamline_params(plot_param_dict)

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
        elif self.GetPlotParam('cnorm_type') == "Log":
            return  mcolors.LogNorm(vmin, vmax)
        else:
            return PowerNormWithNeg(self.GetPlotParam('cpow_num'), vmin, vmax, div_cmap = self.GetPlotParam('UseDivCmap'),midpoint = self.GetPlotParam('div_midpoint'), stretch_colors = self.GetPlotParam('stretch_colors'))

    def ChangePlotType(self, str_arg):
        self.FigWrap.ChangeGraph(str_arg)

    def set_plot_keys(self):
        '''A helper function that will insure that each hdf5 file will only be
        opened once per time step'''
        # First make sure that omega_plasma & xi is loaded so we can fix the
        # x & y distances.

        self.arrs_needed = ['c_omp', 'istep', 'dens']

        # Load ppc if we are normalizing the density
        self.arrs_needed.append('ppc0')
        # To plot rho we need both dens and densi
        if self.GetPlotParam('dens_type') >0: # Load the ion density
            self.arrs_needed.append('densi')

        if self.GetPlotParam('show_streamlines') and self.GetPlotParam('twoD'):
            streamlines.add_streamline_plot_keys(self)

        return self.arrs_needed
    def LoadData(self):
        ''' A Helper function that loads the data for the plot'''
        if self.GetPlotParam('cmap') == 'None':
            if self.GetPlotParam('UseDivCmap'):
                self.cmap = self.parent.MainParamDict['DivColorMap']
            else:
                self.cmap = self.parent.MainParamDict['ColorMap']

        else:
            self.cmap = self.GetPlotParam('cmap')


        if self.GetPlotParam('normalize_density'):
            self.ppc0 = self.FigWrap.LoadKey('ppc0')
        elif not np.isnan(self.FigWrap.LoadKey('ppc0')):
            self.ppc0 = 1.0
        else:
            self.SetPlotParam('normalize_density', False, update_plot = False)

        self.dens_color = new_cmaps.cmaps[self.parent.MainParamDict['ColorMap']](0.5)
        # get c_omp and istep to convert cells to physical units
        self.c_omp = self.FigWrap.LoadKey('c_omp')
        self.istep = self.FigWrap.LoadKey('istep')
        self.dens = self.FigWrap.LoadKey('dens')[:,:,:]

        # Now calculate rho if needed.
        if self.GetPlotParam('dens_type') == 1:
            if 'densi' in self.parent.DataDict.keys():
                self.densi = self.parent.DataDict['densi']
            else:
                self.densi = self.FigWrap.LoadKey('densi')[:,:,:]
                self.parent.DataDict['densi'] = self.densi
        if self.GetPlotParam('dens_type') == 2:
            if 'dense' in self.parent.DataDict.keys():
                self.dense = self.parent.DataDict['dense']
            else:
                self.dense = self.FigWrap.LoadKey('dens')[:,:,:]-self.FigWrap.LoadKey('densi')[:,:,:]
                self.parent.DataDict['dense'] = self.dense

        if self.GetPlotParam('dens_type')== 3: # rho
            if 'rho' in self.parent.DataDict.keys():
                self.rho = self.parent.DataDict['rho']
            else:
                self.rho = 2*self.FigWrap.LoadKey('densi')[:,:,:] - self.dens
                self.parent.DataDict['rho'] = self.rho

        if 'xaxis_values' in self.parent.DataDict.keys():
            # Generate the x and y axes
            self.xaxis_values = self.parent.DataDict['xaxis_values']
        else:
            self.xaxis_values = np.arange(self.dens.shape[2])/self.c_omp*self.istep
            self.parent.DataDict['xaxis_values'] = self.xaxis_values
        # y values not needed so commenting out
        # self.y_values =  np.arange(self.zval.shape[0])/self.c_omp*self.istep

    def draw(self):
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
        self.gs = gridspec.GridSpecFromSubplotSpec(100,100, subplot_spec = self.parent.gs0[self.FigWrap.pos])#, bottom=0.2,left=0.1,right=0.95, top = 0.95)

        # Now that the data is loaded, start making the plots
        if self.GetPlotParam('twoD'):
            # Link up the spatial axes if desired
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

            # First choose the 'zval' to plot, we can only do one because it is 2-d.
            if self.FigWrap.GetPlotParam('dens_type') == 0:
                if self.FigWrap.GetPlotParam('normalize_density'):
                    self.zval = self.dens*self.ppc0**(-1.0)
                    self.two_d_label = r'$n/n_0$'
                else:
                    self.zval = self.dens
                    self.two_d_label = r'$n$'

            if self.FigWrap.GetPlotParam('dens_type') == 1:
                if self.FigWrap.GetPlotParam('normalize_density'):
                    self.zval = self.densi*self.ppc0**(-1.0)
                    self.two_d_label = r'$n_i/n_0$'
                else:
                    self.zval = self.densi
                    self.two_d_label = r'$n_i$'
            if self.FigWrap.GetPlotParam('dens_type') == 2:
                if self.FigWrap.GetPlotParam('normalize_density'):
                    self.zval = self.dense*self.ppc0**(-1.0)
                    self.two_d_label = r'$n_e/n_0$'
                else:
                    self.zval = self.dense
                    self.two_d_label = r'$n_e$'
            if self.FigWrap.GetPlotParam('dens_type') == 3:
                if self.FigWrap.GetPlotParam('normalize_density'):
                    self.zval = self.rho*self.ppc0**(-1.0)
                    self.two_d_label = r'$\rho/n_0$'
                else:
                    self.zval = self.rho
                    self.two_d_label = r'$\rho$'


            if self.parent.MainParamDict['2DSlicePlane'] ==0: # x-y plane
                if self.parent.MainParamDict['ImageAspect']:
                    self.cax = self.axes.imshow(self.zval[self.parent.zSlice,:,:], norm = self.norm(), origin = 'lower')
                else:
                    self.cax = self.axes.imshow(self.zval[self.parent.zSlice,:,:], norm = self.norm(), origin = 'lower',
                                                aspect = 'auto')
            elif self.parent.MainParamDict['2DSlicePlane'] ==1: # x-z plane
                if self.parent.MainParamDict['ImageAspect']:
                    self.cax = self.axes.imshow(self.zval[:,self.parent.ySlice,:], norm = self.norm(), origin = 'lower')
                else:
                    self.cax = self.axes.imshow(self.zval[:,self.parent.ySlice,:], norm = self.norm(), origin = 'lower',
                                                aspect = 'auto')
            elif self.parent.MainParamDict['2DSlicePlane'] ==2: # y-z plane
                if self.parent.MainParamDict['ImageAspect']:
                    self.cax = self.axes.imshow(self.zval[:,:,self.parent.xSlice], norm = self.norm(), origin = 'lower')
                else:
                    self.cax = self.axes.imshow(self.zval[:,:,self.parent.xSlice], norm = self.norm(), origin = 'lower',
                                                aspect = 'auto')


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
            if self.GetPlotParam('UseDivCmap') and self.GetPlotParam('stretch_colors'):
                self.vmax = max(np.abs(self.vmin), self.vmax)
                self.vmin = -self.vmax
            self.cax.norm.vmin = self.vmin
            self.cax.norm.vmax = self.vmax
            self.cax.set_interpolation(self.GetPlotParam('interpolation'))
            self.cax.set_cmap(new_cmaps.cmaps[self.cmap])
            self.cax.set_extent([self.xmin, self.xmax, self.ymin, self.ymax])

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


            self.axC = self.figure.add_subplot(self.gs[self.parent.cbar_extent[0]:self.parent.cbar_extent[1], self.parent.cbar_extent[2]:self.parent.cbar_extent[3]])
            self.parent.cbarList.append(self.axC)
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
            else: #Cbar is on the vertical
                self.cbar = self.axC.imshow(np.transpose(self.gradient)[::-1], aspect='auto',
                                            cmap=new_cmaps.cmaps[self.cmap])
                # Make the colobar axis more like the real colorbar
                self.cbar.set_extent([0, 1.0, 0, 1.0])
                self.axC.tick_params(axis='x',
                                which = 'both', # bothe major and minor ticks
                                top = False, # turn off top ticks
                                bottom = False, # turn off top ticks
                                labelbottom = False, # turn off top ticks
                                labelsize=self.parent.MainParamDict['NumFontSize'])

                self.axC.tick_params(axis='y',          # changes apply to the y-axis
                                which='both',      # both major and minor ticks are affected
                                left=False,      # ticks along the bottom edge are off
                                right=True,         # ticks along the top edge are off
                                labelleft=False,
                                labelright = True,
                                labelsize = self.parent.MainParamDict['NumFontSize'])

            if not self.GetPlotParam('show_cbar'):
                self.axC.set_visible(False)
            else:
                self.CbarTickFormatter()

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
                self.axes.set_xlim(self.xmin,self.xmax)

            if self.parent.MainParamDict['SetyLim']:
                self.axes.set_ylim(self.parent.MainParamDict['yBottom']*self.c_omp/self.istep, self.parent.MainParamDict['yTop']*self.c_omp/self.istep)
            else:
                self.axes.set_ylim(self.ymin, self.ymax)
            self.axes.set_xlabel(r'$x\ [c/\omega_{\rm pe}]$', labelpad = self.parent.MainParamDict['xLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
            if self.parent.MainParamDict['2DSlicePlane'] == 0:
                self.axes.set_ylabel(r'$y\ [c/\omega_{\rm pe}]$', labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
            if self.parent.MainParamDict['2DSlicePlane'] == 1:
                self.axes.set_ylabel(r'$z\ [c/\omega_{\rm pe}]$', labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
            if self.parent.MainParamDict['2DSlicePlane'] == 2:
                self.axes.set_xlabel(r'$y\ [c/\omega_{\rm pe}]$', labelpad = self.parent.MainParamDict['xLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
                self.axes.set_ylabel(r'$z\ [c/\omega_{\rm pe}]$', labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])


        else:
            # Do the 1D Plots
            if self.parent.MainParamDict['LinkSpatial'] != 0 and self.parent.MainParamDict['LinkSpatial'] != 3:
                if self.FigWrap.pos == self.parent.first_x:
                    self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]])
                else:
                    self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]],
                    sharex = self.parent.SubPlotList[self.parent.first_x[0]][self.parent.first_x[1]].graph.axes)
            else:
                self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]])

            # Make the 1-D plots
            if self.parent.MainParamDict['Average1D']:
                self.linedens = self.axes.plot(self.xaxis_values, np.average(self.dens.reshape(-1,self.dens.shape[-1]), axis = 0), color = self.dens_color)
            else:
                self.linedens = self.axes.plot(self.xaxis_values, self.dens[self.parent.zSlice,self.parent.ySlice,:], color = self.dens_color)

            if self.GetPlotParam('dens_type')==1:
                if self.parent.MainParamDict['Average1D']:
                    self.linedens[0].set_data(self.xaxis_values, np.average(self.densi.reshape(-1,self.densi.shape[-1]), axis = 0))
                else: # x-y plane
                    self.linedens[0].set_data(self.xaxis_values, self.densi[self.parent.zSlice,self.parent.ySlice,:])
            if self.GetPlotParam('dens_type')==2:
                if self.parent.MainParamDict['Average1D']:
                    self.linedens[0].set_data(self.xaxis_values, np.average(self.dense.reshape(-1,self.dense.shape[-1]), axis = 0))
                else: # x-y plane
                    self.linedens[0].set_data(self.xaxis_values, self.dense[self.parent.zSlice,self.parent.ySlice,:])
            if self.GetPlotParam('dens_type')==2:
                if self.parent.MainParamDict['Average1D']:
                    self.linedens[0].set_data(self.xaxis_values, np.average(self.rho.reshape(-1,self.rho.shape[-1]), axis = 0))
                else: # x-y plane
                    self.linedens[0].set_data(self.xaxis_values, self.rho[self.parent.zSlice,self.parent.ySlice,:])
            if self.GetPlotParam('normalize_density'):
                self.linedens[0].set_data(self.linedens[0].get_data()[0], self.linedens[0].get_data()[1]*self.ppc0**(-1))

            #### Set the ylims... there is a problem where it scales the ylims for the invisible lines:
            min_max = [self.linedens[0].get_data()[1].min(), self.linedens[0].get_data()[1].max()]
            dist = min_max[1]-min_max[0]
            min_max[0] -= 0.04*dist
            min_max[1] += 0.04*dist
            if self.GetPlotParam('stretch_colors'):
                tmp = max(abs(min_max[0]), abs(min_max[1]))
                min_max = [-tmp, tmp]
            self.axes.set_ylim(min_max)
            self.shock_line =self.axes.axvline(self.parent.shock_loc, linewidth = 1.5, linestyle = '--', color = self.parent.shock_color, path_effects=[PathEffects.Stroke(linewidth=2, foreground='k'),
                    PathEffects.Normal()])
            self.shock_line.set_visible(self.GetPlotParam('show_shock'))

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
                self.axes.set_xlim(self.xaxis_values[0],self.xaxis_values[-1])

            if self.GetPlotParam('set_v_min'):
                self.axes.set_ylim(bottom = self.GetPlotParam('v_min'))
            if self.GetPlotParam('set_v_max'):
                self.axes.set_ylim(top = self.GetPlotParam('v_max'))

            # Handle the axes labeling
            tmp_str = r'$\rm density$'
            if self.GetPlotParam('dens_type') == 1:
                tmp_str = r'$n_i$'
            if self.GetPlotParam('dens_type') == 2:
                tmp_str = r'$n_e$'
            if self.GetPlotParam('dens_type') == 3:
                tmp_str = r'$\rho$'

            if self.GetPlotParam('normalize_density'):
                tmp_str += r'$\ [n_0]$'
            self.axes.set_xlabel(r'$x\ [c/\omega_{\rm pe}]$', labelpad = self.parent.MainParamDict['xLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
            self.axes.set_ylabel(tmp_str, labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])

        if self.GetPlotParam('show_streamlines') and self.GetPlotParam('twoD'):
            streamlines.draw_streamlines(self)

        if self.GetPlotParam('show_cpu_domains'):
            self.FigWrap.SetCpuDomainLines()

    def refresh(self):
        '''This is a function that will be called only if self.axes already
        holds a density type plot. We only update things that have shown.  If
        hasn't changed, or isn't viewed, don't touch it. The difference between this and last
        time, is that we won't actually do any drawing in the plot. The plot
        will be redrawn after all subplots data is changed. '''

        # Main goal, only change what is showing..
        # First do the 1D plots, because it is simpler
        if self.GetPlotParam('twoD') == 0:
            if self.GetPlotParam('dens_type') == 0:
                if self.parent.MainParamDict['Average1D']:
                    self.linedens[0].set_data(self.xaxis_values, np.average(self.dens.reshape(-1,self.dens.shape[-1]), axis = 0))
                else: # x-y plane
                    self.linedens[0].set_data(self.xaxis_values, self.dens[self.parent.zSlice,self.parent.ySlice,:])

            elif self.GetPlotParam('dens_type')==1:
                if self.parent.MainParamDict['Average1D']:
                    self.linedens[0].set_data(self.xaxis_values, np.average(self.densi.reshape(-1,self.densi.shape[-1]), axis = 0))
                else:
                    self.linedens[0].set_data(self.xaxis_values, self.densi[self.parent.zSlice,self.parent.ySlice,:])
            elif self.GetPlotParam('dens_type')==2:
                if self.parent.MainParamDict['Average1D']:
                    self.linedens[0].set_data(self.xaxis_values, np.average(self.dense.reshape(-1,self.dense.shape[-1]), axis = 0))
                else:
                    self.linedens[0].set_data(self.xaxis_values, self.dense[self.parent.zSlice,self.parent.ySlice,:])
            elif self.GetPlotParam('dens_type')==3:
                if self.parent.MainParamDict['Average1D']:
                    self.linedens[0].set_data(self.xaxis_values, np.average(self.rho.reshape(-1,self.rho.shape[-1]), axis = 0))
                else:
                    self.linedens[0].set_data(self.xaxis_values, self.rho[self.parent.zSlice,self.parent.ySlice,:])
            if self.GetPlotParam('normalize_density'):
                self.linedens[0].set_data(self.linedens[0].get_data()[0], self.linedens[0].get_data()[1]*self.ppc0**(-1))

            #### Set the ylims...
            min_max = [self.linedens[0].get_data()[1].min(),self.linedens[0].get_data()[1].max()]
            dist = min_max[1]-min_max[0]
            min_max[0] -= 0.04*dist
            min_max[1] += 0.04*dist
            if self.GetPlotParam('stretch_colors'):
                tmp = max(abs(min_max[0]), abs(min_max[1]))
                min_max = [-tmp, tmp]
            self.axes.set_ylim(min_max)
            if self.GetPlotParam('set_v_min'):
                self.axes.set_ylim(bottom = self.GetPlotParam('v_min'))
            if self.GetPlotParam('set_v_max'):
                self.axes.set_ylim(top = self.GetPlotParam('v_max'))
            if self.GetPlotParam('show_shock'):
                self.shock_line.set_xdata([self.parent.shock_loc,self.parent.shock_loc])

            if self.parent.MainParamDict['SetxLim']:
                if self.parent.MainParamDict['xLimsRelative']:
                    self.axes.set_xlim(self.parent.MainParamDict['xLeft'] + self.parent.shock_loc,
                                       self.parent.MainParamDict['xRight'] + self.parent.shock_loc)
                else:
                    self.axes.set_xlim(self.parent.MainParamDict['xLeft'], self.parent.MainParamDict['xRight'])
            else:
                self.axes.set_xlim(self.xaxis_values[0], self.xaxis_values[-1])


        else: # Now refresh the plot if it is 2D
            if self.GetPlotParam('dens_type') == 0:
                if self.parent.MainParamDict['2DSlicePlane'] == 0: # x-y plane
                    self.cax.set_data(self.dens[self.parent.zSlice,:,:])
                elif self.parent.MainParamDict['2DSlicePlane'] == 1: # x-z plane
                    self.cax.set_data(self.dens[:,self.parent.ySlice,:])
                elif self.parent.MainParamDict['2DSlicePlane'] == 2: # y-z plane
                    self.cax.set_data(self.dens[:,:,self.parent.xSlice])



            elif self.GetPlotParam('dens_type')==1:
                if self.parent.MainParamDict['2DSlicePlane'] == 0: # x-y plane
                    self.cax.set_data(self.densi[self.parent.zSlice,:,:])
                elif self.parent.MainParamDict['2DSlicePlane'] == 1: # x-z plane
                    self.cax.set_data(self.densi[:,self.parent.ySlice,:])
                elif self.parent.MainParamDict['2DSlicePlane'] == 2: # y-z plane
                    self.cax.set_data(self.densi[:,:,self.parent.xSlice])
            elif self.GetPlotParam('dens_type')==2:
                if self.parent.MainParamDict['2DSlicePlane'] == 0: # x-y plane
                    self.cax.set_data(self.dense[self.parent.zSlice,:,:])
                elif self.parent.MainParamDict['2DSlicePlane'] == 1: # x-z plane
                    self.cax.set_data(self.dense[:,self.parent.ySlice,:])
                elif self.parent.MainParamDict['2DSlicePlane'] == 2: # y-z plane
                    self.cax.set_data(self.dense[:,:,self.parent.xSlice])
            elif self.GetPlotParam('dens_type')==3:
                if self.parent.MainParamDict['2DSlicePlane'] == 0: # x-y plane
                    self.cax.set_data(self.rho[self.parent.zSlice,:,:])
                elif self.parent.MainParamDict['2DSlicePlane'] == 1: # x-z plane
                    self.cax.set_data(self.rho[:,self.parent.ySlice,:])
                elif self.parent.MainParamDict['2DSlicePlane'] == 2: # y-z plane
                    self.cax.set_data(self.rho[:,:,self.parent.xSlice])


            if self.GetPlotParam('normalize_density'):
                self.cax.set_data(self.cax.get_array()/self.ppc0)


            self.ymin = 0
            self.ymax =  self.cax.get_array().shape[0]/self.c_omp*self.istep
            self.xmin = 0
            self.xmax =  self.cax.get_array().shape[1]/self.c_omp*self.istep
            self.cax.set_extent([self.xmin,self.xmax, self.ymin, self.ymax])
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

            if self.parent.MainParamDict['2DSlicePlane'] == 0:
                self.axes.set_ylabel(r'$y\ [c/\omega_{\rm pe}]$', labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
            if self.parent.MainParamDict['2DSlicePlane'] == 1:
                self.axes.set_ylabel(r'$z\ [c/\omega_{\rm pe}]$', labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
            if self.parent.MainParamDict['2DSlicePlane'] == 2:
                self.axes.set_xlabel(r'$y\ [c/\omega_{\rm pe}]$', labelpad = self.parent.MainParamDict['xLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
                self.axes.set_ylabel(r'$z\ [c/\omega_{\rm pe}]$', labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
            self.vmin = self.cax.get_array().min()
            if self.GetPlotParam('set_v_min'):
                self.vmin = self.GetPlotParam('v_min')
            self.vmax = self.cax.get_array().max()
            if self.GetPlotParam('set_v_max'):
                self.vmax = self.GetPlotParam('v_max')
            if self.GetPlotParam('UseDivCmap') and self.GetPlotParam('stretch_colors'):
                self.vmax = max(np.abs(self.vmin), self.vmax)
                self.vmin = -self.vmax
            self.cax.norm.vmin = self.vmin
            self.cax.norm.vmax = self.vmax


            if self.GetPlotParam('show_cbar'):
                self.CbarTickFormatter()
            if self.GetPlotParam('show_shock'):
                self.shockline_2d.set_xdata([self.parent.shock_loc,self.parent.shock_loc])

        if self.GetPlotParam('show_streamlines') and self.GetPlotParam('twoD'):
            streamlines.refresh_streamlines(self)

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
#                    self.cbar.set_data(self.gradient)
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
        self.frm = ttk.Frame(self)
        self.frm.pack(fill=Tk.BOTH, expand=True)
        self.protocol('WM_DELETE_WINDOW', self.OnClosing)
        self.bind('<Return>', self.TxtEnter)


        # Create the OptionMenu to chooses the Chart Type:
        self.InterpolVar = Tk.StringVar(self)
        self.InterpolVar.set(self.parent.GetPlotParam('interpolation')) # default value
        self.InterpolVar.trace('w', self.InterpolChanged)

        ttk.Label(self.frm, text="Interpolation Method:").grid(row=0, column = 2)
        InterplChooser = ttk.OptionMenu(self.frm, self.InterpolVar, self.parent.GetPlotParam('interpolation'), *tuple(self.parent.InterpolationMethods))
        InterplChooser.grid(row =0, column = 3, sticky = Tk.W + Tk.E)

        # Create the OptionMenu to chooses the Chart Type:
        self.ctypevar = Tk.StringVar(self)
        self.ctypevar.set(self.parent.chartType) # default value
        self.ctypevar.trace('w', self.ctypeChanged)

        ttk.Label(self.frm, text="Choose Chart Type:").grid(row=0, column = 0)
        ctypeChooser = ttk.OptionMenu(self.frm, self.ctypevar, self.parent.chartType, *tuple(self.parent.ChartTypes))
        ctypeChooser.grid(row =0, column = 1, sticky = Tk.W + Tk.E)


        self.TwoDVar = Tk.IntVar(self) # Create a var to track whether or not to plot in 2-D
        self.TwoDVar.set(self.parent.GetPlotParam('twoD'))
        cb = ttk.Checkbutton(self.frm, text = "Show in 2-D",
                variable = self.TwoDVar,
                command = self.Change2d)
        cb.grid(row = 1, sticky = Tk.W)

        # the Radiobox Control to choose the Field Type
        self.DensList = ['dens', 'dens_i', 'dens_e', 'rho']
        self.DensTypeVar  = Tk.IntVar()
        self.DensTypeVar.set(self.parent.GetPlotParam('dens_type'))

        ttk.Label(self.frm, text='Choose Density:').grid(row = 2, sticky = Tk.W)

        for i in range(len(self.DensList)):
            ttk.Radiobutton(self.frm,
                text=self.DensList[i],
                variable=self.DensTypeVar,
                command = self.RadioField,
                value=i).grid(row = 3+i//2, column = i-2*(i//2), sticky =Tk.W)


        # Control whether or not Cbar is shown
        self.CbarVar = Tk.IntVar()
        self.CbarVar.set(self.parent.GetPlotParam('show_cbar'))
        cb = ttk.Checkbutton(self.frm, text = "Show Color bar",
                        variable = self.CbarVar,
                        command = self.CbarHandler)
        cb.grid(row = 6, sticky = Tk.W)

        # show shock
        self.ShockVar = Tk.IntVar()
        self.ShockVar.set(self.parent.GetPlotParam('show_shock'))
        cb = ttk.Checkbutton(self.frm, text = "Show Shock",
                        variable = self.ShockVar,
                        command = self.ShockVarHandler)
        cb.grid(row = 6, column = 1, sticky = Tk.W)

        # Normalize Density Var
        self.NormDVar = Tk.IntVar()
        self.NormDVar.set(self.parent.GetPlotParam('normalize_density'))
        cb = ttk.Checkbutton(self.frm, text = "Normalize to ppc0",
                        variable = self.NormDVar,
                        command = self.NormPPCHandler)
        cb.grid(row = 7, sticky = Tk.W)

        # show labels
        self.ShowLabels = Tk.IntVar()
        self.ShowLabels.set(self.parent.GetPlotParam('show_labels'))
        cb = ttk.Checkbutton(self.frm, text = "Show Labels 2D",
                        variable = self.ShowLabels,
                        command = self.LabelHandler)
        cb.grid(row = 7, column = 1, sticky = Tk.W)
        # Control whether or not diverging cmap is used
        self.DivVar = Tk.IntVar()
        self.DivVar.set(self.parent.GetPlotParam('UseDivCmap'))
        cb = ttk.Checkbutton(self.frm, text = "Use Diverging Cmap",
                        variable = self.DivVar,
                        command = self.DivHandler)
        cb.grid(row = 8, sticky = Tk.W)

        # Use full div cmap
        self.StretchVar = Tk.IntVar()
        self.StretchVar.set(self.parent.GetPlotParam('stretch_colors'))
        cb = ttk.Checkbutton(self.frm, text = "Symmetric about zero",
                        variable = self.StretchVar,
                        command = self.StretchHandler)
        cb.grid(row = 9, column = 0, columnspan =2, sticky = Tk.W)

        self.CPUVar = Tk.IntVar()
        self.CPUVar.set(self.parent.GetPlotParam('show_cpu_domains'))
        cb = ttk.Checkbutton(self.frm, text = "Show CPU domains",
                        variable = self.CPUVar,
                        command = self.CPUVarHandler)
        cb.grid(row = 10, column = 0, sticky = Tk.W)

        # Create the OptionMenu to chooses the cnorm_type:
        self.cnormvar = Tk.StringVar(self)
        self.cnormvar.set(self.parent.chartType) # default value
        self.cnormvar.trace('w', self.cnormChanged)

        ttk.Label(self.frm, text="Choose Color Norm:").grid(row=6, column = 2)
        cnormChooser = ttk.OptionMenu(self.frm, self.cnormvar, self.parent.GetPlotParam('cnorm_type'), *tuple(['Pow', 'Linear']))
        cnormChooser.grid(row =6, column = 3, sticky = Tk.W + Tk.E)

        # Now the gamma of the pow norm
        self.powGamma = Tk.StringVar()
        self.powGamma.set(str(self.parent.GetPlotParam('cpow_num')))
        ttk.Label(self.frm, text ='gamma =').grid(row = 7, column = 2, sticky =Tk.E)
        ttk.Label(self.frm, text ='If cnorm is Pow =>').grid(row = 8, column = 2,columnspan = 2, sticky =Tk.W)
        ttk.Label(self.frm, text ='sign(data)*|data|**gamma').grid(row = 9, column = 2,columnspan = 2, sticky =Tk.E)

        self.GammaEnter = ttk.Entry(self.frm, textvariable=self.powGamma, width=7)
        self.GammaEnter.grid(row = 7, column = 3)


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


        cb = ttk.Checkbutton(self.frm, text ='Set dens min',
                        variable = self.setZminVar)
        cb.grid(row = 3, column = 2, sticky = Tk.W)
        self.ZminEnter = ttk.Entry(self.frm, textvariable=self.Zmin, width=7)
        self.ZminEnter.grid(row = 3, column = 3)

        cb = ttk.Checkbutton(self.frm, text ='Set dens max',
                        variable = self.setZmaxVar)
        cb.grid(row = 4, column = 2, sticky = Tk.W)

        self.ZmaxEnter = ttk.Entry(self.frm, textvariable=self.Zmax, width=7)
        self.ZmaxEnter.grid(row = 4, column = 3)

        streamlines.add_streamline_buttons(self, self.parent, starting_row=11)

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

    def DivHandler(self, *args):
        if self.parent.GetPlotParam('UseDivCmap')== self.DivVar.get():
            pass
        elif self.parent.GetPlotParam('twoD'):
            self.parent.SetPlotParam('UseDivCmap', self.DivVar.get(),  NeedsRedraw = True)
        else:
            self.parent.SetPlotParam('UseDivCmap', self.DivVar.get(),  update_plot = False)


    def StretchHandler(self, *args):
        if self.parent.GetPlotParam('stretch_colors') == self.StretchVar.get():
            pass
        elif self.parent.GetPlotParam('twoD'):
            self.parent.SetPlotParam('stretch_colors', self.StretchVar.get(), NeedsRedraw = True)
        else:
            self.parent.SetPlotParam('stretch_colors', self.StretchVar.get(), update_plot = True)

    def cnormChanged(self, *args):
        if self.parent.GetPlotParam('cnorm_type')== self.cnormvar.get():
            pass
        elif self.parent.GetPlotParam('twoD'):
            self.parent.SetPlotParam('cnorm_type', self.cnormvar.get(), NeedsRedraw = True)
        else:
            self.parent.SetPlotParam('cnorm_type', self.cnormvar.get(), update_plot = False)


    def NormPPCHandler(self, *args):
        if self.parent.GetPlotParam('normalize_density')== self.NormDVar.get():
            pass
        elif np.isnan(self.parent.ppc0):
            self.NormDVar.set(self.parent.GetPlotParam('normalize_density'))
        else:
            # First see if the plot is 2d and change the labels if necessary
            if self.parent.GetPlotParam('twoD')==1:
                tmp_str = ''
                if self.parent.GetPlotParam('dens_type') == 0:
                    tmp_str += r'$n$'
                if self.parent.GetPlotParam('dens_type') == 1:
                    tmp_str += r'$n_i$'
                if self.parent.GetPlotParam('dens_type') == 2:
                    tmp_str += r'$n_e$'
                if self.parent.GetPlotParam('dens_type') == 3:
                    tmp_str += r'$\rho$'
                if self.NormDVar.get():
                    tmp_str += r'$\ [n_0]$'
                self.parent.an_2d.set_text(tmp_str)
            # Now for the 1d label
            else:
                if self.parent.GetPlotParam('dens_type') == 0:
                    if self.NormDVar.get():
                        self.parent.axes.set_ylabel(r'${\rm density}\ [n_0]$', size = self.parent.parent.MainParamDict['AxLabelSize'])
                    else:
                        self.parent.axes.set_ylabel('density', size = self.parent.parent.MainParamDict['AxLabelSize'])
                if self.parent.GetPlotParam('dens_type') == 1:
                    if self.NormDVar.get():
                        self.parent.axes.set_ylabel(r'$n_i\ [n_0]$', size = self.parent.parent.MainParamDict['AxLabelSize'])
                    else:
                        self.parent.axes.set_ylabel(r'$n_i$', size = self.parent.parent.MainParamDict['AxLabelSize'])
                if self.parent.GetPlotParam('dens_type') == 2:
                    if self.NormDVar.get():
                        self.parent.axes.set_ylabel(r'$n_e\ [n_0]$', size = self.parent.parent.MainParamDict['AxLabelSize'])
                    else:
                        self.parent.axes.set_ylabel(r'$n_e$', size = self.parent.parent.MainParamDict['AxLabelSize'])
                if self.parent.GetPlotParam('dens_type') == 3:
                    if self.NormDVar.get():
                        self.parent.axes.set_ylabel(r'$\rho\ [n_0]$', size = self.parent.parent.MainParamDict['AxLabelSize'])
                    else:
                        self.parent.axes.set_ylabel(r'$\rho$', size = self.parent.parent.MainParamDict['AxLabelSize'])

            self.parent.SetPlotParam('normalize_density', self.NormDVar.get())#, update_plot = self.parent.GetPlotParam('dens_type')==0)

    def LabelHandler(self, *args):
        if self.parent.GetPlotParam('show_labels')== self.ShowLabels.get():
            pass
        else:
            if self.parent.GetPlotParam('twoD'):
                self.parent.an_2d.set_visible(self.ShowLabels.get())

            self.parent.SetPlotParam('show_labels', self.ShowLabels.get(), update_plot =self.parent.GetPlotParam('twoD'))

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

    def TxtEnter(self, e):
        streamlines.streamlines_callback(self, update_plot=True)
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


    def RadioField(self):
        if self.DensTypeVar.get() == self.parent.GetPlotParam('dens_type'):
            pass
        else:
            if self.parent.GetPlotParam('twoD'):
                if self.DensTypeVar.get() == 0:
                    if self.parent.GetPlotParam('normalize_density'):
                        self.parent.an_2d.set_text(r'$n/n_0$')
                    else:
                        self.parent.an_2d.set_text(r'$n$')

                elif self.DensTypeVar.get() == 1:
                    if self.parent.GetPlotParam('normalize_density'):
                        self.parent.an_2d.set_text(r'$n_i/n_0$')
                    else:
                        self.parent.an_2d.set_text(r'$n_i$')

                elif self.DensTypeVar.get() == 2:
                    if self.parent.GetPlotParam('normalize_density'):
                        self.parent.an_2d.set_text(r'$n_e/n_0$')
                    else:
                        self.parent.an_2d.set_text(r'$n_e$')
                elif self.DensTypeVar.get() == 3:
                    if self.parent.GetPlotParam('normalize_density'):
                        self.parent.an_2d.set_text(r'$\rho/n_0$')
                    else:
                        self.parent.an_2d.set_text(r'$\rho$')

            else:
                if self.DensTypeVar.get() == 0:
                    if self.parent.GetPlotParam('normalize_density'):
                        self.parent.axes.set_ylabel(r'${\rm density}\ [n_0]$', size = self.parent.parent.MainParamDict['AxLabelSize'])
                    else:
                        self.parent.axes.set_ylabel(r'$\rm density$', size = self.parent.parent.MainParamDict['AxLabelSize'])


                elif self.DensTypeVar.get() ==1:
                    if self.parent.GetPlotParam('normalize_density'):
                        self.parent.axes.set_ylabel(r'$n_i\ [n_0]$', size = self.parent.parent.MainParamDict['AxLabelSize'])
                    else:
                        self.parent.axes.set_ylabel(r'$n_i$', size = self.parent.parent.MainParamDict['AxLabelSize'])

                elif self.DensTypeVar.get() ==2:
                    if self.parent.GetPlotParam('normalize_density'):
                        self.parent.axes.set_ylabel(r'$n_e\ [n_0]$', size = self.parent.parent.MainParamDict['AxLabelSize'])
                    else:
                        self.parent.axes.set_ylabel(r'$n_e$', size = self.parent.parent.MainParamDict['AxLabelSize'])
                elif self.DensTypeVar.get() ==3:
                    if self.parent.GetPlotParam('normalize_density'):
                        self.parent.axes.set_ylabel(r'$\rho\ [n_0]$', size = self.parent.parent.MainParamDict['AxLabelSize'])
                    else:
                        self.parent.axes.set_ylabel(r'$\rho$', size = self.parent.parent.MainParamDict['AxLabelSize'])


            self.parent.SetPlotParam('dens_type', self.DensTypeVar.get())


    def OnClosing(self):
        self.parent.settings_window = None
        self.destroy()
