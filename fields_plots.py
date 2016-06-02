#!/usr/bin/env python
import Tkinter as Tk
import ttk as ttk
import matplotlib
import numpy as np
import numpy.ma as ma
import new_cmaps
from new_cnorms import PowerNormWithNeg
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as PathEffects

class FieldsPanel:
    # A dictionary of all of the parameters for this plot with the default parameters

    plot_param_dict = {'twoD': 0,
                       'field_type': 0, #0 = B-Field, 1 = E-field
                       'show_x' : 1,
                       'show_y' : 1,
                       'show_z' : 1,
                       'show_cbar': True,
                       'v_min': 0,
                       'v_max' : 10,
                       'set_v_min': False,
                       'set_v_max': False,
                       'show_shock' : False,
                       'OutlineText': True,
                       'spatial_x': True,
                       'spatial_y': None,
                       'normalize_fields': True, # Normalize fields to their upstream values
                       'cnorm_type': 'Linear', # Colormap norm;  options are Log, Pow or Linear
                       'cpow_num': 0.6, # Used in the PowerNorm
                       'interpolation': 'nearest',
                       'cmap': None # If cmap is none, the plot will inherit the parent's cmap
                       }

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
        self.InterpolationMethods = ['nearest', 'bilinear', 'bicubic', 'spline16',
            'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric',
            'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']

    def ChangePlotType(self, str_arg):
        self.FigWrap.ChangeGraph(str_arg)

    def norm(self, vmin=None, vmax=None):
        if self.GetPlotParam('cnorm_type') =="Linear":
            return mcolors.Normalize(vmin, vmax)
        elif self.GetPlotParam('cnorm_type') == "Log":
            return  mcolors.LogNorm(vmin, vmax)
        else:
            return  PowerNormWithNeg(self.FigWrap.GetPlotParam('cpow_num'), vmin, vmax)

    def set_plot_keys(self):
        '''A helper function that will insure that each hdf5 file will only be
        opened once per time step'''
        # First make sure that omega_plasma & xi is loaded so we can fix the
        # x & y distances.


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

        return self.arrs_needed

    def LoadData(self):
        ''' A Helper function that loads the data for the plot'''
        # First see of the x_axis and y_axis values have already been calculated
        # and stored in the DataDict for this time step
        self.c_omp = self.FigWrap.LoadKey('c_omp')[0]
        self.istep = self.FigWrap.LoadKey('istep')[0]
        if self.GetPlotParam('cmap') is None:
            self.cmap = self.parent.cmap
        else:
            self.cmap = self.GetPlotParam('cmap')
        self.xcolor = new_cmaps.cmaps[self.parent.cmap](0.2)
        self.ycolor = new_cmaps.cmaps[self.parent.cmap](0.5)
        self.zcolor = new_cmaps.cmaps[self.parent.cmap](0.8)

        if np.isnan(self.parent.btheta):
            self.SetPlotParam('normalize_fields', 0, update_plot = False)
        # see if the axis values are saved in the data dict
        if 'xaxis_values' in self.parent.DataDict.keys():
            self.xaxis_values = self.parent.DataDict['xaxis_values']
        else:
            # x-values haven't been calculated yet, generate them then save them to the dictionary for later.
            self.xaxis_values = np.arange(self.FigWrap.LoadKey('bx')[0,:,:].shape[1])/self.c_omp*self.istep
        #            print self.xaxis_values
            self.parent.DataDict['xaxis_values'] = np.copy(self.xaxis_values)

        # Now the y_values Not needed so commenting out
        '''
        if 'yaxis_values' in self.parent.DataDict.keys():
            self.yaxis_values = self.parent.DataDict('y_values')
        else:
            self.yaxis_values = np.arange(self.FigWrap.LoadKey('bx')[0,:,:].shape[0])/self.c_omp*self.istep
            self.parent.DataDict('yaxis_values') = self.yaxis_values
        '''


        if self.GetPlotParam('field_type') == 0: # Load the B-Field
            #
            self.fx = self.FigWrap.LoadKey('bx')[0,:,:]
            self.fy = self.FigWrap.LoadKey('by')[0,:,:]
            self.fz = self.FigWrap.LoadKey('bz')[0,:,:]

            kxwarg = 'bxmin_max'
            kywarg = 'bymin_max'
            kzwarg = 'bzmin_max'

            if self.GetPlotParam('normalize_fields'):
                self.fx = self.fx/self.parent.b0
                self.fy = self.fy/self.parent.b0
                self.fz = self.fz/self.parent.b0
                kxwarg += '_normalized'
                kywarg += '_normalized'
                kzwarg += '_normalized'
            self.oneDslice = self.fx.shape[0]/2

            # Have we already calculated min/max?

            if kxwarg in self.parent.DataDict.keys():
                self.bxmin_max = self.parent.DataDict[kxwarg]

            else:
                self.bxmin_max =  self.min_max_finder(self.fx)
                self.parent.DataDict[kxwarg] = list(self.bxmin_max)

            if kywarg in self.parent.DataDict.keys():
                self.bymin_max = self.parent.DataDict[kywarg]

            else:
                self.bymin_max = self.min_max_finder(self.fy)
                self.parent.DataDict[kywarg] = list(self.bymin_max)

            if kzwarg in self.parent.DataDict.keys():
                self.bzmin_max = self.parent.DataDict[kzwarg]

            else:
                self.bzmin_max = self.min_max_finder(self.fz)
                self.parent.DataDict[kzwarg] = list(self.bzmin_max)


        if self.GetPlotParam('field_type') == 1: # Load the e-Field
            self.fx = self.FigWrap.LoadKey('ex')[0,:,:]
            self.fy = self.FigWrap.LoadKey('ey')[0,:,:]
            self.fz = self.FigWrap.LoadKey('ez')[0,:,:]

            kxwarg = 'exmin_max'
            kywarg = 'eymin_max'
            kzwarg = 'ezmin_max'

            if self.GetPlotParam('normalize_fields'):
                self.fx = self.fx/self.parent.e0
                self.fy = self.fy/self.parent.e0
                self.fz = self.fz/self.parent.e0
                kxwarg += '_normalized'
                kywarg += '_normalized'
                kzwarg += '_normalized'

            self.oneDslice = self.fx.shape[0]/2



            # Have we already calculated min/max?
            if kxwarg in self.parent.DataDict.keys():
                self.exmin_max = self.parent.DataDict[kxwarg]

            else:
                self.exmin_max = self.min_max_finder(self.fx)
                self.parent.DataDict[kxwarg] = list(self.exmin_max)

            if kywarg in self.parent.DataDict.keys():
                self.eymin_max = self.parent.DataDict[kywarg]

            else:
                self.eymin_max = self.min_max_finder(self.fy)
                self.parent.DataDict[kywarg] = list(self.eymin_max)

            if kzwarg in self.parent.DataDict.keys():
                self.ezmin_max = self.parent.DataDict[kzwarg]

            else:
                self.ezmin_max = self.min_max_finder(self.fz)
                self.parent.DataDict[kzwarg] = list(self.ezmin_max)

    def min_max_finder(self, arr):
        # find 1d lims
        oneD_lims = [arr[self.oneDslice,:].min(), arr[self.oneDslice,:].max()]
        # now give it a bit of spacing, a 4% percent difference of the distance
        dist = oneD_lims[1]-oneD_lims[0]
        oneD_lims[0] -= 0.04*dist
        oneD_lims[1] += 0.04*dist
        twoD_lims = [arr.min(), arr.max()]
        return [oneD_lims, twoD_lims]

    def draw(self):

        ''' A function that draws the data. In the interest in speeding up the
        code, draw should only be called when you want to recreate the whole
        figure, i.e. it  will be slow. Most times you will only want to update
        what has changed in the figure. This will be done in a function called
        refresh, that should be much much faster.'''

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
        self.gs = gridspec.GridSpecFromSubplotSpec(100,100, subplot_spec = self.parent.gs0[self.FigWrap.pos])

        # Now that the data is loaded, start making the plots
        if self.GetPlotParam('twoD'):
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
            if self.GetPlotParam('show_x'):
                self.zval = self.fx
                if self.GetPlotParam('normalize_fields'):
                    self.two_d_labels = (r'$B_x/B_0$', r'$E_x/E_0$')
                else:
                    self.two_d_labels = (r'$B_x$', r'$E_x$')
                # set the other plot values to zero in the PlotParams
                self.SetPlotParam('show_y', 0, update_plot = False)
                self.SetPlotParam('show_z', 0, update_plot = False)

            elif self.GetPlotParam('show_y'):
                self.zval = self.fy
                if self.GetPlotParam('normalize_fields'):
                    self.two_d_labels = (r'$B_y/B_0$', r'$E_y/E_0$')
                else:
                    self.two_d_labels = (r'$B_y$', r'$E_y$')


                # set the other plot values to zero in the PlotParams
                self.SetPlotParam('show_x', 0, update_plot = False)
                self.SetPlotParam('show_z', 0, update_plot = False)

            else:
                # make sure z is loaded, (something has to be)
                # set the other plot values to zero in the PlotParams
                self.SetPlotParam('show_x', 0, update_plot = False)
                self.SetPlotParam('show_y', 0, update_plot = False)

                # set show_z to 1 to be consistent
                self.SetPlotParam('show_z', 1, update_plot = False)

                self.zval = self.fz
                if self.GetPlotParam('normalize_fields'):
                    self.two_d_labels = (r'$B_z/B_0$', r'$E_z/E_0$')
                else:
                    self.two_d_labels = (r'$B_z$', r'$E_z$')


            self.ymin = 0
            self.ymax =  self.zval.shape[0]/self.c_omp*self.istep
            self.xmin = 0
            self.xmax =  self.zval.shape[1]/self.c_omp*self.istep

            self.vmin = self.zval.min()
            if self.GetPlotParam('set_v_min'):
                self.vmin = self.GetPlotParam('v_min')
            self.vmax = self.zval.max()
            if self.GetPlotParam('set_v_max'):
                self.vmax = self.GetPlotParam('v_max')

            if self.parent.plot_aspect:
                self.cax = self.axes.imshow(self.zval, norm = self.norm(), origin = 'lower')
            else:
                self.cax = self.axes.imshow(self.zval, origin = 'lower', norm = self.norm(),
                                            aspect= 'auto')
            self.cax.set_cmap(new_cmaps.cmaps[self.cmap])
            self.cax.set_extent([self.xmin,self.xmax, self.ymin, self.ymax])
            self.cax.norm.vmin = self.vmin
            self.cax.norm.vmax = self.vmax

            self.axes.add_artist(self.cax)


            self.TwoDan = self.axes.annotate(self.two_d_labels[self.GetPlotParam('field_type')],
                            xy = (0.9,.9),
                            xycoords= 'axes fraction',
                            color = 'white',
                            **self.annotate_kwargs)

            self.axes.set_axis_bgcolor('lightgrey')

            self.axC = self.figure.add_subplot(self.gs[:4,:])
            self.cbar = self.axC.imshow(self.gradient, aspect='auto',
                                        cmap=new_cmaps.cmaps[self.cmap])

            # Make the colobar axis more like the real colorbar
            self.cbar.set_extent([0, 1.0, 0, 1.0])
            self.axC.tick_params(axis='x',
                            which = 'both', # bothe major and minor ticks
                            top = 'off', # turn off top ticks
                            labelsize=self.parent.num_font_size)

            self.axC.tick_params(axis='y',          # changes apply to the y-axis
                            which='both',      # both major and minor ticks are affected
                            left='off',      # ticks along the bottom edge are off
                            right='off',         # ticks along the top edge are off
                            labelleft='off')

            if self.GetPlotParam('show_cbar') == 0:
                self.axC.set_visible = False
            else:
                self.CbarTickFormatter()


            self.shockline_2d = self.axes.axvline(self.parent.shock_loc, linewidth = 1.5, linestyle = '--', color = self.parent.shock_color, path_effects=[PathEffects.Stroke(linewidth=2, foreground='k'),
                                    PathEffects.Normal()])
            self.shockline_2d.set_visible(self.GetPlotParam('show_shock'))

            self.axes.set_axis_bgcolor('lightgrey')
            self.axes.tick_params(labelsize = self.parent.num_font_size, color=tick_color)

            if self.parent.xlim[0]:
                self.axes.set_xlim(self.parent.xlim[1],self.parent.xlim[2])
            else:
                self.axes.set_xlim(self.xmin, self.xmax)
            if self.parent.ylim[0]:
                self.axes.set_ylim(self.parent.ylim[1],self.parent.ylim[2])
            else:
                self.axes.set_ylim(self.ymin, self.ymax)
            self.axes.set_xlabel(r'$x\ [c/\omega_{\rm pe}]$', labelpad = self.parent.xlabel_pad, color = 'black')
            self.axes.set_ylabel(r'$y\ [c/\omega_{\rm pe}]$', labelpad = self.parent.ylabel_pad, color = 'black')

        else:
            if self.parent.LinkSpatial != 0 and self.parent.LinkSpatial != 3:
                if self.FigWrap.pos == self.parent.first_x:
                    self.axes = self.figure.add_subplot(self.gs[18:92,:])
                else:
                    self.axes = self.figure.add_subplot(self.gs[18:92,:],
                    sharex = self.parent.SubPlotList[self.parent.first_x[0]][self.parent.first_x[1]].graph.axes)
            else:
                self.axes = self.figure.add_subplot(self.gs[18:92,:])

            self.annotate_pos = [0.8,0.9]
            self.linex = self.axes.plot(self.xaxis_values, self.fx[self.oneDslice,:], color = self.xcolor)
            self.linex[0].set_visible(self.GetPlotParam('show_x'))
            if self.GetPlotParam('field_type') == 0:
                tmp_str = r'$B_x$'
            else:
                tmp_str = r'$E_x$'
            self.anx = self.axes.annotate(tmp_str, xy = self.annotate_pos,
                            xycoords = 'axes fraction',
                            color = self.xcolor,
                            **self.annotate_kwargs)
            self.anx.set_visible(self.GetPlotParam('show_x'))

            self.annotate_pos[0] += .08
            self.liney = self.axes.plot(self.xaxis_values, self.fy[self.oneDslice,:], color = self.ycolor)
            self.liney[0].set_visible(self.GetPlotParam('show_y'))
            if self.GetPlotParam('field_type') == 0:
                tmp_str = r'$B_y$'
            else:
                tmp_str = r'$E_y$'
            self.any =self.axes.annotate(tmp_str, xy = self.annotate_pos,
                            xycoords= 'axes fraction',
                            color = self.ycolor,
                            **self.annotate_kwargs)

            self.any.set_visible(self.GetPlotParam('show_y'))
            self.annotate_pos[0] += .08


            self.linez = self.axes.plot(self.xaxis_values, self.fz[self.oneDslice,:], color = self.zcolor)
            self.linez[0].set_visible(self.GetPlotParam('show_z'))

            if self.GetPlotParam('field_type') == 0:
                tmp_str = r'$B_z$'
            else:
                tmp_str = r'$E_z$'
            self.anz = self.axes.annotate(tmp_str, xy = self.annotate_pos,
                                xycoords= 'axes fraction',
                                color = self.zcolor,
                                **self.annotate_kwargs
                                )
            self.anz.set_visible(self.GetPlotParam('show_z'))

            self.shock_line = self.axes.axvline(self.parent.shock_loc, linewidth = 1.5, linestyle = '--', color = self.parent.shock_color, path_effects=[PathEffects.Stroke(linewidth=2, foreground='k'),
                        PathEffects.Normal()])

            self.shock_line.set_visible(self.GetPlotParam('show_shock'))

            self.axes.set_axis_bgcolor('lightgrey')
            self.axes.tick_params(labelsize = self.parent.num_font_size, color=tick_color)#, tick1On= False, tick2On= False)


            self.axes.relim()
            self.axes.autoscale_view(scaley = True)

            if self.parent.xlim[0]:
                self.axes.set_xlim(self.parent.xlim[1],self.parent.xlim[2])
            else:
                self.axes.set_xlim(self.xaxis_values[0],self.xaxis_values[-1])


            if self.GetPlotParam('set_v_min'):
                self.axes.set_ylim(ymin = self.GetPlotParam('v_min'))
            if self.GetPlotParam('set_v_max'):
                self.axes.set_ylim(ymax = self.GetPlotParam('v_max'))

            self.axes.set_xlabel(r'$x\ [c/\omega_{\rm pe}]$', labelpad = self.parent.xlabel_pad, color = 'black')
            if self.GetPlotParam('field_type') == 0:
                if self.GetPlotParam('normalize_fields'):
                    self.axes.set_ylabel(r'$B/B_0$', labelpad = self.parent.ylabel_pad, color = 'black')
                else:
                    self.axes.set_ylabel(r'$B$', labelpad = self.parent.ylabel_pad, color = 'black')
            else:
                if self.GetPlotParam('normalize_fields'):
                    self.axes.set_ylabel(r'$E/E_0$', labelpad = self.parent.ylabel_pad, color = 'black')
                else:
                    self.axes.set_ylabel(r'$E$', labelpad = self.parent.ylabel_pad, color = 'black')

    def refresh(self):

        '''This is a function that will be called only if self.axes already
        holds a fields type plot. We only update things that have changed & are
        shown.  If hasn't changed or isn't shown, don't touch it. The difference
        between this and last time, is that we won't actually do any drawing in
        the plot. The plot will be redrawn after all subplots are refreshed. '''


        # Main goal, only change what is showing..
        # First do the 1D plots, because it is simpler
        if self.GetPlotParam('twoD') == 0:
            self.line_ymin_max = [np.inf, -np.inf]
            if self.GetPlotParam('show_x'):
                self.linex[0].set_data(self.xaxis_values, self.fx[self.oneDslice,:])
                if self.GetPlotParam('field_type') == 0:
                    self.line_ymin_max[0] = min(self.bxmin_max[0][0], self.line_ymin_max[0])
                    self.line_ymin_max[1] = max(self.bxmin_max[0][1], self.line_ymin_max[1])
                else:
                    self.line_ymin_max[0] = min(self.exmin_max[0][0], self.line_ymin_max[0])
                    self.line_ymin_max[1] = max(self.exmin_max[0][1], self.line_ymin_max[1])

            if self.GetPlotParam('show_y'):
                self.liney[0].set_data(self.xaxis_values, self.fy[self.oneDslice,:])
                if self.GetPlotParam('field_type') == 0:
                    self.line_ymin_max[0] = min(self.bymin_max[0][0], self.line_ymin_max[0])
                    self.line_ymin_max[1] = max(self.bymin_max[0][1], self.line_ymin_max[1])
                else:
                    self.line_ymin_max[0] = min(self.eymin_max[0][0], self.line_ymin_max[0])
                    self.line_ymin_max[1] = max(self.eymin_max[0][1], self.line_ymin_max[1])


            if self.GetPlotParam('show_z'):
                self.linez[0].set_data(self.xaxis_values, self.fz[self.oneDslice,:])
                if self.GetPlotParam('field_type') == 0:
                    self.line_ymin_max[0] = min(self.bzmin_max[0][0], self.line_ymin_max[0])
                    self.line_ymin_max[1] = max(self.bzmin_max[0][1], self.line_ymin_max[1])
                else:
                    self.line_ymin_max[0] = min(self.ezmin_max[0][0], self.line_ymin_max[0])
                    self.line_ymin_max[1] = max(self.ezmin_max[0][1], self.line_ymin_max[1])

            if not self.GetPlotParam('show_x') and not self.GetPlotParam('show_y') and not self.GetPlotParam('show_z'):
                self.line_ymin_max = [0,1]

            self.axes.set_ylim(self.line_ymin_max)
            if self.GetPlotParam('show_shock'):
                self.shock_line.set_xdata([self.parent.shock_loc,self.parent.shock_loc])

            if self.parent.xlim[0]:
                self.axes.set_xlim(self.parent.xlim[1],self.parent.xlim[2])
            else:
                self.axes.set_xlim(self.xaxis_values[0], self.xaxis_values[-1])

            if self.GetPlotParam('set_v_min'):
                self.axes.set_ylim(ymin = self.GetPlotParam('v_min'))
            if self.GetPlotParam('set_v_max'):
                self.axes.set_ylim(ymax = self.GetPlotParam('v_max'))

            # Code for trying blitting, not worth it.
#            self.axes.draw_artist(self.axes.patch)
#            self.axes.draw_artist(self.linex[0])
#            self.axes.draw_artist(self.liney[0])
#            self.axes.draw_artist(self.linez[0])
#            self.axes.draw_artist(self.anx)
#            self.axes.draw_artist(self.any)
#            self.axes.draw_artist(self.anz)

        else: # Now refresh the plot if it is 2D
            if self.GetPlotParam('show_x'):
                self.cax.set_data(self.fx)
                self.ymin = 0
                self.ymax =  self.fx.shape[0]/self.c_omp*self.istep
                self.xmin = 0
                self.xmax = self.xaxis_values[-1]
                if self.GetPlotParam('field_type') == 0:
                    self.clims = np.copy(self.bxmin_max[1])
                    if self.GetPlotParam('normalize_fields'):
                        self.TwoDan.set_text(r'$B_x/B_0$')
                    else:
                        self.TwoDan.set_text(r'$B_x$')
                else:
                    self.clims = np.copy(self.exmin_max[1])
                    if self.GetPlotParam('normalize_fields'):
                        self.TwoDan.set_text(r'$E_x/E_0$')
                    else:
                        self.TwoDan.set_text(r'$E_x$')

            if self.GetPlotParam('show_y'):
                self.cax.set_data(self.fy)
                self.ymin = 0
                self.ymax =  self.fy.shape[0]/self.c_omp*self.istep
                self.xmin = 0
                self.xmax = self.xaxis_values[-1]
                if self.GetPlotParam('field_type') == 0:
                    self.clims = np.copy(self.bymin_max[1])
                    if self.GetPlotParam('normalize_fields'):
                        self.TwoDan.set_text(r'$B_y/B_0$')
                    else:
                        self.TwoDan.set_text(r'$B_y$')
                else:
                    self.clims = np.copy(self.eymin_max[1])
                    if self.GetPlotParam('normalize_fields'):
                        self.TwoDan.set_text(r'$E_y/E_0$')
                    else:
                        self.TwoDan.set_text(r'$E_y$')
            if self.GetPlotParam('show_z'):
                self.cax.set_data(self.fz)
                self.ymin = 0
                self.ymax =  self.fz.shape[0]/self.c_omp*self.istep
                self.xmin = 0
                self.xmax =  self.xaxis_values[-1]
                if self.GetPlotParam('field_type') == 0:
                    self.clims = np.copy(self.bzmin_max[1])
                    if self.GetPlotParam('normalize_fields'):
                        self.TwoDan.set_text(r'$B_z/B_0$')
                    else:
                        self.TwoDan.set_text(r'$B_z$')
                else:
                    self.clims = np.copy(self.ezmin_max[1])
                    if self.GetPlotParam('normalize_fields'):
                        self.TwoDan.set_text(r'$E_z/E_0$')
                    else:
                        self.TwoDan.set_text(r'$E_z$')
            if self.parent.xlim[0]:
                self.axes.set_xlim(self.parent.xlim[1],self.parent.xlim[2])
            else:
                self.axes.set_xlim(self.xmin,self.xmax)
            if self.parent.ylim[0]:
                self.axes.set_ylim(self.parent.ylim[1],self.parent.ylim[2])
            else:
                self.axes.set_ylim(self.ymin,self.ymax)

            self.cax.set_extent([self.xmin, self.xmax, self.ymin, self.ymax])


            if self.GetPlotParam('set_v_min'):
                self.clims[0] =  self.GetPlotParam('v_min')
            if self.GetPlotParam('set_v_max'):
                self.clims[1] =  self.GetPlotParam('v_max')
            self.cax.set_clim(self.clims)

            self.CbarTickFormatter()

            if self.GetPlotParam('show_shock'):
                self.shockline_2d.set_xdata([self.parent.shock_loc,self.parent.shock_loc])
            self.axes.set_ylabel(r'$y\ [c/\omega_{\rm pe}]$', labelpad = self.parent.ylabel_pad, color = 'black')
            #self.axes.draw_artist(self.axes.patch)
            #self.axes.draw_artist(self.cax)
            #self.axes.draw_artist(self.axes.xaxis)

    def CbarTickFormatter(self):
        ''' A helper function that sets the cbar ticks & labels. This used to be
        easier, but because I am no longer using the colorbar class i have to do
        stuff manually.'''
        self.axC.set_xlim(0,1)
        self.axC.set_ylim(0,1)
        clim = np.copy(self.cax.get_clim())
        if self.GetPlotParam('show_cbar'):
            if self.GetPlotParam('cnorm_type') == "Log":
                ctick_range = np.linspace(0, 1, 6)
                data_val = np.logspace(np.log10(clim[0]),np.log10(clim[1]), 6)
                self.axC.set_xticks(ctick_range)
                ctick_labels = []
                for elm in data_val:
                    tmp_s = '%.2f' % elm
                    ctick_labels.append(tmp_s)
#                print ctick_labels
                self.axC.set_xticklabels(ctick_labels)

            if self.GetPlotParam('cnorm_type') == "Pow":
                ctick_range = np.linspace(0, 1, 6)
                pow_min = np.sign(clim[0])*np.abs(clim[0])**self.FigWrap.GetPlotParam('cpow_num')
                pow_max = np.sign(clim[1])*np.abs(clim[1])**self.FigWrap.GetPlotParam('cpow_num')
                data_val = np.sign(np.linspace(pow_min, pow_max, 6))*np.abs(np.linspace(pow_min, pow_max, 6))**(1.0/self.FigWrap.GetPlotParam('cpow_num'))
                self.axC.set_xticks(ctick_range)

                ctick_labels = []
                for elm in data_val:
                    tmp_s = '%.2f' % elm
                    ctick_labels.append(tmp_s)
                self.axC.set_xticklabels(ctick_labels)

            if self.GetPlotParam('cnorm_type') == "Linear":
                ctick_range = np.linspace(0, 1, 6)
                data_val = np.linspace(clim[0],clim[1], 6)
                self.axC.set_xticks(ctick_range)

                ctick_labels = []
                for elm in data_val:
                    if np.abs(elm)<1E-2:
                        tmp_s = '0'
                    else:
                        tmp_s = '%.2f' % elm
                    ctick_labels.append(tmp_s)
                self.axC.set_xticklabels(ctick_labels)

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

        self.wm_title('Field Plot (%d,%d) Settings' % self.parent.FigWrap.pos)
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
        self.FieldList = ['B Field', 'E field']
        self.FieldTypeVar  = Tk.IntVar()
        self.FieldTypeVar.set(self.parent.GetPlotParam('field_type'))

        ttk.Label(frm, text='Choose Field:').grid(row = 2, sticky = Tk.W)

        for i in range(len(self.FieldList)):
            ttk.Radiobutton(frm,
                text=self.FieldList[i],
                variable=self.FieldTypeVar,
                command = self.RadioField,
                value=i).grid(row = 3+i, sticky =Tk.W)

        # the Check boxes for the dimension
        ttk.Label(frm, text='Dimenison:').grid(row = 1, column = 1, sticky = Tk.W)

        self.ShowXVar = Tk.IntVar(self) # Create a var to track whether or not to plot in 2-D
        self.ShowXVar.set(self.parent.GetPlotParam('show_x'))
        cb = ttk.Checkbutton(frm, text = "Show x",
            variable = self.ShowXVar,
            command = self.Selector)
        cb.grid(row = 2, column = 1, sticky = Tk.W)

        self.ShowYVar = Tk.IntVar(self) # Create a var to track whether or not to plot in 2-D
        self.ShowYVar.set(self.parent.GetPlotParam('show_y'))
        cb = ttk.Checkbutton(frm, text = "Show y",
            variable = self.ShowYVar,
            command = self.Selector)
        cb.grid(row = 3, column = 1, sticky = Tk.W)

        self.ShowZVar = Tk.IntVar(self) # Create a var to track whether or not to plot in 2-D
        self.ShowZVar.set(self.parent.GetPlotParam('show_z'))
        cb = ttk.Checkbutton(frm, text = "Show z",
            variable = self.ShowZVar,
            command = self.Selector)
        cb.grid(row = 4, column = 1, sticky = Tk.W)


        # Control whether or not Cbar is shown
        self.CbarVar = Tk.IntVar()
        self.CbarVar.set(self.parent.GetPlotParam('show_cbar'))
        cb = ttk.Checkbutton(frm, text = "Show Color bar",
                        variable = self.CbarVar,
                        command = self.CbarHandler)
        cb.grid(row = 6, sticky = Tk.W)

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


        cb = ttk.Checkbutton(frm, text ='Set B or E min',
                        variable = self.setZminVar)
        cb.grid(row = 3, column = 2, sticky = Tk.W)
        self.ZminEnter = ttk.Entry(frm, textvariable=self.Zmin, width=7)
        self.ZminEnter.grid(row = 3, column = 3)

        cb = ttk.Checkbutton(frm, text ='Set B or E max',
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

        self.NormFieldVar = Tk.IntVar()
        self.NormFieldVar.set(self.parent.GetPlotParam('normalize_fields'))
        cb = ttk.Checkbutton(frm, text = "Normalize Fields",
                        variable = self.NormFieldVar,
                        command = self.NormFieldHandler)
        cb.grid(row = 6, column = 1, sticky = Tk.W)


    def CbarHandler(self, *args):
        if self.parent.GetPlotParam('show_cbar')== self.CbarVar.get():
            pass
        else:
            if self.parent.GetPlotParam('twoD'):
                self.parent.axC.set_visible(self.CbarVar.get())

            self.parent.SetPlotParam('show_cbar', self.CbarVar.get(), update_plot = self.parent.GetPlotParam('twoD'))

    def ShockVarHandler(self, *args):
        if self.parent.GetPlotParam('show_shock')== self.ShockVar.get():
            pass
        else:
            if self.parent.GetPlotParam('twoD'):
                self.parent.shockline_2d.set_visible(self.ShockVar.get())
            else:
                self.parent.shock_line.set_visible(self.ShockVar.get())

            self.parent.SetPlotParam('show_shock', self.ShockVar.get())

    def NormFieldHandler(self, *args):
        if self.parent.GetPlotParam('normalize_fields') == self.NormFieldVar.get():
            pass
        else:
            if ~self.parent.GetPlotParam('twoD'):
                if self.parent.GetPlotParam('field_type') == 0:
                    if self.NormFieldVar.get():
                        self.parent.axes.set_ylabel(r'$B/B_0$', labelpad = self.parent.parent.ylabel_pad, color = 'black')
                    else:
                        self.parent.axes.set_ylabel(r'$B$', labelpad = self.parent.parent.ylabel_pad, color = 'black')
                else:
                    if self.NormFieldVar.get():
                        self.parent.axes.set_ylabel(r'$E/E_0$', labelpad = self.parent.parent.ylabel_pad, color = 'black')
                    else:
                        self.parent.axes.set_ylabel(r'$E$', labelpad = self.parent.parent.ylabel_pad, color = 'black')

            self.parent.SetPlotParam('normalize_fields', self.NormFieldVar.get())


    def Change2d(self):

        if self.TwoDVar.get() == self.parent.GetPlotParam('twoD'):
            pass
        else:
            if self.TwoDVar.get():
                # Make sure only one dimension checked
                if self.parent.GetPlotParam('show_x'):
                    self.ShowYVar.set(0)
                    self.ShowZVar.set(0)

                elif self.parent.GetPlotParam('show_y'):
                    self.ShowZVar.set(0)

                elif ~self.parent.GetPlotParam('show_z'):
                    self.ShowXVar.set(1)

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
        if self.setZminVar.get() == self.parent.GetPlotParam('set_v_min'):
            pass
        else:
            self.parent.SetPlotParam('set_v_min', self.setZminVar.get())

    def setZmaxChanged(self, *args):
        if self.setZmaxVar.get() == self.parent.GetPlotParam('set_v_max'):
            pass
        else:
            self.parent.SetPlotParam('set_v_max', self.setZmaxVar.get())

    def RadioField(self):
        if self.FieldTypeVar.get() == self.parent.GetPlotParam('field_type'):
            pass
        else:
            if self.FieldTypeVar.get() == 0:
                if self.parent.GetPlotParam('normalize_fields'):
                    self.parent.axes.set_ylabel(r'$B/B_0$', labelpad = self.parent.parent.ylabel_pad, color = 'black')
                else:
                    self.parent.axes.set_ylabel(r'$B$', labelpad = self.parent.parent.ylabel_pad, color = 'black')
                self.parent.anx.set_text(r'$B_x$')
                self.parent.any.set_text(r'$B_y$')
                self.parent.anz.set_text(r'$B_z$')
            else:
                if self.parent.GetPlotParam('normalize_fields'):
                    self.parent.axes.set_ylabel(r'$E/E_0$', labelpad = self.parent.parent.ylabel_pad, color = 'black')
                else:
                    self.parent.axes.set_ylabel(r'$E$', labelpad = self.parent.parent.ylabel_pad, color = 'black')
                self.parent.anx.set_text(r'$E_x$')
                self.parent.any.set_text(r'$E_y$')
                self.parent.anz.set_text(r'$E_z$')

            self.parent.SetPlotParam('field_type', self.FieldTypeVar.get())

    def Selector(self):
        # First check if it is 2-D:
        if self.parent.GetPlotParam('twoD'):

            if self.ShowXVar.get() == 0 and self.ShowYVar.get() == 0 and self.ShowZVar.get() == 0:
                # All are zero, something must be selected for this plot
                self.ShowXVar.set(1)


            if self.parent.GetPlotParam('show_x') != self.ShowXVar.get():
                # set the other plot values to zero in the PlotParams
                self.parent.SetPlotParam('show_y', 0, update_plot = False)

                self.parent.SetPlotParam('show_z', 0, update_plot = False)

                # Uncheck the boxes
                self.ShowYVar.set(self.parent.GetPlotParam('show_y'))
                self.ShowZVar.set(self.parent.GetPlotParam('show_z'))

                self.parent.SetPlotParam('show_x', self.ShowXVar.get())


            elif self.parent.GetPlotParam('show_y') != self.ShowYVar.get():
                # set the other plot values to zero in the PlotParams
                self.parent.SetPlotParam('show_x', 0 ,update_plot = False)
                self.parent.SetPlotParam('show_z', 0 ,update_plot = False)

                # Uncheck the boxes
                self.ShowXVar.set(self.parent.GetPlotParam('show_x'))
                self.ShowZVar.set(self.parent.GetPlotParam('show_z'))


                self.parent.SetPlotParam('show_y', self.ShowYVar.get())

            elif self.parent.GetPlotParam('show_z') != self.ShowZVar.get():
                # set the other plot values to zero in the PlotParams
                self.parent.SetPlotParam('show_x', 0 ,update_plot = False)
                self.parent.SetPlotParam('show_y', 0 ,update_plot = False)

                # Uncheck the boxes

                self.ShowXVar.set(self.parent.GetPlotParam('show_x'))
                self.ShowYVar.set(self.parent.GetPlotParam('show_y'))

                self.parent.SetPlotParam('show_z', self.ShowZVar.get())


        else:

            if self.parent.GetPlotParam('show_x') != self.ShowXVar.get():
                self.parent.linex[0].set_visible(self.ShowXVar.get())
                self.parent.anx.set_visible(self.ShowXVar.get())
                self.parent.SetPlotParam('show_x', self.ShowXVar.get())

            elif self.parent.GetPlotParam('show_y') != self.ShowYVar.get():
                self.parent.liney[0].set_visible(self.ShowYVar.get())
                self.parent.any.set_visible(self.ShowYVar.get())
                self.parent.SetPlotParam('show_y', self.ShowYVar.get())

            elif self.parent.GetPlotParam('show_z') != self.ShowZVar.get():
                self.parent.linez[0].set_visible(self.ShowZVar.get())
                self.parent.anz.set_visible(self.ShowZVar.get())
                self.parent.SetPlotParam('show_z', self.ShowZVar.get())

    def TxtEnter(self, e):
        self.FieldsCallback()

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


    def OnClosing(self):
        self.parent.settings_window = None
        self.destroy()
