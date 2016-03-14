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
import time

class FieldsPanel:
    # A dictionary of all of the parameters for this plot with the default parameters

    plot_param_dict = {'twoD': 0,
                       'field_type': 0, #0 = B-Field, 1 = E-field
                       'show_x' : 1,
                       'show_y' : 1,
                       'show_z' : 1,
                       'show_cbar': True,
                       'z_min': 0,
                       'z_max' : 10,
                       'set_z_min': False,
                       'set_z_max': False,
                       'show_shock' : False,
                       'OutlineText': True,
                       'spatial_x': True,
                       'spatial_y': None,
                       'interpolation': 'hermite'}

    b_key_names = ['bx', 'by', 'bz']
    e_key_names = ['ex', 'ey', 'ez']

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
        self.arrs_needed = ['c_omp', 'istep', 'bx', 'by', 'bz', 'ex', 'ey', 'ez']
        return self.arrs_needed

    def LoadData(self):
        ''' A Helper function that loads the data for the plot'''
        # First see of the x_axis and y_axis values have already been calculated
        # and stored in the DataDict for this time step
        self.c_omp = self.FigWrap.LoadKey('c_omp')[0]
        self.istep = self.FigWrap.LoadKey('istep')[0]
        self.xcolor = new_cmaps.cmaps[self.parent.cmap](0.2)
        self.ycolor = new_cmaps.cmaps[self.parent.cmap](0.5)
        self.zcolor = new_cmaps.cmaps[self.parent.cmap](0.8)

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
            self.fx = self.FigWrap.LoadKey('bx')[0,:,:]
            self.fy = self.FigWrap.LoadKey('by')[0,:,:]
            self.fz = self.FigWrap.LoadKey('bz')[0,:,:]


            # Have we already calculated min/max?
            if 'bxmin_max' in self.parent.DataDict.keys():
                self.bxmin_max = self.parent.DataDict['bxmin_max']

            else:
                self.bxmin_max =  [self.fx.min(), self.fx.max()]
                self.parent.DataDict['bxmin_max'] = list(self.bxmin_max)

            if 'bymin_max' in self.parent.DataDict.keys():
                self.bymin_max = self.parent.DataDict['bymin_max']

            else:
                self.bymin_max =  [self.fy.min(), self.fy.max()]
                self.parent.DataDict['bymin_max'] = list(self.bymin_max)

            if 'bzmin_max' in self.parent.DataDict.keys():
                self.bzmin_max = self.parent.DataDict['bzmin_max']

            else:
                self.bzmin_max =  [self.fz.min(), self.fz.max()]
                self.parent.DataDict['bzmin_max'] = list(self.bzmin_max)

            self.annotate_text_list = [r'$B_x$', r'$B_y$',r'$B_z$']
        if self.GetPlotParam('field_type') == 1: # Load the e-Field
            self.fx = self.FigWrap.LoadKey('ex')[0,:,:]
            self.fy = self.FigWrap.LoadKey('ey')[0,:,:]
            self.fz = self.FigWrap.LoadKey('ez')[0,:,:]

            # Have we already calculated min/max?


            if 'exmin_max' in self.parent.DataDict.keys():
                self.exmin_max = self.parent.DataDict['exmin_max']

            else:
                self.exmin_max =  [self.fx.min(), self.fx.max()]
                self.parent.DataDict['exmin_max'] = list(self.exmin_max)

            if 'eymin_max' in self.parent.DataDict.keys():
                self.eymin_max = self.parent.DataDict['eymin_max']

            else:
                self.eymin_max =  [self.fy.min(), self.fy.max()]
                self.parent.DataDict['eymin_max'] = list(self.eymin_max)

            if 'ezmin_max' in self.parent.DataDict.keys():
                self.ezmin_max = self.parent.DataDict['ezmin_max']

            else:
                self.ezmin_max =  [self.fz.min(), self.fz.max()]
                self.parent.DataDict['ezmin_max'] = list(self.ezmin_max)
            self.annotate_text_list = [r'$E_x$',r'$E_y$',r'$E_z$']
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

        #load the data
        self.LoadData()
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
                self.two_d_labels = (r'$B_x$', r'$E_x$')

                # set the other plot values to zero in the PlotParams
                self.SetPlotParam('show_y', 0, update_plot = False)
                self.SetPlotParam('show_z', 0, update_plot = False)

            elif self.GetPlotParam('show_y'):
                self.zval = self.fy
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
                self.two_d_labels = (r'$B_z$', r'$E_z$')


            self.ymin = 0
            self.ymax =  self.zval.shape[0]/self.c_omp*self.istep
            self.xmin = 0
            self.xmax =  self.zval.shape[1]/self.c_omp*self.istep

            self.vmin = None
            if self.GetPlotParam('set_z_min'):
                self.vmin = self.GetPlotParam('z_min')
            self.vmax = None
            if self.GetPlotParam('set_z_max'):
                self.vmax = self.GetPlotParam('z_max')
            if self.parent.plot_aspect:
                self.cax = self.axes.imshow(self.zval,
                    cmap = new_cmaps.cmaps[self.parent.cmap],
                    origin = 'lower',
                    extent = (self.xmin,self.xmax, self.ymin, self.ymax),
                    vmin = self.vmin,
                    vmax = self.vmax,
                    interpolation=self.GetPlotParam('interpolation'))
            else:
                self.cax = self.axes.imshow(self.zval,
                    cmap = new_cmaps.cmaps[self.parent.cmap],
                    origin = 'lower',
                    aspect = 'auto',
                    extent = (self.xmin,self.xmax, self.ymin, self.ymax),
                    vmin = self.vmin,
                    vmax = self.vmax,
                    interpolation=self.GetPlotParam('interpolation'))



            self.TwoDan = self.axes.annotate(self.two_d_labels[self.GetPlotParam('field_type')],
                            xy = (0.9,.9),
                            xycoords= 'axes fraction',
                            color = 'white',
                            **self.annotate_kwargs)

            self.axes.set_axis_bgcolor('lightgrey')

            self.axC = self.figure.add_subplot(self.gs[:4,:])
            self.cbar = self.figure.colorbar(self.cax, ax = self.axes, cax = self.axC, orientation = 'horizontal')

            cmin = self.zval.min()
            if self.vmin:
                cmin = self.vmin
            cmax = self.zval.max()
            if self.vmax:
                cmax = self.vmax
            self.cbar.set_ticks(np.linspace(cmin, cmax, 5))
            self.cbar.ax.tick_params(labelsize=self.parent.num_font_size)

            if self.GetPlotParam('show_cbar') == 0:
                self.cbar.set_visible = False



            self.shockline_2d = self.axes.axvline(self.parent.shock_loc, linewidth = 1.5, linestyle = '--', color = self.parent.shock_color, path_effects=[PathEffects.Stroke(linewidth=2, foreground='k'),
                                    PathEffects.Normal()])
            self.shockline_2d.set_visible(self.GetPlotParam('show_shock'))

            self.axes.set_axis_bgcolor('lightgrey')
            self.axes.tick_params(labelsize = self.parent.num_font_size, color=tick_color)
#        self.axes.set_xlim(self.xmin,self.xmax)
            if self.parent.xlim[0]:
                self.axes.set_xlim(self.parent.xlim[1],self.parent.xlim[2])
            if self.parent.ylim[0]:
                self.axes.set_ylim(self.parent.ylim[1],self.parent.ylim[2])
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
            self.linex = self.axes.plot(self.xaxis_values, self.fx[self.fx.shape[0]/2,:], color = self.xcolor)
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
            self.liney = self.axes.plot(self.xaxis_values, self.fy[self.fy.shape[0]/2,:], color = self.ycolor)
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


            self.linez = self.axes.plot(self.xaxis_values, self.fz[self.fz.shape[0]/2,:], color = self.zcolor)
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

            if self.parent.xlim[0]:
                self.axes.set_xlim(self.parent.xlim[1],self.parent.xlim[2])
            else:
                self.axes.set_xlim(self.xaxis_values[0],self.xaxis_values[-1])


            self.axes.relim(visible_only=True)
            self.axes.autoscale_view(scaley = True)
            if self.GetPlotParam('set_z_min'):
                self.axes.set_ylim(ymin = self.GetPlotParam('z_min'))
            if self.GetPlotParam('set_z_max'):
                self.axes.set_ylim(ymax = self.GetPlotParam('z_max'))

            self.axes.set_xlabel(r'$x\ [c/\omega_{\rm pe}]$', labelpad = self.parent.xlabel_pad, color = 'black')

    def refresh(self):

        tic = time.time()
        '''This is a function that will be called only if self.axes already
        holds a fields type plot. We only update things that have changed.  If
        hasn't changed, don't touch it. The difference between this and last
        time, is that we won't actually do any drawing in the plot. The plot
        will be redrawn after all subplots are refreshed. '''
        self.LoadData()

        # Main goal, only change what is showing..
        # First do the 1D plots, because it is simpler
        if self.GetPlotParam('twoD') == 0:
            if self.GetPlotParam('show_x'):
                self.linex[0].set_data(self.xaxis_values, self.fx[self.fx.shape[0]/2,:])
                self.linex[0].set_color(self.xcolor)
                self.anx.set_color(self.xcolor)


            if self.GetPlotParam('show_y'):
                self.liney[0].set_data(self.xaxis_values, self.fy[self.fy.shape[0]/2,:])
                self.liney[0].set_color(self.ycolor)
                self.any.set_color(self.ycolor)


            if self.GetPlotParam('show_z'):
                self.linez[0].set_data(self.xaxis_values, self.fz[self.fz.shape[0]/2,:])
                self.linez[0].set_color(self.zcolor)
                self.anz.set_color(self.zcolor)


            if self.GetPlotParam('show_shock'):
                self.shock_line.set_xdata([self.parent.shock_loc,self.parent.shock_loc])

            self.axes.relim(visible_only=True)
            self.axes.autoscale_view(scaley = True)
            if self.parent.xlim[0]:
                self.axes.set_xlim(self.parent.xlim[1],self.parent.xlim[2])
            else:
                self.axes.set_xlim(self.xaxis_values[0], self.xaxis_values[-1])

            if self.GetPlotParam('set_z_min'):
                self.axes.set_ylim(ymin = self.GetPlotParam('z_min'))
            if self.GetPlotParam('set_z_max'):
                self.axes.set_ylim(ymax = self.GetPlotParam('z_max'))
            self.axes.draw_artist(self.axes.patch)
            self.axes.draw_artist(self.linex[0])
            self.axes.draw_artist(self.liney[0])
            self.axes.draw_artist(self.linez[0])
            self.axes.draw_artist(self.anx)
            self.axes.draw_artist(self.any)
            self.axes.draw_artist(self.anz)

        else: # Now refresh the plot if it is 2D
            if self.GetPlotParam('show_x'):
                self.cax.set_data(self.fx)
                self.ymin = 0
                self.ymax =  self.fx.shape[0]/self.c_omp*self.istep
                self.xmin = 0
                self.xmax =  self.fx.shape[1]/self.c_omp*self.istep
                if self.GetPlotParam('field_type') == 0:
                    self.clims = self.bxmin_max
                    self.TwoDan.set_text(r'$B_x$')
                else:
                    self.clims = self.exmin_max
                    self.TwoDan.set_text(r'$E_x$')

            if self.GetPlotParam('show_y'):
                self.cax.set_data(self.fy)
                self.ymin = 0
                self.ymax =  self.fy.shape[0]/self.c_omp*self.istep
                self.xmin = 0
                self.xmax =  self.fy.shape[1]/self.c_omp*self.istep
                if self.GetPlotParam('field_type') == 0:
                    self.clims = self.bymin_max
                    self.TwoDan.set_text(r'$B_y$')
                else:
                    self.clims = self.eymin_max
                    self.TwoDan.set_text(r'$E_y$')

            if self.GetPlotParam('show_z'):
                self.cax.set_data(self.fz)
                self.ymin = 0
                self.ymax =  self.fz.shape[0]/self.c_omp*self.istep
                self.xmin = 0
                self.xmax =  self.fz.shape[1]/self.c_omp*self.istep
                if self.GetPlotParam('field_type') == 0:
                    self.clims = self.bzmin_max
                    self.TwoDan.set_text(r'$B_z$')
                else:
                    self.clims = self.ezmin_max
                    self.TwoDan.set_text(r'$E_z$')

            self.cax.set_extent([self.xmin,self.xmax, self.ymin, self.ymax])


            self.cax.set_clim(self.clims)

            self.climArgs = {}
            if self.GetPlotParam('set_z_min'):
                self.climArgs['vmin'] =  self.GetPlotParam('z_min')
            if self.GetPlotParam('set_z_max'):
                self.climArgs['vmax'] =  self.GetPlotParam('z_max')
            if len(self.climArgs)>0:
                self.cax.set_clim(**self.climArgs)

            self.cbar.set_ticks(np.linspace(self.cax.get_clim()[0],self.cax.get_clim()[1], 5))

            if self.GetPlotParam('show_shock'):
                self.shockline_2d.set_xdata([self.parent.shock_loc,self.parent.shock_loc])
            self.axes.draw_artist(self.axes.patch)
            self.axes.draw_artist(self.cax)

        toc = time.time()
        print toc-tic


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
        cmapChooser = apply(ttk.OptionMenu, (frm, self.ctypevar, self.parent.chartType) + tuple(self.parent.ChartTypes))
        cmapChooser.grid(row =0, column = 1, sticky = Tk.W + Tk.E)


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
                        command = lambda:
                        self.parent.SetPlotParam('show_cbar', self.CbarVar.get()))
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

    def RadioField(self):
        if self.FieldTypeVar.get() == self.parent.GetPlotParam('field_type'):
            pass
        else:
            if self.FieldTypeVar.get() == 0:
                self.parent.anx.set_text(r'$B_x$')
                self.parent.any.set_text(r'$B_y$')
                self.parent.anz.set_text(r'$B_z$')
            else:
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
