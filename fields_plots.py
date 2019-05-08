#!/usr/bin/env python
import tkinter as Tk
from tkinter import ttk, messagebox
import matplotlib
import numpy as np
import numpy.ma as ma
import new_cmaps
import sys, traceback
from new_cnorms import PowerNormWithNeg, PowerNormFunc
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from modest_image import imshow
import matplotlib.patheffects as PathEffects
import matplotlib.transforms as mtransforms

class FieldsPanel:
    # A dictionary of all of the parameters for this plot with the default parameters
    example = """## WHATEVER IS TYPED HERE IS EVALUATED AS PURE PYTHON. THERE IS NO ERROR CHECKING
## OR ANY SANITIZATION OF USER INPUT. YOU WILL INHERIT THE NAMESPACE OF THE MAIN
## PROGRAM, BUT YOU CAN IMPORT OTHER LIBRARIES, DEFINE HELPER FUNCTIONS, WHATEVER.
## JUST BE SURE AT SOME POINT YOU DEFINE A FUNCTION NAMED FieldFunc THAT RETURNS
## SOMETHING THE SAME SHAPE AS YOUR FIELD ARRAYS. SIMULATION DATA CAN ONLY BE
## ACCESSED INSIDE OF FieldFunc.
#
## IT'S EASY TO DO BAD THINGS HERE... TYPE CAREFULLY :)
#
#def FieldFunc(bx, by, bz):
#    # Be sure to include all the neccesary data you need to calculate your
#    # derived field quantity as arguments to the 'FieldFunc' function.
#    # The only valid arguments to field function are things saved in the Tristan
#    # HDF5 files: e.g., ui, bx, jz...etc. The argumes return the raw tristan arrays.
#
#    # You must return an array the same shape as the fields array, or an array that
#    # is the same length as the x axis of the simulation (and then checking 1D)
#
#    return bx**2+by**2+bz**2
#    """
    plot_param_dict = {'twoD': 0,
                       'field_type': 0, #0 = B-Field, 1 = E-field, 2 Currents, 3 = UserDefined quantity
                       'cmdstr1': example,
                       'cmdstr2': example,
                       'cmdstr3': example,
                       'OneDOnly': [False, False, False],
                       'yaxis_label': ['$B$','$E$','$J$','$B$'],
                       '2D_label': [['$B_x$','$B_y$','$B_z$'],
                                    ['$E_x$','$E_y$','$E_z$'],
                                    ['$J_x$','$J_y$','$J_z$'],
                                    ['$B_\mathrm{tot}$','$B_\mathrm{tot}$','$B_\mathrm{tot}$']],
                       '1D_label': [['$B_x$','$B_y$','$B_z$'],
                                    ['$E_x$','$E_y$','$E_z$'],
                                    ['$J_x$','$J_y$','$J_z$'],
                                    ['$B_\mathrm{tot}$','$B_\mathrm{tot}$','$B_\mathrm{tot}$']],
                       'show_x' : 1,
                       'show_y' : 1,
                       'show_z' : 1,
                       'show_cbar': True,
                       'v_min': 0,
                       'v_max' : 10,
                       'set_v_min': False,
                       'set_v_max': False,
                       'show_shock' : False,
                       'show_FFT_region': False,
                       'OutlineText': True,
                       'spatial_x': True,
                       'spatial_y': False,
                       'normalize_fields': True, # Normalize fields to their upstream values
                       'cnorm_type': 'Linear', # Colormap norm;  options are Log, Pow or Linear
                       'cpow_num': 1.0, # Used in the PowerNorm
                       'div_midpoint': 0.0, # The cpow color norm normalizes data to [0,1] using np.sign(x-midpoint)*np.abs(x-midpoint)**(-cpow_num) -> [0,midpoint,1] if it is a divering cmap or [0,1] if it is not a divering cmap
                       'interpolation': 'none',
                       'cmap': 'None', # If cmap is none, the plot will inherit the parent's cmap
                       'UseDivCmap': True, # Use a diverging cmap for the 2d plots
                       'stretch_colors': False, # If stretch colors is false, then for a diverging cmap the plot ensures -b and b are the same distance from the midpoint of the cmap.
                       'show_cpu_domains': False # plots lines showing how the CPUs are divvying up the computational region
                       }
    BoolList = ['twoD', 'show_cbar', 'set_v_min', 'set_v_max',
             'show_shock', 'OutlineText', 'spatial_x', 'spatial_y',
             'Show_FFT_region', 'UseDivCmap', 'stretch_colors', 'normalize_fields',
             'show_x', 'show_y', 'show_z', 'show_cpu_domains']
    IntList = ['field_type']
    FloatList = ['v_min', 'v_max', 'cpow_num', 'div_midpoint']
    #StrList = ['interpolation', 'cnorm_type', 'cmap']
    StrList = ['cnorm_type', 'cmap', 'interpolation', 'cmdstr1', 'cmdstr2', 'cmdstr3']
    SpecialList = ['yaxis_label', '2D_label', '2D_label', 'OneDOnly']
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


    def ChangePlotType(self, str_arg):
        self.FigWrap.ChangeGraph(str_arg)

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

    def set_plot_keys(self):
        '''A helper function that will insure that each hdf5 file will only be
        opened once per time step'''
        # First make sure that omega_plasma & xi is loaded so we can fix the
        # x & y distances.

        # Then see if we are plotting E-field or B-Field
        if self.GetPlotParam('field_type') == 0: # Load the B-Field
            self.arrs_needed = ['c_omp', 'istep', 'bx']#, 'by', 'bz']
            if self.GetPlotParam('show_y'):
                self.arrs_needed.append('by')
            if self.GetPlotParam('show_z'):
                self.arrs_needed.append('bz')

        if self.GetPlotParam('field_type') == 1: # Load the E-Field
            self.arrs_needed = ['c_omp', 'istep', 'ex']
            if self.GetPlotParam('show_y'):
                self.arrs_needed.append('ey')
            if self.GetPlotParam('show_z'):
                self.arrs_needed.append('ez')

        if self.GetPlotParam('field_type') == 2: # Load the currents
            self.arrs_needed = ['c_omp', 'istep', 'jx']
            if self.GetPlotParam('show_y'):
                self.arrs_needed.append('jy')
            if self.GetPlotParam('show_z'):
                self.arrs_needed.append('jz')

        if self.GetPlotParam('field_type') == 3: # Check what the user wants.
            self.arrs_needed = ['c_omp', 'istep', 'bx']
            if self.GetPlotParam('show_x'):
                for line in self.GetPlotParam('cmdstr1').splitlines():
                    if line[1:15] == 'def FieldFunc(':
                        self.f1args = [elm.strip() for elm in line[15:-2].split(',')]
                        self.arrs_needed += self.f1args
            if self.GetPlotParam('show_y'):
                for line in self.GetPlotParam('cmdstr2').splitlines():
                    if line[1:15] == 'def FieldFunc(':
                        self.f2args = [elm.strip() for elm in line[15:-2].split(',')]
                        self.arrs_needed += self.f2args
            if self.GetPlotParam('show_z'):
                for line in self.GetPlotParam('cmdstr3').splitlines():
                    if line[1:15] == 'def FieldFunc(':
                        self.f3args = [elm.strip() for elm in line[15:-2].split(',')]
                        self.arrs_needed += self.f3args
        return self.arrs_needed

    def LoadData(self):
        ''' A Helper function that loads the data for the plot'''
        # First see of the x_axis and y_axis values have already been calculated
        # and stored in the DataDict for this time step
        self.c_omp = self.FigWrap.LoadKey('c_omp')[0]
        self.istep = self.FigWrap.LoadKey('istep')[0]
        if self.GetPlotParam('cmap') == 'None':
            if self.GetPlotParam('UseDivCmap'):
                self.cmap = self.parent.MainParamDict['DivColorMap']
            else:
                self.cmap = self.parent.MainParamDict['ColorMap']
        else:
            self.cmap = self.GetPlotParam('cmap')
        self.xcolor = new_cmaps.cmaps[self.parent.MainParamDict['ColorMap']](0.2)
        self.ycolor = new_cmaps.cmaps[self.parent.MainParamDict['ColorMap']](0.5)
        self.zcolor = new_cmaps.cmaps[self.parent.MainParamDict['ColorMap']](0.8)

        if np.isnan(self.parent.btheta):
            # Maybe B_0 is 0????
            self.SetPlotParam('normalize_fields', 0, update_plot = False)

        # see if the axis values are saved in the data dict
        if 'xaxis_values' in self.parent.DataDict.keys():
            self.xaxis_values = self.parent.DataDict['xaxis_values']
        else:
            # x-values haven't been calculated yet, generate them then save them to the dictionary for later.
            if self.GetPlotParam('field_type') ==0 or self.GetPlotParam('field_type') == 3:
                self.xaxis_values = np.arange(self.FigWrap.LoadKey('bx').shape[2])/self.c_omp*self.istep
            elif self.GetPlotParam('field_type') ==1:
                self.xaxis_values = np.arange(self.FigWrap.LoadKey('ex').shape[2])/self.c_omp*self.istep
            elif self.GetPlotParam('field_type') ==2:
                self.xaxis_values = np.arange(self.FigWrap.LoadKey('jx').shape[2])/self.c_omp*self.istep
            self.parent.DataDict['xaxis_values'] = np.copy(self.xaxis_values)

        self.flagx = 0 # 0 means it didn't plot, 1 means it is 1D only, 2 means it returned a 3d object
        self.flagy = 0
        self.flagz = 0

        if self.GetPlotParam('field_type') == 0: # Load the B-Field
            if self.GetPlotParam('show_x'):
                self.flagx = 2
                if self.GetPlotParam('normalize_fields'):
                    self.fx = self.FigWrap.LoadKey('bx')*self.parent.b0**-1
                else:
                    self.fx = self.FigWrap.LoadKey('bx')
            if self.GetPlotParam('show_y'):
                self.flagy = 2
                if self.GetPlotParam('normalize_fields'):
                    self.fy = self.FigWrap.LoadKey('by')*self.parent.b0**-1
                else:
                    self.fy = self.FigWrap.LoadKey('by')

            if self.GetPlotParam('show_z'):
                self.flagz = 2
                if self.GetPlotParam('normalize_fields'):
                    self.fz =self.FigWrap.LoadKey('bz')*self.parent.b0**-1
                else:
                    self.fz =self.FigWrap.LoadKey('bz')

        if self.GetPlotParam('field_type') == 1: # Load the E-Field
            if self.GetPlotParam('show_x'):
                self.flagx = 2
                if self.GetPlotParam('normalize_fields'):
                    self.fx = self.FigWrap.LoadKey('ex')*self.parent.e0**-1
                else:
                    self.fx = self.FigWrap.LoadKey('ex')
            if self.GetPlotParam('show_y'):
                self.flagy = 2
                if self.GetPlotParam('normalize_fields'):
                    self.fy = self.FigWrap.LoadKey('ey')*self.parent.e0**-1
                else:
                    self.fy = self.FigWrap.LoadKey('ey')

            if self.GetPlotParam('show_z'):
                self.flagz = 2
                if self.GetPlotParam('normalize_fields'):
                    self.fz =self.FigWrap.LoadKey('ez')*self.parent.e0**-1
                else:
                    self.fz =self.FigWrap.LoadKey('ez')


        elif self.GetPlotParam('field_type') == 2: # Load the currents

            if self.GetPlotParam('show_x'):
                self.fx = self.FigWrap.LoadKey('jx')
                self.flagx = 2

            if self.GetPlotParam('show_y'):
                self.fy = self.FigWrap.LoadKey('jy')
                self.flagy = 2

            if self.GetPlotParam('show_z'):
                self.fz = self.FigWrap.LoadKey('jz')
                self.flagz = 2

        elif self.GetPlotParam('field_type') == 3: # User Defined fields
            if self.GetPlotParam('show_x'):
                if not set(self.f1args).isdisjoint(self.parent.prtl_keys):
                    keyx = hash(self.GetPlotParam('cmdstr1')+str(self.parent.stride)+str(self.GetPlotParam('OneDOnly')[0]))
                else:
                    keyx = hash(self.GetPlotParam('cmdstr1')+str(self.GetPlotParam('OneDOnly')[0]) )
                if keyx in self.parent.DataDict.keys():
                    self.fx = self.parent.DataDict[keyx]
                    if self.GetPlotParam('OneDOnly')[0]:
                        self.flagx = 1
                    else:
                        self.flagx = 2
                else:
                    try:
                        tmpcstr = ''
                        for line in self.GetPlotParam('cmdstr1').splitlines():
                            tmpcstr += line[1:] +'\n'
                        tmpcstr += 'self.fx = FieldFunc(*[self.FigWrap.LoadKey(k) for k in self.f1args])'
                        eval(compile(tmpcstr, '<string>', 'exec'))
                        self.parent.DataDict[keyx] = self.fx
                        if self.GetPlotParam('OneDOnly')[0]:
                            self.flagx = 1
                        else:
                            self.flagx = 2
                    except:
                        tb_lines = traceback.format_exc(sys.exc_info()[2]).splitlines()
                        tb_lines.pop(1)
                        tb_lines[1] = ''
                        err_msg = ''
                        for l in tb_lines:
                            if l[0:17] == '  File "<string>"':
                                err_msg += '  User Defined Function,'
                                err_msg += l[18:] +'\n'
                            else:
                                err_msg += l+'\n'
                        messagebox.showinfo('Error when evaluating user defined function 1:', err_msg)

                        self.fx = np.NAN
                        self.flagx = 0

            if self.GetPlotParam('show_y'):
                if not set(self.f2args).isdisjoint(self.parent.prtl_keys):
                    keyy = hash(self.GetPlotParam('cmdstr2')+str(self.parent.stride)+str(self.GetPlotParam('OneDOnly')[1]))
                else:
                    keyy = hash(self.GetPlotParam('cmdstr2')+str(self.GetPlotParam('OneDOnly')[1]))
                if keyy in self.parent.DataDict.keys():
                    self.fy = self.parent.DataDict[keyy]
                    if self.GetPlotParam('OneDOnly')[0]:
                        self.flagy = 1
                    else:
                        self.flagy = 2
                else:
                    try:
                        tmpcstr = ''
                        for line in self.GetPlotParam('cmdstr2').splitlines():
                            tmpcstr += line[1:] +'\n'
                        tmpcstr += 'self.fy = FieldFunc(*[self.FigWrap.LoadKey(k) for k in self.f2args])'
                        eval(compile(tmpcstr, '<string>', 'exec'))
                        self.parent.DataDict[keyy] = self.fy
                        if self.GetPlotParam('OneDOnly')[1]:
                            self.flagy = 1
                        else:
                            self.flagy = 2
                    except:
                        tb_lines = traceback.format_exc(sys.exc_info()[2]).splitlines()
                        tb_lines.pop(1)
                        tb_lines[1] = ''
                        err_msg = ''
                        for l in tb_lines:
                            if l[0:17] == '  File "<string>"':
                                err_msg += '  User Defined Function,'
                                err_msg += l[18:] +'\n'
                            else:
                                err_msg += l+'\n'
                        messagebox.showinfo('Error when evaluating user defined function 2:', err_msg)
                        self.fy = np.NAN
                        self.flagy = 0

            if self.GetPlotParam('show_z'):
                if not set(self.f3args).isdisjoint(self.parent.prtl_keys):
                    keyz = hash(self.GetPlotParam('cmdstr3')+str(self.parent.stride)+str(self.GetPlotParam('OneDOnly')[2]))
                else:
                    keyz = hash(self.GetPlotParam('cmdstr3')+str(self.GetPlotParam('OneDOnly')[2]))
                if keyz in self.parent.DataDict.keys():
                    self.fz = self.parent.DataDict[keyz]
                    if self.GetPlotParam('OneDOnly')[2]:
                        self.flagz = 1
                    else:
                        self.flagz = 2
                else:
                    try:
                        tmpcstr = ''
                        for line in self.GetPlotParam('cmdstr3').splitlines():
                            tmpcstr += line[1:] +'\n'
                        tmpcstr += 'self.fz = FieldFunc(*[self.FigWrap.LoadKey(k) for k in self.f3args])'
                        eval(compile(tmpcstr, '<string>', 'exec'))
                        self.parent.DataDict[keyz] = self.fz
                        if self.GetPlotParam('OneDOnly')[2]:
                            self.flagz = 1
                        else:
                            self.flagz = 2

                    except:
                        tb_lines = traceback.format_exc(sys.exc_info()[2]).splitlines()
                        tb_lines.pop(1)
                        tb_lines[1] = ''
                        err_msg = ''
                        for l in tb_lines:
                            if l[0:17] == '  File "<string>"':
                                err_msg += '  User Defined Function,'
                                err_msg += l[18:] +'\n'
                            else:
                                err_msg += l+'\n'
                        messagebox.showinfo('Error when evaluating user defined function 3:', err_msg)

                        self.fz = np.NAN
                        self.flagz = 0

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

            # First choose the 'zval' to plot, we can only do one because it is 2-d.
            self.plotFlag = -1
            if self.GetPlotParam('show_x') and self.flagx == 2:
                if self.parent.MainParamDict['2DSlicePlane'] == 0: # Show the x-y plane
                    if self.parent.MainParamDict['ImageAspect']:
                        self.cax = imshow(self.axes, self.fx[self.parent.zSlice,:,:], norm = self.norm(), origin = 'lower')
                    else:
                        self.cax = imshow(self.axes,self.fx[self.parent.zSlice,:,:], origin = 'lower', norm = self.norm(),
                                                    aspect= 'auto')

                elif self.parent.MainParamDict['2DSlicePlane'] == 1: # Show the x-z plane
                    if self.parent.MainParamDict['ImageAspect']:
                        self.cax = imshow(self.axes, self.fx[:,self.parent.ySlice,:], norm = self.norm(), origin = 'lower')
                    else:
                        self.cax = imshow(self.axes, self.fx[:,self.parent.ySlice,:], origin = 'lower', norm = self.norm(),
                                                    aspect= 'auto')

                self.plotFlag = 0
                self.SetPlotParam('show_y', 0, update_plot = False)
                self.SetPlotParam('show_z', 0, update_plot = False)

            elif self.GetPlotParam('show_y') and self.flagy == 2:
                if self.parent.MainParamDict['2DSlicePlane'] == 0: # Show the x-y plane
                    if self.parent.MainParamDict['ImageAspect']:
                        self.cax = imshow(self.axes, self.fy[self.parent.zSlice,:,:], norm = self.norm(), origin = 'lower')
                    else:
                        self.cax = imshow(self.axes, self.fy[self.parent.zSlice,:,:], origin = 'lower', norm = self.norm(),
                                                    aspect= 'auto')
                elif self.parent.MainParamDict['2DSlicePlane'] == 1: # Show the x-z plane
                    if self.parent.MainParamDict['ImageAspect']:
                        self.cax = imshow(self.axes, self.fy[:,self.parent.ySlice,:], norm = self.norm(), origin = 'lower')
                    else:
                        self.cax = imshow(self.axes,self.fy[:,self.parent.ySlice,:], origin = 'lower', norm = self.norm(),
                                                    aspect= 'auto')
                self.plotFlag = 1
                self.SetPlotParam('show_x', 0, update_plot = False)
                self.SetPlotParam('show_z', 0, update_plot = False)


            elif self.GetPlotParam('show_z') and self.flagz == 2:
                # make sure z is loaded, (something has to be)
                # set the other plot values to zero in the PlotParams
                if self.parent.MainParamDict['2DSlicePlane'] == 0: # Show the x-y plane
                    if self.parent.MainParamDict['ImageAspect']:
                        self.cax = imshow(self.axes,self.fz[self.parent.zSlice,:,:], norm = self.norm(), origin = 'lower')
                    else:
                        self.cax = imshow(self.axes, self.fz[self.parent.zSlice,:,:], origin = 'lower', norm = self.norm(),
                                                    aspect= 'auto')
                elif self.parent.MainParamDict['2DSlicePlane'] == 1: # Show the x-z plane
                    if self.parent.MainParamDict['ImageAspect']:
                        self.cax = imshow(self.axes,self.fz[:,self.parent.ySlice,:], norm = self.norm(), origin = 'lower')
                    else:
                        self.cax = imshow(self.axes,self.fz[:,self.parent.ySlice,:], origin = 'lower', norm = self.norm(),
                                                    aspect= 'auto')

                self.plotFlag = 2
                self.SetPlotParam('show_x', 0, update_plot = False)
                self.SetPlotParam('show_y', 0, update_plot = False)

            else:
                if self.parent.MainParamDict['ImageAspect']:
                    self.cax = self.axes.imshow(np.ones([2,2]), norm = self.norm(), origin = 'lower')
                else:
                    self.cax = self.axes.imshow(np.ones([2,2]), norm = self.norm(), origin = 'lower', aspect = 'auto')
                self.cax.set_data(np.ma.masked_array(np.empty([2,2]), mask = np.ones([2,2])))
            self.ymin = 0
            self.ymax =  self.cax.get_array().shape[0]/self.c_omp*self.istep
            self.xmin = 0
            self.xmax =  self.cax.get_array().shape[1]/self.c_omp*self.istep
            self.cax.set_cmap(new_cmaps.cmaps[self.cmap])

            if self.plotFlag>=0:

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
                self.cax.set_extent([self.xmin,self.xmax, self.ymin, self.ymax])



            self.axes.add_artist(self.cax)
            self.anntext =''
            if self.plotFlag >= 0:
                self.anntext = self.GetPlotParam('2D_label')[self.GetPlotParam('field_type')][self.plotFlag]
                if self.GetPlotParam('field_type') ==0  and self.GetPlotParam('normalize_fields'):
                    self.anntext +=r'$/B_0$'
                if self.GetPlotParam('field_type') ==1  and self.GetPlotParam('normalize_fields'):
                    self.anntext +=r'$/E_0$'

            self.TwoDan = self.axes.annotate(self.anntext,
                        xy = (0.9,.9),
                        xycoords= 'axes fraction',
                        color = 'white',
                        **self.annotate_kwargs)

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
            else:
                self.cbar = self.axC.imshow(np.transpose(self.gradient)[::-1], aspect='auto',
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
                                right=True,         # ticks along the top edge are off
                                labelleft = False,
                                labelright  = True,
                                labelsize=self.parent.MainParamDict['NumFontSize'])

            if self.GetPlotParam('show_cbar') == 0 or self.plotFlag == -1:
                self.axC.set_visible(False)
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
            self.cax.set_interpolation(self.GetPlotParam('interpolation'))
            if self.parent.MainParamDict['SetyLim']:
                self.axes.set_ylim(self.parent.MainParamDict['yBottom'],self.parent.MainParamDict['yTop'])
            else:
                self.axes.set_ylim(self.ymin, self.ymax)
            self.axes.set_xlabel(r'$x\ [c/\omega_{\rm pe}]$', labelpad = self.parent.MainParamDict['xLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
            if self.parent.MainParamDict['2DSlicePlane'] == 0:
                self.axes.set_ylabel(r'$y\ [c/\omega_{\rm pe}]$', labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
            if self.parent.MainParamDict['2DSlicePlane'] == 1:
                self.axes.set_ylabel(r'$z\ [c/\omega_{\rm pe}]$', labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])

        else: # It's 1D
            if self.parent.MainParamDict['LinkSpatial'] != 0 and self.parent.MainParamDict['LinkSpatial'] != 3:
                if self.FigWrap.pos == self.parent.first_x:
                    self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]])
                else:
                    self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]],
                    sharex = self.parent.SubPlotList[self.parent.first_x[0]][self.parent.first_x[1]].graph.axes)
            else:
                self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]])

            self.annotate_pos = [0.8,0.9]
            self.xmin, self.xmax = self.xaxis_values[0], self.xaxis_values[-1]

            min_max = [np.inf, -np.inf]
            if self.flagx > 0 and self.GetPlotParam('show_x'):
                if self.flagx == 1 and len(self.fx.shape) == 1:
                    self.linex = self.axes.plot(self.xaxis_values, self.fx, color = self.xcolor)
                elif self.parent.MainParamDict['Average1D']:
                    self.linex = self.axes.plot(self.xaxis_values, np.average(self.fx.reshape(-1,self.fx.shape[-1]), axis =0), color = self.xcolor)
                else:
                    self.linex = self.axes.plot(self.xaxis_values, self.fx[self.parent.zSlice,self.parent.ySlice,:], color = self.xcolor)
                min_max[0]=min(min_max[0],self.linex[0].get_data()[1].min())
                min_max[1]=max(min_max[1],self.linex[0].get_data()[1].max())
            else:
                self.linex = self.axes.plot(np.arange(10), np.arange(10), color = self.xcolor)
                self.linex[0].set_visible(False)

            self.anx = self.axes.annotate(self.GetPlotParam('1D_label')[self.GetPlotParam('field_type')][0], xy = self.annotate_pos,
                            xycoords = 'axes fraction',
                            color = self.xcolor,
                            **self.annotate_kwargs)
            self.anx.set_visible(self.GetPlotParam('show_x'))

            self.annotate_pos[0] += .08
            if self.flagy >0 and self.GetPlotParam('show_y'):
                if self.flagy == 1 and len(self.flagy.shape) == 1:
                    self.liney = self.axes.plot(self.xaxis_values, self.fy, color = self.ycolor)
                elif self.parent.MainParamDict['Average1D']:
                    self.liney = self.axes.plot(self.xaxis_values, np.average(self.fy.reshape(-1,self.fy.shape[-1]), axis = 0), color = self.ycolor)
                else:
                    self.liney = self.axes.plot(self.xaxis_values, self.fy[self.parent.zSlice,self.parent.ySlice,:], color = self.ycolor)

                min_max[0]=min(min_max[0],self.liney[0].get_data()[1].min())
                min_max[1]=max(min_max[1],self.liney[0].get_data()[1].max())
            else:
                self.liney = self.axes.plot(np.arange(10), np.arange(10), color = self.ycolor)
                self.liney[0].set_visible(False)

            self.any =self.axes.annotate(self.GetPlotParam('1D_label')[self.GetPlotParam('field_type')][1], xy = self.annotate_pos,
                            xycoords= 'axes fraction',
                            color = self.ycolor,
                            **self.annotate_kwargs)

            self.any.set_visible(self.GetPlotParam('show_y'))
            self.annotate_pos[0] += .08

            if self.flagz and self.GetPlotParam('show_z'):
                if self.flagx == 1 and len(self.fz.shape) == 1:
                    self.linez = self.axes.plot(self.xaxis_values, self.fz, color = self.zcolor)
                if self.parent.MainParamDict['Average1D']:
                    self.linez = self.axes.plot(self.xaxis_values, np.average(self.fz.reshape(-1,self.fz.shape[-1]), axis = 0), color = self.zcolor)
                else: # In the x-y plane
                    self.linez = self.axes.plot(self.xaxis_values, self.fz[self.parent.zSlice,self.parent.ySlice,:], color = self.zcolor)
                min_max[0]=min(min_max[0],self.linez[0].get_data()[1].min())
                min_max[1]=max(min_max[1],self.linez[0].get_data()[1].max())

            else:
                self.linez = self.axes.plot(np.arange(10), np.arange(10), color = self.zcolor)
                self.linez[0].set_visible(False)

            if np.isinf(min_max[0]):
                min_max[0]=None
                min_max[1]=None
            else:
                dist = min_max[1]-min_max[0]
                min_max[0] -= 0.04*dist
                min_max[1] += 0.04*dist
            self.axes.set_ylim(min_max)

            self.anz = self.axes.annotate(self.GetPlotParam('1D_label')[self.GetPlotParam('field_type')][2], xy = self.annotate_pos,
                                xycoords= 'axes fraction',
                                color = self.zcolor,
                                **self.annotate_kwargs
                                )
            self.anz.set_visible(self.GetPlotParam('show_z'))

            self.shock_line = self.axes.axvline(self.parent.shock_loc, linewidth = 1.5, linestyle = '--', color = self.parent.shock_color, path_effects=[PathEffects.Stroke(linewidth=2, foreground='k'),
                        PathEffects.Normal()])

            self.shock_line.set_visible(self.GetPlotParam('show_shock'))

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
                self.axes.set_xlim(self.xaxis_values[0],self.xaxis_values[-1])


            if self.GetPlotParam('set_v_min'):
                self.axes.set_ylim(ymin = self.GetPlotParam('v_min'))
            if self.GetPlotParam('set_v_max'):
                self.axes.set_ylim(ymax = self.GetPlotParam('v_max'))

            self.axes.set_xlabel(r'$x\ [c/\omega_{\rm pe}]$', labelpad = self.parent.MainParamDict['xLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
            tmplblstr = self.GetPlotParam('yaxis_label')[self.GetPlotParam('field_type')]

            if self.GetPlotParam('normalize_fields'):
                if self.GetPlotParam('field_type') ==0:
                    tmplblstr +=r'$/B_0$'
                elif self.GetPlotParam('field_type') ==1:
                    tmplblstr +=r'$/E_0$'

            self.axes.set_ylabel(tmplblstr, labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
        ####
        # FFT REGION PLOTTING CODE
        ####
        self.lineleft = self.axes.axvline(0, linewidth = 1.5, linestyle = ':', color = self.parent.FFT_color)
        self.lineright = self.axes.axvline(0, linewidth = 1.5, linestyle = ':', color = self.parent.FFT_color)
        self.lineleft.set_visible(self.GetPlotParam('show_FFT_region'))
        self.lineright.set_visible(self.GetPlotParam('show_FFT_region'))

        if self.GetPlotParam('show_FFT_region'):
            self.left_loc = self.parent.MainParamDict['FFTLeft'] + self.parent.shock_loc*self.parent.MainParamDict['FFTRelative']
            self.left_loc = max(self.left_loc, self.xmin)
            self.lineleft.set_xdata([self.left_loc,self.left_loc])

            self.right_loc = self.parent.MainParamDict['FFTRight'] + self.parent.shock_loc*self.parent.MainParamDict['FFTRelative']
            self.right_loc = min(self.right_loc, self.xmax)
            self.lineright.set_xdata([self.right_loc,self.right_loc])

        ####
        #
        # Code to show the CPU domains
        #
        ####

        if self.GetPlotParam('show_cpu_domains'):
            self.FigWrap.SetCpuDomainLines()
    def refresh(self):

        '''This is a function that will be called only if self.axes already
        holds a fields type plot. We only update things that have changed & are
        shown.  If hasn't changed or isn't shown, don't touch it. The difference
        between this and last time, is that we won't actually do any drawing in
        the plot. The plot will be redrawn after all subplots are refreshed. '''


        # Main goal, only change what is showing..

        self.xmin, self.xmax = self.xaxis_values[0], self.xaxis_values[-1]
        self.lineleft.set_visible(self.GetPlotParam('show_FFT_region'))
        self.lineright.set_visible(self.GetPlotParam('show_FFT_region'))

        if self.GetPlotParam('show_FFT_region'):
            # Update the position of the FFT region
            self.left_loc = self.parent.MainParamDict['FFTLeft'] + self.parent.shock_loc*self.parent.MainParamDict['FFTRelative']
            self.left_loc = max(self.left_loc, self.xmin)
            self.lineleft.set_xdata([self.left_loc,self.left_loc])

            self.right_loc = self.parent.MainParamDict['FFTRight'] + self.parent.shock_loc*self.parent.MainParamDict['FFTRelative']
            self.right_loc = min(self.right_loc, self.xmax)
            self.lineright.set_xdata([self.right_loc,self.right_loc])

        # Now do the 1D plots, because it is simpler
        if self.GetPlotParam('twoD') == 0:
            min_max = [np.inf, -np.inf]
            if self.GetPlotParam('show_x') and self.flagx:
                if self.flagx == 1 and len(self.fx.shape) == 1:
                    self.linex[0].set_data(self.xaxis_values, self.fx)
                elif self.parent.MainParamDict['Average1D']:
                    self.linex[0].set_data(self.xaxis_values, np.average(self.fx.reshape(-1,self.fx.shape[-1]), axis =0))
                else: # In the x-y plane
                    self.linex[0].set_data(self.xaxis_values, self.fx[self.parent.zSlice,self.parent.ySlice,:])
                self.linex[0].set_visible(True)
                self.anx.set_visible(True)
                min_max[0]=min(min_max[0],self.linex[0].get_data()[1].min())
                min_max[1]=max(min_max[1],self.linex[0].get_data()[1].max())

            if self.GetPlotParam('show_y') and self.flagy:
                if self.flagy == 1 and len(self.fy.shape) == 1:
                    self.liney[0].set_data(self.xaxis_values, self.fy)
                elif self.parent.MainParamDict['Average1D']:
                    self.liney[0].set_data(self.xaxis_values, np.average(self.fy.reshape(-1,self.fy.shape[-1]), axis =0))
                else:
                    self.liney[0].set_data(self.xaxis_values, self.fy[self.parent.zSlice,self.parent.ySlice,:])
                self.liney[0].set_visible(True)
                self.any.set_visible(True)
                min_max[0]=min(min_max[0],self.liney[0].get_data()[1].min())
                min_max[1]=max(min_max[1],self.liney[0].get_data()[1].max())

            if self.GetPlotParam('show_z'):
                if self.flagz ==1 and len(self.fz.shape) == 1:
                    self.linez[0].set_data(self.xaxis_values, self.fz)
                elif self.parent.MainParamDict['Average1D']:
                    self.linez[0].set_data(self.xaxis_values, np.average(self.fz.reshape(-1,self.fz.shape[-1]), axis =0))
                else:
                    self.linez[0].set_data(self.xaxis_values, self.fz[self.parent.zSlice,self.parent.ySlice,:])
                self.linez[0].set_visible(True)
                self.anz.set_visible(True)
                min_max[0]=min(min_max[0],self.linez[0].get_data()[1].min())
                min_max[1]=max(min_max[1],self.linez[0].get_data()[1].max())

            if np.isinf(min_max[0]):
                min_max[0]=None
                min_max[1]=None
            else:
                dist = min_max[1]-min_max[0]
                min_max[0] -= 0.04*dist
                min_max[1] += 0.04*dist
            self.axes.set_ylim(min_max)

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

            if self.GetPlotParam('set_v_min'):
                self.axes.set_ylim(ymin = self.GetPlotParam('v_min'))
            if self.GetPlotParam('set_v_max'):
                self.axes.set_ylim(ymax = self.GetPlotParam('v_max'))


        else: # Now refresh the plot if it is 2D
            self.plotFlag = -1

            if self.GetPlotParam('show_x') and self.flagx >1:
                self.plotFlag = 0
                if self.parent.MainParamDict['2DSlicePlane'] == 0: #x-y plane
                    self.cax.set_data(self.fx[self.parent.zSlice,:,:])
                elif self.parent.MainParamDict['2DSlicePlane'] == 1: #x-z plane
                    self.cax.set_data(self.fx[:,self.parent.ySlice,:])

            elif self.GetPlotParam('show_y') and self.flagy >1:
                self.plotFlag = 1
                if self.parent.MainParamDict['2DSlicePlane'] == 0: #x-y plane
                    self.cax.set_data(self.fy[self.parent.zSlice,:,:])
                elif self.parent.MainParamDict['2DSlicePlane'] == 1: #x-z plane
                    self.cax.set_data(self.fy[:,self.parent.ySlice,:])

            elif self.GetPlotParam('show_z') and self.flagz>1:
                self.plotFlag = 2
                if self.parent.MainParamDict['2DSlicePlane'] == 0: #x-y plane
                    self.cax.set_data(self.fz[self.parent.zSlice,:,:])
                elif self.parent.MainParamDict['2DSlicePlane'] == 1: #x-z plane
                    self.cax.set_data(self.fz[:,self.parent.ySlice,:])
            else:
                self.cax.set_data(np.ma.masked_array(np.empty([2,2]), mask = np.ones([2,2])))
                self.clims = [None, None]

            self.axC.set_visible(self.plotFlag !=-1)
            if self.plotFlag != -1:
                self.ymin = 0
                self.ymax =  self.cax.get_array().shape[0]/self.c_omp*self.istep
                self.xmin = 0
                self.xmax =  self.xaxis_values[-1]
                self.clims = [self.cax.get_array().min(), self.cax.get_array().max()]


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
            if self.plotFlag >= 0:
                self.anntext = self.GetPlotParam('2D_label')[self.GetPlotParam('field_type')][self.plotFlag]
                if self.GetPlotParam('field_type') ==0  and self.GetPlotParam('normalize_fields'):
                    self.anntext +=r'$/B_0$'
                if self.GetPlotParam('field_type') ==1  and self.GetPlotParam('normalize_fields'):
                    self.anntext +=r'$/E_0$'
                self.TwoDan.set_text(self.anntext)
            else:
                self.TwoDan.set_text('')

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

            if self.GetPlotParam('show_shock'):
                self.shockline_2d.set_xdata([self.parent.shock_loc,self.parent.shock_loc])
            if self.parent.MainParamDict['2DSlicePlane'] == 0:
                self.axes.set_ylabel(r'$y\ [c/\omega_{\rm pe}]$', labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
            if self.parent.MainParamDict['2DSlicePlane'] == 1:
                self.axes.set_ylabel(r'$z\ [c/\omega_{\rm pe}]$', labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])

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
            self.settings_window = FieldSettings(self)
        else:
            self.settings_window.destroy()
            self.settings_window = FieldSettings(self)


class FieldSettings(Tk.Toplevel):
    def __init__(self, parent):
        self.parent = parent
        Tk.Toplevel.__init__(self)
        self.def1_window = None
        self.def2_window = None
        self.def3_window = None
        self.wm_title('Fields & Currents Plot (%d,%d) Settings' % self.parent.FigWrap.pos)
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
        self.FieldList = ['B Field', 'E field', 'J [current]', 'User Defined']
        self.FieldTypeVar  = Tk.IntVar()
        self.FieldTypeVar.set(self.parent.GetPlotParam('field_type'))

        ttk.Label(self.frm, text='Choose Field:').grid(row = 2, sticky = Tk.W)

        for i in range(len(self.FieldList)):
            ttk.Radiobutton(self.frm,
                text=self.FieldList[i],
                variable=self.FieldTypeVar,
                command = self.RadioField,
                value=i).grid(row = 3+i, sticky =Tk.W)

        # the Check boxes for the dimension
        self.label = ttk.Label(self.frm, text='Dimension:')
        self.label.grid(row = 2, column = 1, sticky = Tk.W)

        self.ShowXVar = Tk.IntVar(self) # Create a var to track whether or not to show X
        self.ShowXVar.set(self.parent.GetPlotParam('show_x'))
        self.cbx = ttk.Checkbutton(self.frm, text = "Show x",
            variable = self.ShowXVar,
            command = self.Selector)
        self.cbx.grid(row = 3, column = 1, sticky = Tk.W)

        self.ShowYVar = Tk.IntVar(self) # Create a var to track whether or not to plot Y
        self.ShowYVar.set(self.parent.GetPlotParam('show_y'))
        self.cby = ttk.Checkbutton(self.frm, text = "Show y",
            variable = self.ShowYVar,
            command = self.Selector)
        self.cby.grid(row = 4, column = 1, sticky = Tk.W)

        self.ShowZVar = Tk.IntVar(self) # Create a var to track whether or not to plot Z
        self.ShowZVar.set(self.parent.GetPlotParam('show_z'))
        self.cbz = ttk.Checkbutton(self.frm, text = "Show z",
            variable = self.ShowZVar,
            command = self.Selector)
        self.cbz.grid(row = 5, column = 1, sticky = Tk.W)
        if self.FieldTypeVar.get()==3:
            # ADD BUTTONS TO DEFINE THE FUNCTIONS
            self.df1button = ttk.Button(self.frm, text = 'Def F1', command = self.OpenDef1)
            self.df1button.grid(row =3, column =2)
            self.df2button = ttk.Button(self.frm, text = 'Def F2', command = self.OpenDef2)
            self.df2button.grid(row =4, column =2)
            self.df3button = ttk.Button(self.frm, text = 'Def F3', command = self.OpenDef3)
            self.df3button.grid(row =5, column =2)
            # CHANGE LABELS
            self.cbx.config(text='Show F1')
            self.cby.config(text='Show F2')
            self.cbz.config(text='Show F3')
            self.label.config(text='Choose Function:')
        # Control whether or not Cbar is shown
        self.CbarVar = Tk.IntVar()
        self.CbarVar.set(self.parent.GetPlotParam('show_cbar'))
        cb = ttk.Checkbutton(self.frm, text = "Show Color bar",
                        variable = self.CbarVar,
                        command = self.CbarHandler)
        cb.grid(row = 7, sticky = Tk.W)

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
        cb = ttk.Checkbutton(self.frm, text = "Asymmetric Color Space",
                        variable = self.StretchVar,
                        command = self.StretchHandler)
        cb.grid(row = 8, column = 1, sticky = Tk.W)

        # Create the OptionMenu to chooses the cnorm_type:
        self.cnormvar = Tk.StringVar(self)
        self.cnormvar.set(self.parent.chartType) # default value
        self.cnormvar.trace('w', self.cnormChanged)

        ttk.Label(self.frm, text="Choose Color Norm:").grid(row=6, column = 3)
        cnormChooser = ttk.OptionMenu(self.frm, self.cnormvar, self.parent.GetPlotParam('cnorm_type'), *tuple(['Pow', 'Linear']))
        cnormChooser.grid(row =6, column = 4, sticky = Tk.W + Tk.E)

        # Now the gamma of the pow norm
        self.powGamma = Tk.StringVar()
        self.powGamma.set(str(self.parent.GetPlotParam('cpow_num')))
        ttk.Label(self.frm, text ='gamma =').grid(row = 7, column = 3, sticky =Tk.E)
        ttk.Label(self.frm, text ='If cnorm is Pow =>').grid(row = 8, column = 3,columnspan = 2, sticky =Tk.N)
        ttk.Label(self.frm, text ='sign(data)*|data|**gamma').grid(row = 9, column = 3,columnspan = 2, sticky =Tk.E)

        self.GammaEnter = ttk.Entry(self.frm, textvariable=self.powGamma, width=7)
        self.GammaEnter.grid(row = 7, column = 4)



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


        cb = ttk.Checkbutton(self.frm, text ='Set B or E min',
                        variable = self.setZminVar)
        cb.grid(row = 3, column = 3, sticky = Tk.W)
        self.ZminEnter = ttk.Entry(self.frm, textvariable=self.Zmin, width=7)
        self.ZminEnter.grid(row = 3, column = 4)

        cb = ttk.Checkbutton(self.frm, text ='Set B or E max',
                        variable = self.setZmaxVar)
        cb.grid(row = 4, column = 3, sticky = Tk.W)

        self.ZmaxEnter = ttk.Entry(self.frm, textvariable=self.Zmax, width=7)
        self.ZmaxEnter.grid(row = 4, column = 4)

        self.ShockVar = Tk.IntVar()
        self.ShockVar.set(self.parent.GetPlotParam('show_shock'))
        cb = ttk.Checkbutton(self.frm, text = "Show Shock",
                        variable = self.ShockVar,
                        command = self.ShockVarHandler)
        cb.grid(row = 9, column = 1, sticky = Tk.W)

        self.FFTVar = Tk.IntVar()
        self.FFTVar.set(self.parent.GetPlotParam('show_FFT_region'))
        cb = ttk.Checkbutton(self.frm, text = "Show FFT Region",
                        variable = self.FFTVar,
                        command = self.FFTVarHandler)
        cb.grid(row = 9, column = 0, sticky = Tk.W)

        self.CPUVar = Tk.IntVar()
        self.CPUVar.set(self.parent.GetPlotParam('show_cpu_domains'))
        cb = ttk.Checkbutton(self.frm, text = "Show CPU domains",
                        variable = self.CPUVar,
                        command = self.CPUVarHandler)
        cb.grid(row = 10, column = 0, sticky = Tk.W)

        self.NormFieldVar = Tk.IntVar()
        self.NormFieldVar.set(self.parent.GetPlotParam('normalize_fields'))
        cb = ttk.Checkbutton(self.frm, text = "Normalize Fields",
                        variable = self.NormFieldVar,
                        command = self.NormFieldHandler)
        cb.grid(row = 7, column = 1, sticky = Tk.W)


    def CbarHandler(self, *args):
        if self.parent.GetPlotParam('show_cbar')== self.CbarVar.get():
            pass
        else:
            if self.parent.GetPlotParam('twoD'):
                self.parent.axC.set_visible(self.CbarVar.get())

            self.parent.SetPlotParam('show_cbar', self.CbarVar.get(), update_plot = self.parent.GetPlotParam('twoD'))

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

    def OpenDef1(self):
        if self.def1_window is None:
            self.def1_window = UserDefSettings(self, self.parent,1)
        else:
            self.def1_window.destroy()
            self.def1_window = UserDefSettings(self, self.parent,1)
    def OpenDef2(self):
        if self.def2_window is None:
            self.def2_window = UserDefSettings(self, self.parent,2)
        else:
            self.def2_window.destroy()
            self.def2_window = UserDefSettings(self, self.parent,2)
    def OpenDef3(self):
        if self.def3_window is None:
            self.def3_window = UserDefSettings(self, self.parent,3)
        else:
            self.def3_window.destroy()
            self.def3_window = UserDefSettings(self, self.parent,3)

    def ShockVarHandler(self, *args):
        if self.parent.GetPlotParam('show_shock')== self.ShockVar.get():
            pass
        else:
            if self.parent.GetPlotParam('twoD'):
                self.parent.shockline_2d.set_visible(self.ShockVar.get())
            else:
                self.parent.shock_line.set_visible(self.ShockVar.get())

            self.parent.SetPlotParam('show_shock', self.ShockVar.get())

    def FFTVarHandler(self, *args):
        if self.parent.GetPlotParam('show_FFT_region')== self.FFTVar.get():
            pass
        else:
            self.parent.SetPlotParam('show_FFT_region', self.FFTVar.get(), update_plot = False)
            self.parent.lineleft.set_visible(self.parent.GetPlotParam('show_FFT_region'))
            self.parent.lineright.set_visible(self.parent.GetPlotParam('show_FFT_region'))

            ### The .parent.parent is less than ideal.... consider re-writing.
            if self.parent.GetPlotParam('show_FFT_region'):
                self.parent.left_loc = self.parent.parent.MainParamDict['FFTLeft'] + self.parent.parent.shock_loc*self.parent.parent.MainParamDict['FFTRelative']
                self.parent.left_loc = max(self.parent.left_loc, self.parent.xmin)
                self.parent.lineleft.set_xdata([self.parent.left_loc,self.parent.left_loc])

                self.parent.right_loc = self.parent.parent.MainParamDict['FFTRight'] + self.parent.parent.shock_loc*self.parent.parent.MainParamDict['FFTRelative']
                self.parent.right_loc = min(self.parent.right_loc, self.parent.xmax)
                self.parent.lineright.set_xdata([self.parent.right_loc,self.parent.right_loc])
            self.parent.parent.canvas.draw()
            self.parent.parent.canvas.get_tk_widget().update_idletasks()

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

    def NormFieldHandler(self, *args):
        if self.parent.GetPlotParam('normalize_fields') == self.NormFieldVar.get():
            pass
        else:
            if ~self.parent.GetPlotParam('twoD'):
                tmplblstr = self.parent.GetPlotParam('yaxis_label')[self.FieldTypeVar.get()]
                if self.NormFieldVar.get():
                    if self.parent.GetPlotParam('field_type') ==0:
                        tmplblstr +=r'$/B_0$'
                    elif self.parent.GetPlotParam('field_type') ==1:
                        tmplblstr +=r'$/E_0$'

                self.parent.axes.set_ylabel(tmplblstr, labelpad = self.parent.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.parent.MainParamDict['AxLabelSize'])


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

    def RadioField(self):
        if self.FieldTypeVar.get() == self.parent.GetPlotParam('field_type'):
            pass
        else:
            if not self.parent.GetPlotParam('twoD'):
                tmplblstr = self.parent.GetPlotParam('yaxis_label')[self.FieldTypeVar.get()]
                if self.parent.GetPlotParam('normalize_fields'):
                    if self.FieldTypeVar.get() ==0:
                        tmplblstr +=r'$/B_0$'
                    elif self.FieldTypeVar.get() ==1:
                        tmplblstr +=r'$/E_0$'

                self.parent.axes.set_ylabel(tmplblstr, labelpad = self.parent.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.parent.MainParamDict['AxLabelSize'])


                self.parent.anx.set_text(self.parent.GetPlotParam('1D_label')[self.FieldTypeVar.get()][0])
                self.parent.any.set_text(self.parent.GetPlotParam('1D_label')[self.FieldTypeVar.get()][1])
                self.parent.anz.set_text(self.parent.GetPlotParam('1D_label')[self.FieldTypeVar.get()][2])

            ####
            #
            # Plot to add UserDef fields
            #
            #####
            if self.FieldTypeVar.get()==3:
                # ADD BUTTONS TO DEFINE THE FUNCTIONS
                self.df1button = ttk.Button(self.frm, text = 'Def F1', command = self.OpenDef1)
                self.df1button.grid(row =3, column =2)
                self.df2button = ttk.Button(self.frm, text = 'Def F2', command = self.OpenDef2)
                self.df2button.grid(row =4, column =2)
                self.df3button = ttk.Button(self.frm, text = 'Def F3', command = self.OpenDef3)
                self.df3button.grid(row =5, column =2)
                # CHANGE THE LABELS OF ALL THE CHECKBUTTONS
                self.label.config(text='Choose Function:')
                self.cbx.config(text='Show F1')
                self.cby.config(text='Show F2')
                self.cbz.config(text='Show F3')
                # TURN OFF ALL THE LINES
                self.ShowXVar.set(False)
                self.ShowYVar.set(False)
                self.ShowZVar.set(False)
                self.parent.SetPlotParam('show_x', False, update_plot = False)
                self.parent.SetPlotParam('show_y', False, update_plot = False)
                self.parent.SetPlotParam('show_z', False, update_plot = False)
                if ~self.parent.GetPlotParam('twoD'):
                    self.parent.linex[0].set_visible(False)
                    self.parent.anx.set_visible(False)
                    self.parent.liney[0].set_visible(False)
                    self.parent.any.set_visible(False)
                    self.parent.linez[0].set_visible(False)
                    self.parent.anz.set_visible(False)

            elif self.parent.GetPlotParam('field_type') ==3 :
                ### DESTROY THE buttons
                self.df1button.destroy()
                self.df2button.destroy()
                self.df3button.destroy()
                self.label.config(text='Dimension:')
                self.cbx.config(text='Show x')
                self.cby.config(text='Show y')
                self.cbz.config(text='Show z')
                self.ShowXVar.set(False)
                self.ShowYVar.set(False)
                self.ShowZVar.set(False)
                self.parent.SetPlotParam('show_x', False, update_plot = False)
                self.parent.SetPlotParam('show_y', False, update_plot = False)
                self.parent.SetPlotParam('show_z', False, update_plot = False)
                if ~self.parent.GetPlotParam('twoD'):
                    self.parent.linex[0].set_visible(False)
                    self.parent.anx.set_visible(False)
                    self.parent.liney[0].set_visible(False)
                    self.parent.any.set_visible(False)
                    self.parent.linez[0].set_visible(False)
                    self.parent.anz.set_visible(False)

            self.parent.SetPlotParam('field_type', self.FieldTypeVar.get())

    def Selector(self):
        # First check if it is 2-D:
        if self.parent.GetPlotParam('twoD'):

            #if self.ShowXVar.get() == 0 and self.ShowYVar.get() == 0 and self.ShowZVar.get() == 0:
            #    # All are zero, something must be selected for this plot
            #    self.ShowXVar.set(1)


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
        try:
            #make sure the user types in a float, no longer check to see if it has changed, because of precision issues.
            self.parent.SetPlotParam('v_min', float(self.Zmin.get()), update_plot = False)
            if self.parent.GetPlotParam('set_v_min'):
                if self.parent.GetPlotParam('twoD'):
                    self.parent.cax.norm.vmin = self.parent.GetPlotParam('v_min')
                else:
                    self.parent.axes.set_ylim(ymin = self.parent.GetPlotParam('v_min') )
        except ValueError:
            #if they type in random stuff, just set it ot the param value
            self.Zmin.set(str(self.parent.GetPlotParam('v_min')))

        try:
            #make sure the user types in a float, no longer check to see if it has changed, because of precision issues.
            self.parent.SetPlotParam('v_max', float(self.Zmax.get()), update_plot = False)
            if self.parent.GetPlotParam('set_v_max'):
                if self.parent.GetPlotParam('twoD'):
                    self.parent.cax.norm.vmax = self.parent.GetPlotParam('v_max')
                else:
                    self.parent.axes.set_ylim(ymax = self.parent.GetPlotParam('v_max') )
        except ValueError:
            #if they type in random stuff, just set it ot the param value
            self.Zmax.set(str(self.parent.GetPlotParam('v_max')))
        self.parent.SetPlotParam('v_max', self.parent.GetPlotParam('v_max'))


    def OnClosing(self):
        self.parent.settings_window = None
        self.destroy()

class UserDefSettings(Tk.Toplevel):
    def __init__(self, parent, subplot, fnum):
        self.parent = parent
        self.subplot = subplot
        self.fnum = fnum
        Tk.Toplevel.__init__(self)

        self.wm_title('Define Fuction %d' % fnum)
        self.parent = parent

        S = Tk.Scrollbar(self)
        self.T = Tk.Text(self, height=25, width=100)
        S.pack(side=Tk.RIGHT, fill=Tk.Y)
        self.T.pack(side=Tk.TOP, fill=Tk.Y)
        S.config(command=self.T.yview)
        self.T.config(yscrollcommand=S.set)

        tmpstr = ''
        for line in self.subplot.GetPlotParam('cmdstr'+str(self.fnum)).splitlines():
            tmpstr += line[1:] +'\n'
        self.T.insert(Tk.END, tmpstr)
        miniframe = ttk.Frame(self)
        ttk.Label(miniframe, text ="1D y-label:").grid(row=0, column =2)
        self.ylabel = Tk.StringVar()
        self.ylabel.set(self.subplot.GetPlotParam('yaxis_label')[self.subplot.GetPlotParam('field_type')])
        ttk.Entry(miniframe, textvariable=self.ylabel, width=15).grid(row = 0, column = 3)

        ttk.Label(miniframe, text ="1D label:").grid(row=0, column =0)
        self.oneDlabel = Tk.StringVar()
        self.oneDlabel.set(self.subplot.GetPlotParam('1D_label')[3][self.fnum-1])
        ttk.Entry(miniframe, textvariable=self.oneDlabel, width=15).grid(row = 0, column = 1)

        ttk.Label(miniframe, text ="2D label:").grid(row=1, column =0)
        self.twoDlabel = Tk.StringVar()
        self.twoDlabel.set(self.subplot.GetPlotParam('2D_label')[3][self.fnum-1])
        ttk.Entry(miniframe, textvariable=self.twoDlabel, width=15).grid(row = 1, column = 1)

        self.OneDVar = Tk.IntVar()
        self.OneDVar.set(self.subplot.GetPlotParam('OneDOnly')[self.fnum-1])
        ttk.Checkbutton(miniframe, text = 'Returns a 1D array along x', variable = self.OneDVar).grid(row = 1, column = 2, columnspan= 2, sticky = Tk.W)
        miniframe.pack(side=Tk.TOP)
        ttk.Button(self, text = 'Save F'+str(self.fnum), command = self.SaveStr).pack(side =Tk.TOP)
    def SaveStr(self):
        tmpstr = ''
        for line in self.T.get(1.0, Tk.END).splitlines():
            tmpstr += '#' + line + '\n'

        self.subplot.SetPlotParam('cmdstr'+str(self.fnum), tmpstr, update_plot=False)

        ### THIS IS SLOPPY!
        self.subplot.SetPlotParam('yaxis_label',self.subplot.GetPlotParam('yaxis_label')[0:3]+ [self.ylabel.get()], update_plot =False)
        tmplist = list(self.subplot.GetPlotParam('2D_label')[3])
        tmplist[self.fnum-1] = self.twoDlabel.get()
        tmplist2 = list(self.subplot.GetPlotParam('2D_label')[0:3])
        tmplist2.append(tmplist)
        self.subplot.SetPlotParam('2D_label',tmplist2, update_plot =False)

        tmplist = self.subplot.GetPlotParam('1D_label')[3]
        tmplist[self.fnum-1] = self.oneDlabel.get()
        tmplist2 = list(self.subplot.GetPlotParam('1D_label')[0:3])
        tmplist2.append(tmplist)

        self.subplot.SetPlotParam('1D_label',tmplist2, update_plot =False)
        if ~self.subplot.GetPlotParam('twoD'):
            self.subplot.axes.set_ylabel(self.subplot.GetPlotParam('yaxis_label')[3])
            self.subplot.anx.set_text(self.subplot.GetPlotParam('1D_label')[3][0])
            self.subplot.any.set_text(self.subplot.GetPlotParam('1D_label')[3][1])
            self.subplot.anz.set_text(self.subplot.GetPlotParam('1D_label')[3][2])

        #self.subplot.GetPlotParam('1D_label')[self.subplot.GetPlotParam('field_type')][self.fnum-1] = self.oneDlabel.get()
        #self.subplot.GetPlotParam('2D_label')[self.subplot.GetPlotParam('field_type')][self.fnum-1] = self.twoDlabel.get()
        self.subplot.GetPlotParam('OneDOnly')[self.fnum -1] = self.OneDVar.get()
        if self.fnum ==1:
            self.subplot.SetPlotParam('show_x', True)
            self.parent.ShowXVar.set(True)
        if self.fnum ==2:
            self.subplot.SetPlotParam('show_y', True)
            self.parent.ShowYVar.set(True)
        if self.fnum ==3:
            self.subplot.SetPlotParam('show_z', True)
            self.parent.ShowZVar.set(True)
        self.OnClosing()

    def OnClosing(self):
        if self.fnum ==1:
            self.parent.def1_window = None
        if self.fnum ==2:
            self.parent.def2_window = None
        else:
            self.parent.def3_window = None
        self.destroy()
