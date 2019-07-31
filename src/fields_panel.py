import matplotlib
import numpy as np
import numpy.ma as ma

import sys, traceback
sys.path.append('../')
import new_cmaps
from new_cnorms import PowerNormWithNeg, PowerNormFunc
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
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
                       'show_cpu_domains': False, # plots lines showing how the CPUs are divvying up the computational region
                       'face_color': 'gainsboro' }

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
        self.chartType = 'FieldsPlot'
        self.figure = self.parent.figure
        #self.SetPlotParam('spatial_y', self.GetPlotParam('twoD'), update_plot = False)
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




    def draw(self, output):
        ''' A Helper function that loads the data for the plot'''
        # First see of the x_axis and y_axis values have already been calculated
        # and stored in the DataDict for this time step
        self.c_omp = getattr(output,'c_omp')
        self.istep = getattr(output,'istep')
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
            self.param_dict['normalize_fields'] = 0


        # x-values haven't been calculated yet, generate them then save them to the dictionary for later.
        if self.GetPlotParam('field_type') ==0 or self.GetPlotParam('field_type') == 3:
            self.xaxis_values = np.arange(getattr(output,'bx').shape[2])/self.c_omp*self.istep
        elif self.GetPlotParam('field_type') ==1:
            self.xaxis_values = np.arange(getattr(output,'ex').shape[2])/self.c_omp*self.istep
        elif self.GetPlotParam('field_type') ==2:
            self.xaxis_values = np.arange(getattr(output,'jx').shape[2])/self.c_omp*self.istep

        self.flagx = 0 # 0 means it didn't plot, 1 means it is 1D only, 2 means it returned a 3d object
        self.flagy = 0
        self.flagz = 0

        if self.GetPlotParam('field_type') == 0: # Load the B-Field
            if self.GetPlotParam('show_x'):
                self.flagx = 2
                if self.GetPlotParam('normalize_fields'):
                    self.fx = getattr(output,'bx')*self.parent.b0**-1
                else:
                    self.fx = getattr(output,'bx')
            if self.GetPlotParam('show_y'):
                self.flagy = 2
                if self.GetPlotParam('normalize_fields'):
                    self.fy = getattr(output,'by')*self.parent.b0**-1
                else:
                    self.fy = getattr(output,'by')

            if self.GetPlotParam('show_z'):
                self.flagz = 2
                if self.GetPlotParam('normalize_fields'):
                    self.fz =getattr(output,'bz')*self.parent.b0**-1
                else:
                    self.fz =getattr(output,'bz')

        if self.GetPlotParam('field_type') == 1: # Load the E-Field
            if self.GetPlotParam('show_x'):
                self.flagx = 2
                if self.GetPlotParam('normalize_fields'):
                    self.fx = getattr(output,'ex')*self.parent.e0**-1
                else:
                    self.fx = getattr(output,'ex')
            if self.GetPlotParam('show_y'):
                self.flagy = 2
                if self.GetPlotParam('normalize_fields'):
                    self.fy = getattr(output,'ey')*self.parent.e0**-1
                else:
                    self.fy = getattr(output,'ey')

            if self.GetPlotParam('show_z'):
                self.flagz = 2
                if self.GetPlotParam('normalize_fields'):
                    self.fz =getattr(output,'ez')*self.parent.e0**-1
                else:
                    self.fz =getattr(output,'ez')


        elif self.GetPlotParam('field_type') == 2: # Load the currents

            if self.GetPlotParam('show_x'):
                self.fx = getattr(output,'jx')
                self.flagx = 2

            if self.GetPlotParam('show_y'):
                self.fy = getattr(output,'jy')
                self.flagy = 2

            if self.GetPlotParam('show_z'):
                self.fz = getattr(output,'jz')
                self.flagz = 2

        elif self.GetPlotParam('field_type') == 3: # User Defined fields
            if self.GetPlotParam('show_x'):
                try:
                    tmpcstr = ''
                    for line in self.GetPlotParam('cmdstr1').splitlines():
                        tmpcstr += line[1:] +'\n'
                    tmpcstr += 'self.fx = FieldFunc(*[getattr(output,k) for k in self.f1args])'
                    exec(compile(tmpcstr,'<string>', 'exec'), locals(), locals())#, '<string>', 'exec'), **{'self':self})
                    #print(FieldFunc)
                    if self.GetPlotParam('OneDOnly')[0]:
                        self.flagx = 1
                    else:
                        self.flagx = 2
                except:
                    print(sys.exc_info())
                    """
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

                    """
                    messagebox.showinfo('Error when evaluating user defined function 1:', print(sys.exc_info()))#(err_msg)

                    self.fx = np.NAN
                    self.flagx = 0

            if self.GetPlotParam('show_y'):
                try:
                    tmpcstr = ''
                    for line in self.GetPlotParam('cmdstr2').splitlines():
                        tmpcstr += line[1:] +'\n'
                    tmpcstr += 'self.fy = FieldFunc(*[getattr(output,k) for k in self.f2args])'
                    eval(compile(tmpcstr, '<string>', 'exec'), locals(), locals())
                    if self.GetPlotParam('OneDOnly')[1]:
                        self.flagy = 1
                    else:
                        self.flagy = 2
                except:
                    print(sys.exc_info())
                    """
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

                        """
                    messagebox.showinfo('Error when evaluating user defined function 2:', print(sys.exc_info()))#(err_msg)
                    self.fy = np.NAN
                    self.flagy = 0

            if self.GetPlotParam('show_z'):
                try:
                    tmpcstr = ''
                    for line in self.GetPlotParam('cmdstr3').splitlines():
                        tmpcstr += line[1:] +'\n'
                    tmpcstr += 'self.fz = FieldFunc(*[getattr(output,k) for k in self.f3args])'
                    eval(compile(tmpcstr, '<string>', 'exec'), locals(), locals())
                    if self.GetPlotParam('OneDOnly')[2]:
                        self.flagz = 1
                    else:
                        self.flagz = 2

                except:
                    print(sys.exc_info())
                    """
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

                    """
                    messagebox.showinfo('Error when evaluating user defined function 3:', print(sys.exc_info()))#(err_msg)

                    self.fz = np.NAN
                    self.flagz = 0
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

            # First choose the 'zval' to plot, we can only do one because it is 2-d.
            self.plotFlag = -1
            if self.GetPlotParam('show_x') and self.flagx == 2:
                if self.parent.MainParamDict['2DSlicePlane'] == 0: # Show the x-y plane
                    if self.parent.MainParamDict['ImageAspect']:
                        self.cax = self.axes.imshow(self.fx[self.parent.zSlice,:,:], norm = self.norm(), origin = 'lower')
                    else:
                        self.cax = self.axes.imshow(self.fx[self.parent.zSlice,:,:], origin = 'lower', norm = self.norm(),
                                                    aspect= 'auto')

                elif self.parent.MainParamDict['2DSlicePlane'] == 1: # Show the x-z plane
                    if self.parent.MainParamDict['ImageAspect']:
                        self.cax = self.axes.imshow(self.fx[:,self.parent.ySlice,:], norm = self.norm(), origin = 'lower')
                    else:
                        self.cax = self.axes.imshow(self.fx[:,self.parent.ySlice,:], origin = 'lower', norm = self.norm(),
                                                    aspect= 'auto')

                self.plotFlag = 0

            elif self.GetPlotParam('show_y') and self.flagy == 2:
                if self.parent.MainParamDict['2DSlicePlane'] == 0: # Show the x-y plane
                    if self.parent.MainParamDict['ImageAspect']:
                        self.cax = self.axes.imshow(self.fy[self.parent.zSlice,:,:], norm = self.norm(), origin = 'lower')
                    else:
                        self.cax = self.axes.imshow(self.fy[self.parent.zSlice,:,:], origin = 'lower', norm = self.norm(),
                                                    aspect= 'auto')
                elif self.parent.MainParamDict['2DSlicePlane'] == 1: # Show the x-z plane
                    if self.parent.MainParamDict['ImageAspect']:
                        self.cax = self.axes.imshow(self.fy[:,self.parent.ySlice,:], norm = self.norm(), origin = 'lower')
                    else:
                        self.cax = self.axes.imshow(self.fy[:,self.parent.ySlice,:], origin = 'lower', norm = self.norm(),
                                                    aspect= 'auto')
                self.plotFlag = 1


            elif self.GetPlotParam('show_z') and self.flagz == 2:
                # make sure z is loaded, (something has to be)
                # set the other plot values to zero in the PlotParams
                if self.parent.MainParamDict['2DSlicePlane'] == 0: # Show the x-y plane
                    if self.parent.MainParamDict['ImageAspect']:
                        self.cax = self.axes.imshow(self.fz[self.parent.zSlice,:,:], norm = self.norm(), origin = 'lower')
                    else:
                        self.cax = self.axes.imshow(self.fz[self.parent.zSlice,:,:], origin = 'lower', norm = self.norm(),
                                                    aspect= 'auto')
                elif self.parent.MainParamDict['2DSlicePlane'] == 1: # Show the x-z plane
                    if self.parent.MainParamDict['ImageAspect']:
                        self.cax = self.axes.imshow(self.fz[:,self.parent.ySlice,:], norm = self.norm(), origin = 'lower')
                    else:
                        self.cax = self.axes.imshow(self.fz[:,self.parent.ySlice,:], origin = 'lower', norm = self.norm(),
                                                    aspect= 'auto')

                self.plotFlag = 2


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
                                            cmap=new_cmaps.cmaps[self.cmap], origin='upper')

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
            self.parent.SetCpuDomainLines()
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
        return self.param_dict[keyname]
