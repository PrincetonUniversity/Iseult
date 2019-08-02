import matplotlib
import numpy as np
import numpy.ma as ma
import sys
sys.path.append('../')

import new_cmaps
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as PathEffects

class EnergyPanel:
    # A dictionary of all of the parameters for this plot with the default parameters

    plot_param_dict = {'twoD' : 1,
                       'masked': 1,
                       'cnorm_type': 'Log',
                       'prtl_type': 0,
                       'show_cbar': True,
                       'weighted': False,
                       'show_shock': False,
                       'show_int_region': True,
                       'set_color_limits': False,
                       'xbins' : 200,
                       'ebins' : 200,
                       'v_min': -2.0,
                       'v_max' : 0,
                       'set_v_min': False,
                       'set_v_max': False,
                       'set_y_min' : False,
                       'y_min': 1.0,
                       'set_y_max': False,
                       'y_max': 200.0,
                       'spatial_x': True,
                       'spatial_y': False,
                       'interpolation': 'nearest',
                       'face_color': 'gainsboro'}


    prtl_opts = ['proton', 'electron']
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
        self.figure = self.parent.figure
        self.chartType = 'EnergyPlot'
        self.InterpolationMethods = ['none','nearest', 'bilinear', 'bicubic', 'spline16',
            'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric',
            'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']
        # A variable that controls whether the energy integration region
        # is shown
        # Figure out the energy color the intergration region
        if self.GetPlotParam('prtl_type') == 1: #electons
            self.energy_color = self.parent.electron_color
        else:
            self.energy_color = self.parent.ion_color
        # A list that will hold any lines for the integration region


    def norm(self, vmin=None,vmax=None):
        if self.GetPlotParam('cnorm_type') == 'Log':
            return  mcolors.LogNorm(vmin, vmax)
        else:
            return mcolors.Normalize(vmin, vmax)

    def UpdateLabelsandColors(self):
        self.x_label = r'$x\ [c/\omega_{\rm pe}]$'

        if self.GetPlotParam('prtl_type') == 0: #protons
            self.energy_color = self.parent.ion_color
            self.y_label  = r'$E_p\ [m_i c^2]$'
        for line in self.IntRegionLines:
            line.set_color(self.energy_color)

        if self.GetPlotParam('prtl_type') == 1: #electons
            self.energy_color = self.parent.electron_color
            self.y_label  = r'$E_{e}\ [m_i c^2]$'
    def update_data(self, output):
        # Generate the X-axis values
        self.c_omp = getattr(output,'c_omp')
        self.istep = getattr(output,'istep')
        self.weights = None
        self.x_values = None
        self.y_values = None

        # Choose the particle type and px, py, or pz
        if self.GetPlotParam('prtl_type') == 0: #protons
            self.energy_color = self.parent.ion_color
            self.x_values = getattr(output,'xi')/self.c_omp
            if self.GetPlotParam('weighted'):
                self.weights = getattr(output,'chi')

            u = getattr(output,'ui')
            v = getattr(output,'vi')
            w = getattr(output,'wi')

        if self.GetPlotParam('prtl_type') == 1: #electons
            self.energy_color = self.parent.electron_color
            self.x_values = getattr(output,'xe')/self.c_omp

            if self.GetPlotParam('weighted'):
                self.weights = getattr(output,'che')
            u = getattr(output,'ue')
            v = getattr(output,'ve')
            w = getattr(output,'we')

        self.y_values = np.sqrt(u**2+v**2+w**2+1)-1
        if self.GetPlotParam('prtl_type') == 1:
            self.y_values *= getattr(output,'me')/getattr(output,'mi')
        self.Ymin = min(self.y_values)
        self.Ymax = max(self.y_values)
        self.Ymax = self.Ymax if ( self.Ymin != self.Ymax ) else self.Ymin+1

        self.xmin = 0
        self.xmax = getattr(output,'bx').shape[2]/self.c_omp*self.istep
        self.xmax = self.xmax if ( self.xmin != self.xmax ) else self.xmin+1

        self.hist2d = np.histogram2d(self.y_values, self.x_values,
                        bins = [self.GetPlotParam('ebins'), self.GetPlotParam('xbins')],
                        range = [[self.Ymin,self.Ymax],[0,self.xmax]],
                        weights = self.weights)

        if self.GetPlotParam('masked'):
            zval = ma.masked_array(self.hist2d[0])
            zval[zval == 0] = ma.masked
            zval *= float(zval.max())**(-1)
            tmplist = [zval[~zval.mask].min(), zval.max()]
        else:
            zval = np.copy(self.hist2d[0])
            zval[zval==0] = 0.5
            zval *= float(zval.max())**(-1)
            tmplist = [zval.min(), zval.max()]

        self.hist2d = zval, self.hist2d[1], self.hist2d[2], tmplist

    def draw(self):
        self.IntRegionLines = []

        self.UpdateLabelsandColors()

        self.xmin = self.hist2d[2][0]
        self.xmax = self.hist2d[2][-1]

        self.ymin = self.hist2d[1][0]
        self.ymax = self.hist2d[1][-1]


        if self.GetPlotParam('masked'):
            self.tick_color = 'k'
        else:
            self.tick_color = 'white'


        self.clim = list(self.hist2d[3])
        if self.GetPlotParam('set_v_min'):
            self.clim[0] = 10**self.GetPlotParam('v_min')
        if self.GetPlotParam('set_v_max'):
            self.clim[1] = 10**self.GetPlotParam('v_max')


        self.gs = gridspec.GridSpecFromSubplotSpec(100,100, subplot_spec = self.parent.gs0[self.pos])#, bottom=0.2,left=0.1,right=0.95, top = 0.95)

        self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]])

        self.cax = self.axes.imshow(self.hist2d[0],
                                    cmap = new_cmaps.cmaps[self.parent.MainParamDict['ColorMap']],
                                    norm = self.norm(), origin = 'lower',
                                    aspect = 'auto',
                                    interpolation=self.GetPlotParam('interpolation'))

        self.cax.set_extent([self.xmin, self.xmax, self.ymin, self.ymax])

        self.cax.set_clim(self.clim)

        self.shock_line = self.axes.axvline(self.parent.shock_loc, linewidth = 1.5, linestyle = '--', color = self.parent.shock_color, path_effects=[PathEffects.Stroke(linewidth=2, foreground='k'),
                   PathEffects.Normal()])
        if not self.GetPlotParam('show_shock'):
            self.shock_line.set_visible(False)

        self.axC = self.figure.add_subplot(self.gs[self.parent.cbar_extent[0]:self.parent.cbar_extent[1], self.parent.cbar_extent[2]:self.parent.cbar_extent[3]])

        # Technically I should use the colorbar class here,
        # but I found it annoying in some of it's limitations.
        if self.parent.MainParamDict['HorizontalCbars']:
            self.cbar = self.axC.imshow(self.gradient, aspect='auto',
                                    cmap=new_cmaps.cmaps[self.parent.MainParamDict['ColorMap']])
            # Make the colobar axis more like the real colorbar
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
                                    cmap=new_cmaps.cmaps[self.parent.MainParamDict['ColorMap']], origin='upper')
            # Make the colobar axis more like the real colorbar
            self.axC.tick_params(axis='x',
                                which = 'both', # bothe major and minor ticks
                                top = False, # turn off top ticks
                                bottom = False,
                                labelbottom = False,
                                labelsize=self.parent.MainParamDict['NumFontSize'])

            self.axC.tick_params(axis='y',          # changes apply to the y-axis
                                which='both',      # both major and minor ticks are affected
                                left= False,      # ticks along the bottom edge are off
                                right= True,         # ticks along the top edge are off
                                labelleft= False,
                                labelright=True,
                                labelsize=self.parent.MainParamDict['NumFontSize'])


        if not self.GetPlotParam('show_cbar'):
            self.axC.set_visible(False)


        if int(matplotlib.__version__[0]) < 2:
            self.axes.set_axis_bgcolor(self.GetPlotParam('face_color'))
        else:
            self.axes.set_facecolor(self.GetPlotParam('face_color'))

        self.axes.tick_params(labelsize = self.parent.MainParamDict['NumFontSize'], color=self.tick_color)

        self.axes.set_xlabel(self.x_label, labelpad = self.parent.MainParamDict['xLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
        self.axes.set_ylabel(self.y_label, labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])

        self.refresh()

    def refresh(self):
        '''This is a function that will be called only if self.axes already
        holds a density type plot. We only update things that have shown. If
        hasn't changed, or isn't viewed, don't touch it. The difference between this and last
        time, is that we won't actually do any drawing in the plot. The plot
        will be redrawn after all subplots data is changed. '''


        # Main goal, only change what is showing..


        self.xmin = self.hist2d[2][0]
        self.xmax = self.hist2d[2][-1]
        self.ymin = self.hist2d[1][0]
        self.ymax = self.hist2d[1][-1]
        self.clim = list(self.hist2d[3])


        self.cax.set_data(self.hist2d[0])

        self.cax.set_extent([self.xmin,self.xmax, self.ymin, self.ymax])



        if self.GetPlotParam('set_v_min'):
            self.clim[0] =  10**self.GetPlotParam('v_min')
        if self.GetPlotParam('set_v_max'):
            self.clim[1] =  10**self.GetPlotParam('v_max')

        self.cax.set_clim(self.clim)

        if self.GetPlotParam('show_cbar'):
            self.CbarTickFormatter()
        if self.GetPlotParam('show_shock'):
            self.shock_line.set_xdata([self.parent.shock_loc,self.parent.shock_loc])

        self.UpdateLabelsandColors()
        self.axes.set_xlabel(self.x_label, labelpad = self.parent.MainParamDict['xLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])
        self.axes.set_ylabel(self.y_label, labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])

        if self.GetPlotParam('set_y_min'):
            self.ymin = self.GetPlotParam('y_min')
        if self.GetPlotParam('set_y_max'):
            self.ymax = self.GetPlotParam('y_max')
        self.axes.set_ylim(self.ymin, self.ymax)

        if self.parent.MainParamDict['SetxLim'] and self.parent.MainParamDict['LinkSpatial'] == 1:
            if self.parent.MainParamDict['xLimsRelative']:
                self.axes.set_xlim(self.parent.MainParamDict['xLeft'] + self.parent.shock_loc,
                                   self.parent.MainParamDict['xRight'] + self.parent.shock_loc)
            else:
                self.axes.set_xlim(self.parent.MainParamDict['xLeft'], self.parent.MainParamDict['xRight'])
        else:
            self.axes.set_xlim(self.xmin,self.xmax)

    def CbarTickFormatter(self):
        ''' A helper function that sets the cbar ticks & labels. This used to be
        easier, but because I am no longer using the colorbar class i have to do
        stuff manually.'''
        clim = np.copy(self.cax.get_clim())
        if self.GetPlotParam('show_cbar'):
            if self.GetPlotParam('cnorm_type') == "Log":
                if self.parent.MainParamDict['HorizontalCbars']:
                    self.cbar.set_extent([np.log10(clim[0]),np.log10(clim[1]),0,1])
                    self.axC.set_xlim(np.log10(clim[0]),np.log10(clim[1]))
                    self.axC.xaxis.set_label_position("top")
                    if self.GetPlotParam('prtl_type') ==0:
                        self.axC.set_xlabel(r'$\log{\ \ f_i(p)}$', size = self.parent.MainParamDict['AxLabelSize'])
                    else:
                        self.axC.set_xlabel(r'$\log{\ \ f_e(p)}$', size = self.parent.MainParamDict['AxLabelSize'])

                else:
                    self.cbar.set_extent([0,1,np.log10(clim[0]),np.log10(clim[1])])
                    self.axC.set_ylim(np.log10(clim[0]),np.log10(clim[1]))
                    self.axC.locator_params(axis='y',nbins=6)
                    self.axC.yaxis.set_label_position("right")
                    if self.GetPlotParam('prtl_type') ==0:
                        self.axC.set_ylabel(r'$\log{\ \ f_i(p)}$', labelpad =self.parent.MainParamDict['cbarLabelPad'], rotation = -90, size = self.parent.MainParamDict['AxLabelSize'])
                    else:
                        self.axC.set_ylabel(r'$\log{\ \ f_e(p)}$', labelpad =self.parent.MainParamDict['cbarLabelPad'], rotation = -90, size = self.parent.MainParamDict['AxLabelSize'])

            else:# self.GetPlotParam('cnorm_type') == "Linear":
                if self.parent.MainParamDict['HorizontalCbars']:
                    self.cbar.set_extent([clim[0], clim[1], 0, 1])
                    self.axC.set_xlim(clim[0], clim[1])
                    self.axC.xaxis.set_label_position("top")
                    if self.GetPlotParam('prtl_type') ==0:
                        self.axC.set_xlabel(r'$f_i(p)$', size = self.parent.MainParamDict['AxLabelSize'])
                    else:
                        self.axC.set_xlabel(r'$f_e(p)$', size = self.parent.MainParamDict['AxLabelSize'])

                else:
                    self.cbar.set_extent([0, 1, clim[0], clim[1]])
                    self.axC.set_ylim(clim[0], clim[1])
                    self.axC.locator_params(axis='y', nbins=6)
                    self.axC.yaxis.set_label_position("right")
                    if self.GetPlotParam('prtl_type') ==0:
                        self.axC.set_ylabel(r'$f_i(p)$', labelpad =self.parent.MainParamDict['cbarLabelPad'], rotation = -90, size = self.parent.MainParamDict['AxLabelSize'])
                    else:
                        self.axC.set_ylabel(r'$f_e(p)$', labelpad =self.parent.MainParamDict['cbarLabelPad'], rotation = -90, size = self.parent.MainParamDict['AxLabelSize'])


    def GetPlotParam(self, keyname):
        return self.param_dict[keyname]
