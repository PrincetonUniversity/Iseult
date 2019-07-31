#!/usr/bin/env python
import matplotlib
import numpy as np
import numpy.ma as ma
import sys
sys.path.append('../')
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
                       'spatial_y': False,
                       'face_color': 'gainsboro'}

    # We need the types of all the parameters for the config file

    def __init__(self, parent, pos, param_dict):
        self.param_dict = {}
        for key, val in self.plot_param_dict.items():
            self.param_dict[key] = val
        for key, val in param_dict.items():
            self.param_dict[key] = val
        self.pos = pos
        self.parent = parent
        self.chartType = 'FFTPlots'
        self.figure = self.parent.figure
        self.InterpolationMethods = ['none','nearest', 'bilinear', 'bicubic', 'spline16',
            'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric',
            'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']
        self.ylabel_list = [r'$|\delta B_z(k)|/B_0$', r'$|\delta B_y(k)|/B_0$', r'$|E_x(k)|/E_0$', r'$\chi(k)\ [\ \! ^\circ]$']



    def LimFinder(self, arr):
        oneD_lims = [arr.min(), arr.max()]
        # now give it a bit of spacing, a 4% percent difference of the distance
        dist = oneD_lims[1]-oneD_lims[0]
        oneD_lims[0] -= 0.04*dist
        oneD_lims[1] += 0.04*dist
        return oneD_lims

    def draw(self, output):

        ''' A function that draws the data. In the interest in speeding up the
        code, draw should only be called when you want to recreate the whole
        figure, i.e. it  will be slow. Most times you will only want to update
        what has changed in the figure. This will be done in a function called
        refresh, that should be much much faster.'''
        ''' A Helper function that loads the data for the plot'''
        # Load c_omp, istep, and the line color
        self.c_omp = getattr(output,'c_omp')
        self.istep = getattr(output,'istep')
        self.line_color = new_cmaps.cmaps[self.parent.MainParamDict['ColorMap']](0.5)

        # Load the x-axis
        # x-values haven't been calculated yet, generate them then save them to the dictionary for later.
        self.xaxis_values = np.arange(getattr(output,'bx')[0,:,:].shape[1])/self.c_omp*self.istep

        # Load the left & right locations
        self.left_loc = self.parent.MainParamDict['FFTLeft'] + self.parent.shock_loc*self.parent.MainParamDict['FFTRelative']
        self.left_loc = max(self.left_loc, self.xaxis_values[0])

        self.right_loc = self.parent.MainParamDict['FFTRight'] + self.parent.shock_loc*self.parent.MainParamDict['FFTRelative']
        self.right_loc = min(self.right_loc, self.xaxis_values[-1])


        # A list that will make sure that the data has the same int region
        self.region_args = list([self.left_loc, self.right_loc])
        # Check if the region is the same in the DataDict

        ####
        # First find the left & right indices
        ####
        iL = self.xaxis_values.searchsorted(self.left_loc)
        iR = self.xaxis_values.searchsorted(self.right_loc)
        self.all_min_max = []
        # Calculate K_axis
        self.k_axis = np.arange(iR-iL)*(2*np.pi/(self.xaxis_values[1]-self.xaxis_values[0]))/(iR-iL)-(2*np.pi/(self.xaxis_values[1]-self.xaxis_values[0]))*.5
        # Calculate all of the FFTs, just simpler to do it this way...
        bz = getattr(output,'bz')[0,:,:]

        self.BzFFT = np.fft.fft(bz[self.parent.ySlice,iL:iR]*self.parent.b0**(-1.0))
        # Shift the fft so it is centered
        self.BzFFT = np.fft.fftshift(self.BzFFT)
        self.all_min_max.append(self.LimFinder(np.abs(self.BzFFT)))

        bx = getattr(output,'bx')[0,:,:]
        by = getattr(output,'by')[0,:,:]
        b_perp_in_plane = by*self.parent.b0**(-1.0)
        self.BperpFFT = np.fft.fft(b_perp_in_plane[self.parent.ySlice,iL:iR])
        # Shift the fft so it is centered
        self.BperpFFT = np.fft.fftshift(self.BperpFFT)

        self.all_min_max.append(self.LimFinder(np.abs(self.BperpFFT)))
        ex = getattr(output,'ex')[0,:,:]
        self.ExFFT = np.fft.fft(ex[self.parent.ySlice,iL:iR]*self.parent.e0**(-1.0))
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

        self.klims = [[0, 5*kmaxval], [kmaxval/5., 50*kmaxval]]

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

        # Set the tick color
        tick_color = 'black'

        # Create a gridspec to handle spacing better
        self.gs = gridspec.GridSpecFromSubplotSpec(100,100, subplot_spec = self.parent.gs0[self.pos])


        # See if there are other k-axes, and share with them
        self.axes = self.figure.add_subplot(self.gs[self.parent.axes_extent[0]:self.parent.axes_extent[1], self.parent.axes_extent[2]:self.parent.axes_extent[3]])


        self.line = self.axes.plot(self.k_axis, self.y, color = self.line_color)

        if int(matplotlib.__version__[0]) < 2:
            self.axes.set_axis_bgcolor(self.GetPlotParam('face_color'))
        else:
            self.axes.set_facecolor(self.GetPlotParam('face_color'))

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
            self.axes.set_ylim(bottom = self.GetPlotParam('y_min'))
        if self.GetPlotParam('set_y_max'):
            self.axes.set_ylim(top = self.GetPlotParam('y_max'))

        self.axes.set_xlabel(r'$k \ [\omega_{pe}/c]$', labelpad = self.parent.MainParamDict['xLabelPad'], color = 'black', size  = self.parent.MainParamDict['AxLabelSize'])
        self.axes.set_ylabel(self.ylabel, labelpad = self.parent.MainParamDict['yLabelPad'], color = 'black', size = self.parent.MainParamDict['AxLabelSize'])

    def GetPlotParam(self, keyname):
        return self.param_dict[keyname]
