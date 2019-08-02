import os,sys, subprocess, yaml, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as cm
import new_cmaps
import numpy as np
from collections import deque
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from tristan_sim import TristanSim
from phase_panel import PhasePanel
from fields_panel import FieldsPanel
from density_panel import DensPanel
from spectra_panel import SpectralPanel
from mag_panel import BPanel
from energy_panel import EnergyPanel
from fft_panel import FFTPanel
from total_energy_panel import TotEnergyPanel
from moments_panel import MomentsPanel
from matplotlib.backends.backend_agg import FigureCanvasAgg
from PIL import Image

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['image.origin'] = 'upper'
import argparse


"""
### NEED THIS!

"""


class Oengus():
    """ We simply derive a new class of Frame as the man frame of our app"""
    def __init__(self, preset_view='Default', sim=None, name =''):
        self.sim_name = name
        self.sim = sim
        self.IseultDir = os.path.join(os.path.dirname(__file__), '..')
        self.dirname = sim.dir

        # Read Config File
        try:
            with open(os.path.join(self.IseultDir, '.iseult_configs', preset_view.strip().replace(' ', '_') +'.yml')) as f:
                self.cfgDict = yaml.safe_load(f)
        except:
            print('Cannot find/load ' +  preset_view.strip().replace(' ', '_') +'.yml in .iseult_configs. If the name of view contains whitespace,')
            print('either it must be enclosed in quotation marks or given with whitespace replaced by _.')
            print('Name is case sensitive. Reverting to Default view')
            with open(os.path.join(self.IseultDir, '.iseult_configs', 'Default.yml')) as f:
                sefl.cfgDict = yaml.safe_load(f)
        self.GenMainParamDict()

        # Create the figure

        if self.MainParamDict['HorizontalCbars']:
            self.axes_extent = self.MainParamDict['HAxesExtent']
            self.cbar_extent = self.MainParamDict['HCbarExtent']
            self.SubPlotParams = self.MainParamDict['HSubPlotParams']

        else:
            self.axes_extent = self.MainParamDict['VAxesExtent']
            self.cbar_extent = self.MainParamDict['VCbarExtent']
            self.SubPlotParams = self.MainParamDict['VSubPlotParams']
        self.figure = plt.figure(figsize = self.MainParamDict['FigSize'], dpi = self.MainParamDict['dpi'], edgecolor = 'none', facecolor = 'w')

        self.figure.subplots_adjust( **self.SubPlotParams)
        self.canvas = FigureCanvasAgg(self.figure)


        # Make the object hold the timestep info
        # Some options to set the way the spectral lines are dashed
        self.spect_plot_counter = 0
        self.dashes_options = [[],[3,1],[5,1],[1,1]]

        # divy up the figure into a bunch of subplots using GridSpec.
        self.gs0 = gridspec.GridSpec(self.MainParamDict['NumOfRows'],self.MainParamDict['NumOfCols'])

        # Create the list of all of subplot wrappers
        self.SubPlotList = [[] for i in range(self.MainParamDict['MaxRows'])]
        self.showingCPUs = False
        self.showingTotEnergy = False
        PlotTypeDict = {
            'PhasePlot': PhasePanel,
            'EnergyPlot': EnergyPanel,
            'FieldsPlot': FieldsPanel,
            'DensityPlot': DensPanel,
            'SpectraPlot': SpectralPanel,
            'MagPlots': BPanel,
            'FFTPlots': FFTPanel,
            'TotalEnergyPlot': TotEnergyPanel,
            'Moments': MomentsPanel
        }
        for i in range(self.MainParamDict['NumOfRows']):
            for j in range(self.MainParamDict['NumOfCols']):
                tmp_str = f"Chart{i}_{j}"
                if tmp_str in self.cfgDict.keys():
                    tmpchart_type = self.cfgDict[tmp_str]['ChartType']
                    self.SubPlotList[i].append(PlotTypeDict[tmpchart_type](self, (i,j), self.cfgDict[tmp_str]))
                    self.showingTotEnergy += tmpchart_type == 'TotalEnergyPlot'
                    try:
                        self.showingCPUs += self.cfgDict[tmp_str]['show_cpu_domains']
                    except KeyError:
                        pass
                else:
                    # The graph isn't specified in the config file, just set it equal to phase plots
                    self.SubPlotList[i].append(PlotTypeDict['PhasePlot'](self, (i,j), {}))



        ##
        #
        # Open TristanSim
        #
        ##


        # previous objects
        if self.showingTotEnergy:
            self.calc_total_energy()


        ### CALC THE SHOCK speed
        self.btheta = np.nan
        if 'sigma' in self.sim._h5Key2FileDict.keys():
            if 'btheta' in self.sim._h5Key2FileDict.keys():
                if np.abs(self.sim[0].sigma)!=0:
                    self.btheta = self.sim[0].btheta
        o = self.sim[0]
        nxf0 = o.by.shape[1]
        if np.isnan(self.btheta):
            self.b0 = 1.0
            self.e0 = 1.0
        else:
            # Normalize by b0
            self.bx0 = o.bx[0,-1,-10]
            self.by0 = o.by[0,-1,-10]
            self.bz0 = o.bz[0,-1,-10]
            self.b0 = np.sqrt(self.bx0**2+self.by0**2+self.bz0**2)
            self.ex0 = o.ex[0,-1,-2]
            self.ey0 = o.ey[0,-1,-2]
            self.ez0 = o.ez[0,-1,-2]
            self.e0 = np.sqrt(self.ex0**2+self.ey0**2+self.ez0**2)

        dens_arr =np.copy(self.sim[-1].dens[0,:,:])
        final_time = self.sim[-1].time
        istep = self.sim[-1].istep
        c_omp = self.sim[-1].c_omp

        # Find out where the shock is at the last time step.
        jstart = int(min(10*c_omp/istep, nxf0))
        # build the final x_axis of the plot

        xaxis_final = np.arange(dens_arr.shape[1])/c_omp*istep
        # Find the shock by seeing where the density is 1/2 of it's
        # max value.

        dens_half_max = max(dens_arr[dens_arr.shape[0]//2,jstart:])*.5

        # Find the farthest location where the average density is greater
        # than half max
        ishock_final = np.where(dens_arr[dens_arr.shape[0]//2,jstart:]>=dens_half_max)[0][-1]
        xshock_final = xaxis_final[ishock_final]
        xshock_final -= np.min(self.sim[-1].xe[self.sim[-1].xe!=0])/self.sim[-1].c_omp
        self.shock_speed = xshock_final/final_time
        self.create_graphs()
    def GenMainParamDict(self, config_file = None):
        ''' The function that reads in a config file and then makes MainParamDict to hold all of the main iseult parameters.
            It also sets all of the plots parameters.'''

        #config = configparser.RawConfigParser()



        # Since configparser reads in strings we have to format the data.
        # First create MainParamDict with the default parameters,
        # the dictionary that will hold the parameters for the program.
        # See ./iseult_configs/Default.cfg for a description of what each parameter does.
        self.MainParamDict = {'zSlice': 0.0, # THIS IS A float WHICH IS THE RELATIVE POSITION OF THE 2D SLICE 0->1
                              '2DSlicePlane': 0, # 0 = x-y plane, 1 == x-z plane
                              'Average1D': 0,
                              'ySlice': 0.5, # THIS IS A FLOAT WHICH IS THE RELATIVE POSITION OF THE 1D SLICE 0->1
                              'WindowSize': '1200x700',
                              'yTop': 100.0,
                              'yBottom': 0.0,
                              'Reload2End': True,
                              'ColorMap': 'viridis',
                              'FFTLeft': 0.0,
                              'ShowTitle': True,
                              'ImageAspect': 0,
                              'WaitTime': 0.01,
                              'MaxCols': 8,
                              'VAxesExtent': [4, 90, 0, 92],
                              'kRight': 1.0,
                              'DoLorentzBoost': False,
                              'NumOfRows': 3,
                              'MaxRows': 8,
                              'SetkLim': False,
                              'VCbarExtent': [4, 90, 94, 97],
                              'SkipSize': 5,
                              'xLeft': 0.0,
                              'NumFontSize': 11,
                              'AxLabelSize': 11,
                              'FFTRelative': True,
                              'NumOfCols': 2,
                              'VSubPlotParams': {'right': 0.95,
                                                 'bottom': 0.06,
                                                 'top': 0.93,
                                                 'wspace': 0.23,
                                                 'hspace': 0.15,
                                                 'left': 0.06},
                              'HAxesExtent': [18, 92, 0, -1],
                              'SetyLim': False,
                              'HSubPlotParams': {'right': 0.95,
                                                 'bottom': 0.06,
                                                 'top': 0.91,
                                                 'wspace': 0.15,
                                                 'hspace': 0.3,
                                                 'left': 0.06},
                              'yLabelPad': 0,
                              'cbarLabelPad': 15,
                              'SetxLim': False,
                              'xLimsRelative': False,
                              'ConstantShockVel': True,
                              'xRight': 100.0,
                              'LinkSpatial': 2,
                              'HCbarExtent': [0, 4, 0, -1],
                              'Recording': False,
                              'xLabelPad': 0,
                              'annotateTextSize': 18,
                              'FFTRight': 200.0,
                              'ClearFig': True,
                              'HorizontalCbars': False,
                              'DivColorMap': 'BuYlRd',
                              'LinkK': True,
                              'GammaBoost': 0.0,
                              'kLeft': 0.1,
                              'LoopPlayback': True,
                              'PrtlStride': 5,
                              'electron_color': '#fca636',
                              'electron_fit_color': 'yellow',
                              'ion_color': '#d6556d',
                              'ion_fit_color': 'r',
                              'shock_color': 'w',
                              'FigSize':  [12.0, 6.22],
                              'dpi': 100,
                              'FFT_color': 'k',
                              'legendLabelSize':11}
        for key, val in self.cfgDict['MainParamDict'].items():
            self.MainParamDict[key] = val
        self.electron_color = self.MainParamDict['electron_color']
        self.ion_color = self.MainParamDict['ion_color']
        self.shock_color = self.MainParamDict['shock_color']
        self.ion_fit_color = self.MainParamDict['ion_fit_color']
        self.electron_fit_color = self.MainParamDict['electron_fit_color']
        self.FFT_color = self.MainParamDict['FFT_color']


    def calc_total_energy(self):
        self.TotalEnergyTimes = []
        self.TotalElectronEnergy = []
        self.TotalIonEnergy = []#
        self.TotalMagEnergy = []
        self.TotalElectricEnergy =[]
        self.TotalBzEnergy = []
        for o in self.sim:
            self.TotalEnergyTimes.append(o.time)
            self.TotalElectronEnergy.append(np.sum(np.sqrt(o.ue*o.ue + o.ve*o.ve + o.we*o.we +1)-1)*o.stride*abs(o.qi)*o.c**2)
            self.TotalIonEnergy.append(np.sum(np.sqrt(o.ui*o.ui + o.vi*o.vi + o.wi*o.wi +1)-1)*o.stride*abs(o.qi)*o.mi/o.me*o.c**2)
            self.TotalMagEnergy.append(np.sum(o.bz*o.bz+o.bx*o.bx+o.by*o.by)*o.istep**2*.5)
            self.TotalElectricEnergy.append(np.sum(o.ez*o.ez+o.ex*o.ex+o.ey*o.ey)*o.istep**2*.5)
            self.TotalBzEnergy.append(np.sum(o.bz*o.bz)*o.istep**2*.5)
            o.clear()
        self.TotalElectronEnergy = np.array(self.TotalElectronEnergy)
        self.TotalIonEnergy = np.array(self.TotalIonEnergy)
        self.TotalMagEnergy = np.array(self.TotalMagEnergy)
        self.TotalElectricEnergy = np.array(self.TotalElectricEnergy)
        self.TotalBzEnergy = np.array(self.TotalBzEnergy)

    def create_graphs(self):
        o = self.sim[0]
        # FIND THE SLICE
        self.MaxZInd = o.bx.shape[0]-1
        self.MaxYInd = o.bx.shape[1]-1
        self.MaxXInd = o.bx.shape[2]-1

        self.ySlice = int(np.around(self.MainParamDict['ySlice']*self.MaxYInd))
        self.zSlice = int(np.around(self.MainParamDict['zSlice']*self.MaxZInd))
        if self.MainParamDict['ConstantShockVel']:
            self.shock_loc = o.time*self.shock_speed
            self.shock_loc += np.min(o.xe[o.xe!=0])/o.c_omp
        else:
            jstart = int(min(10*o.c_omp/o.istep, o.dens[0,:,:].shape[1]))
            cur_xaxis = np.arange(o.dens[0,:,:].shape[1])/o.c_omp*o.istep
            # Find the shock by seeing where the density is 1/2 of it's
            # max value.

            dens_half_max = max(o.dens[0,:,:][o.dens[0,:,:].shape[0]//2,jstart:])*.5

            # Find the farthest location where the average density is greater
            # than half max
            ishock = np.where(o.dens[0,:,:][o.dens[0,:,:].shape[0]/2,jstart:]>=dens_half_max)[0][-1]
            self.shock_loc = cur_xaxis[ishock]

            #self.cpu_x_locs = np.cumsum(self.DataDict['mx']-5)/self.DataDict['c_omp'][0]
            #self.cpu_y_locs = np.cumsum(self.DataDict['my']-5)/self.DataDict['c_omp'][0]

        # Now that the DataDict is created, iterate over all the subplots and
        # load the data into them:
        for i in range(self.MainParamDict['NumOfRows']):
            for j in range(self.MainParamDict['NumOfCols']):
                self.SubPlotList[i][j].update_data(o)

        for i in range(self.MainParamDict['NumOfRows']):
            for j in range(self.MainParamDict['NumOfCols']):
                self.SubPlotList[i][j].draw()
        if self.showingCPUs:
            if 'my' in self.sim._h5Key2FileDict.keys():
                cpu_y_locs = np.cumsum(o.my-5)/o.c_omp
            else:
                tmpSize = ((self.MaxYInd+1)*o.istep)//(o.my0-5)
                cpu_y_locs = np.cumsum(np.ones(tmpSize)*(o.my0)-5)/o.c_omp
            if 'mx' in self.sim._h5Key2FileDict.keys():
                cpu_x_locs = np.cumsum(o.mx-5)/o.c_omp
            else:
                tmpSize = ((self.MaxXInd+1)*o.istep)//(o.mx0-5)
                cpu_x_locs = np.cumsum(np.ones(tmpSize)*(o.mx0)-5)/o.c_omp


            for i in range(self.MainParamDict['NumOfRows']):
                for j in range(self.MainParamDict['NumOfCols']):
                    try:
                        if self.SubPlotList[i][j].param_dict['show_cpu_domains']:
                            for k in range(len(self.parent.cpu_x_locs)):
                                self.SubPlotList[i][j].axes.axvline(cpu_x_locs[k], linewidth = 1, linestyle = ':',color = 'w')
                            for k in range(len(self.parent.cpu_y_locs)):
                                self.SubPlotList[i][j].axes.axvline(cpu_y_locs[k], linewidth = 1, linestyle = ':',color = 'w')

                    except KeyError:
                        pass

        if self.MainParamDict['ShowTitle']:
            if len(self.sim_name) == 0:
                self.figure.suptitle(os.path.abspath(self.dirname)+ '/*.'+o.fnum+' at time t = %d $\omega_{pe}^{-1}$'  % round(o.time), size = 15)
            else:
                self.figure.suptitle(self.sim_name +', t = %d $\omega_{pe}^{-1}$'  % round(o.time), size = 15)
        ####
        #
        # Write the lines to the phase plots
        #
        ####

        # first find all the phase plots that need writing to
        self.phase_plot_list = []
        self.spectral_plot_list = []

        for i in range(self.MainParamDict['NumOfRows']):
            for j in range(self.MainParamDict['NumOfCols']):
                if self.SubPlotList[i][j].chartType =='PhasePlot' or self.SubPlotList[i][j].chartType =='EnergyPlot':
                    if self.SubPlotList[i][j].GetPlotParam('show_int_region'):
                        self.phase_plot_list.append([i,j])
                if self.SubPlotList[i][j].chartType =='SpectraPlot':
                    self.spectral_plot_list.append([i,j])

        for pos in self.phase_plot_list:
            if self.SubPlotList[pos[0]][pos[1]].GetPlotParam('prtl_type') == 0:
                for spos in self.spectral_plot_list:
                    if self.SubPlotList[spos[0]][spos[1]].GetPlotParam('show_ions'):
                        k = min(self.SubPlotList[spos[0]][spos[1]].spect_num, len(self.dashes_options)-1)
                        # Append the left line to the list
                        self.SubPlotList[pos[0]][pos[1]].IntRegionLines.append(self.SubPlotList[pos[0]][pos[1]].axes.axvline(
                        max(self.SubPlotList[spos[0]][spos[1]].i_left_loc, self.SubPlotList[pos[0]][pos[1]].xmin+1),
                        linewidth = 1.5, linestyle = '-', color = self.ion_color))
                        # Choose the left dashes pattern
                        self.SubPlotList[pos[0]][pos[1]].IntRegionLines[-1].set_dashes(self.dashes_options[k])

                        # Append the right line to the list
                        self.SubPlotList[pos[0]][pos[1]].IntRegionLines.append(self.SubPlotList[pos[0]][pos[1]].axes.axvline(
                        min(self.SubPlotList[spos[0]][spos[1]].i_right_loc, self.SubPlotList[pos[0]][pos[1]].xmax+1),
                        linewidth = 1.5, linestyle = '-', color = self.ion_color))
                        # Choose the right dashes pattern
                        self.SubPlotList[pos[0]][pos[1]].IntRegionLines[-1].set_dashes(self.dashes_options[k])
            else:
                for spos in self.spectral_plot_list:
                    if self.SubPlotList[spos[0]][spos[1]].GetPlotParam('show_electrons'):
                        k = min(self.SubPlotList[spos[0]][spos[1]].spect_num, len(self.dashes_options)-1)
                        # Append the left line to the list
                        self.SubPlotList[pos[0]][pos[1]].IntRegionLines.append(self.SubPlotList[pos[0]][pos[1]].axes.axvline(
                        max(self.SubPlotList[spos[0]][spos[1]].e_left_loc, self.SubPlotList[pos[0]][pos[1]].xmin+1),
                        linewidth = 1.5, linestyle = '-', color = self.electron_color))
                        # Choose the left dashes pattern
                        self.SubPlotList[pos[0]][pos[1]].IntRegionLines[-1].set_dashes(self.dashes_options[k])

                        # Append the right line to the list
                        self.SubPlotList[pos[0]][pos[1]].IntRegionLines.append(self.SubPlotList[pos[0]][pos[1]].axes.axvline(
                        min(self.SubPlotList[spos[0]][spos[1]].e_right_loc, self.SubPlotList[pos[0]][pos[1]].xmax+1),
                        linewidth = 1.5, linestyle = '-', color = self.electron_color))
                        # Choose the right dashes pattern
                        self.SubPlotList[pos[0]][pos[1]].IntRegionLines[-1].set_dashes(self.dashes_options[k])
        self.canvas.draw()

    def draw_output(self, n):
        o = self.sim[n]
        # FIND THE SLICE
        self.MaxZInd = o.bx.shape[0]-1
        self.MaxYInd = o.bx.shape[1]-1
        self.MaxXInd = o.bx.shape[2]-1

        self.ySlice = int(np.around(self.MainParamDict['ySlice']*self.MaxYInd))
        self.zSlice = int(np.around(self.MainParamDict['zSlice']*self.MaxZInd))
        if self.MainParamDict['ConstantShockVel']:
            self.shock_loc = o.time*self.shock_speed
            self.shock_loc += np.min(o.xe[o.xe!=0])/o.c_omp
        else:
            jstart = int(min(10*o.c_omp/o.istep, o.dens[0,:,:].shape[1]))
            cur_xaxis = np.arange(o.dens[0,:,:].shape[1])/o.c_omp*o.istep
            # Find the shock by seeing where the density is 1/2 of it's
            # max value.

            dens_half_max = max(o.dens[0,:,:][o.dens[0,:,:].shape[0]//2,jstart:])*.5

            # Find the farthest location where the average density is greater
            # than half max
            ishock = np.where(o.dens[0,:,:][o.dens[0,:,:].shape[0]/2,jstart:]>=dens_half_max)[0][-1]
            self.shock_loc = cur_xaxis[ishock]

            #self.cpu_x_locs = np.cumsum(self.DataDict['mx']-5)/self.DataDict['c_omp'][0]
            #self.cpu_y_locs = np.cumsum(self.DataDict['my']-5)/self.DataDict['c_omp'][0]

        # Now that the DataDict is created, iterate over all the subplots and
        # load the data into them:
        for i in range(self.MainParamDict['NumOfRows']):
            for j in range(self.MainParamDict['NumOfCols']):
                self.SubPlotList[i][j].update_data(o)

        for i in range(self.MainParamDict['NumOfRows']):
            for j in range(self.MainParamDict['NumOfCols']):
                self.SubPlotList[i][j].refresh()
        if self.showingCPUs:
            if 'my' in self.sim._h5Key2FileDict.keys():
                cpu_y_locs = np.cumsum(o.my-5)/o.c_omp
            else:
                tmpSize = ((self.MaxYInd+1)*o.istep)//(o.my0-5)
                cpu_y_locs = np.cumsum(np.ones(tmpSize)*(o.my0)-5)/o.c_omp
            if 'mx' in self.sim._h5Key2FileDict.keys():
                cpu_x_locs = np.cumsum(o.mx-5)/o.c_omp
            else:
                tmpSize = ((self.MaxXInd+1)*o.istep)//(o.mx0-5)
                cpu_x_locs = np.cumsum(np.ones(tmpSize)*(o.mx0)-5)/o.c_omp


            for i in range(self.MainParamDict['NumOfRows']):
                for j in range(self.MainParamDict['NumOfCols']):
                    try:
                        if self.SubPlotList[i][j].param_dict['show_cpu_domains']:
                            for k in range(len(self.parent.cpu_x_locs)):
                                self.SubPlotList[i][j].axes.axvline(cpu_x_locs[k], linewidth = 1, linestyle = ':',color = 'w')
                            for k in range(len(self.parent.cpu_y_locs)):
                                self.SubPlotList[i][j].axes.axvline(cpu_y_locs[k], linewidth = 1, linestyle = ':',color = 'w')

                    except KeyError:
                        pass

        if self.MainParamDict['ShowTitle']:
            if len(self.sim_name) == 0:
                self.figure.suptitle(os.path.abspath(self.dirname)+ '/*.'+o.fnum+' at time t = %d $\omega_{pe}^{-1}$'  % round(o.time), size = 15)
            else:
                self.figure.suptitle(self.sim_name +', t = %d $\omega_{pe}^{-1}$'  % round(o.time), size = 15)
        ####
        #
        # Write the lines to the phase plots
        #
        ####

        # first find all the phase plots that need writing to
        self.phase_plot_list = []
        self.spectral_plot_list = []

        for i in range(self.MainParamDict['NumOfRows']):
            for j in range(self.MainParamDict['NumOfCols']):
                if self.SubPlotList[i][j].chartType =='PhasePlot' or self.SubPlotList[i][j].chartType =='EnergyPlot':
                    if self.SubPlotList[i][j].GetPlotParam('show_int_region'):
                        self.phase_plot_list.append([i,j])
                if self.SubPlotList[i][j].chartType =='SpectraPlot':
                    self.spectral_plot_list.append([i,j])

        for pos in self.phase_plot_list:
            if self.SubPlotList[pos[0]][pos[1]].GetPlotParam('prtl_type') == 0:
                for spos in self.spectral_plot_list:
                    if self.SubPlotList[spos[0]][spos[1]].GetPlotParam('show_ions'):
                        k = min(self.SubPlotList[spos[0]][spos[1]].spect_num, len(self.dashes_options)-1)
                        # Append the left line to the list
                        self.SubPlotList[pos[0]][pos[1]].IntRegionLines.append(self.SubPlotList[pos[0]][pos[1]].axes.axvline(
                        max(self.SubPlotList[spos[0]][spos[1]].i_left_loc, self.SubPlotList[pos[0]][pos[1]].xmin+1),
                        linewidth = 1.5, linestyle = '-', color = self.ion_color))
                        # Choose the left dashes pattern
                        self.SubPlotList[pos[0]][pos[1]].IntRegionLines[-1].set_dashes(self.dashes_options[k])

                        # Append the right line to the list
                        self.SubPlotList[pos[0]][pos[1]].IntRegionLines.append(self.SubPlotList[pos[0]][pos[1]].axes.axvline(
                        min(self.SubPlotList[spos[0]][spos[1]].i_right_loc, self.SubPlotList[pos[0]][pos[1]].xmax+1),
                        linewidth = 1.5, linestyle = '-', color = self.ion_color))
                        # Choose the right dashes pattern
                        self.SubPlotList[pos[0]][pos[1]].IntRegionLines[-1].set_dashes(self.dashes_options[k])
            else:
                for spos in self.spectral_plot_list:
                    if self.SubPlotList[spos[0]][spos[1]].GetPlotParam('show_electrons'):
                        k = min(self.SubPlotList[spos[0]][spos[1]].spect_num, len(self.dashes_options)-1)
                        # Append the left line to the list
                        self.SubPlotList[pos[0]][pos[1]].IntRegionLines.append(self.SubPlotList[pos[0]][pos[1]].axes.axvline(
                        max(self.SubPlotList[spos[0]][spos[1]].e_left_loc, self.SubPlotList[pos[0]][pos[1]].xmin+1),
                        linewidth = 1.5, linestyle = '-', color = self.electron_color))
                        # Choose the left dashes pattern
                        self.SubPlotList[pos[0]][pos[1]].IntRegionLines[-1].set_dashes(self.dashes_options[k])

                        # Append the right line to the list
                        self.SubPlotList[pos[0]][pos[1]].IntRegionLines.append(self.SubPlotList[pos[0]][pos[1]].axes.axvline(
                        min(self.SubPlotList[spos[0]][spos[1]].e_right_loc, self.SubPlotList[pos[0]][pos[1]].xmax+1),
                        linewidth = 1.5, linestyle = '-', color = self.electron_color))
                        # Choose the right dashes pattern
                        self.SubPlotList[pos[0]][pos[1]].IntRegionLines[-1].set_dashes(self.dashes_options[k])

        self.canvas.draw()

        s, (width, height) = self.canvas.print_to_buffer()
        return Image.frombytes('RGBA', (width, height), s)
def runMe(cmd_args):
    tic = time.time()
    cmdout = ['ffmpeg',
                '-y', '-f', 'image2pipe', # overwrite, image2 is a colorspace thing.
                '-framerate', str(int(cmd_args.framerate)), # Set framerate to the the user selected option
                 '-pattern_type', 'glob', '-i', '-', # Not sure what this does... I am going to get rid of it
                 '-codec', 'copy',  # save as a *.mov
                 cmd_args.outmovie]#, '&']#, # output name,

    pipe = subprocess.Popen(cmdout, stdin=subprocess.PIPE)
    sims = []
    iseult_figs = []

    for i in range(len(cmd_args.O)):
        dirname= os.curdir
        dirlist = os.listdir(dirname)
        if len(cmd_args.O[i])>0:
            dirname = os.path.join(dirname, cmd_args.O[i])
        elif 'output' in dirlist:
            dirname = os.path.join(dirname, 'output')
        curname = ''
        if i<len(cmd_args.name):
            curname = cmd_args.name[i]
        curSim = TristanSim(dirname)

        cntxt = {'preset_view':cmd_args.p,
            'sim': curSim,
            'name':curname
            }
        sims.append(curSim)
        iseult_figs.append(Oengus(**cntxt))
    for s in sims:
        s.tlist = np.array([o.time for o in s])
    tSteps= []
    for s in sims:
        if len(tSteps) < len(s.tlist):
            # list(s.tlist) instead of s.tlist here is to force a deep copy. quirk of python.
            tSteps = list(s.tlist)


    for t in tSteps:
        imgs = []
        for s, ifig in zip(sims, iseult_figs):
            n = np.where(np.min(np.abs(s.tlist-t)) == np.abs(s.tlist-t))[0][0]
            imgs.append(ifig.draw_output(n))

        imgs_comb = np.vstack(list(np.asarray(i) for i in imgs))

        # save that beautiful picture
        imgs_comb = Image.fromarray(imgs_comb)
        imgs_comb.save(pipe.stdin, 'PNG')
        print(f"saving image {n} to pipe")
    pipe.stdin.close()
    pipe.wait()

    # Make sure all went well
    if pipe.returncode != 0:
        raise sp.CalledProcessError(pipe.returncode, cmd_out)
"""
# list of tuples containing directories we want to analyze, and the
# name of the simulations
dirList = []

sims = []
for outdir, name in dirList:
    tmp = TristanSim(outdir)
    tmp.name = name
    sims.append(tmp)

for s in sims:
    s.tlist = np.array([o.time for o in s])


# Calculate the shock speed for each sim
for s in sims:
    ysize = s[-1].dens.shape[1]
    dens1D = s[-1].dens[0, ysize//2,:]
    ishock = -1
    while dens1D[ishock] < 2*s[-1].ppc0:
        ishock -= 1

    xshock = (len(dens1D)+ishock)*s[-1].istep/s[-1].c_omp

    # Measure relative to the wall
    xshock -= np.min(s[-1].xe[s[-1].xe!=0])/s[-1].c_omp

    s.vshock = xshock/s[-1].time

tSteps= []
for s in sims:
    if len(tSteps) < len(s.tlist):
        # list(s.tlist) instead of s.tlist here is to force a deep copy. quirk of python.
        tSteps = list(s.tlist)



###
#
# CREATE THE SCRATCH DIR WHERE WE WILL SAVE OUR IMAGES
#
###
try:
    os.makedirs('tmp_erase_me_plz')

except:
    raise Exception('Cannot create directory tmp_erase_me_plz. Either it already exists or you do not have write priviledges')


fig = plt.figure(figsize=(11, 7)) # creates a fig of (width, height)

xmin = -1500 # relative to xshock
xmax = 3000
rowNum = len(sims)
k = 1
e_color = '#D93E30'
i_color = '#01261F'
for t in tSteps:
    print(k)
    plt.clf() # clear the figure
    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=.95, wspace=0.3, hspace=0.5)
    rows = []
    for i in range(rowNum):
        row = []
        for j in range(3):
            row.append(plt.subplot(rowNum, 3, i*(3)+j+1))
        rows.append(row)

    # Now iterate over sims and rows
    for row, s in zip(rows, sims):
        # get the right output time
        n = np.where(np.min(np.abs(s.tlist-t)) == np.abs(s.tlist-t))[0][0]
        o = s[n]

        # convert our relative lims to current xlims
        xshock = o.time*s.vshock + np.min(o.xe[o.xe!=0])/o.c_omp

        # Make the phase diagram
        hist2D(o.xe/o.c_omp, o.ue, bins=[400,200], # [xbins, ybins]
               ax = row[0], normhist=False, vmin=1, vmax =2E4, cmap = BlGrRdYl, clabel=r'$f_i(p)$')
        row[0].axvline(xshock -1500, color=e_color)
        row[0].axvline(xshock -500, color=e_color)
        row[0].set_ylim(-2, 12)
        row[0].set_ylabel(r'$\gamma_e\beta_e$')

        # Make the line plot
        ymid = o.by.shape[1]//2
        x_arr = np.arange(o.by.shape[2])*o.istep/o.c_omp
        row[1].plot(x_arr, o.by[0, 0,:], label = r'$B_y$')
        row[1].plot(x_arr, o.bz[0, 0,:], label = r'$B_z$')
        row[1].legend(loc='upper right')
        row[1].set_ylim(-.001, .0058)
        row[1].set_title(s.name + r", $\omega_{pe}t=$"f" {int(o.time)}")
        row[1].set_ylabel(r'$B$')
        tristanSpect(o, species='ion', ax =row[2], xLeft = xshock -1500, xRight = xshock-500, color = i_color)
        tristanSpect(o, species='lec', ax =row[2], xLeft = xshock -1500, xRight = xshock-500, color = e_color)
        row[2].set_yscale('log')
        row[2].set_xscale('log')
        row[2].set_ylim(10,8E5)
        row[2].set_xlim(1E-6, 2E2)
        for ax in row[0:2]:
            ax.set_xlim(xshock+xmin, xshock+xmax)
            ax.set_xlabel(r'$x\ [c/\omega_{pe}]$')
        # You don't have to do this, but it is nice to free up some memory here
        o.clear()

    #plt.tight_layout()
    plt.savefig(os.path.join('tmp_erase_me_plz', 'out'+f"{k}".zfill(3) + '.png'))
    k+=1


#cmdstring = ['xterm', '-e','ffmpeg',
#             '-y', '-f', 'image2', # overwrite, image2 is a colorspace thing.
#             '-framerate', str(int(FPS)), # Set framerate to the the user selected option
#             '-pattern_type', 'glob', '-i', os.path.join('tmp_erase_me_plz','*.png'), # Not sure what this does... I am going to get rid of it
#             '-codec', 'copy',  # save as a *.mov
#             MOVIE_NAME]#, '&']#, # output name,
#subprocess.call(cmdstring)


"""
"""
for name in os.listdir(os.path.join('tmp_erase_me_plz')):
    os.remove(os.path.join('tmp_erase_me_plz', name))
os.rmdir('tmp_erase_me_plz')

if __name__ == "__main__":
    print('hi')
"""
