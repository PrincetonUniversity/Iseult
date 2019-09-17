#! /usr/bin/env python
import os, sys # Used to make the code portable
import yaml, configparser
import numpy as np
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))
import new_cmaps
def convertOldConfig(config_file):
    ''' The function that reads in a config file and then makes MainParamDict to hold all of the main iseult parameters.
        It also sets all of the plots parameters.'''
    cfgDict = {}
    config = configparser.RawConfigParser()

    config.read(config_file)
    cfgDict['general'] =  {'ConfigName': config.get('general', 'ConfigName')}
    ###
    #
    # FIRST LOAD THE  main param dictionary
    #
    ###

    # Since configparser reads in strings we have to format the data.
    # First create MainParamDict with the default parameters,
    # the dictionary that will hold the parameters for the program.
    # See ./iseult_configs/Default.cfg for a description of what each parameter does.
    cfgDict['MainParamDict'] = {'zSlice': 0.0, # THIS IS A float WHICH IS THE RELATIVE POSITION OF THE 2D SLICE 0->1
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
                          'FFT_color': 'k',
                          'legendLabelSize':11}

    # The list of things that should be formatted as booleans.
    BoolList = ['Reload2End', 'ClearFig', 'ShowTitle', 'DoLorentzBoost',
                'HorizontalCbars', 'SetxLim', 'SetyLim', 'SetkLim', 'LinkK',
                'LoopPlayback', 'Recording', 'FFTRelative', 'xLimsRelative',
                'ConstantShockVel', 'Average1D']


    for elm in BoolList:
        if elm.lower() in config.options('main'):
            cfgDict['MainParamDict'][elm] = config.getboolean('main', elm)

    # The list of things that should be formatted as ints.
    IntList = ['NumFontSize', 'AxLabelSize', 'xLabelPad', 'yLabelPad', #'MaxRows','MaxCols',  No longer saving this to config files.
               'NumOfRows', 'NumOfCols', 'ImageAspect', 'legendLabelSize'
               'LinkSpatial', 'SkipSize', 'PrtlStride','cbarLabelPad', 'annotateTextSize']

    for elm in IntList:
        if elm.lower() in config.options('main'):
            cfgDict['MainParamDict'][elm] = int(config.getfloat('main', elm))


    # The list of things that should be formatted as Floats.
    FloatList = [                'GammaBoost',
                #'ElectronLeft', 'ElectronRight', 'IonLeft', 'IonRight',
                'FFTLeft', 'FFTRight',
                'xLeft', 'xRight', 'yBottom', 'yTop', 'kLeft', 'kRight',
                'WaitTime', 'ySlice', 'zSlice']

    for elm in FloatList:
        if elm.lower() in config.options('main'):
            cfgDict['MainParamDict'][elm] = config.getfloat('main', elm)

    StrList = ['ColorMap', 'DivColorMap', 'WindowSize']

    for elm in StrList:
        if elm.lower() in config.options('main'):
            cfgDict['MainParamDict'][elm] = config.get('main', elm)

    # Some special parsing rules. First some lists
    IntListsList = ['HAxesExtent', 'HCbarExtent', 'VAxesExtent', 'VCbarExtent']
    for elm in IntListsList:
        if elm.lower() in config.options('main'):
            tmplist = config.get('main', elm).split(',')
            for i in range(len(tmplist)):
                tmplist[i] = int(tmplist[i])
            cfgDict['MainParamDict'][elm] = list(tmplist)

    # Now some dicts
    DictList = ['HSubPlotParams', 'VSubPlotParams']
    for elm in DictList:
        if elm.lower() in config.options('main'):
            tmplist = config.get('main', elm).split(',')
            cfgDict['MainParamDict'][elm]['left'] = float(tmplist[0])
            cfgDict['MainParamDict'][elm]['right'] = float(tmplist[1])
            cfgDict['MainParamDict'][elm]['top'] = float(tmplist[2])
            cfgDict['MainParamDict'][elm]['bottom'] = float(tmplist[3])
            cfgDict['MainParamDict'][elm]['wspace'] = float(tmplist[4])
            cfgDict['MainParamDict'][elm]['hspace'] = float(tmplist[5])
    cmaps_with_green = ['viridis', 'Rainbow + White', 'Blue/Green/Red/Yellow', 'Cube YF', 'Linear_L']
    if cfgDict['MainParamDict']['ColorMap'] in cmaps_with_green:
        cfgDict['MainParamDict']['ion_color'] = "#{0:02x}{1:02x}{2:02x}".format(int(np.round(new_cmaps.cmaps['plasma'](0.55)[0]*255)), int(np.round(new_cmaps.cmaps['plasma'](0.55)[1]*255)), int(np.round(new_cmaps.cmaps['plasma'](0.55)[2]*255)))
        cfgDict['MainParamDict']['electron_color'] ="#{0:02x}{1:02x}{2:02x}".format(int(np.round(new_cmaps.cmaps['plasma'](0.8)[0]*255)), int(np.round(new_cmaps.cmaps['plasma'](0.8)[1]*255)), int(np.round(new_cmaps.cmaps['plasma'](0.8)[2]*255)))
        cfgDict['MainParamDict']['ion_fit_color'] = 'r'
        cfgDict['MainParamDict']['electron_fit_color'] = 'yellow'

    else:
        cfgDict['MainParamDict']['ion_color'] = "#{0:02x}{1:02x}{2:02x}".format(int(np.round(new_cmaps.cmaps['viridis'](0.45)[0]*255)), int(np.round(new_cmaps.cmaps['viridis'](0.45)[1]*255)), int(np.round(new_cmaps.cmaps['viridis'](0.45)[2]*255)))
        cfgDict['MainParamDict']['electron_color'] ="#{0:02x}{1:02x}{2:02x}".format(int(np.round(new_cmaps.cmaps['viridis'](0.75)[0]*255)), int(np.round(new_cmaps.cmaps['viridis'](0.75)[1]*255)), int(np.round(new_cmaps.cmaps['viridis'](0.75)[2]*255)))
        cfgDict['MainParamDict']['ion_fit_color'] = 'mediumturquoise'
        cfgDict['MainParamDict']['electron_fit_color'] = 'lime'

    ###
    #
    # NOW LOAD ALL THE SUBPLOTS
    #
    ###
    subplotCfg = {'PhasePlot': {},
                  'EnergyPlot':{},
                  'FieldsPlot': {},
                  'DensityPlot': {},
                  'SpectraPlot': {},
                  'MagPlots': {},
                  'FFTPlots': {},
                  'TotalEnergyPlot': {},
                  'Moments':{},
                         }
    subplotCfg['DensityPlot']['plot_param_dict'] = {'twoD': 0,
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
    subplotCfg['DensityPlot']['BoolList'] = ['twoD', 'show_cbar', 'set_color_limits', 'set_v_min', 'set_v_max',
                   'show_labels', 'show_shock', 'OutlineText', 'spatial_x', 'spatial_y',
                   'normalize_density', 'UseDivCmap', 'stretch_colors', 'show_cpu_domains']
    subplotCfg['DensityPlot']['IntList'] = ['dens_type']
    subplotCfg['DensityPlot']['FloatList'] = ['v_min', 'v_max', 'cpow_num', 'div_midpoint']
    #StrList = ['interpolation', 'cnorm_type', 'cmap']
    subplotCfg['DensityPlot']['StrList'] = ['interpolation', 'cnorm_type', 'cmap']

    subplotCfg['EnergyPlot']['plot_param_dict'] = {'twoD' : 1,
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

    subplotCfg['EnergyPlot']['BoolList'] = ['twoD', 'masked', 'weighted', 'show_cbar', 'show_shock', 'show_int_region', 'set_color_limits',
                'set_v_min', 'set_v_max', 'set_y_min', 'set_y_max', 'spatial_x', 'spatial_y']
    subplotCfg['EnergyPlot']['IntList']= ['prtl_type', 'xbins', 'ebins']
    subplotCfg['EnergyPlot']['FloatList'] = ['v_min', 'v_max', 'y_min', 'y_max', 'cpow_num']
    #StrList = ['interpolation', 'cnorm_type']
    subplotCfg['EnergyPlot']['StrList'] = ['cnorm_type'] # No longer loading interpolation from config files

    subplotCfg['FFTPlots']['plot_param_dict'] = {'twoD': 0,
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
    subplotCfg['FFTPlots']['BoolList'] = ['twoD', 'set_y_min', 'set_y_max',
                'xLog', 'yLog', 'spatial_x', 'spatial_y']
    subplotCfg['FFTPlots']['IntList'] = ['FFT_type']
    subplotCfg['FFTPlots']['FloatList'] = ['y_min', 'y_max']
    subplotCfg['FFTPlots']['StrList'] = []
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
    subplotCfg['FieldsPlot']['plot_param_dict'] = {'twoD': 0,
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

    subplotCfg['FieldsPlot']['BoolList'] = ['twoD', 'show_cbar', 'set_v_min', 'set_v_max',
             'show_shock', 'OutlineText', 'spatial_x', 'spatial_y',
             'Show_FFT_region', 'UseDivCmap', 'stretch_colors', 'normalize_fields',
             'show_x', 'show_y', 'show_z', 'show_cpu_domains']
    subplotCfg['FieldsPlot']['IntList'] = ['field_type']
    subplotCfg['FieldsPlot']['FloatList'] = ['v_min', 'v_max', 'cpow_num', 'div_midpoint']
    #StrList = ['interpolation', 'cnorm_type', 'cmap']
    subplotCfg['FieldsPlot']['StrList'] = ['cnorm_type', 'cmap', 'interpolation', 'cmdstr1', 'cmdstr2', 'cmdstr3']
    subplotCfg['FieldsPlot']['SpecialList'] = ['yaxis_label', '2D_label', '2D_label', 'OneDOnly']

    subplotCfg['MagPlots']['plot_param_dict'] = {'twoD': 0,
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
    subplotCfg['MagPlots']['BoolList'] = ['twoD', 'show_cbar', 'set_v_min', 'set_v_max',
             'show_shock', 'OutlineText', 'spatial_x', 'spatial_y', 'Show_FFT_region',
             'UseDivCmap', 'stretch_colors']
    subplotCfg['MagPlots']['IntList'] = ['mag_plot_type']
    subplotCfg['MagPlots']['FloatList'] = ['v_min', 'v_max', 'cpow_num', 'div_midpoint']
    #StrList = ['interpolation', 'cnorm_type', 'cmap'] # Not loading interpolation anymore
    subplotCfg['MagPlots']['StrList'] = ['interpolation','cnorm_type', 'cmap']

    subplotCfg['Moments']['plot_param_dict'] = {'twoD': 0,
                       'm_type': 0, # 0 = average_velocity, 1 = average_momentum, 2 = Energy
                       'v_min': 0,
                       'v_max' : 10,
                       'set_v_min': False,
                       'set_v_max': False,
                       'show_x': True,
                       'show_y': False,
                       'show_z': False,
                       'show_ions': True,
                       'show_electrons': True,
                       'show_total': False,
                       'UpstreamFrame': False,
                       'weighted': False,
                       'xbins': 100,
                       'spatial_x': True,
                       'show_legend': True,
                       'spatial_y': False,
                       'logy': False,
                       'legend_loc': 'N/A',
                       'face_color': 'gainsboro'
                       }
                         # We need the types of all the parameters for the config file
    subplotCfg['Moments']['BoolList'] = ['twoD', 'set_v_min', 'set_v_max',
                'UpstreamFrame', 'show_ions', 'show_electrons', 'show_x',
                'show_y', 'show_z', 'spatial_x', 'spatial_y', 'weighted', 'show_legend']
    subplotCfg['Moments']['IntList'] = ['m_type','xbins']
    subplotCfg['Moments']['FloatList'] = ['v_min', 'v_max']
    subplotCfg['Moments']['StrList'] = ['legend_loc']

    subplotCfg['PhasePlot']['plot_param_dict'] = {'twoD' : 1,
                   'mom_dim': 0,
                   'masked': 1,
                   'cnorm_type': 'Log', #Colormap normalization. Opts are Log or Linear
                   'prtl_type': 0,
                   'cpow_num': 0.6,
                   'show_cbar': True,
                   'weighted': False,
                   'show_shock': False,
                   'show_int_region': True,
                   'xbins' : 200,
                   'pbins' : 200,
                   'v_min': -2.0,
                   'v_max' : 0,
                   'set_v_min': False,
                   'set_v_max': False,
                   'p_min': -2.0,
                   'p_max' : 2,
                   'set_E_min' : False,
                   'E_min': 1.0,
                   'set_E_max': False,
                   'E_max': 200.0,
                   'set_p_min': False,
                   'set_p_max': False,
                   'spatial_x': True,
                   'spatial_y': False,
                   'interpolation': 'nearest',
                   'face_color': 'gainsboro'}

    subplotCfg['PhasePlot']['BoolList'] = ['twoD', 'masked', 'weighted', 'show_cbar', 'show_shock', 'show_int_region',
                'set_v_min', 'set_v_max', 'set_p_min', 'set_p_max', 'set_E_min', 'Set_E_max', 'spatial_x', 'spatial_y']
    subplotCfg['PhasePlot']['IntList'] = ['prtl_type','mom_dim', 'xbins', 'pbins']
    subplotCfg['PhasePlot']['FloatList'] = ['v_min', 'v_max', 'p_min', 'p_max', 'cpow_num', 'E_min', 'E_max']
    #StrList = ['interpolation', 'cnorm_type']
    subplotCfg['PhasePlot']['StrList'] = ['cnorm_type'] # No longer loading interpolation from config files

    subplotCfg['SpectraPlot']['plot_param_dict'] = {'spectral_type': 0, #0 dn/dp, 1 = dn/dE
                       'twoD': 0, # Plot is not 2D
                       'show_ions': 1, # Include ion spectrum?
                       'show_electrons': 1, # Include electron Spectrum?
                       'rest_frame': False, # Show the data in the rest frame?
                       'set_ylim': True,
                       'set_xlim': True,
                       'x_min': 0.05,
                       'x_max': 200,
                       'spatial_x': False, # Important info for the main program
                       'spatial_y': False,
                       'show_legend': True, # Include legend?
                       'normalize_spectra': True, # To normalze the spectrum so it's integral = 1
                       'y_min': -6, # in logspace
                       'y_max': 0, # in logspace
                       'MeasureEpsP': False, #
                       'IonLeft': -10000.0,
                       'PowerLawIonMax': 10.0,
                       'GammaIonInjection': 1.0,
                       'DoPowerLawFitIon': False,
                       'DoPowerLawFitElectron': False,
                       'DelGami': 0.06,
                       'DelGame': 0.03,
                       'SetTi': False,
                       'PowerLawIonMin': 1.0,
                       'SetTe': False,
                       'GammaElectronInjection': 30.0,
                       'PowerLawElectronMin': 1.0,
                       'ElectronRight': 0.0,
                       'PrtlIntegrationRelative': True,
                       'PowerLawElectronMax': 10.0,
                       'ElectronLeft': -10000.0,
                       'IonRight': 0.0,
                       'MeasureEpsE': False,
                       'eNormalizer' : 0.,
                       'iNormalizer' : 0.,
                       'BoostedIons' : False, #Show a energy boosted ion plot to compare to electrons
                       'T_legend_loc':  'N/A', # location of temperature legend.
                       'PL_legend_loc': 'N/A',
                       'face_color': 'gainsboro'} # lcation of power law legend
    subplotCfg['SpectraPlot']['BoolList'] = ['twoD', 'set_ylim', 'set_xlim',
                'spatial_x', 'spatial_y', 'normalize_spectra',
                'show_ions', 'show_electrons', 'rest_frame','show_legend',
                'PrtlIntegrationRelative', 'BoostedIons'
                'SetTe', 'SetTi','MeasureEpsP', 'MeasureEpsE',
                'DoPowerLawFitElectron', 'DoPowerLawFitIon']
    subplotCfg['SpectraPlot']['IntList'] = ['spectral_type']
    subplotCfg['SpectraPlot']['FloatList'] = ['y_min', 'y_max', 'x_min', 'x_max',
                 'DelGame', 'DelGami', 'GammaIonInjection', 'GammaElectronInjection',
                 'PowerLawElectronMin', 'PowerLawElectronMax',
                 'PowerLawIonMin', 'PowerLawIonMax', 'eNormalizer', 'iNormalizer',
                 'ElectronLeft', 'ElectronRight', 'IonLeft', 'IonRight',]
    subplotCfg['SpectraPlot']['StrList'] = ['T_legend_loc', 'PL_legend_loc']
    subplotCfg['TotalEnergyPlot']['plot_param_dict'] = {'twoD': 0,
                   'show_prtl_KE': True,
                   'show_field_E': True,
                   'show_ion_E': False,
                   'show_electron_E': False,
                   'show_total_E': True,
                   'show_Bz_energy': False,
                   'show_B_E': False,
                   'show_E_E': False,
                   'y_min': 0,
                   'y_max' : 10,
                   'set_y_min': False,
                   'set_y_max': False,
                   'show_legend': True,
                   'show_current_time': True,
                   'x_min': 0,
                   'x_max' : 10,
                   'set_x_min': False,
                   'set_x_max': False,
                   'yLog': True,
                   'spatial_x': False,
                   'spatial_y': False,
                   'legend_loc': 'N/A',
                   'face_color': 'gainsboro'}

    # We need the types of all the parameters for the config file
    subplotCfg['TotalEnergyPlot']['BoolList'] = ['twoD', 'set_y_min', 'set_y_max','show_prtl_KE', 'show_field_E','show_B_E', 'show_E_E',
                'show_ion_E', 'show_electron_E', 'show_total_E', 'yLog', 'set_x_min',
                'set_x_max', 'show_legend','show_current_time', 'show_Bz_energy'] # spatial_x and spatial_y should never be true, remove from boollist
    subplotCfg['TotalEnergyPlot']['IntList'] = ['E_type']
    subplotCfg['TotalEnergyPlot']['FloatList'] = ['y_min', 'y_max','x_min', 'x_max']
    subplotCfg['TotalEnergyPlot']['StrList'] = ['legend_loc']

    for i in range(cfgDict['MainParamDict']['NumOfRows']):
        for j in range(cfgDict['MainParamDict']['NumOfCols']):
            tmp_str = 'Chart' + str(i) + ',' + str(j)
            new_str = 'Chart' + str(i) + '_' + str(j)
            cfgDict[new_str] = {}
            if tmp_str in config.sections():
                tmpchart_type = config.get(tmp_str, 'ChartType')
                cfgDict[new_str]['ChartType'] = tmpchart_type
                #Now load in all the parameters from the config file
                for param in subplotCfg[tmpchart_type]['BoolList']:
                    if param.lower() in config.options(tmp_str):
                        cfgDict[new_str][param] = config.getboolean(tmp_str, param)

                for param in subplotCfg[tmpchart_type]['IntList']:
                    if param.lower() in config.options(tmp_str):
                        cfgDict[new_str][param] = int(config.getfloat(tmp_str, param))

                for param in subplotCfg[tmpchart_type]['FloatList']:
                    if param.lower() in config.options(tmp_str):
                        cfgDict[new_str][param] = config.getfloat(tmp_str, param)

                for param in subplotCfg[tmpchart_type]['StrList']:
                    if param.lower() in config.options(tmp_str):
                        cfgDict[new_str][param] = config.get(tmp_str, param)
                try:
                    for param in subplotCfg[tmpchart_type]['SpecialList']:
                        if param.lower() in config.options(tmp_str):
                            if tmpchart_type == 'FieldsPlot':
                                if param == 'yaxis_label' or 'OneDOnly':
                                    tstr = config.get(tmp_str, param)
                                    cfgDict[new_str][param] = tstr.strip('[').strip(']').strip().split(',')
                            if tmpchart_type == 'FieldsPlot':
                                if param == '2D_label' or param =='1D_label':
                                    tstr = config.get(tmp_str, param)
                                    print(tstr)
                                    tstr = tstr.replace('[', '')
                                    tstr = tstr.replace(']', '')
                                    tstr = tstr.replace(' ', '')
                                    tstr = tstr.replace("'", '')
                                    tstr = tstr.replace('\\\\', '\\')
                                    flattened_list = tstr.split(',')
                                    print(flattened_list)
                                    # NOW MAKE IT A LIST OF LISTs
                                    cfgDict[new_str][param] = [*map(list, zip(*[iter(flattened_list)]*3))]
                                    print(cfgDict[new_str][param])
                except KeyError:
                    pass
            else:
                # The graph isn't specified in the config file, just set it equal to a phase plot
                cfgDict[new_str]['ChartType'] = 'PhasePlot'


    # WRITE THE NEW YAML FILE
    cfgfile = config_file.rstrip('.cfg') + '.yml'
    with open(cfgfile, 'w') as cfgFile:
        cfgFile.write(yaml.safe_dump(cfgDict))

if __name__ == '__main__':
    cfg_files = os.listdir(os.path.join(os.curdir, '.iseult_configs'))

    for cfile in cfg_files:
        if cfile[-3:] == 'cfg':
             convertOldConfig(os.path.join(os.path.join(os.curdir,'.iseult_configs'), cfile))
