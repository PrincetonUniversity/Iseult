#! /usr/bin/env python
import re # regular expressions
import os, sys # Used to make the code portable
import h5py # Allows us the read the data files
import time, string, io
from PIL import Image
import matplotlib
matplotlib.use('TkAgg')
import new_cmaps
import numpy as np
from collections import deque
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.animation as manimation
from matplotlib.figure import Figure
from phase_plots import PhasePanel
from fields_plots import FieldsPanel
from density_plots import DensPanel
from spectra import SpectralPanel
from mag_plots import BPanel
from energy_plots import EnergyPanel
from fft_plots import FFTPanel
from total_energy_plots import TotEnergyPanel
from moments import MomentsPanel
from functools import partial
import subprocess, yaml
#import datetime
#from ThreeD_mag_plots import ThreeDBPanel STILL TESTING

# I don't think that matplotlib allows multi-threading, in the interactive mode.
# This is a flag that i have so I can mess around trying to get it to work.
Use_MultiProcess = False # DO NOT SET TO TRUE!
import time
import tkinter as Tk
from tkinter import ttk, filedialog, messagebox

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['image.resample'] = False
matplotlib.rcParams['image.origin'] = 'upper'

import argparse

def destroy(e):
    sys.exit()

class MyCustomToolbar(NavigationToolbar2Tk):
    def __init__(self, plotCanvas, parent):
        # create the default toolbar
        # plotCanvas is the tk Canvas we want to link to the toolbar,
        # parent is the iseult main app
        NavigationToolbar2Tk.__init__(self, plotCanvas, parent)
        #print(self._nav_stack)
        self.parent = parent

class Spinbox(ttk.Entry):
    def __init__(self, master=None, **kw):
        ttk.Entry.__init__(self, master, "ttk::spinbox", **kw)

    def current(self, newindex=None):
        return self.tk.call(self._w, 'current', index)

    def set(self, value):
        return self.tk.call(self._w, 'set', value)

class SubPlotWrapper:
    """A simple class that will eventually hold all of the information
    about each sub_plot in the Figure"""

    def __init__(self, parent, figure=None, pos = None, subplot_spec = None, ctype=None, graph = None):
        # What happens when a SubPlotWrapper is initiated
        self.parent = parent # Define the parent. MUST BE MainApp OF ISEULT
        self.chartType = 'PhasePlot' # The default chartype
        # A dictionary that contains all of the plot types.
        # The panels must follow certain conventions. TODO: Make a class and
        # Have the panels extend this class.
        self.PlotTypeDict = {'PhasePlot': PhasePanel,
                             'EnergyPlot': EnergyPanel,
                             'FieldsPlot': FieldsPanel,
                             'DensityPlot': DensPanel,
                             'SpectraPlot': SpectralPanel,
                             'MagPlots': BPanel,
                             'FFTPlots': FFTPanel,
                             'TotalEnergyPlot': TotEnergyPanel,
                             'Moments': MomentsPanel
                             }
        #####
        #
        # First we create dictionary of dictionarys that will store all of the
        # plot params values at self.PlotParamsDict['ctype']['param_name'],
        #
        # We also need a dictionary that returns a list of all of the params of
        # a certain type, used in loading the config files.
        # self.ParamsTypeDict['ctype']['BoolList'] returns all of the Booleans
        # stored in the self.PlotParamsDict['ctype'] dictionary.
        # BoolList, IntList, FloatList and StrList are the options, and the only
        # types that are allowed self.PlotParamsDict['ctype'] dictionary.
        #
        ####
        self.GenParamDict()
        self.figure = figure
        self.subplot_spec = subplot_spec
        self.pos = pos
        self.graph = graph # The panel class-- e.g. PhasesPanel, FieldsPanel...etc
        self.Changedto1D = False # needed to keep track of color bars and views
        self.Changedto2D = False # needed to keep track of color bars and views
        #
    def GetKeys(self):
        ''' A function that returns a list of all of the keys required to plot
        the subplot contained within SubPlotWrapper. the set_plot_keys function
        must be defined in each of the subplot panel classes.'''
        return self.graph.set_plot_keys()

    def LoadKey(self, h5key):
        '''This is a function the graph should call to load a particular key
        stored in the Tristan outputfiles'''
        return self.parent.DataDict[h5key]
    def LoadData(self):

        ''' LoadData is called by MainApp, it is defined by the subplot panel
        class, but basically it is a function that should load all of the raw
        output data of the current time slice required to make the plot and
        then perform all the necessary calculations to calculate the quantity
        of interest.'''
        self.graph.LoadData()
    def ChangeGraph(self, str_arg):
        '''ChangeGraph changes the plotted graph to the one given by str_arg.
        str_arg must be a key in self.PlotTypeDict'''

        # First check if the current plot has a color bar
        tmpIs2D = self.PlotParamsDict[self.chartType]['twoD']

        # Change the graph type
        self.chartType = str_arg
        # put a list of the previous chart types in iseult

        self.graph = self.PlotTypeDict[self.chartType](self.parent, self)

        # If the graph changes from 1D or 2D, we need to save this
        if tmpIs2D != self.PlotParamsDict[self.chartType]['twoD']:
            self.Changedto1D = not self.PlotParamsDict[self.chartType]['twoD']
            self.Changedto2D = self.PlotParamsDict[self.chartType]['twoD']

        self.parent.RenewCanvas(ForceRedraw = True)

    def GenParamDict(self):
        '''First we create dictionary of dictionarys that will store all of the
        plot params values at self.PlotParamsDict['ctype']['param_name'],

        We also need a dictionary that returns a list of all of the params of
        a certain type, used in loading the config files.
        self.ParamsTypeDict['ctype']['BoolList'] returns all of the Booleans
        stored in the self.PlotParamsDict['ctype'] dictionary.
        BoolList, IntList, FloatList and StrList are the options, and the only
        types that are allowed self.PlotParamsDict['ctype'] dictionary.
        Generate a dictionary that will store all of the params at dict['ctype']['param_name']
        '''

        self.PlotParamsDict = {plot_type: '' for plot_type in self.PlotTypeDict.keys()}
        for elm in self.PlotTypeDict.keys():

            self.PlotParamsDict[elm] = {key: self.PlotTypeDict[elm].plot_param_dict[key] for key in self.PlotTypeDict[elm].plot_param_dict.keys()}




    def RestoreDefaultPlotParams(self, ctype = None, RestoreAll = False):
        ''' Restore the PlotParamsDictionary to the default values contained in
        the class definition.'''
        if ctype is None: # restore the currently shown plot
            ctype = self.chartType
        if RestoreAll: # Restore all of the plot types for this subplot
            self.GenParamDict()
        else: # Only restore the values of the ctype listed.
            self.PlotParamsDict[ctype] = {key: self.PlotTypeDict[ctype].plot_param_dict[key] for key in self.PlotTypeDict[ctype].plot_param_dict.keys()}

    def SetPlotParam(self, pname, val, ctype = None, update_plot = True, NeedsRedraw = False):
        ''' A function that changes plot param 'pname' for the ctype chart to val
        in the dictionary held by SubPlotWrapper. If update_plot is true,
        the figure will be refreshed. If NeedsRedraw is true, the figure will
        be cleared then re-drawn.'''


        if ctype is None:
            ctype = self.chartType
        # Check to see if a plot is changed from 1d to 2d
        if pname =='twoD':
            if self.PlotParamsDict[ctype][pname] == 1 and val == 0:
                self.Changedto1D = True
                NeedsRedraw = True
            if self.PlotParamsDict[ctype][pname] == 0 and val == 1:
                self.Changedto2D = True
                NeedsRedraw = True

        self.PlotParamsDict[ctype][pname] = val
        if update_plot or NeedsRedraw:
            self.parent.RenewCanvas(ForceRedraw = NeedsRedraw)



    def GetPlotParam(self, pname, ctype = None):
        ''' A function that returns the value of the plot param 'pname' for the
        ctype chart. If ctype is None, the currently shown subplot is used'''

        if ctype is None:
            ctype = self.chartType
        return self.PlotParamsDict[ctype][pname]

    def SetGraph(self, ctype = None):
        ''' SetGraph is useful if you want to change the plot type without redrawing the figure
        (e.g. if you want to set many graphs at once, like when loading in a config file.)'''
        if ctype:
            self.chartType = ctype
        self.graph = self.PlotTypeDict[self.chartType](self.parent, self)

    def DrawGraph(self):
        ''' This function calls a function that must be defined in the subplotpanel
         class, e.g. FieldsPanel.... It creates all of the axes used in the panel,
         and writes the data. It should be called after clearing a figure, the
         chartype is changed or when first initializing the figure. It will be
         called if RenewCanvas(ForceRedraw = True)'''

        self.graph.draw()

    def RefreshGraph(self):
        ''' This function calls a function that must be defined in the subplotpanel
         class, e.g. FieldsPanel.... It only updates things held by  the panel to the data in output files with the current timestep.
         It should be called when stepping through the times, or possibly when a plot param changes.
         It will be called if RenewCanvas(ForceRedraw = False)'''
        self.graph.refresh()

    def OpenSubplotSettings(self):

        ''' A function that that must be defined in the subplotpanel class, e.g.
        FieldsPanel.... Opens up the pop-up that allows one to change parameters
        of the plot, change chart type, etc.'''

        self.graph.OpenSettings()

    def SetCpuDomainLines(self):
        '''This function sets the Cpu lines up. It should only be called when
        redrawing the axes and after the axes is creates as it creates all of
        the line objects.'''

        # regardless if it is 1D or 2D we'll show the x_domains...
        # This could change if we decide to add the ability to show transverse 1D slices
        self.cpu_x_lines = []
        self.cpu_y_lines = []
        for i in range(len(self.parent.cpu_x_locs)):
            self.cpu_x_lines.append(self.graph.axes.axvline(self.parent.cpu_x_locs[i], linewidth = 1, linestyle = ':',color = 'w') )
        if self.GetPlotParam('twoD'):
            for i in range(len(self.parent.cpu_y_locs)):
                self.cpu_y_lines.append(self.graph.axes.axhline(self.parent.cpu_y_locs[i], linewidth = 1, linestyle = ':',color = 'w'))

    def UpdateCpuDomainLines(self):
        '''This updates the location of the Cpu lines. It should only be called
        when refreshing the axes as it requires the line objects to already be
        created.'''
        # regardless if it is 1D or 2D we'll show the x_domains...
        # This could change if we decide to add the ability to show transverse 1D slices
        for i in range(len(self.cpu_x_lines)):
            self.cpu_x_lines[i].set_xdata([self.parent.cpu_x_locs[i],self.parent.cpu_x_locs[i]])

        if self.GetPlotParam('twoD'):
            for i in range(len(self.parent.cpu_y_locs)):
                self.cpu_y_lines[i].set_ydata([self.parent.cpu_y_locs[i],self.parent.cpu_y_locs[i]])


    def RemoveCpuDomainLines(self):
        '''This removes the Cpu lines. It should only be called
        when the user unselects show CPU domains.'''
        # iterate over the line list and destroy the objects
        for i in range(len(self.cpu_x_lines)):
            self.cpu_x_lines.pop().remove()
        for i in range(len(self.cpu_y_lines)):
            self.cpu_y_lines.pop().remove()
class Knob:
    """
    ---- Taken from the Matplotlib gallery
    Knob - simple class with a "setKnob" method.
    A Knob instance is attached to a Param instance, e.g., param.attach(knob)
    Base class is for documentation purposes.
    """
    def setKnob(self, value):
        pass

class Param:
    """
    ---- Taken from the Matplotlib gallery
    The idea of the "Param" class is that some parameter in the GUI may have
    several knobs that both control it and reflect the parameter's state, e.g.
    a slider, text, and dragging can all change the value of the frequency in
    the waveform of this example.
    The class allows a cleaner way to update/"feedback" to the other knobs when
    one is being changed.  Also, this class handles min/max constraints for all
    the knobs.
    Idea - knob list - in "set" method, knob object is passed as well
      - the other knobs in the knob list have a "set" method which gets
        called for the others.
    """
    def __init__(self, initialValue=None, minimum=0., maximum=1.):
        self.minimum = minimum
        self.maximum = maximum
        if initialValue != self.constrain(initialValue):
            raise ValueError('illegal initial value')
        self.value = initialValue
        self.knobs = []

    def attach(self, knob):
        self.knobs += [knob]

    def set(self, value, knob=None):
        if self.value != self.constrain(value):
            self.value = self.constrain(value)
            for feedbackKnob in self.knobs:
                if feedbackKnob != knob:
                    feedbackKnob.setKnob(self.value)
        # Adding a new feature that allows one to loop backwards or forwards:
        elif self.maximum != self.minimum:

            if self.value == self.maximum:
                self.value = self.minimum
                for feedbackKnob in self.knobs:
                    if feedbackKnob != knob:
                        feedbackKnob.setKnob(self.value)

            elif self.value == self.minimum:
                self.value = self.maximum
                for feedbackKnob in self.knobs:
                    if feedbackKnob != knob:
                        feedbackKnob.setKnob(self.value)
        return self.value

    def setMax(self, max_arg, knob=None):
        self.maximum = max_arg
        self.value = self.constrain(self.value)
        for feedbackKnob in self.knobs:
            if feedbackKnob != knob:
                feedbackKnob.setKnob(self.value)
        return self.value

    def constrain(self, value):
        if value <= self.minimum:
            value = self.minimum
        if value >= self.maximum:
            value = self.maximum
        return value


class PlaybackBar(Tk.Frame):

    """
    A Class that will handle the time-stepping in Iseult, and has the
    following, a step left button, a play/pause button, a step right button, a
    playbar, and a settings button.
    """

    def __init__(self, parent, param, canvas = None):
        Tk.Frame.__init__(self)
        self.parent = parent
        self.playPressed = False

        # This param should be the time-step of the simulation
        self.param = param

        # make a button that skips left
        self.skipLB = ttk.Button(self, text = '<', command = self.SkipLeft)
        self.skipLB.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=0)

        # make the play button
        self.playB = ttk.Button(self, text = 'Play', command = self.PlayHandler)
        self.playB.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=0)

        # a button that skips right
        self.skipRB = ttk.Button(self, text = '>', command = self.SkipRight)
        self.skipRB.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=0)

        # An entry box that will let us choose the time-step
        ttk.Label(self, text='n= ').pack(side=Tk.LEFT, fill=Tk.BOTH, expand=0)

        # A StringVar for a box to type in a frame num, linked to self.param
        self.tstep = Tk.StringVar()
        # set it to the param value
        self.tstep.set(str(self.param.value))

        # the entry box
        self.txtEnter = ttk.Entry(self, textvariable=self.tstep, width=6)
        self.txtEnter.pack(side=Tk.LEFT, fill = Tk.BOTH, expand = 0)

        # A slider that will show the progress in the simulation as well as
        # allow us to select a time. Now the slider just changes the tstep box
        self.slider = ttk.Scale(self, from_=self.param.minimum, to=self.param.maximum, command = self.ScaleHandler)
        self.slider.set(self.param.value)
        self.slider.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)
        # bind releasing the moust button to updating the plots.
        self.slider.bind("<ButtonRelease-1>", self.UpdateValue)

        new_frame = ttk.Frame(self)
        self.LoopVar = Tk.IntVar()
        self.LoopVar.set(self.parent.MainParamDict['LoopPlayback'])
        self.LoopVar.trace('w', self.LoopChanged)
        self.RecordFrames = ttk.Checkbutton(new_frame, text = 'Loop',
                                            variable = self.LoopVar)
        self.RecordFrames.pack(side=Tk.TOP, fill=Tk.BOTH, expand=0)


        self.RecVar = Tk.IntVar()
        self.RecVar.set(self.parent.MainParamDict['Recording'])
        self.RecVar.trace('w', self.RecChanged)
        ttk.Checkbutton(new_frame, text = 'Record',
                        variable = self.RecVar).pack(side=Tk.TOP, fill=Tk.BOTH, expand=0)
        new_frame.pack(side= Tk.LEFT, fill = Tk.BOTH, expand =0)

        # a measurement button that should lauch a window to take measurements.
        self.MeasuresB= ttk.Button(self, text='FFT', command=self.OpenMeasures)
        self.MeasuresB.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=0)


        # a settings button that should lauch some global settings.
        self.SettingsB= ttk.Button(self, text='Settings', command=self.parent.OpenSettings)
        self.SettingsB.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=0)

        # a reload button that reloads the files and then refreshes the plot
        ttk.Button(self, text = 'Reload', command = self.OnReload).pack(side=Tk.LEFT, fill=Tk.BOTH, expand=0)
        # a refresh button that refreshing the current timestep
        ttk.Button(self, text = 'Refresh', command = self.OnRefresh).pack(side=Tk.LEFT, fill=Tk.BOTH, expand=0)
        #attach the parameter to the Playbackbar
        self.param.attach(self)

    def OnReload(self, *args):
        self.parent.ReloadPath()
        self.parent.RenewCanvas()

    def OnRefresh(self, *args):
        self.parent.RefreshTimeStep()
        self.parent.RenewCanvas()

    def RecChanged(self, *args):
        if self.RecVar.get() == self.parent.MainParamDict['Recording']:
            pass
        else:
            self.parent.MainParamDict['Recording'] = self.RecVar.get()
            if self.parent.MainParamDict['Recording'] == 1:
                self.parent.PrintFig()

    def LoopChanged(self, *args):
        if self.LoopVar.get() == self.parent.MainParamDict['LoopPlayback']:
            pass
        else:
            self.parent.MainParamDict['LoopPlayback'] = self.LoopVar.get()

    def SkipLeft(self, e = None):
        self.param.set(self.param.value - self.parent.MainParamDict['SkipSize'])

    def SkipRight(self, e = None):
        self.param.set(self.param.value + self.parent.MainParamDict['SkipSize'])

    def PlayHandler(self, e = None):
        if not self.playPressed:
            # Set the value of play pressed to true, change the button name to
            # pause, turn off clear_fig, and start the play loop.
            self.playPressed = True
            self.parent.RenewCanvas()
            """
            if not self.parent.MainParamDict['Recording']:
                #self.parent.HashIseultState()
                already_saved = False
                if self.parent.TimeStep.value in self.parent.SavedHashes.keys(): # we have already saved an image for this TimeStep
                    # is the current state of Iseult equal to the state when we saved said image?
                    already_saved = self.parent.SavedHashes[self.parent.TimeStep.value] ==  self.parent.StateHash

                if not already_saved:
                    self.parent.RenewCanvas()

                # Prevent the window from being resized
                self.parent.resizable(0,0)
                #            self.parent.MainParamDict['ClearFig'] = False
                tmp_size = self.parent.f.get_size_inches()*self.parent.f.dpi

                # Create the figure
                self.MovieFrame = ttk.Frame(self.parent)
                self.parent.MovieFig = Figure(figsize = self.parent.f.get_size_inches(), dpi = self.parent.f.dpi, edgecolor = 'none')#, facecolor = '0.75')
                self.parent.MovieFig.subplots_adjust(left = 0, right = 1, top = 1, bottom = 0 , wspace = 0, hspace = 0)
                # a tk.DrawingArea
                self.parent.MovieCanvas = FigureCanvasTkAgg(self.parent.MovieFig, master=self.MovieFrame)
                #            self.parent.MovieCanvas = Tk.Canvas(self.parent, width=tmp_size[0], height=tmp_size[1])

                im = Image.frombuffer('RGBA', (int(tmp_size[0]), int(tmp_size[1])), self.parent.SavedImgStr[self.parent.TimeStep.value], 'raw', 'RGBA', 0, 1)
                self.MovieFrame.place(in_=self.parent, relx=0.5, y=0, anchor=Tk.N)#, bordermode="outside")
                self.parent.MovieCanvas._tkcanvas.pack(side=Tk.RIGHT, fill=Tk.BOTH, expand=1)
                self.parent.MovieAx = self.parent.MovieFig.add_subplot(111)
                self.parent.MovieAx.axis('off')
                self.parent.MovieIm = self.parent.MovieAx.imshow(im, interpolation = 'nearest')
                self.parent.MovieCanvas.get_tk_widget().update_idletasks()
            """
            self.playB.config(text='Pause')

            self.after(int(self.parent.MainParamDict['WaitTime']*1E3), self.blink)
        else:
            self.parent.resizable(1,1)
            # pause the play loop, turn clear fig back on, and set the button name back to play
            self.playPressed = False
            try:
                self.MovieFrame.destroy()
            except AttributeError:
                pass
            self.parent.RenewCanvas()
#            self.parent.MainParamDict['ClearFig'] = True
            self.playB.config(text='Play')


    def OpenMeasures(self):
        if self.parent.measure_window is None:
            self.parent.measure_window = MeasureFrame(self.parent)
        else:
            self.parent.measure_window.destroy()
            self.parent.measure_window = MeasureFrame(self.parent)


    def blink(self):
        if self.playPressed:
            # First check to see if the timestep can get larger
            if self.param.value == self.param.maximum and not self.parent.MainParamDict['LoopPlayback']:
                # push pause button
                self.PlayHandler()

            # otherwise skip right by size skip size
            else:
                self.param.set(self.param.value + self.parent.MainParamDict['SkipSize'])

            # start loopin'
            self.after(int(self.parent.MainParamDict['WaitTime']*1E3), self.blink)


    def TextCallback(self):
        try:
            #make sure the user types in a int
            if int(self.tstep.get()) != self.param.value:
                self.param.set(int(float(self.tstep.get())))
        except ValueError:
            #if they type in random stuff, just set it ot the param value
            self.tstep.set(str(self.param.value))

    def ScaleHandler(self, e):
        # if changing the scale will change the value of the parameter, do so
        try:
            if int(self.tstep.get()) != int(self.slider.get()):
                self.tstep.set(str(int(self.slider.get())))
        except ValueError:
            #if they type in random stuff, just set it ot the param value
            self.tstep.set(str(int(self.slider.get())))

    def UpdateValue(self, *args):
        if int(self.slider.get()) != self.param.value:
            self.param.set(int(self.slider.get()))
    def setKnob(self, value):
        pass
#        #set the text entry value
#        self.tstep.set(str(value))
        #set the slider
#        self.slider.set(value)


class SaveDialog(Tk.Toplevel):

    def __init__(self, parent, title = None):

        Tk.Toplevel.__init__(self, parent)
        self.transient(parent)

        if title:
            self.title(title)

        self.parent = parent

        self.result = None

        body = ttk.Frame(self)
        self.initial_focus = self.body(body)
#        body.pack(fill=Tk.BOTH)#, expand=True)
        body.pack(fill = Tk.BOTH, anchor = Tk.CENTER, expand=1)

        self.buttonbox()

        self.grab_set()

        if not self.initial_focus:
            self.initial_focus = self

        self.protocol("WM_DELETE_WINDOW", self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,
                                  parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    #
    # construction hooks

    def body(self, master):
        # create dialog body.  return widget that should have
        # initial focus.  this method should be overridden
        ttk.Label(master, text="Name of View:").grid(row=0)
        self.e1 = ttk.Entry(master, width=17)
        self.e1.grid(row=0, column=1, sticky = Tk.E)

    def buttonbox(self):
        # add standard button box. override if you don't want the
        # standard buttons

        box = ttk.Frame(self)

        w = ttk.Button(box, text="Save", width=10, command=self.ok, default=Tk.ACTIVE)
        w.pack(side=Tk.LEFT, padx=5, pady=5)
        w = ttk.Button(box, text="Cancel", width=10, command=self.cancel)
        w.pack(side=Tk.LEFT, padx=5, pady=5)

        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

        box.pack()

    #
    # standard button semantics

    def ok(self, event=None):

        if not self.validate():
            self.initial_focus.focus_set() # put focus back
            return

        self.withdraw()
        self.update_idletasks()

        self.apply()

        self.cancel()

    def cancel(self, event=None):

        # put focus back to the parent window
        self.parent.focus_set()
        self.destroy()

    #
    # command hooks

    def validate(self):
        ''' Check to make sure the config file doesn't already exist'''
        Name = str(self.e1.get())
#        AlreadyExists = False
#        os.listdir(os.path.join(self.parent.IseultDir, '.iseult_configs'))
        if Name == '':
            messagebox.showwarning(
                "Bad input",
                "Field must contain a name, please try again"
            )
        else:
            return 1 # override

    def apply(self):
        ''' Save the config file'''
        self.parent.SaveIseultState(os.path.join(self.parent.IseultDir, '.iseult_configs', str(self.e1.get()).strip().replace(' ', '_') +'.yml'), str(self.e1.get()).strip())
class MaxNDialog(Tk.Toplevel):

    def __init__(self, parent, title = None):

        Tk.Toplevel.__init__(self, parent)
        self.transient(parent)

        if title:
            self.title(title)

        self.parent = parent

        self.result = None

        body = ttk.Frame(self)
        self.initial_focus = self.body(body)
#        body.pack(fill=Tk.BOTH)#, expand=True)
        body.pack(fill = Tk.BOTH, anchor = Tk.CENTER, expand=1)

        self.buttonbox()

        self.grab_set()

        if not self.initial_focus:
            self.initial_focus = self

        self.protocol("WM_DELETE_WINDOW", self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,
                                  parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    #
    # construction hooks

    def body(self, master):
        # create dialog body.  return widget that should have
        # initial focus.  this method should be overridden
        ttk.Label(master, text="Max Frame (-1 for last frame):").grid(row=0)
        self.e1 = ttk.Entry(master, width=17)
        self.e1.insert(0,str(self.parent.cmd_args.n))
        self.e1.grid(row=0, column=1, sticky = Tk.E)

    def buttonbox(self):
        # add standard button box. override if you don't want the
        # standard buttons

        box = ttk.Frame(self)

        w = ttk.Button(box, text="OK", width=10, command=self.ok, default=Tk.ACTIVE)
        w.pack(side=Tk.LEFT, padx=5, pady=5)

        self.bind("<Return>", self.ok)
        box.pack()

    #
    # standard button semantics

    def ok(self, event=None):

        if not self.validate():
            self.initial_focus.focus_set() # put focus back
            return

        self.withdraw()
        self.update_idletasks()

        self.apply()

        self.cancel()

    def cancel(self, event=None):

        # put focus back to the parent window
        self.parent.focus_set()
        self.destroy()

    #
    # command hooks

    def validate(self):
        ''' Check to make sure the user put a good input in as max file'''
        self.N = str(self.e1.get())
        try:
            self.N = int(self.e1.get())
        except ValueError:
            self.N = ''
        if self.N == '':
            messagebox.showwarning(
                "Bad input",
                "Max N must contain an int, please try again"
            )
        else:
            return 1 # override

    def apply(self):
        '''Update the -n option'''
        self.parent.cmd_args.n = int(self.N)

class MovieDialog(Tk.Toplevel):

    def __init__(self, parent, title = None):

        Tk.Toplevel.__init__(self, parent)
        self.transient(parent)

        if title:
            self.title(title)

        self.parent = parent

        self.result = None

        body = ttk.Frame(self)
        self.initial_focus = self.body(body)
#        body.pack(fill=Tk.BOTH)#, expand=True)
        body.pack(fill = Tk.BOTH, anchor = Tk.CENTER, expand=1)

        self.buttonbox()

        self.grab_set()

        if not self.initial_focus:
            self.initial_focus = self

        self.protocol("WM_DELETE_WINDOW", self.cancel)

        self.geometry("+%d+%d" % (parent.winfo_rootx()+50,
                                  parent.winfo_rooty()+50))

        self.initial_focus.focus_set()

        self.wait_window(self)

    #
    # construction hooks

    def body(self, master):
        # create dialog body.  return widget that should have
        # initial focus.  this method should be overridden
        ttk.Label(master, text="Name of Movie:").grid(row=0)
        self.e1 = ttk.Entry(master, width=17)
        self.e1.grid(row=0, column=1, sticky = Tk.E)

        ttk.Label(master, text="First Frame:").grid(row=1)
        self.e2 = ttk.Entry(master, width=17)
        self.e2.grid(row=1, column=1, sticky = Tk.E)

        ttk.Label(master, text="Last Frame (-1 for final frame):").grid(row=2)
        self.e3 = ttk.Entry(master, width=17)
        self.e3.grid(row=2, column=1, sticky = Tk.E)

        ttk.Label(master, text="Step Size:").grid(row=3)
        self.e4 = ttk.Entry(master, width=17)
        self.e4.grid(row=3, column=1, sticky = Tk.E)

        ttk.Label(master, text="Frames Per Second:").grid(row=4)
        self.e5 = ttk.Entry(master, width=17)
        self.e5.grid(row=4, column=1, sticky = Tk.E)

    def buttonbox(self):
        # add standard button box. override if you don't want the
        # standard buttons

        box = ttk.Frame(self)

        w = ttk.Button(box, text="Save", width=10, command=self.ok, default=Tk.ACTIVE)
        w.pack(side=Tk.LEFT, padx=5, pady=5)
        w = ttk.Button(box, text="Cancel", width=10, command=self.cancel)
        w.pack(side=Tk.LEFT, padx=5, pady=5)

        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

        box.pack()

    #
    # standard button semantics

    def ok(self, event=None):

        if not self.validate():
            self.initial_focus.focus_set() # put focus back
            return

        self.withdraw()
        self.update_idletasks()

        self.apply()

        self.cancel()

    def cancel(self, event=None):

        # put focus back to the parent window
        self.parent.focus_set()
        self.destroy()

    #
    # command hooks

    def validate(self):
        ''' Check to make sure the Movie will work'''
        self.Name = str(self.e1.get())
        try:
            self.StartFrame = int(self.e2.get())
        except ValueError:
            self.StartFrame = ''
        try:
            self.EndFrame = int(self.e3.get())
        except ValueError:
            self.EndFrame = ''
        try:
            self.Step = int(self.e4.get())
        except ValueError:
            self.Step = ''
        try:
            self.FPS = int(self.e5.get())
        except ValueError:
            self.FPS = ''

        if self.Name != '':
            self.Name = str(self.e1.get()).strip().replace(' ', '_') +'.mov'
        if self.StartFrame <0:
            self.StartFrame = len(self.parent.PathDict['Param'])+self.StartFrame + 1
        if self.EndFrame <0:
            self.EndFrame = len(self.parent.PathDict['Param'])+self.EndFrame + 1

        if self.Name == '':
            messagebox.showwarning(
                "Bad input",
                "Field must contain a name, please try again"
            )

        elif self.StartFrame == '':
            messagebox.showwarning(
                "Bad input",
                "StartFrame must contain an int, please try again"
            )
        elif self.EndFrame == '':
            messagebox.showwarning(
                "Bad input",
                "EndFrame must contain an int, please try again"
            )
        elif self.StartFrame == 0:
            messagebox.showwarning(
                "Bad input",
                "Starting frame cannot be zero"
            )
        elif self.EndFrame == 0:
            messagebox.showwarning(
                "Bad input",
                "Ending frame cannot be zero"
            )


        elif self.Step == '':
            messagebox.showwarning(
                "Bad input",
                "Step must contain an int, please try again"
            )
        elif self.Step <=0:
            messagebox.showwarning(
                "Bad input",
                "Step must be an integer >0, please try again"
            )
        elif self.FPS == '':
            messagebox.showwarning(
                "Bad input",
                "FPS must contain an int >0, please try again"
            )
        elif self.FPS <= 0:
            messagebox.showwarning(
                "Bad input",
                "FPS must contain an int >0, please try again"
            )
        else:
            return 1 # override

    def apply(self):
        ''' Save the Movie'''
        self.parent.MakeAMovie(fname = self.Name,
                                start = self.StartFrame,
                                stop = self.EndFrame,
                                step = self.Step,
                                FPS = self.FPS)


class SettingsFrame(Tk.Toplevel):
    def __init__(self, parent):

        Tk.Toplevel.__init__(self)
        self.wm_title('General Settings')
        self.protocol('WM_DELETE_WINDOW', self.OnClosing)

        self.bind('<Return>', self.SettingsCallback)

        self.parent = parent
        frm = ttk.Frame(self)
        frm.pack(fill=Tk.BOTH, expand=True)

        # Make an entry to change the skip size
        self.skipSize = Tk.StringVar(self)
        self.skipSize.set(self.parent.MainParamDict['SkipSize']) # default value
        self.skipSize.trace('w', self.SkipSizeChanged)
        ttk.Label(frm, text="Skip Size:").grid(row=0)
        self.skipEnter = ttk.Entry(frm, textvariable=self.skipSize, width = 6)
        self.skipEnter.grid(row =0, column = 1, sticky = Tk.W + Tk.E)

        # Make an button to change the wait time
        self.waitTime = Tk.StringVar(self)
        self.waitTime.set(self.parent.MainParamDict['WaitTime']) # default value
        self.waitTime.trace('w', self.WaitTimeChanged)
        ttk.Label(frm, text="Playback Wait Time:").grid(row=1)
        self.waitEnter = ttk.Entry(frm, textvariable=self.waitTime, width = 6)
        self.waitEnter.grid(row =1, column = 1, sticky = Tk.W + Tk.E)

        # Have a list of the color maps
        self.cmapvar = Tk.StringVar(self)
        self.cmapvar.set(self.parent.MainParamDict['ColorMap']) # default value
        self.cmapvar.trace('w', self.CmapChanged)

        ttk.Label(frm, text="Color map:").grid(row=2)
        cmapChooser = ttk.OptionMenu(frm, self.cmapvar, self.parent.MainParamDict['ColorMap'], *tuple(new_cmaps.sequential))
        cmapChooser.grid(row =2, column = 1, sticky = Tk.W + Tk.E)

        # Have a list of the color maps
        self.divcmapList = new_cmaps.cmaps.keys()
        self.div_cmapvar = Tk.StringVar(self)
        self.div_cmapvar.set(self.parent.MainParamDict['DivColorMap']) # default value
        self.div_cmapvar.trace('w', self.DivCmapChanged)

        ttk.Label(frm, text="Diverging Cmap:").grid(row=3)
        cmapChooser = ttk.OptionMenu(frm, self.div_cmapvar, self.parent.MainParamDict['DivColorMap'], *tuple(new_cmaps.diverging))
        cmapChooser.grid(row =3, column = 1, sticky = Tk.W + Tk.E)


        # Make an entry to change the number of columns
        self.columnNum = Tk.StringVar(self)
        self.columnNum.set(self.parent.MainParamDict['NumOfCols']) # default value
        self.columnNum.trace('w', self.ColumnNumChanged)
        ttk.Label(frm, text="# of columns:").grid(row=4)
        self.ColumnSpin = Spinbox(frm,  from_=1, to=self.parent.MainParamDict['MaxCols'], textvariable=self.columnNum, width = 6)
        self.ColumnSpin.grid(row =4, column = 1, sticky = Tk.W + Tk.E)

        # Make an entry to change the number of columns
        self.rowNum = Tk.StringVar(self)
        self.rowNum.set(self.parent.MainParamDict['NumOfRows']) # default value
        self.rowNum.trace('w', self.RowNumChanged)
        ttk.Label(frm, text="# of rows:").grid(row=5)
        self.RowSpin = Spinbox(frm, from_=1, to=self.parent.MainParamDict['MaxRows'], textvariable=self.rowNum, width = 6)
        self.RowSpin.grid(row =5, column = 1, sticky = Tk.W + Tk.E)

        self.PrtlStrideVar = Tk.StringVar()
        self.PrtlStrideVar.set(str(self.parent.MainParamDict['PrtlStride']))
        ttk.Entry(frm, textvariable = self.PrtlStrideVar, width =6).grid(row =6, column =1, sticky = Tk.W +Tk.E)
        ttk.Label(frm, text='Particle stride').grid(row= 6,column =0)

        # Control whether or not Title is shown
        self.TitleVar = Tk.IntVar()
        self.TitleVar.set(self.parent.MainParamDict['ShowTitle'])
        self.TitleVar.trace('w', self.TitleChanged)

        self.LimVar = Tk.IntVar()
        self.LimVar.set(self.parent.MainParamDict['SetxLim'])
        self.LimVar.trace('w', self.LimChanged)

        self.xleft = Tk.StringVar()
        self.xleft.set(str(self.parent.MainParamDict['xLeft']))
        self.xright = Tk.StringVar()
        self.xright.set(str(self.parent.MainParamDict['xRight']))


        ttk.Label(frm, text = 'min').grid(row= 7, column = 1, sticky = Tk.N)
        ttk.Label(frm, text = 'max').grid(row= 7, column = 2, sticky = Tk.N)
        cb = ttk.Checkbutton(frm, text ='Set xlim',
                        variable = self.LimVar)
        cb.grid(row = 8, sticky = Tk.N)
        ttk.Entry(frm, textvariable=self.xleft, width = 8).grid(row = 8, column =1, sticky = Tk.N)
        ttk.Entry(frm, textvariable=self.xright, width = 8).grid(row = 8, column =2, sticky = Tk.N)



        self.yLimVar = Tk.IntVar()
        self.yLimVar.set(self.parent.MainParamDict['SetyLim'])
        self.yLimVar.trace('w', self.yLimChanged)



        self.yleft = Tk.StringVar()
        self.yleft.set(str(self.parent.MainParamDict['yBottom']))
        self.yright = Tk.StringVar()
        self.yright.set(str(self.parent.MainParamDict['yTop']))


        ttk.Checkbutton(frm, text ='Set ylim',
                        variable = self.yLimVar).grid(row = 9, sticky = Tk.N)
        ttk.Entry(frm, textvariable=self.yleft, width = 8 ).grid(row = 9, column =1, sticky = Tk.N)
        ttk.Entry(frm, textvariable=self.yright, width =8 ).grid(row = 9, column =2, sticky = Tk.N)

        self.kLimVar = Tk.IntVar()
        self.kLimVar.set(self.parent.MainParamDict['SetkLim'])
        self.kLimVar.trace('w', self.kLimChanged)



        self.kleft = Tk.StringVar()
        self.kleft.set(str(self.parent.MainParamDict['kLeft']))
        self.kright = Tk.StringVar()
        self.kright.set(str(self.parent.MainParamDict['kRight']))


        ttk.Checkbutton(frm, text ='Set klim', variable = self.kLimVar).grid(row = 10, sticky = Tk.N)
        ttk.Entry(frm, textvariable=self.kleft, width = 8 ).grid(row = 10, column =1, sticky = Tk.N)
        ttk.Entry(frm, textvariable=self.kright, width =8 ).grid(row = 10, column =2, sticky = Tk.N)

        self.xRelVar = Tk.IntVar()
        self.xRelVar.set(self.parent.MainParamDict['xLimsRelative'])
        self.xRelVar.trace('w', self.xRelChanged)
        ttk.Checkbutton(frm, text = "x limits & zooms relative to shock",
                        variable = self.xRelVar).grid(row = 11, columnspan = 3, sticky = Tk.W)

        framecb = ttk.Frame(frm)

        ttk.Label(framecb, text='Choose 2D plane:').pack(side = Tk.LEFT, expand = 0)
        self.PlaneVar = Tk.IntVar()
        self.PlaneVar.set(self.parent.MainParamDict['2DSlicePlane'])
        self.xybutton = ttk.Radiobutton(framecb,
                            text='x-y',
                            variable=self.PlaneVar,
                            command = self.RadioPlane,
                            value=0)
        self.xybutton.pack(side = Tk.LEFT, expand = 0)
        self.xzbutton = ttk.Radiobutton(framecb,
                            text='x-z',
                            variable=self.PlaneVar,
                            command = self.RadioPlane,
                            value=1)
        self.xzbutton.pack(side = Tk.LEFT, expand = 0)
        framecb.grid(row = 12, columnspan = 4)
        framey = ttk.Frame(frm)
        self.ySliceVar = Tk.IntVar()
        self.ySliceVar.set(self.parent.ySlice)
        self.units_listy = []
        for i in range(self.parent.MaxYInd+1):
            self.units_listy.append(str(i*self.parent.istep/self.parent.c_omp))

        self.ySliceVarC_omp = Tk.StringVar()
        self.ySliceVarC_omp.set(self.units_listy[self.ySliceVar.get()])

        labely = ttk.Label(framey, text='y-slice')#
        labely.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=0)


        # A slider that will select the 2D slice in the simulation
        self.slidery = ttk.Scale(framey, from_=0, to=self.parent.MaxYInd, command = self.yScaleHandler)
        self.slidery.set(self.ySliceVar.get())
        self.slidery.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)


        self.txtEntery = ttk.Entry(framey, textvariable=self.ySliceVarC_omp, width=6)
        self.txtEntery.pack(side=Tk.LEFT, fill = Tk.BOTH, expand = 0)
        if self.parent.MaxYInd ==0:
            self.txtEntery.state(['disabled'])
            self.slidery.state(['disabled'])
        ttk.Label(framey, text='[c_omp]').pack(side=Tk.LEFT, fill=Tk.BOTH, expand=0)
        # bind releasing the moust button to updating the plots.
        self.slidery.bind("<ButtonRelease-1>", self.yUpdateValue)


        framey.grid(row = 13, columnspan =4)

        framez = ttk.Frame(frm)
        self.zSliceVar = Tk.IntVar()
        self.zSliceVar.set(int(np.around(self.parent.MainParamDict['zSlice']*self.parent.MaxZInd)))

        self.units_listz = []
        for i in range(self.parent.MaxZInd+1):
            self.units_listz.append(str(i*self.parent.istep/self.parent.c_omp))

        self.zSliceVarC_omp = Tk.StringVar()
        self.zSliceVarC_omp.set(self.units_listz[self.zSliceVar.get()])

        # An entry box that will let us choose the time-step
        ttk.Label(framez, text='z-slice').pack(side=Tk.LEFT, fill=Tk.BOTH, expand=0)

        # A slider that will select the 2D slice in the simulation
        self.sliderz = ttk.Scale(framez, from_=0, to=self.parent.MaxZInd, command = self.zScaleHandler)
        self.sliderz.set(self.zSliceVar.get())
        self.sliderz.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)

        self.txtEnterz = ttk.Entry(framez, textvariable=self.zSliceVarC_omp, width=6)
        self.txtEnterz.pack(side=Tk.LEFT, fill = Tk.BOTH, expand = 0)
        ttk.Label(framez, text='[c_omp]').pack(side=Tk.LEFT, fill=Tk.BOTH, expand=0)
        # bind releasing the moust button to updating the plots.
        self.sliderz.bind("<ButtonRelease-1>", self.zUpdateValue)
        if self.parent.MaxZInd ==0:
            self.xzbutton.state(['disabled'])
            self.txtEnterz.state(['disabled'])
            self.sliderz.state(['disabled'])


        framez.grid(row = 14, columnspan =4)

        cb = ttk.Checkbutton(frm, text = "Show Title",
                        variable = self.TitleVar)
        cb.grid(row = 15, sticky = Tk.W)
        # Control whether or not axes are shared with a radio box:
        self.toLinkList = ['None', 'All spatial', 'All non p-x', 'All 2-D spatial']
        self.LinkedVar = Tk.IntVar()
        self.LinkedVar.set(self.parent.MainParamDict['LinkSpatial'])

        ttk.Label(frm, text='Share spatial axes:').grid(row = 0, column = 2, sticky = Tk.W)

        for i in range(len(self.toLinkList)):
            ttk.Radiobutton(frm,
                    text=self.toLinkList[i],
                    variable=self.LinkedVar,
                    command = self.RadioLinked,
                    value=i).grid(row = 1+i, column = 2, sticky =Tk.N)

        self.AspectVar = Tk.IntVar()
        self.AspectVar.set(self.parent.MainParamDict['ImageAspect'])
        self.AspectVar.trace('w', self.AspectVarChanged)

        cb = ttk.Checkbutton(frm, text = "Aspect = 1",
                                variable = self.AspectVar)
        cb.grid(row = 15, column = 1, sticky = Tk.W)

        self.ConstantShockVar = Tk.IntVar()
        self.ConstantShockVar.set(self.parent.MainParamDict['ConstantShockVel'])
        self.ConstantShockVar.trace('w', self.ShockSpeedVarChanged)

        cb = ttk.Checkbutton(frm, text = "Constant Shock v",
                                variable = self.ConstantShockVar)
        cb.grid(row = 15, column = 2, sticky = Tk.W)

        self.Average1DVar = Tk.IntVar()
        self.Average1DVar.set(self.parent.MainParamDict['Average1D'])
        self.Average1DVar.trace('w', self.AverageChanged)
        ttk.Checkbutton(frm, text='1D Average',variable = self.Average1DVar).grid(row = 16, column = 2, sticky = Tk.W)

        self.CbarOrientation = Tk.IntVar()
        self.CbarOrientation.set(self.parent.MainParamDict['HorizontalCbars'])
        self.CbarOrientation.trace('w', self.OrientationChanged)

        cb = ttk.Checkbutton(frm, text = "Horizontal Cbars",
                                variable = self.CbarOrientation)
        cb.grid(row = 16, sticky = Tk.W)


        self.LinkKVar = Tk.IntVar()
        self.LinkKVar.set(self.parent.MainParamDict['LinkK'])
        self.LinkKVar.trace('w', self.LinkKChanged)

        cb = ttk.Checkbutton(frm, text = "Share k-axes",
                                variable = self.LinkKVar)
        cb.grid(row = 16, column =1, sticky = Tk.W)



        self.LorentzBoostVar = Tk.IntVar()
        self.LorentzBoostVar.set(self.parent.MainParamDict['DoLorentzBoost'])
        self.LorentzBoostVar.trace('w', self.LorentzBoostChanged)
        cb = ttk.Checkbutton(frm, text='Boost PhasePlots', variable =  self.LorentzBoostVar).grid(row = 17, sticky = Tk.W)
        ttk.Label(frm, text='Gamma/Beta = \r (- for left boost)').grid(row= 17, rowspan = 2,column =1, sticky = Tk.E)
        self.GammaVar = Tk.StringVar()
        self.GammaVar.set(str(self.parent.MainParamDict['GammaBoost']))
        ttk.Entry(frm, textvariable=self.GammaVar, width = 7).grid(row = 17, column = 2, sticky = Tk.N)

    def yScaleHandler(self, e):
        # if changing the scale will change the value of the parameter, do so
        if self.ySliceVar.get() != int(self.slidery.get()):
            self.ySliceVar.set(int(self.slidery.get()))
            self.ySliceVarC_omp.set(self.units_listy[self.ySliceVar.get()])

    def zScaleHandler(self, e):
        # if changing the scale will change the value of the parameter, do so
        if self.zSliceVar.get() != int(self.sliderz.get()):
            self.zSliceVar.set(int(self.sliderz.get()))
            self.zSliceVarC_omp.set(self.units_listz[self.zSliceVar.get()])

    def zUpdateValue(self, e):
        if self.zSliceVar.get() == self.parent.zSlice:
            pass

        else:
            self.parent.MainParamDict['zSlice'] = float(self.zSliceVar.get())/self.parent.MaxZInd
            self.zSliceVarC_omp.set(self.units_listz[self.zSliceVar.get()])
            self.parent.RenewCanvas()

    def yUpdateValue(self, e):
        if self.ySliceVar.get() == self.parent.ySlice:
            pass

        else:
            self.parent.MainParamDict['ySlice'] = float(self.ySliceVar.get())/self.parent.MaxYInd
            self.ySliceVarC_omp.set(self.units_listy[self.ySliceVar.get()])
            self.parent.RenewCanvas()


    def AspectVarChanged(self, *args):
        if self.AspectVar.get() == self.parent.MainParamDict['ImageAspect']:
            pass

        else:
            self.parent.MainParamDict['ImageAspect'] = self.AspectVar.get()
            self.parent.RenewCanvas(ForceRedraw = True)


    def ShockSpeedVarChanged(self, *args):
        if self.parent.MainParamDict['ConstantShockVel'] != self.ConstantShockVar.get():
            self.parent.MainParamDict['ConstantShockVel'] = self.ConstantShockVar.get()
            self.parent.RenewCanvas(ForceRedraw = True)
    def AverageChanged(self, *args):
        if self.parent.MainParamDict['Average1D'] != self.Average1DVar.get():
            self.parent.MainParamDict['Average1D'] = self.Average1DVar.get()
            self.parent.RenewCanvas()

    def OrientationChanged(self, *args):
        if self.CbarOrientation.get() == self.parent.MainParamDict['HorizontalCbars']:
            pass

        else:
            if self.CbarOrientation.get():
                self.parent.axes_extent = self.parent.MainParamDict['HAxesExtent']
                self.parent.cbar_extent = self.parent.MainParamDict['HCbarExtent']
                self.parent.SubPlotParams = self.parent.MainParamDict['HSubPlotParams']

            else:
                self.parent.axes_extent = self.parent.MainParamDict['VAxesExtent']
                self.parent.cbar_extent = self.parent.MainParamDict['VCbarExtent']
                self.parent.SubPlotParams = self.parent.MainParamDict['VSubPlotParams']
            self.parent.MainParamDict['HorizontalCbars'] = self.CbarOrientation.get()
            self.parent.f.subplots_adjust( **self.parent.SubPlotParams)
            self.parent.RenewCanvas(ForceRedraw=True)

    def LorentzBoostChanged(self, *args):
        if self.LorentzBoostVar.get() == self.parent.MainParamDict['DoLorentzBoost']:
            pass

        else:
            self.parent.MainParamDict['DoLorentzBoost'] = self.LorentzBoostVar.get()
            self.parent.RenewCanvas()

    def TitleChanged(self, *args):
        if self.TitleVar.get()==self.parent.MainParamDict['ShowTitle']:
            pass
        else:
            self.parent.MainParamDict['ShowTitle'] = self.TitleVar.get()
            if self.TitleVar.get() == False:
                self.parent.f.suptitle('')

            self.parent.RenewCanvas()

    def RadioLinked(self, *args):
        # If the shared axes are changed, the whole plot must be redrawn
        if self.LinkedVar.get() == self.parent.MainParamDict['LinkSpatial']:
            pass
        else:
            self.parent.MainParamDict['LinkSpatial'] = self.LinkedVar.get()
            self.parent.RenewCanvas(ForceRedraw = True)
    def RadioPlane(self, *args):
        # If the shared axes are changed, the whole plot must be redrawn
        if self.PlaneVar.get() == self.parent.MainParamDict['2DSlicePlane']:
            pass
        else:
            self.parent.MainParamDict['2DSlicePlane'] = self.PlaneVar.get()
            self.parent.RenewCanvas(    )


    def LinkKChanged(self, *args):
        # If the shared axes are changed, the whole plot must be redrawn
        if self.LinkKVar.get() == self.parent.MainParamDict['LinkK']:
            pass
        else:
            self.parent.MainParamDict['LinkK'] = self.LinkKVar.get()
            self.parent.RenewCanvas(ForceRedraw = True)

    def xRelChanged(self, *args):
        # If the shared axes are changed, the whole plot must be redrawn
        if self.xRelVar.get() == self.parent.MainParamDict['xLimsRelative']:
            pass
        else:
            self.parent.MainParamDict['xLimsRelative'] = self.xRelVar.get()
            self.parent.RenewCanvas()


    def CmapChanged(self, *args):
    # Note here that Tkinter passes an event object to onselect()
        if self.cmapvar.get() == self.parent.MainParamDict['ColorMap']:
            pass
        else:
            self.parent.MainParamDict['ColorMap'] = self.cmapvar.get()
            if self.parent.MainParamDict['ColorMap'] in self.parent.cmaps_with_green:
                self.parent.ion_color = "#{0:02x}{1:02x}{2:02x}".format(int(np.round(new_cmaps.cmaps['plasma'](0.55)[0]*255)), int(np.round(new_cmaps.cmaps['plasma'](0.55)[1]*255)), int(np.round(new_cmaps.cmaps['plasma'](0.55)[2]*255)))
                self.parent.electron_color ="#{0:02x}{1:02x}{2:02x}".format(int(np.round(new_cmaps.cmaps['plasma'](0.8)[0]*255)), int(np.round(new_cmaps.cmaps['plasma'](0.8)[1]*255)), int(np.round(new_cmaps.cmaps['plasma'](0.8)[2]*255)))

                self.parent.ion_fit_color = 'r'
                self.parent.electron_fit_color = 'yellow'

            else:
                self.parent.ion_color = "#{0:02x}{1:02x}{2:02x}".format(int(np.round(new_cmaps.cmaps['viridis'](0.45)[0]*255)), int(np.round(new_cmaps.cmaps['viridis'](0.45)[1]*255)), int(np.round(new_cmaps.cmaps['viridis'](0.45)[2]*255)))
                self.parent.electron_color ="#{0:02x}{1:02x}{2:02x}".format(int(np.round(new_cmaps.cmaps['viridis'](0.75)[0]*255)), int(np.round(new_cmaps.cmaps['viridis'](0.75)[1]*255)), int(np.round(new_cmaps.cmaps['viridis'](0.75)[2]*255)))

                self.parent.ion_fit_color = 'mediumturquoise'
                self.parent.electron_fit_color = 'lime'


            self.parent.RenewCanvas(ForceRedraw = True)

    def DivCmapChanged(self, *args):
    # Note here that Tkinter passes an event object to onselect()
        if self.div_cmapvar.get() == self.parent.MainParamDict['DivColorMap']:
            pass
        else:
            self.parent.MainParamDict['DivColorMap'] = self.div_cmapvar.get()
            self.parent.RenewCanvas(ForceRedraw = True)


    def SkipSizeChanged(self, *args):
    # Note here that Tkinter passes an event object to SkipSizeChange()
        try:
            if self.skipSize.get() == '':
                pass
            else:
                self.parent.MainParamDict['SkipSize'] = int(self.skipSize.get())
        except ValueError:
            self.skipSize.set(self.parent.MainParamDict['SkipSize'])

    def RowNumChanged(self, *args):
        try:
            if self.rowNum.get() == '':
                pass
            if int(self.rowNum.get())<1:
                self.rowNum.set(1)
            if int(self.rowNum.get())>self.parent.MainParamDict['MaxRows']:
                self.rowNum.set(self.parent.MainParamDict['MaxRows'])
            if int(self.rowNum.get()) != self.parent.MainParamDict['NumOfRows']:
                self.parent.MainParamDict['NumOfRows'] = int(self.rowNum.get())
                self.parent.UpdateGridSpec()
        except ValueError:
            self.rowNum.set(self.parent.MainParamDict['NumOfRows'])

    def ColumnNumChanged(self, *args):
        try:
            if self.columnNum.get() == '':
                pass
            if int(self.columnNum.get())<1:
                self.columnNum.set(1)
            if int(self.columnNum.get())>self.parent.MainParamDict['MaxCols']:
                self.columnNum.set(self.parent.MainParamDict['MaxCols'])
            if int(self.columnNum.get()) != self.parent.MainParamDict['NumOfCols']:
                self.parent.MainParamDict['NumOfCols'] = int(self.columnNum.get())
                self.parent.UpdateGridSpec()

        except ValueError:
            self.columnNum.set(self.parent.MainParamDict['NumOfCols'])

    def WaitTimeChanged(self, *args):
    # Note here that Tkinter passes an event object to onselect()
        try:
            if self.waitTime.get() == '':
                pass
            else:
                self.parent.MainParamDict['WaitTime'] = float(self.waitTime.get())
        except ValueError:
            self.waitTime.set(self.parent.MainParamDict['WaitTime'])

    def CheckIfLimsChanged(self):
        to_reload = False
        tmplist = [self.xleft, self.xright, self.yleft, self.yright, self.kleft, self.kright]
        limkeys = ['xLeft', 'xRight', 'yBottom', 'yTop', 'kLeft', 'kRight']
        setKeys = ['SetxLim', 'SetyLim', 'SetkLim']
        for j in range(6):
            setlims = self.parent.MainParamDict[setKeys[j//2]]
            tmpkey = limkeys[j]

            try:
            #make sure the user types in a a number and that it has changed.
                if np.abs(float(tmplist[j].get()) - self.parent.MainParamDict[tmpkey]) > 1E-4:
                    self.parent.MainParamDict[tmpkey] = float(tmplist[j].get())
                    to_reload += setlims

            except ValueError:
                #if they type in random stuff, just set it ot the param value
                tmplist[j].set(str(self.parent.MainParamDict[tmpkey]))
        return to_reload

    def CheckIfGammaChanged(self):
        to_reload = False
        try:
        #make sure the user types in a float
            if np.abs(float(self.GammaVar.get()) - self.parent.MainParamDict['GammaBoost']) > 1E-8:
                self.parent.MainParamDict['GammaBoost'] = float(self.GammaVar.get())
                to_reload += True

        except ValueError:
            #if they type in random stuff, just set it to the param value
            self.GammaVar.set(str(self.parent.MainParamDict['GammaBoost']))
        return to_reload*self.parent.MainParamDict['DoLorentzBoost']

    def CheckIfStrideChanged(self):
        to_reload = False

        try:
            #make sure the user types in a int
            if int(self.PrtlStrideVar.get()) <= 0:
                self.PrtlStrideVar.set(str(self.parent.MainParamDict['PrtlStride']))
            if int(self.PrtlStrideVar.get()) != self.parent.MainParamDict['PrtlStride']:
                self.parent.MainParamDict['PrtlStride'] = int(self.PrtlStrideVar.get())
                self.parent.stride = self.parent.MainParamDict['PrtlStride']
                self.parent.StrideChanged()
                to_reload += True

        except ValueError:
            #if they type in random stuff, just set it to the param value
            self.PrtlStrideVar.set(str(self.parent.MainParamDict['PrtlStride']))
        return to_reload

    def CheckIfSliceChanged(self):
        to_reload = False
        try:
            #make sure the user types in a float
            self.ySliceVar.set(int(np.around(float(self.ySliceVarC_omp.get())*self.parent.c_omp/self.parent.istep)))
            if int(self.ySliceVar.get()) < 0:
                self.ySliceVar.set(0)

            elif int(self.ySliceVar.get()) > self.parent.MaxYInd:
                self.ySliceVar.set(self.parent.MaxYInd)
            self.ySliceVarC_omp.set(self.units_listy[self.ySliceVar.get()])
            if self.ySliceVar.get() != int(np.around(self.parent.MainParamDict['ySlice']*self.parent.MaxYInd)):
                self.parent.MainParamDict['ySlice'] = float(self.ySliceVar.get())/self.parent.MaxYInd
                self.slidery.set(self.ySliceVar.get())
                to_reload += True
        except ValueError:
            #if they type in random stuff, just set it to the param value
            self.ySliceVarC_omp.set(self.units_listy[self.ySliceVar.get()])

        try:
            #make sure the user types in a float
            self.zSliceVar.set(int(np.around(float(self.zSliceVarC_omp.get())*self.parent.c_omp/self.parent.istep)))
            if int(self.zSliceVar.get()) < 0:
                self.zSliceVar.set(0)

            elif int(self.zSliceVar.get()) > self.parent.MaxZInd:
                self.zSliceVar.set(self.parent.MaxZInd)
            self.zSliceVarC_omp.set(self.units_listz[self.zSliceVar.get()])
            if self.zSliceVar.get() != int(np.around(self.parent.MainParamDict['zSlice']*self.parent.MaxZInd)):
                self.parent.MainParamDict['zSlice'] = float(self.OneDSliceVar.get())/self.parent.MaxZInd
                self.sliderz.set(self.zSliceVar.get())
                to_reload += True

        except ValueError:
            #if they type in random stuff, just set it to the param value
            self.TwoDSliceVarC_omp.set(self.units_list2D[self.TwoDSliceVar.get()])
        return to_reload


    def LimChanged(self, *args):
        if self.LimVar.get()==self.parent.MainParamDict['SetxLim']:
            pass
        else:
            self.parent.MainParamDict['SetxLim'] = self.LimVar.get()
            self.parent.RenewCanvas()

    def yLimChanged(self, *args):
        if self.yLimVar.get()==self.parent.MainParamDict['SetyLim']:
            pass
        else:
            self.parent.MainParamDict['SetyLim'] = self.yLimVar.get()
            self.parent.RenewCanvas()

    def kLimChanged(self, *args):
        if self.kLimVar.get()==self.parent.MainParamDict['SetkLim']:
            pass
        else:
            self.parent.MainParamDict['SetkLim'] = self.kLimVar.get()
            self.parent.RenewCanvas()


    def SettingsCallback(self, e):
        to_reload = self.CheckIfLimsChanged()
        to_reload += self.CheckIfGammaChanged()
        to_reload += self.CheckIfStrideChanged()
        to_reload += self.CheckIfSliceChanged()
        if to_reload:
            self.parent.RenewCanvas()



    def OnClosing(self):
        self.parent.settings_window = None
        self.destroy()

class MeasureFrame(Tk.Toplevel):
    def __init__(self, parent):

        Tk.Toplevel.__init__(self)
        self.wm_title('Take Measurements')
        self.protocol('WM_DELETE_WINDOW', self.OnClosing)


        self.parent = parent

        self.bind('<Return>', self.TxtEnter)
        frm = ttk.Frame(self)
        frm.pack(fill=Tk.BOTH, expand=True)


        ttk.Label(frm, text='NOTE: Spectral Measurements have been moved ' +'\r' + 'to the spectral subplot settings window.').grid(row = 0, rowspan = 2,columnspan = 3, sticky = Tk.W)



        # Make an entry to change the integration region
        # A StringVar for a box to type in a value for the left ion region
        self.FFTLVar = Tk.StringVar()
        # set it to the left value
        self.FFTLVar.set(str(self.parent.MainParamDict['FFTLeft']))

        # A StringVar for a box to type in a value for the right ion region
        self.FFTRVar = Tk.StringVar()
        # set it to the right value
        self.FFTRVar.set(str(self.parent.MainParamDict['FFTRight']))

        ttk.Label(frm, text='left').grid(row = 2, column = 1, sticky = Tk.N)
        ttk.Label(frm, text='right').grid(row = 2, column = 2, sticky = Tk.N)

        ttk.Label(frm, text='FFT region:').grid(row = 3, sticky = Tk.W)
        ttk.Entry(frm, textvariable=self.FFTLVar, width=7).grid(row =3, column = 1, sticky = Tk.W + Tk.E)

        ttk.Entry(frm, textvariable=self.FFTRVar, width=7).grid(row = 3, column =2, sticky = Tk.W + Tk.E)

        self.FFTRelVar = Tk.IntVar()
        self.FFTRelVar.set(self.parent.MainParamDict['FFTRelative'])
        self.FFTRelVar.trace('w', self.FFTRelChanged)
        cb = ttk.Checkbutton(frm, text = "FFT Region relative to shock?",
                        variable = self.FFTRelVar)
        cb.grid(row = 4, columnspan = 3, sticky = Tk.W)


    def CheckIfFloatChanged(self, tkVar, paramKey):
        to_reload = False
        try:
            #make sure the user types in a int
            if np.abs(float(tkVar.get())- self.parent.MainParamDict[paramKey])>1E-6:
                self.parent.MainParamDict[paramKey] = float(tkVar.get())
                to_reload = True
            return to_reload

        except ValueError:
            #if they type in random stuff, just set it to the param value
            tkVar.set(str(self.parent.MainParamDict[paramKey]))
            return to_reload

    def TxtEnter(self, e):
        self.MeasuresCallback()


    def FFTRelChanged(self, *args):
        if self.FFTRelVar.get()==self.parent.MainParamDict['FFTRelative']:
            pass
        else:
            self.parent.MainParamDict['FFTRelative'] = self.FFTRelVar.get()
            self.parent.RenewCanvas()



    def MeasuresCallback(self):
        tkvarIntList = [self.FFTLVar, self.FFTRVar]
        IntValList = ['FFTLeft', 'FFTRight']

        to_reload = False

        for j in range(len(tkvarIntList)):
            to_reload += self.CheckIfFloatChanged(tkvarIntList[j], IntValList[j])

        if to_reload:
            self.parent.RenewCanvas()

    def OnClosing(self):
        self.parent.settings_window = None
        self.destroy()


class MainApp(Tk.Tk):
    """ We simply derive a new class of Frame as the man frame of our app"""
    def __init__(self, name,cmd_args):

        Tk.Tk.__init__(self)
        self.update_idletasks()
        menubar = Tk.Menu(self)
        self.wm_title(name)
        self.settings_window = None
        self.measure_window = None


        self.cmd_args = cmd_args
#        if self.cmd_args.r:
#            self.iconify()
        # A variable that keeps track of the first graph with spatial x & y axes
        self.first_x = None
        self.first_y = None

        # An int that stores the current stride
        self.stride = 0

        self.IseultDir = os.path.join(os.path.dirname(__file__),'..')

        # a list of cmaps with orange prtl colors
        self.cmaps_with_green = ['viridis', 'Rainbow + White', 'Blue/Green/Red/Yellow', 'Cube YF', 'Linear_L']



        fileMenu = Tk.Menu(menubar, tearoff=False)
        self.presetMenu = Tk.Menu(menubar, tearoff=False, postcommand=self.ViewUpdate)
        menubar.add_cascade(label="File", underline=0, menu=fileMenu)
        fileMenu.add_command(label= 'Open Directory', command = self.OnOpen, accelerator='Command+o')

        fileMenu.add_command(label="Exit", underline=1,
                             command=quit, accelerator="Ctrl+Q")
        fileMenu.add_command(label= 'Save Current State', command = self.OpenSaveDialog)
        fileMenu.add_command(label= 'Make a Movie', command = self.OpenMovieDialog)
        fileMenu.add_command(label= 'Reset Session', command = self.ResetSession)


        self.bind_all("<Control-q>", self.quit)
        self.bind_all("<Command-o>", self.OnOpen)
        self.bind_all("S", self.OpenSettings)

        # create a bunch of regular expressions used to search for files
        f_re = re.compile('flds.tot.*')
        prtl_re = re.compile('prtl.tot.*')
        s_re = re.compile('spect.*')
        param_re = re.compile('param.*')
        self.re_list = [f_re, prtl_re, s_re, param_re]


        # A list that will keep track of whether a given axes is a colorbar or not:
        self.cbarList = []

        # The dictionary that holdsd the paths
        self.PathDict = {'Flds': [], 'Prtl': [], 'Param': [], 'Spect': []}

        # A dictionary that allows use to see in what HDF5 file each key is stored.
        # i.e. {'ui': 'Prtl', 'ue': 'Flds', etc...},  Originally I generated the
        # key dictionary automatically, but I don't think that is safe anymore.

        self.H5KeyDict = {u'mx0': 'Param',
                          u'teststarti': 'Param',
                          u'teststartl': 'Param',
                          u'sizex': 'Param',
                          u'sizey': 'Param',
                          u'c_omp': 'Param',
                          u'qi': 'Param',
                          u'istep1': 'Param',
                          u'my0': 'Param',
                          u'dlapion': 'Param',
                          u'testendion': 'Param',
                          u'caseinit': 'Param',
                          u'pltstart': 'Param',
                          u'stride': 'Param',
                          u'ntimes': 'Param',
                          u'cooling': 'Param',
                          u'btheta': 'Param',
                          u'c': 'Param',
                          u'acool': 'Param',
                          u'istep': 'Param',
                          u'delgam': 'Param',
                          u'me': 'Param',
                          u'dlaplec': 'Param',
                          u'mi': 'Param',
                          u'torqint': 'Param',
                          u'mx': 'Param',
                          u'mz0': 'Param',
                          u'yi': 'Prtl',
                          u'proci': 'Prtl',
                          u'proce': 'Prtl',
                          u'ye': 'Prtl',
                          u'zi': 'Prtl',
                          u'ze': 'Prtl',
                          u'xsl': 'Spect',
                          u'umean': 'Spect',
                          u'spece': 'Spect',
                          u'v3xi': 'Flds',
                          u'ey': 'Flds',
                          u'ex': 'Flds',
                          u'ez': 'Flds',
                          u'specp': 'Spect',
                          u'densi': 'Flds',
                          u'specprest': 'Spect',
                          u'we': 'Prtl',
                          u'jx': 'Flds',
                          u'jy': 'Flds',
                          u'jz': 'Flds',
                          u'gmax': 'Spect',
                          u'gmin': 'Spect',
                          'spect_dens': 'Spect',
                          u'wi': 'Prtl',
                          u'bx': 'Flds',
                          u'by': 'Flds',
                          u'bz': 'Flds',
                          u'dgam': 'Spect',
                          u'gamma': 'Spect',
                          u'xi': 'Prtl',
                          u'xe': 'Prtl',
                          u'che': 'Prtl',
                          u'chi': 'Prtl',
                          u'ui': 'Prtl',
                          u'ue': 'Prtl',
                          u've': 'Prtl',
                          u'gamma0': 'Param',
                          u'vi': 'Prtl',
                          u'my': 'Param',
                          u'specerest': 'Spect',
                          u'v3yi': 'Flds',
                          u'walloc': 'Param',
                          u'testendlec': 'Param',
                          u'v3x': 'Flds',
                          u'v3y': 'Flds',
                          u'v3z': 'Flds',
                          u'xinject2': 'Param',
                          u'gammae': 'Prtl',
                          u'bphi': 'Param',
                          u'gammai': 'Prtl',
                          u'dummy': 'Param',
                          u'dens': 'Flds',
                          u'sigma': 'Param',
                          u'interval': 'Param',
                          u'inde': 'Prtl',
                          u'v3zi': 'Flds',
                          u'time': 'Param',
                          u'splitratio': 'Param',
                          u'indi': 'Prtl',
                          u'ppc0': 'Param'}
        self.prtl_keys = []
        for k, v in self.H5KeyDict.items():
            if v =='Prtl':
                self.prtl_keys.append(k)

        # Create the figure
        self.f = Figure(figsize = (2,2), dpi = 100, edgecolor = 'none', facecolor = 'w')

        # a tk.DrawingArea

        self.canvas = FigureCanvasTkAgg(self.f, master=self)

        self.GenMainParamDict()

        # now root.geometry() returns valid size/placement
        self.minsize(self.winfo_width(), self.winfo_height())
        self.geometry(self.MainParamDict['WindowSize'])

        if self.MainParamDict['HorizontalCbars']:
            self.axes_extent = self.MainParamDict['HAxesExtent']
            self.cbar_extent = self.MainParamDict['HCbarExtent']
            self.SubPlotParams = self.MainParamDict['HSubPlotParams']

        else:
            self.axes_extent = self.MainParamDict['VAxesExtent']
            self.cbar_extent = self.MainParamDict['VCbarExtent']
            self.SubPlotParams = self.MainParamDict['VSubPlotParams']
        self.f.subplots_adjust( **self.SubPlotParams)

        # Make the object hold the timestep info
        self.TimeStep = Param(1, minimum=1, maximum=1000)
        self.playbackbar = PlaybackBar(self, self.TimeStep, canvas = self.canvas)

        # Add the toolbar
        self.toolbar =  MyCustomToolbar(self.canvas, self)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=Tk.RIGHT, fill=Tk.BOTH, expand=1)

        # Some options to set the way the spectral lines are dashed
        self.dashes_options = [[],[3,1],[5,1],[1,1]]
        # Look for the tristan output files and load the file paths into
        # previous objects
        self.dirname = os.curdir
        if len(self.cmd_args.O[0])>0:
            self.dirname = os.path.join(self.dirname, self.cmd_args.O[0])

        self.findDir()


        self.TimeStep.attach(self)
        self.InitializeCanvas()

        menubar.add_cascade(label='Preset Views', underline=0, menu = self.presetMenu)
        self.playbackbar.pack(side=Tk.TOP, fill=Tk.BOTH, expand=0)
        self.update()


        self.config(menu=menubar)

        self.bind('<Return>', self.TxtEnter)
        self.bind('<Left>', self.playbackbar.SkipLeft)
        self.bind('<Right>', self.playbackbar.SkipRight)
        self.bind('r', self.playbackbar.OnReload)
        self.bind('<space>', self.playbackbar.PlayHandler)
        if self.cmd_args.b :
            self.after(0,self.MakeAMovie('out.mov', 1, -1, 1, 10))
            self.after(1, self.quit())
        self.update()
    def ViewUpdate(self):
        tmpdir = list(os.listdir(os.path.join(self.IseultDir, '.iseult_configs')))
        tmpdir.sort()
        for cfile in tmpdir:
            if cfile.split('.')[-1]=='yml':
                with open(os.path.join(os.path.join(self.IseultDir, '.iseult_configs'), cfile), 'r') as f:
                    cfgDict=yaml.safe_load(f)
                try:
                    if 'general' in cfgDict.keys():
                        if 'ConfigName'  in cfgDict['general'].keys():
                            tmpstr = cfgDict['general']['ConfigName']
                            try:
                                self.presetMenu.delete(tmpstr)
                            except:
                                pass
                            self.presetMenu.add_command(label = tmpstr, command = partial(self.LoadConfig, str(os.path.join(self.IseultDir,'.iseult_configs', cfile))))
                except:
                    pass
    def StrideChanged(self):
        # first we have to remove the calculated energy time steps
        self.TotalEnergyTimeSteps = []
        self.TotalEnergyTimes = np.array([])
        self.TotalIonEnergy = np.array([])
        self.TotalElectronEnergy = np.array([])

        self.TotalMagEnergy = np.array([])
        self.TotalBxEnergy = np.array([])
        self.TotalByEnergy = np.array([])
        self.TotalBzEnergy = np.array([])

        self.TotalExEnergy = np.array([])
        self.TotalEyEnergy = np.array([])
        self.TotalEzEnergy = np.array([])

        self.TotalElectricEnergy = np.array([])

        # figure out all keys that have 'Prtl'
        # now we have to go through the data dictionary and remove the particle info
        for DataDict in self.ListOfDataDict:
            for k in self.prtl_keys:
                DataDict.pop(k, None)


    def quit(self, event):
        print("quitting...")
        sys.exit(0)

    def GenMainParamDict(self, config_file = None):
        ''' The function that reads in a config file and then makes MainParamDict to hold all of the main iseult parameters.
            It also sets all of the plots parameters.'''

        #config = configparser.RawConfigParser()

        if config_file is None:
            try:
                with open(os.path.join(self.IseultDir, '.iseult_configs', self.cmd_args.p.strip().replace(' ', '_') +'.yml')) as f:
                    cfgDict = yaml.safe_load(f)
            except:
                print('Cannot find/load ' +  self.cmd_args.p.strip().replace(' ', '_') +'.yml in .iseult_configs. If the name of view contains whitespace,')
                print('either it must be enclosed in quotation marks or given with whitespace replaced by _.')
                print('Name is case sensitive. Reverting to Default view')
                with open(os.path.join(self.IseultDir, '.iseult_configs', 'Default.yml')) as f:
                    cfgDict = yaml.safe_load(f)
        else:
            with open(config_file) as f:
                cfgDict = yaml.safe_load(f)

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
                              'FFT_color': 'k',
                              'legendLabelSize':11}
        for key, val in cfgDict['MainParamDict'].items():
            self.MainParamDict[key] = val
        self.electron_color = self.MainParamDict['electron_color']
        self.ion_color = self.MainParamDict['ion_color']
        self.shock_color = self.MainParamDict['shock_color']
        self.ion_fit_color = self.MainParamDict['ion_fit_color']
        self.electron_fit_color = self.MainParamDict['electron_fit_color']
        self.FFT_color = self.MainParamDict['FFT_color']

        # if stride is 0 that means it has only been initialized... set to default
        if self.stride == 0:
            self.stride = self.MainParamDict['PrtlStride']

    def SaveIseultState(self, cfgfile, cfgname):
        #config = configparser.RawConfigParser()

        # When adding sections or items, add them in the reverse order of
        # how you want them to be displayed in the actual file.
        # In addition, please note that using RawConfigParser's and the raw
        # mode of ConfigParser's respective set functions, you can assign
        # non-string values to keys internally, but will receive an error
        # when attempting to write to a file or when you get it in non-raw
        # mode. SafeConfigParser does not allow such assignments to take place.
        #config.add_section('general')
        cfgDict = {}
        #config.set('general', 'ConfigName', cfgname)

        #config.add_section('main')
        cfgDict['general'] = {'ConfigName': cfgname}
        #DictList = ['HSubPlotParams', 'VSubPlotParams']
        #IntListsList = ['HAxesExtent', 'HCbarExtent', 'VAxesExtent', 'VCbarExtent']

        # Update the 'WindowSize' attribute to the current window size
        self.MainParamDict['WindowSize'] = str(self.winfo_width())+'x'+str(self.winfo_height())
        # Get figsize and dpi

        self.MainParamDict['FigSize'] = [float(self.f.get_size_inches()[0]),float(self.f.get_size_inches()[1])]

        #print(self.MainParamDict['FigSize'], type(self.MainParamDict['FigSize']))
        self.MainParamDict['dpi'] = self.f.dpi
        #print(self.MainParamDict['FigSize'], self.MainParamDict['dpi'] )
        # Update the current subplot params


        tmp_param_str = 'HSubPlotParams' if self.MainParamDict['HorizontalCbars'] else 'VSubPlotParams'
        try:
            self.MainParamDict[tmp_param_str]['left']=float(self.f.subplotpars.left)
            self.MainParamDict[tmp_param_str]['right']=float(self.f.subplotpars.right)
            self.MainParamDict[tmp_param_str]['top']=float(self.f.subplotpars.top)
            self.MainParamDict[tmp_param_str]['bottom']=float(self.f.subplotpars.bottom)
            self.MainParamDict[tmp_param_str]['wspace']=float(self.f.subplotpars.wspace)
            self.MainParamDict[tmp_param_str]['hspace']=float(self.f.subplotpars.hspace)

        except:
            pass
        cfgDict['MainParamDict'] = self.MainParamDict
        cfgDict['MainParamDict']['electron_color'] = self.electron_color
        cfgDict['MainParamDict']['ion_color'] = self.ion_color
        cfgDict['MainParamDict']['shock_color'] = self.shock_color
        cfgDict['MainParamDict']['ion_fit_color'] = self.ion_fit_color
        cfgDict['MainParamDict']['electron_fit_color'] = self.electron_fit_color
        cfgDict['MainParamDict']['FFT_color'] = self.FFT_color
        self.SaveLLoc()
        for i in range(self.MainParamDict['NumOfRows']):
            for j in range(self.MainParamDict['NumOfCols']):
                tmp_str = 'Chart' + str(i) + '_' + str(j)
                #config.add_section(tmp_str)
                tmp_ctype = self.SubPlotList[i][j].chartType
                #$config.set(tmp_str, 'ChartType', tmp_ctype)
                #for key in self.SubPlotList[i][j].PlotParamsDict[tmp_ctype].keys():
                    #config.set(tmp_str, key, str(s[key]))
                cfgDict[tmp_str] = self.SubPlotList[i][j].PlotParamsDict[tmp_ctype]
                cfgDict[tmp_str]['ChartType'] = tmp_ctype

        #print(yaml.dump(cfgDict))
        # Writing our configuration file to 'example.cfg'

        with open(cfgfile, 'w') as cfgFile:
            cfgFile.write(yaml.safe_dump(cfgDict))

    def GenH5Dict(self):
        '''Loads all of the files and then finds all of the keys in
        the file to load data. Deprecated'''
        for pkey in self.PathDict.keys():
            with h5py.File(os.path.join(self.dirname,self.PathDict[pkey][0]), 'r') as f:
                # Because dens is in both spect* files and flds* files,
                for h5key in f.keys():
                    if h5key == 'dens' and pkey == 'Spect':
                        self.H5KeyDict['spect_dens'] = pkey
                    else:
                        self.H5KeyDict[h5key] = pkey

        print(self.H5KeyDict)

    def ReloadPath(self):
        """ This function updates the current pathdictionary"""
        dirlist = os.listdir(self.dirname)

        if int(self.cmd_args.n)!=-1:
            self.CheckMaxNPopUp()

        # Create a dictionary of all the paths to the files
        self.PathDict = {'Flds': [], 'Prtl': [], 'Param': [], 'Spect': []}

        # create a bunch of regular expressions used to search for files
        f_re = re.compile('flds.tot.*')
        prtl_re = re.compile('prtl.tot.*')
        s_re = re.compile('spect.*')
        param_re = re.compile('param.*')
        self.PathDict['Flds']= list(filter(f_re.match, os.listdir(self.dirname)))
        self.PathDict['Flds'].sort()
        self.PathDict['Prtl']= list(filter(prtl_re.match, os.listdir(self.dirname)))
        self.PathDict['Prtl'].sort()
        self.PathDict['Spect']= list(filter(s_re.match, os.listdir(self.dirname)))
        self.PathDict['Spect'].sort()
        self.PathDict['Param']= list(filter(param_re.match, os.listdir(self.dirname)))
        self.PathDict['Param'].sort()

        ### iterate through the Paths and just get the .nnn number
        if len(self.PathDict['Param']) > 0:
            self.length_of_outfiles = len(self.PathDict['Param'][0].split('.')[-1])
        else:
            self.length_of_outfiles = 3

        for key in self.PathDict.keys():
            for i in range(len(self.PathDict[key])):
                try:
                    self.PathDict[key][i] = int(self.PathDict[key][i].split('.')[-1])
                except ValueError:
                    self.PathDict[key].pop(i)
                except IndexError:
                    pass

        ### GET THE NUMBERS THAT HAVE ALL 4 FILES:

        allFour = set(self.PathDict['Param'])
        for key in self.PathDict.keys():
            allFour &= set(self.PathDict[key])
        allFour = sorted(allFour)

        if int(self.cmd_args.n) != -1 and len(allFour)>0:
            while allFour[-1] > int(self.cmd_args.n) and len(allFour)>0:
                allFour.pop(-1)

        # Rebuild the pathdict only with files that have all 4 things
        self.PathDict['Flds'] = ['flds.tot.'+str(elm).zfill(self.length_of_outfiles) for elm in allFour]
        self.PathDict['Prtl'] = ['prtl.tot.'+str(elm).zfill(self.length_of_outfiles) for elm in allFour]
        self.PathDict['Spect'] = ['spect.'+str(elm).zfill(self.length_of_outfiles) for elm in allFour]
        self.PathDict['Param'] = ['param.'+str(elm).zfill(self.length_of_outfiles) for elm in allFour]
        self.TimeStep.setMax(len(self.PathDict['Flds']))
        self.playbackbar.slider.config(to =(len(self.PathDict['Flds'])))
        if self.MainParamDict['Reload2End']:
            self.TimeStep.value = len(self.PathDict['Flds'])
            self.playbackbar.slider.set(self.TimeStep.value)
        self.shock_finder()

    def CheckMaxNPopUp(self):
        MaxNDialog(self)
    def pathOK(self):
        """ Test to see if the current path contains tristan files
        using regular expressions, then generate the lists of files
        to iterate over"""
        dirlist = os.listdir(self.dirname)
        if 'output' in dirlist:
            self.dirname = os.path.join(self.dirname, 'output')

        is_okay = True

        # Create a dictionary of all the paths to the files
        self.PathDict = {'Flds': [], 'Prtl': [], 'Param': [], 'Spect': []}

        # create a bunch of regular expressions used to search for files
        f_re = re.compile('flds.tot.*')
        prtl_re = re.compile('prtl.tot.*')
        s_re = re.compile('spect.*')
        param_re = re.compile('param.*')



        self.PathDict['Flds']= list(filter(f_re.match, os.listdir(self.dirname)))
        self.PathDict['Flds'].sort()
        self.PathDict['Prtl']= list(filter(prtl_re.match, os.listdir(self.dirname)))
        self.PathDict['Prtl'].sort()
        self.PathDict['Spect']= list(filter(s_re.match, os.listdir(self.dirname)))
        self.PathDict['Spect'].sort()
        self.PathDict['Param']= list(filter(param_re.match, os.listdir(self.dirname)))
        self.PathDict['Param'].sort()

        ### iterate through the Paths and just get the .nnn number
        if len(self.PathDict['Param']) > 0:
            self.length_of_outfiles = len(self.PathDict['Param'][0].split('.')[-1])
        else:
            self.length_of_outfiles = 3

        for key in self.PathDict.keys():
            for i in range(len(self.PathDict[key])):
                try:
                    self.PathDict[key][i] = int(self.PathDict[key][i].split('.')[-1])
                except ValueError:
                    self.PathDict[key].pop(i)
                except IndexError:
                    pass

        ### GET THE NUMBERS THAT HAVE ALL 4 FILES:

        allFour = set(self.PathDict['Param'])
        for key in self.PathDict.keys():
            allFour &= set(self.PathDict[key])
        allFour = sorted(allFour)

        if int(self.cmd_args.n) != -1 and len(allFour)>0:
            while allFour[-1] > int(self.cmd_args.n) and len(allFour)>0:
                allFour.pop(-1)

        is_okay = len(allFour)>0
        # Rebuild the pathdict only with files that have all 4 things

        self.PathDict['Flds'] = ['flds.tot.'+str(elm).zfill(self.length_of_outfiles) for elm in allFour]
        self.PathDict['Prtl'] = ['prtl.tot.'+str(elm).zfill(self.length_of_outfiles) for elm in allFour]
        self.PathDict['Spect'] = ['spect.'+str(elm).zfill(self.length_of_outfiles) for elm in allFour]
        self.PathDict['Param'] = ['param.'+str(elm).zfill(self.length_of_outfiles) for elm in allFour]
        if is_okay:
            self.NewDirectory = True
            self.TimeStep.setMax(len(self.PathDict['Flds']))
            self.playbackbar.slider.config(to =(len(self.PathDict['Flds'])))
            if self.MainParamDict['Reload2End']:
                self.TimeStep.value = len(self.PathDict['Flds'])
                self.playbackbar.slider.set(self.TimeStep.value)
            self.shock_finder()
            self.movie_dir = ''

        return is_okay


    def OnOpen(self, e = None):
        """open a file"""

        if self.cmd_args.n != -1:
            self.CheckMaxNPopUp()

        tmpdir = filedialog.askdirectory(title = 'Choose the directory of the output files', **self.dir_opt)
        if tmpdir == '':
            self.findDir()

        else:
            self.dirname = tmpdir
        if not self.pathOK():
#            p = MyDalog(self, 'Directory must contain either the output directory or all of the following: \n flds.tot.*, ptrl.tot.*, params.*, spect.*', title = 'Cannot find output files')
#            self.wait_window(p.top)
            self.findDir()
        else:
            self.ReDrawCanvas()

    def ResetSession(self, e = None):
        """open a file"""
        if int(self.cmd_args.n) != -1:
            self.CheckMaxNPopUp()
        self.pathOK()
        self.ReDrawCanvas()

    def findDir(self, dlgstr = 'Choose the directory of the output files.'):
        """Look for /ouput folder, where the simulation results are
        stored. If output files are already in the path, they are
        automatically loaded"""
        # defining options for opening a directory
        self.dir_opt = {}
        self.dir_opt['initialdir'] = os.curdir
        self.dir_opt['mustexist'] = True
        self.dir_opt['parent'] = self

        if not self.pathOK():
            tmpdir = filedialog.askdirectory(title = dlgstr, **self.dir_opt)
            if tmpdir != '':
                self.dirname = tmpdir
            if not self.pathOK():
#                p = MyDialog(self, 'Directory must contain either the output directory or all of the following: \n flds.tot.*, ptrl.tot.*, params.*, spect.*', title = 'Cannot find output files')
#                self.wait_window(p.top)
                self.findDir()


    def InitializeCanvas(self, config_file = None):
        '''Initializes the figure, and then packs it into the main window.
        Should only be called once.'''
        if config_file is None:
            try:
                with open(os.path.join(self.IseultDir, '.iseult_configs', self.cmd_args.p.strip().replace(' ', '_') +'.yml')) as f:
                    cfgDict = yaml.safe_load(f)
            except:
                print('Cannot find/load ' +  self.cmd_args.p.strip().replace(' ', '_') +'.yml in .iseult_configs. If the name of view contains whitespace,')
                print('either it must be enclosed in quotation marks or given with whitespace removed.')
                print('Name is case sensitive. Reverting to Default view')
                with open(os.path.join(self.IseultDir, '.iseult_configs', 'Default.yml')) as f:
                    cfgDict = yaml.safe_load(f)
        else:
            with open(config_file) as f:
                cfgDict = yaml.safe_load(f)


        # divy up the figure into a bunch of subplots using GridSpec.
        self.gs0 = gridspec.GridSpec(self.MainParamDict['NumOfRows'],self.MainParamDict['NumOfCols'])

        # Create the list of all of subplot wrappers
        self.SubPlotList = []
        for i in range(self.MainParamDict['MaxRows']):
            tmplist = [SubPlotWrapper(self, figure = self.f, pos=(i,j)) for j in range(self.MainParamDict['MaxCols'])]
            self.SubPlotList.append(tmplist)
        for i in range(self.MainParamDict['MaxRows']):
            for j in range(self.MainParamDict['MaxCols']):
                tmp_str = f"Chart{i}_{j}"
                if tmp_str in cfgDict.keys():
                    tmpchart_type = cfgDict[tmp_str]['ChartType']
                    self.SubPlotList[i][j].SetGraph(tmpchart_type)
                    for key, val in cfgDict[tmp_str].items():
                        self.SubPlotList[i][j].PlotParamsDict[tmpchart_type][key] = val

                else:
                    # The graph isn't specifiedin the config file, just set it equal to phase plots
                    self.SubPlotList[i][j].SetGraph('PhasePlot')


        # Make a list that will hold the previous ctype
        self.MakePrevCtypeList()
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        self.ReDrawCanvas()
        self.f.canvas.mpl_connect('button_press_event', self.onclick)

    def LoadConfig(self, config_file):
        # First get rid of any & all pop up windows:
        if self.settings_window is not None:
            self.settings_window.destroy()
        if self.measure_window is not None:
            self.measure_window.destroy()
        # Go through each sub-plot destroying any pop-up and
        # restoring to default params
        for i in range(self.MainParamDict['NumOfRows']):
            for j in range(self.MainParamDict['NumOfCols']):
                self.SubPlotList[i][j].RestoreDefaultPlotParams()
                try:
                    self.SubPlotList[i][j].graph.settings_window.destroy()
                except:
                    pass
        # Read in the config file
        #config = configparser.RawConfigParser()
        #config.read(config_file)
        cfgDict = {}
        with open(config_file, 'r') as f:
            cfgDict = yaml.safe_load(f)
        # Generate the Main Param Dict
        self.GenMainParamDict(config_file)

        #Loading a config file may change the stride... watch out!
        if self.stride != self.MainParamDict['PrtlStride']:
            self.stride = self.MainParamDict['PrtlStride']
            self.StrideChanged()
        # Load in all the subplot params
        for i in range(self.MainParamDict['NumOfRows']):
            for j in range(self.MainParamDict['NumOfCols']):
                tmp_str = 'Chart' + str(i) + '_' + str(j)

                if tmp_str in cfgDict.keys():

                    tmpchart_type = cfgDict[tmp_str]['ChartType']
                    self.SubPlotList[i][j].SetGraph(tmpchart_type)
                    for key, val in cfgDict[tmp_str].items():
                        self.SubPlotList[i][j].PlotParamsDict[tmpchart_type][key] = val
                else:
                    # The graph isn't specified in the config file, just set it equal to a phase plot
                    self.SubPlotList[i][j].SetGraph('PhasePlot')
        # There are a few parameters that need to be loaded separately, mainly in the playbackbar.
        self.playbackbar.RecVar.set(self.MainParamDict['Recording'])
        self.playbackbar.LoopVar.set(self.MainParamDict['LoopPlayback'])

        # refresh the geometry
        print(self.MainParamDict['WindowSize'])
        self.geometry(self.MainParamDict['WindowSize'])
        if self.MainParamDict['HorizontalCbars']:
            self.axes_extent = self.MainParamDict['HAxesExtent']
            self.cbar_extent = self.MainParamDict['HCbarExtent']
            self.SubPlotParams = self.MainParamDict['HSubPlotParams']

        else:
            self.axes_extent = self.MainParamDict['VAxesExtent']
            self.cbar_extent = self.MainParamDict['VCbarExtent']
            self.SubPlotParams = self.MainParamDict['VSubPlotParams']
        self.f.subplots_adjust( **self.SubPlotParams)
        # refresh the gridspec and re-draw all of the subplots
        self.UpdateGridSpec()

    def UpdateGridSpec(self, *args):
        '''A function that handles updates the gridspec that divides up of the
        plot into X x Y subplots'''
        # To prevent orphaned windows, we have to kill all of the windows of the
        # subplots that are no longer shown.

        for i in range(self.MainParamDict['MaxRows']):
            for j in range(self.MainParamDict['MaxCols']):
                if i < self.MainParamDict['NumOfRows'] and j < self.MainParamDict['NumOfCols']:
                    pass
                elif self.SubPlotList[i][j].graph.settings_window is not None:
                    self.SubPlotList[i][j].graph.settings_window.destroy()

        self.gs0 = gridspec.GridSpec(self.MainParamDict['NumOfRows'],self.MainParamDict['NumOfCols'])
        self.RenewCanvas(keep_view = False, ForceRedraw = True)

    def LoadAllKeys(self):
        ''' A function that will find out will arrays need to be loaded for
        to draw the graphs. Then it will save all the data necessaru to
        If the time hasn't changed, it will only load new keys.'''
        # Make a dictionary that stores all of the keys we will need to load
        # to draw the graphs.
        self.ToLoad = {'Flds': [], 'Prtl': [], 'Param': [], 'Spect': []}
        # we always load time because it is needed to calculate the shock location
        self.ToLoad[self.H5KeyDict['time']].append('time')
        # We always load enough to calculate xmin, xmax, ymin, ymax & the cpu domains:
        self.ToLoad[self.H5KeyDict['c_omp']].append('c_omp')
        self.ToLoad[self.H5KeyDict['istep']].append('istep')
        self.ToLoad[self.H5KeyDict['dens']].append('dens')
        self.ToLoad[self.H5KeyDict['mx']].append('mx')
        self.ToLoad[self.H5KeyDict['my']].append('my')
        # look at each subplot and see what is needed
        for i in range(self.MainParamDict['NumOfRows']):
            for j in range(self.MainParamDict['NumOfCols']):
                # for each subplot, see what keys are needed
                tmpList = self.SubPlotList[i][j].GetKeys()

                for elm in tmpList:
                    # find out what type of file the key is stored in
                    ftype = self.H5KeyDict[elm]
                    # add the key to the list of that file type
                    self.ToLoad[ftype].append(elm)

        # Check to make sure the 2DSlice is OK...
        # Grab c_omp & istep
        with h5py.File(os.path.join(self.dirname,self.PathDict['Param'][self.TimeStep.value-1]), 'r') as f:
            self.c_omp = f['c_omp'][0]
            self.istep = f['istep'][0]

        # FIND THE SLICE
        with h5py.File(os.path.join(self.dirname,self.PathDict['Flds'][self.TimeStep.value-1]), 'r') as f:
            self.MaxZInd = f['bx'].shape[0]-1
            self.MaxYInd = f['bx'].shape[1]-1
            self.MaxXInd = f['bx'].shape[2]-1

            self.ySlice = int(np.around(self.MainParamDict['ySlice']*self.MaxYInd))
            self.zSlice = int(np.around(self.MainParamDict['zSlice']*self.MaxZInd))


        # See if we are in a new Directory
        if self.NewDirectory:
            # Create a new Dictionary that will have StateHashes of visited steps
            self.SavedHashes = {}
            self.SavedImgStr = {}
            self.SavedImgSize = {}
            self.diff_from_home = []


            self.TotalEnergyTimeSteps = []
            self.TotalEnergyTimes = np.array([])
            self.TotalIonEnergy = np.array([])
            self.TotalElectronEnergy = np.array([])

            self.TotalMagEnergy = np.array([])
            self.TotalBxEnergy = np.array([])
            self.TotalByEnergy = np.array([])
            self.TotalBzEnergy = np.array([])

            self.TotalExEnergy = np.array([])
            self.TotalEyEnergy = np.array([])
            self.TotalEzEnergy = np.array([])
            self.TotalElectricEnergy = np.array([])


            # Make a list of timesteps we have already loaded.
            self.timestep_visited = []

            # Timestep queue that ensures that we delete the oldest viewed
            # timestep if memory gets too large
            self.timestep_queue = deque()

            # For each timestep we visit, we will load a dictionary and place it in a list
            self.ListOfDataDict = []

            self.NewDirectory = False
        # see if one of the plots is the total energy panel
        self.showing_total_energy_plt = False
        for i in range(self.MainParamDict['NumOfRows']):
            for j in range(self.MainParamDict['NumOfCols']):
                self.showing_total_energy_plt += str(self.SubPlotList[i][j].chartType) == 'TotalEnergyPlot'

        if not self.TimeStep.value in self.TotalEnergyTimeSteps:
            if self.showing_total_energy_plt:
                tmp_L = ['ui', 'vi', 'wi', 'ue', 've', 'we', 'mi', 'me', 'stride', 'bx', 'by', 'bz', 'ex', 'ey', 'ez','qi','c']
                for elm in tmp_L:
                    self.ToLoad[self.H5KeyDict[elm]].append(elm)

        if self.TimeStep.value in self.timestep_visited:
            cur_ind = self.timestep_visited.index(self.TimeStep.value)
            self.timestep_queue.remove(self.TimeStep.value)
            self.DataDict = self.ListOfDataDict[cur_ind]
            for pkey in self.ToLoad.keys():
                tmplist = list(set(self.ToLoad[pkey])) # get rid of duplicate keys
                tmplist2 = np.copy(tmplist)

                # get rid of keys that are already loaded
                for i in range(len(tmplist2)):
                    if tmplist2[i] in self.DataDict.keys():
                        tmplist.remove(tmplist2[i])
                # Now iterate over each path key and create a datadictionary
                if len(tmplist)> 0:
                    if pkey == 'Prtl': # we load particle arrays with a stride because they are very expensive
                        with h5py.File(os.path.join(self.dirname,self.PathDict[pkey][self.TimeStep.value-1]), 'r') as f:
                            for elm in tmplist:
                                try:
                                    # Load all the keys
                                    self.DataDict[elm] = f[elm][::self.MainParamDict['PrtlStride']]

                                except KeyError:
                                    raise
                    else:
                        with h5py.File(os.path.join(self.dirname,self.PathDict[pkey][self.TimeStep.value-1]), 'r') as f:
                            for elm in tmplist:
                                try:
                                    # Load all the keys
                                    if elm == 'spect_dens':
                                        self.DataDict[elm] = f['dens'][:]
                                    else:
                                        self.DataDict[elm] = f[elm][:]

                                except KeyError:
                                    if elm == 'sizex':
                                        self.DataDict[elm] = 1
                                    elif elm == 'c':
                                        self.DataDict[elm]= 0.45
                                    elif elm == 'ppc0':
                                        self.DataDict[elm] = np.NaN
                                    elif elm == 'my':
                                        tmpSize = ((self.MaxYInd+1)*f['istep'][0])//(f['my0'][0]-5)
                                        self.DataDict[elm] = np.ones(tmpSize)*(f['my0'][0])
                                    elif elm == 'mx':
                                        tmpSize = ((self.MaxXInd+1)*f['istep'][0])//(f['mx0'][0]-5)
                                        self.DataDict[elm] = np.ones(tmpSize)*(f['mx0'][0])


                                    else:
                                        raise

            self.timestep_queue.append(self.TimeStep.value)

        else:
            # The time has not already been visited so we have to reload everything
            self.DataDict = {}
            for pkey in self.ToLoad.keys():
                tmplist = list(set(self.ToLoad[pkey])) # get rid of duplicate keys
                # Load the file
                if len(tmplist)> 0:
                    if pkey =='Prtl': # we load particle arrays with a stride because they are expensive
                        with h5py.File(os.path.join(self.dirname,self.PathDict[pkey][self.TimeStep.value-1]), 'r') as f:
                            for elm in tmplist:
                                try:
                                    # Load all the key
                                    self.DataDict[elm] = f[elm][::self.MainParamDict['PrtlStride']]

                                except KeyError:
                                    raise
                    else:
                        with h5py.File(os.path.join(self.dirname,self.PathDict[pkey][self.TimeStep.value-1]), 'r') as f:
                            for elm in tmplist:
                                try:
                                    # Load all the keys
                                    if elm == 'spect_dens':
                                        self.DataDict[elm] = f['dens'][:]
                                    else:
                                        self.DataDict[elm] = f[elm][:]

                                except KeyError:
                                    if elm == 'sizex':
                                        self.DataDict[elm] = 1
                                    elif elm == 'c':
                                        self.DataDict[elm]= 0.45
                                    elif elm == 'ppc0':
                                        self.DataDict[elm] = np.NaN
                                    elif elm == 'my':
                                        tmpSize = ((self.MaxYInd+1)*f['istep'][0])//(f['my0'][0]-5)
                                        self.DataDict[elm] = np.ones(tmpSize)*(f['my0'][0])
                                    elif elm == 'mx':
                                        tmpSize = ((self.MaxXInd+1)*f['istep'][0])//(f['mx0'][0]-5)
                                        self.DataDict[elm] = np.ones(tmpSize)*(f['mx0'][0])

                                    else:
                                        raise

            # don't keep more than 30 time steps in memory because of RAM issues
            if len(self.timestep_visited)>30:
                oldest_time = self.timestep_queue.popleft()
                oldest_ind = self.timestep_visited.index(oldest_time)
                self.timestep_visited.remove(oldest_time)
                self.ListOfDataDict.pop(oldest_ind)
            self.timestep_visited.append(self.TimeStep.value)
            self.ListOfDataDict.append(self.DataDict)
            self.timestep_queue.append(self.TimeStep.value)

        if not self.TimeStep.value in self.TotalEnergyTimeSteps:
            if self.showing_total_energy_plt:
                self.TotalEnergyTimeSteps.append(self.TimeStep.value)
                self.TotalEnergyTimeSteps.sort()
                ind = self.TotalEnergyTimes.searchsorted(self.DataDict['time'][0])
                self.TotalEnergyTimes = np.append(np.append(self.TotalEnergyTimes[0:ind],self.DataDict['time'][0]),self.TotalEnergyTimes[ind:])

                TotalElectronKE = self.DataDict['ue']*self.DataDict['ue']
                TotalElectronKE += self.DataDict['ve']*self.DataDict['ve']
                TotalElectronKE += self.DataDict['we']*self.DataDict['we']+1
                TotalElectronKE = np.sum(np.sqrt(TotalElectronKE)-1)
                #TotalElectronKE += -len(self.DataDict['we'])

                TotalElectronKE *= self.DataDict['stride'][0]*self.MainParamDict['PrtlStride'] # multiply by the stride.
                TotalElectronKE *= np.abs(self.DataDict['qi'][0])*self.DataDict['c'][0]**2 # * m_e c^2, mass of particle is its charge, qe/me=1



                TotalIonKE = self.DataDict['ui']*self.DataDict['ui']
                TotalIonKE += self.DataDict['vi']*self.DataDict['vi']
                TotalIonKE += self.DataDict['wi']*self.DataDict['wi']+1
                TotalIonKE = np.sum(np.sqrt(TotalIonKE)-1)
                #TotalIonKE += -len(self.DataDict['we'])

                TotalIonKE *= self.DataDict['stride'][0]*self.MainParamDict['PrtlStride'] # multiply by the stride
                TotalIonKE *= self.DataDict['mi'][0]/self.DataDict['me']*np.abs(self.DataDict['qi'][0])*self.DataDict['c'][0]**2 #mass of particle is its charge, qe/me=1

                TotalKE = (TotalElectronKE +TotalIonKE)
                # Divide by x size
#                TotalKEDensity *= (self.DataDict['dens'][0,:,:].shape[1]/self.DataDict['c_omp'][0]*self.DataDict['istep'][0])**-1
                # Divide by y size
#                TotalKEDensity *= (self.DataDict['dens'][0,:,:].shape[0]/self.DataDict['c_omp'][0]*self.DataDict['istep'][0])**-1

                self.TotalElectronEnergy = np.append(np.append(self.TotalElectronEnergy[0:ind],TotalElectronKE),self.TotalElectronEnergy[ind:])
                self.TotalIonEnergy = np.append(np.append(self.TotalIonEnergy[0:ind],TotalIonKE),self.TotalIonEnergy[ind:])


                BxEnergy = np.sum(self.DataDict['bx'][:,:,:]*self.DataDict['bx'][:,:,:]) * self.DataDict['istep'][0]**2*.5
                ByEnergy = np.sum(self.DataDict['by'][:,:,:]*self.DataDict['by'][:,:,:]) * self.DataDict['istep'][0]**2*.5
                BzEnergy = np.sum(self.DataDict['bz'][:,:,:]*self.DataDict['bz'][:,:,:]) * self.DataDict['istep'][0]**2*.5
                #TotalBEnergy = BxEnergy + ByEnergy + BzEnergy

                ExEnergy = np.sum(self.DataDict['ex'][:,:,:]*self.DataDict['ex'][:,:,:]) * self.DataDict['istep'][0]**2*.5
                EyEnergy = np.sum(self.DataDict['ey'][:,:,:]*self.DataDict['ey'][:,:,:]) * self.DataDict['istep'][0]**2*.5
                EzEnergy = np.sum(self.DataDict['ez'][:,:,:]*self.DataDict['ez'][:,:,:]) * self.DataDict['istep'][0]**2*.5

                #TotalEEnergy = ExEnergy + EyEnergy + EzEnergy

                # sum over the array and then divide by the number of points len(x)*len(y)
                self.TotalBxEnergy = np.append(np.append(self.TotalBxEnergy[0:ind],BxEnergy), self.TotalBxEnergy[ind:])
                self.TotalByEnergy = np.append(np.append(self.TotalByEnergy[0:ind],ByEnergy), self.TotalByEnergy[ind:])
                self.TotalBzEnergy = np.append(np.append(self.TotalBzEnergy[0:ind],BzEnergy), self.TotalBzEnergy[ind:])
                self.TotalMagEnergy = self.TotalBxEnergy + self.TotalByEnergy + self.TotalBzEnergy

                self.TotalExEnergy = np.append(np.append(self.TotalExEnergy[0:ind],ExEnergy), self.TotalExEnergy[ind:])
                self.TotalEyEnergy = np.append(np.append(self.TotalEyEnergy[0:ind],EyEnergy), self.TotalEyEnergy[ind:])
                self.TotalEzEnergy = np.append(np.append(self.TotalEzEnergy[0:ind],EzEnergy), self.TotalEzEnergy[ind:])
                self.TotalElectricEnergy = self.TotalExEnergy + self.TotalEyEnergy + self.TotalEzEnergy

        if self.MainParamDict['ConstantShockVel']:
            # We can just calculate the time * self.shock_speed
            if np.isnan(self.prev_shock_loc):
                # If self.prev_shock_loc is NaN, that means this is the first time
                # the shock has been found, and the previous and current shock_loc
                # should be the same.

                # First calculate the new shock location
                self.shock_loc = self.DataDict['time'][0]*self.shock_speed
                # Set previous shock loc to current location
                self.prev_shock_loc = np.copy(self.shock_loc)
            else:
                # First save the previous shock location,
                self.prev_shock_loc = np.copy(self.shock_loc)
                # Now calculate the new shock location
                self.shock_loc = self.DataDict['time'][0]*self.shock_speed

        else:
            # Let's see if the shock_loc is in the DataDict
            if not 'shock_loc' in self.DataDict.keys():
                # Have to figure out where the shock is

                jstart = int(min(10*self.DataDict['c_omp']/self.DataDict['istep'], self.DataDict['dens'][0,:,:].shape[1]))
                cur_xaxis = np.arange(self.DataDict['dens'][0,:,:].shape[1])/self.DataDict['c_omp']*self.DataDict['istep']
                # Find the shock by seeing where the density is 1/2 of it's
                # max value.

                dens_half_max = max(self.DataDict['dens'][0,:,:][self.DataDict['dens'][0,:,:].shape[0]//2,jstart:])*.5

                # Find the farthest location where the average density is greater
                # than half max
                ishock = np.where(self.DataDict['dens'][0,:,:][self.DataDict['dens'][0,:,:].shape[0]//2,jstart:]>=dens_half_max)[0][-1]
                self.DataDict['shock_loc'] = cur_xaxis[ishock]

            if np.isnan(self.prev_shock_loc):
                # If self.prev_shock_loc is NaN, that means this is the first time
                # the shock has been found, and the previous and current shock_loc
                # should be the same.

                # First calculate the new shock location
                self.shock_loc = self.DataDict['shock_loc']
                # Set previous shock loc to current location
                self.prev_shock_loc = np.copy(self.shock_loc)
            else:
                # First save the previous shock location,
                self.prev_shock_loc = np.copy(self.shock_loc)
                # Now calculate the new shock location
                self.shock_loc = self.DataDict['shock_loc']

        self.cpu_x_locs = np.cumsum(self.DataDict['mx']-5)/self.DataDict['c_omp'][0]
        self.cpu_y_locs = np.cumsum(self.DataDict['my']-5)/self.DataDict['c_omp'][0]
        # Now that the DataDict is created, iterate over all the subplots and
        # load the data into them:
        for i in range(self.MainParamDict['NumOfRows']):
            for j in range(self.MainParamDict['NumOfCols']):
                self.SubPlotList[i][j].LoadData()

    def RefreshTimeStep(self):
        ''' A function that will find out will arrays need to be loaded for
        to draw the graphs. Then it will save all the data necessaru to
        If the time hasn't changed, it will only load new keys.'''
        if self.TimeStep.value in self.timestep_visited:
            cur_ind = self.timestep_visited.index(self.TimeStep.value)
            self.timestep_visited.pop(cur_ind)
            self.ListOfDataDict.pop(cur_ind)
            self.timestep_queue.remove(self.TimeStep.value)

        if self.TimeStep.value in self.TotalEnergyTimeSteps:
            self.TotalEnergyTimeSteps.remove(self.TimeStep.value)
            ind = self.TotalEnergyTimes.searchsorted(self.DataDict['time'][0])
            if ind < len(self.TotalEnergyTimes)-1:
                self.TotalEnergyTimes = np.append(self.TotalEnergyTimes[0:ind],self.TotalEnergyTimes[ind+1:])
                self.TotalElectronEnergy = np.append(self.TotalElectronEnergy[0:ind],self.TotalElectronEnergy[ind+1:])
                self.TotalIonEnergy = np.append(self.TotalIonEnergy[0:ind],self.TotalIonEnergy[ind+1:])
                self.TotalMagEnergy = np.append(self.TotalMagEnergy[0:ind], self.TotalMagEnergy[ind+1:])
                self.TotalBzEnergy = np.append(self.TotalBzEnergy[0:ind], self.TotalBzEnergy[ind+1:])
                self.TotalByEnergy = np.append(self.TotalByEnergy[0:ind], self.TotalByEnergy[ind+1:])
                self.TotalBxEnergy = np.append(self.TotalBxEnergy[0:ind], self.TotalBxEnergy[ind+1:])
                self.TotalEzEnergy = np.append(self.TotalEzEnergy[0:ind], self.TotalEzEnergy[ind+1:])
                self.TotalEyEnergy = np.append(self.TotalEyEnergy[0:ind], self.TotalEyEnergy[ind+1:])
                self.TotalExEnergy = np.append(self.TotalExEnergy[0:ind], self.TotalExEnergy[ind+1:])
                self.TotalElectricEnergy = np.append(self.TotalElectricEnergy[0:ind], self.TotalElectricEnergy[ind+1:])
            else:
                self.TotalEnergyTimes = self.TotalEnergyTimes[0:ind]
                self.TotalElectronEnergy = self.TotalElectronEnergy[0:ind]
                self.TotalIonEnergy = self.TotalIonEnergy[0:ind]
                self.TotalMagEnergy = self.TotalMagEnergy[0:ind]
                self.TotalBxEnergy = self.TotalBxEnergy[0:ind]
                self.TotalByEnergy = self.TotalByEnergy[0:ind]
                self.TotalBzEnergy = self.TotalBzEnergy[0:ind]
                self.TotalExEnergy = self.TotalExEnergy[0:ind]
                self.TotalEyEnergy = self.TotalEyEnergy[0:ind]
                self.TotalEzEnergy = self.TotalEzEnergy[0:ind]
                self.TotalElectricEnergy = self.TotalElectricEnergy[0:ind]

    def MakePrevCtypeList(self):
        self.prev_ctype_list = []
        for i in range(self.MainParamDict['NumOfRows']):
            tmp_ctype_l = []
            for j in range(self.MainParamDict['NumOfCols']):
                tmp_ctype_l.append(str(self.SubPlotList[i][j].chartType))
            self.prev_ctype_list.append(tmp_ctype_l)

    def SaveLLoc(self):
        for i in range(self.MainParamDict['NumOfRows']):
            for j in range(self.MainParamDict['NumOfCols']):
                subplot = self.SubPlotList[i][j]
                if subplot.chartType == 'Moments' or subplot.chartType == 'TotalEnergyPlot':
                    try:
                        if subplot.graph.legend._get_loc() != 1:
                            subplot.SetPlotParam('legend_loc', ' '.join(str(x) for x in subplot.graph.legend._get_loc()), update_plot = False)
                    except:
                        pass
                if subplot.chartType == 'SpectraPlot':
                    try:
                        if subplot.graph.legDelta._get_loc() != 1:
                            subplot.SetPlotParam('PL_legend_loc', ' '.join(str(x) for x in subplot.graph.legDelta._get_loc()), update_plot = False)
                    except:
                        pass
                    try:
                        if subplot.graph.legT._get_loc() != 2:
                            subplot.SetPlotParam('T_legend_loc', ' '.join(str(x) for x in subplot.graph.legT._get_loc()), update_plot = False)
                    except:
                        pass
    def SetLLoc(self):
        for i in range(self.MainParamDict['NumOfRows']):
            for j in range(self.MainParamDict['NumOfCols']):
                subplot = self.SubPlotList[i][j]
                if subplot.chartType == 'Moments' or subplot.chartType == 'TotalEnergyPlot':
                    if subplot.GetPlotParam('legend_loc') != 'N/A':
                        tmp_tup = float(subplot.GetPlotParam('legend_loc').split()[0]),float(subplot.GetPlotParam('legend_loc').split()[1])
                        subplot.graph.legend._set_loc(tmp_tup)
                if subplot.chartType == 'SpectraPlot':
                    if subplot.GetPlotParam('T_legend_loc') != 'N/A':
                        tmp_tup = float(subplot.GetPlotParam('T_legend_loc').split()[0]),float(subplot.GetPlotParam('T_legend_loc').split()[1])
                        try:
                            subplot.graph.legT._set_loc(tmp_tup)
                        except:
                            pass
                    if subplot.GetPlotParam('PL_legend_loc') != 'N/A':
                        tmp_tup = float(subplot.GetPlotParam('PL_legend_loc').split()[0]),float(subplot.GetPlotParam('PL_legend_loc').split()[1])
                        try:
                            subplot.graph.legDelta._set_loc(tmp_tup)
                        except:
                            pass

    """
    def FindCbars(self, prev = False):
        ''' A function that will find where all the cbars are in the current view '''
        self.cbarList = []
        for i in range(self.MainParamDict['NumOfRows']):
            for j in range(self.MainParamDict['NumOfCols']):
                self.IsCbarList.append(False)
                if prev ==True:
                    if self.SubPlotList[i][j].GetPlotParam('twoD') == 1 and not self.SubPlotList[i][j].Changedto2D:
                    # Note the axes still show up in the view if they are set to zero so we have to do it this way.
                        self.IsCbarList.append(True)
                    elif self.SubPlotList[i][j].Changedto1D:
                        self.IsCbarList.append(True)
                elif self.SubPlotList[i][j].GetPlotParam('twoD') == 1:
                    self.IsCbarList.append(True)

    """
    def SaveView(self):
        # A function that will make sure our view will stay the same as the
        # plot updates.

        cur_view = []
        for ax, (view, (pos_orig, pos_active)) in self.toolbar._nav_stack().items():
            #print(type(ax))
            if ax in self.cbarList:
                continue
            cur_view.append(view)
        #    Go to the home view
        self.toolbar._nav_stack.home()
        #self.toolbar.home()
        home_view = []
        for ax, (view, (pos_orig, pos_active)) in self.toolbar._nav_stack().items():
            #print(type(ax))
            if ax in self.cbarList:
                continue
            home_view.append(view)
        #home_view =  list(self.toolbar._nav_stack.__call__())

            #print(view, pos_orig)
            # Find cbars
            #self.FindCbars(prev=True)
            # Filter out the colorbar axes
            #print(self.IsCbarList)
        try:
            self.is_changed_list = []
            self.diff_from_home = []
            self.old_views = []
            if cur_view is not None:
                for i in range(len(cur_view)):
                    is_changed =[]
                    diff_list = []
                    for j in range(4):
                        num_changed = home_view[i][j]-cur_view[i][j] != 0.0
                        is_changed.append(num_changed)
                        if num_changed:
                            if self.MainParamDict['xLimsRelative'] and j < 2:
                                #define the difference relative to the shock loc
                                diff_list.append(cur_view[i][j]-self.shock_loc)
                            else:
                                # define the difference relative to the home loc
                                diff_list.append(cur_view[i][j])

                        else:
                            # They haven't zoomed in, diff should be zero,
                            # but I'm making it a string so cur_view-shock_loc can be
                            # equal to zero and the hash still distinguish between the two cases.
                            diff_list.append('n/a')


                    self.is_changed_list.append(is_changed)
                    self.old_views.append(cur_view[i])
                    self.diff_from_home.append(diff_list)

        except IndexError:
            pass

    def LoadView(self):

        #self.toolbar._nav_stack.clear()
        #self.toolbar.home()
        #self.set_history_buttons()
        #self.toolbar._nav_stack.home()
        #self.toolbar.set_history_buttons()
        #self.toolbar._update_view()

        #self._update_view()

        self.toolbar.push_current()

        cur_view = []
        tmpList = []
        for ax, (view, (pos_orig, pos_active)) in self.toolbar._nav_stack().items():

            tmpList.append(view)
            if ax in self.cbarList:
                continue

            cur_view.append(view)
            #ax._set_view(view)
            #self.toolbar._update_view()
        # Find the cbars in the current plot
        #self.FindCbars()
        try:
            # put the parts that have changed from the old view
            # into the proper place in the next view
            m = 0 # a counter that allows us to go from labeling the plots in [i][j] to 1d
            k = 0 # a counter that skips over the colorbars
            for i in range(self.MainParamDict['NumOfRows']):
                for j in range(self.MainParamDict['NumOfCols']):
                    tmp_old_view = list(self.old_views.pop(0))
                    tmp_new_view = list(cur_view[k])
                    self.SubPlotList[i][j].graph.axes._set_view(cur_view[k])
                    if self.prev_ctype_list[i][j] == self.SubPlotList[i][j].chartType:
                        # see if the view has changed from the home view
                        is_changed = self.is_changed_list[m]
                        if self.SubPlotList[i][j].Changedto2D or self.SubPlotList[i][j].Changedto1D:
                            # only keep the x values if they have changed
                            for n in range(2):
                                if is_changed[n]:
                                    if self.SubPlotList[i][j].PlotParamsDict[self.SubPlotList[i][j].chartType]['spatial_x']:
                                        tmp_new_view[n] = tmp_old_view[n]+self.MainParamDict['xLimsRelative']*(self.shock_loc-self.prev_shock_loc)
                                    else:
                                        tmp_new_view[n] = tmp_old_view[n]
                        else:
                            # Keep any y or x that is changed
                            for n in range(4):
                                if is_changed[n]:
                                    tmp_new_view[n] = tmp_old_view[n]
                                    if n < 2:
                                        if self.SubPlotList[i][j].PlotParamsDict[self.SubPlotList[i][j].chartType]['spatial_x']:
                                            tmp_new_view[n] = tmp_old_view[n]+self.MainParamDict['xLimsRelative']*(self.shock_loc-self.prev_shock_loc)
                                        else:
                                            tmp_new_view[n] = tmp_old_view[n]

                    cur_view[k] = tmp_new_view

                    self.SubPlotList[i][j].graph.axes._set_view(cur_view[k])

                    # Handle the counting of the 'views' array in matplotlib
                    #skip over colorbar axes
                    m += 1
                    k += 1
                    self.SubPlotList[i][j].Changedto1D = False
                    self.SubPlotList[i][j].Changedto2D = False
            self.toolbar.push_current()
            #self.toolbar._nav_stack.push(cur_view)
            #print(len(self.toolbar._nav_stack._elements))
            #self.toolbar.set_history_buttons()
            #self.toolbar._update_view()
        except IndexError:
            pass

    def RenewCanvas(self, keep_view = True, ForceRedraw = False):

        '''We have two way of updated the graphs: 1) by refreshing them using
        self.RefreshCanvas, we don't recreate all of the artists that matplotlib
        needs to make the plot work. self.RefreshCanvas should be fast. Two we
        can ReDraw the canvas using self.ReDrawCanvas. This recreates all the
        artists and will be slow. Sometimes the graph must be redrawn however,
        if the GridSpec changed, more plots are added, the chartype changed, if
        the plot went from 2d to 1D, etc.. If any change occurs that requires a
        redraw, renewcanvas must be called with ForceRedraw = True. '''

        self.SaveLLoc()
        if ForceRedraw:
            self.ReDrawCanvas(keep_view = keep_view)
        else:
            self.RefreshCanvas(keep_view = keep_view)
        # Record the current ctypes for later
        self.MakePrevCtypeList()


        # remove some unnecessary data
        tmp_list = ['ui', 'vi', 'wi', 'ue', 've', 'we', 'che', 'chi']
        for elm in tmp_list:
            self.DataDict.pop(elm, None)
        self.SetLLoc()
        # Save the image for quick playback later
        #self.SaveTmpFig()




    def HashIseultState(self):
        ''' A function that saves a hash of the current state of Iseult. Used to
        determine if we can just show a saved image of a previous timeslice,
        or if we must reload it.'''

        #First update the main param dict to save the current window size:
        self.MainParamDict['WindowSize'] = str(self.winfo_width())+'x'+str(self.winfo_height())
        # a tuple that will eventually be hashed.
        state_tuple = self.freeze(self.diff_from_home)
        # keys we should skip over when making the hash.
        SkipList = ['Reload2End', 'WaitTime', 'MaxCols', 'MaxRows', 'SkipSize', 'Recording', 'ClearFig']
        for key in self.MainParamDict.keys():
            if key in SkipList:
                pass
            else:
                state_tuple += key, self.freeze(self.MainParamDict[key])
        for i in range(self.MainParamDict['NumOfRows']):
            for j in range(self.MainParamDict['NumOfCols']): #add every chart's param dictionary
                tmp_str = 'Chart' + str(i) + ',' + str(j)
                tmp_ctype = self.SubPlotList[i][j].chartType
                state_tuple += tmp_str, tmp_ctype, self.freeze(self.SubPlotList[i][j].PlotParamsDict[tmp_ctype])

        # Now save the difference of the zoom from the home view of the current plot


        if self.showing_total_energy_plt:
            state_tuple += self.freeze(self.TotalEnergyTimeSteps)
        # add to the state_tuple the last modification time of all the output files:
        for key in self.PathDict.keys():
            state_tuple += os.path.getmtime(os.path.join(self.dirname,self.PathDict[key][self.TimeStep.value-1])),
#        fname = 'iseult_img_'+ str(self.TimeStep.value).zfill(3)+'.png'
        self.StateHash = hash(state_tuple)
#        print self.freeze(self.MainParamDict)

    def SaveTmpFig(self):
        self.HashIseultState()
        already_saved = False
        if self.TimeStep.value in self.SavedHashes.keys(): # we have already saved an image for this TimeStep
            # is the current state of Iseult equal to the state when we saved said image?
            already_saved = self.SavedHashes[self.TimeStep.value] ==  self.StateHash


        if not already_saved: # nope, better save it again!
            self.SavedHashes[self.TimeStep.value] =  self.StateHash
            self.SavedImgSize[self.TimeStep.value] =int(self.f.get_size_inches()[0]*self.f.dpi), int(self.f.get_size_inches()[1]*self.f.dpi)

            ram = io.BytesIO()
            self.f.savefig(ram, format='raw', dpi=self.f.dpi, facecolor=self.f.get_facecolor())
            ram.seek(0)
            self.SavedImgStr[self.TimeStep.value] = ram.read() # Save the image into SavedImgs
            ram.close()

    def freeze(self, d):
        ''' This function takes in a dictionary or list, which are not hashable,
        and returns a frozen set, which can be used to hash the dictionary'''
        if isinstance(d, dict):
            return frozenset((key, self.freeze(value)) for key, value in d.items())
        elif isinstance(d, list):
            return tuple(self.freeze(value) for value in d)
        return d

    def ReDrawCanvas(self, keep_view = True):
        #  We need to see if the user has moved around the zoom level in python.
        # First we see if there are any views in the toolbar
        cur_view =  self.toolbar._nav_stack.__call__()
        if cur_view is None:
            keep_view = False
        if self.NewDirectory:
            keep_view = False
        if keep_view:
            self.SaveView()

        self.f.clf()
        #
        if self.MainParamDict['ClearFig']:
            self.canvas.draw()

        self.LoadAllKeys()


        # Calculate the new xmin, and xmax

        # Find the first position with a physical x,y & k axis:
        self.first_x = None
        self.first_y = None
        self.first_k = None
        k = 0
        # find the first spatial x and y
        for i in range(self.MainParamDict['NumOfRows']):
            for j in range(self.MainParamDict['NumOfCols']):

                # First handle the axes sharing
                if self.SubPlotList[i][j].chartType == 'FFTPlots':
                    # The plot type is a spectral plot, which has no spatial dim
                    if self.first_k is None:
                        self.first_k = (i,j)
                elif self.SubPlotList[i][j].chartType == 'SpectraPlot':
                    # The plot type is a spectral plot, which has no spatial dim
                    pass
                elif self.MainParamDict['LinkSpatial'] != 1 and self.SubPlotList[i][j].chartType == 'PhasePlot':
                    # If this is the case we don't care about the phase plots
                    # as we don't want to share the axes
                    pass
                elif self.MainParamDict['LinkSpatial'] != 1 and self.SubPlotList[i][j].chartType == 'EnergyPlot':
                    # If this is the case we don't care about the phase plots
                    # as we don't want to share the axes
                    pass
                elif self.MainParamDict['LinkSpatial'] == 3 and self.SubPlotList[i][j].GetPlotParam('twoD'):
                    # If the plot is twoD share the axes
                    if self.first_x is None and self.SubPlotList[i][j].GetPlotParam('spatial_x'):
                        self.first_x = (i,j)
                    if self.first_y is None and self.SubPlotList[i][j].GetPlotParam('spatial_y'):
                        self.first_y = (i,j)

                else:
                    # Just find the first spatial x and y direction.
                    if self.first_x is None and self.SubPlotList[i][j].GetPlotParam('spatial_x'):
                        self.first_x = (i,j)
                    if self.first_y is None and self.SubPlotList[i][j].GetPlotParam('spatial_y'):
                        self.first_y = (i,j)

                # Now... We can draw the graph.
                self.SubPlotList[i][j].DrawGraph()

        if self.MainParamDict['ShowTitle']:
            tmpstr = self.PathDict['Prtl'][self.TimeStep.value-1].split('.')[-1]
            self.f.suptitle(os.path.abspath(self.dirname)+ '/*.'+tmpstr+' at time t = %d $\omega_{pe}^{-1}$'  % round(self.DataDict['time'][0]), size = 15)
        if keep_view:
            self.LoadView()


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
                        k = min(self.SubPlotList[spos[0]][spos[1]].graph.spect_num, len(self.dashes_options)-1)
                        # Append the left line to the list
                        self.SubPlotList[pos[0]][pos[1]].graph.IntRegionLines.append(self.SubPlotList[pos[0]][pos[1]].graph.axes.axvline(
                        max(self.SubPlotList[spos[0]][spos[1]].graph.i_left_loc, self.SubPlotList[pos[0]][pos[1]].graph.xmin+1),
                        linewidth = 1.5, linestyle = '-', color = self.ion_color))
                        # Choose the left dashes pattern
                        self.SubPlotList[pos[0]][pos[1]].graph.IntRegionLines[-1].set_dashes(self.dashes_options[k])

                        # Append the right line to the list
                        self.SubPlotList[pos[0]][pos[1]].graph.IntRegionLines.append(self.SubPlotList[pos[0]][pos[1]].graph.axes.axvline(
                        min(self.SubPlotList[spos[0]][spos[1]].graph.i_right_loc, self.SubPlotList[pos[0]][pos[1]].graph.xmax+1),
                        linewidth = 1.5, linestyle = '-', color = self.ion_color))
                        # Choose the right dashes pattern
                        self.SubPlotList[pos[0]][pos[1]].graph.IntRegionLines[-1].set_dashes(self.dashes_options[k])
            else:
                for spos in self.spectral_plot_list:
                    if self.SubPlotList[spos[0]][spos[1]].GetPlotParam('show_electrons'):
                        k = min(self.SubPlotList[spos[0]][spos[1]].graph.spect_num, len(self.dashes_options)-1)
                        # Append the left line to the list
                        self.SubPlotList[pos[0]][pos[1]].graph.IntRegionLines.append(self.SubPlotList[pos[0]][pos[1]].graph.axes.axvline(
                        max(self.SubPlotList[spos[0]][spos[1]].graph.e_left_loc, self.SubPlotList[pos[0]][pos[1]].graph.xmin+1),
                        linewidth = 1.5, linestyle = '-', color = self.electron_color))
                        # Choose the left dashes pattern
                        self.SubPlotList[pos[0]][pos[1]].graph.IntRegionLines[-1].set_dashes(self.dashes_options[k])

                        # Append the right line to the list
                        self.SubPlotList[pos[0]][pos[1]].graph.IntRegionLines.append(self.SubPlotList[pos[0]][pos[1]].graph.axes.axvline(
                        min(self.SubPlotList[spos[0]][spos[1]].graph.e_right_loc, self.SubPlotList[pos[0]][pos[1]].graph.xmax+1),
                        linewidth = 1.5, linestyle = '-', color = self.electron_color))
                        # Choose the right dashes pattern
                        self.SubPlotList[pos[0]][pos[1]].graph.IntRegionLines[-1].set_dashes(self.dashes_options[k])

        self.canvas.draw()
        self.canvas.get_tk_widget().update_idletasks()


        if self.MainParamDict['Recording']:
            self.PrintFig()


    def RefreshCanvas(self, keep_view = True):
        #  We need to see if the user has moved around the zoom level in python.
        # First we see if there are any views in the toolbar
        cur_view =  self.toolbar._nav_stack.__call__()
        if cur_view is None:

            keep_view = False

            self.diff_from_home = []
            for i in range(self.MainParamDict['NumOfRows']*self.MainParamDict['NumOfCols']):
                self.diff_from_home.append(['n/a', 'n/a', 'n/a', 'n/a'])


        if self.NewDirectory:
            keep_view = False
        if keep_view:
            self.SaveView()


        self.toolbar._nav_stack.clear()

        self.LoadAllKeys()


        # By design, the first_x and first_y cannot change if the graph is
        # being refreshed. Any call that would require this needs a redraw
        # Now we refresh the graph.
        for i in range(self.MainParamDict['NumOfRows']):
            for j in range(self.MainParamDict['NumOfCols']):

                self.SubPlotList[i][j].RefreshGraph()

        if self.MainParamDict['ShowTitle']:
            tmpstr = self.PathDict['Prtl'][self.TimeStep.value-1].split('.')[-1]
            self.f.suptitle(os.path.abspath(self.dirname)+ '/*.'+tmpstr+' at time t = %d $\omega_{pe}^{-1}$'  % round(self.DataDict['time'][0]), size = 15)

        if keep_view:
            self.LoadView()


        for pos in self.phase_plot_list:
            i = 0
            if self.SubPlotList[pos[0]][pos[1]].GetPlotParam('prtl_type') == 0:
                for spos in self.spectral_plot_list:
                    if self.SubPlotList[spos[0]][spos[1]].GetPlotParam('show_ions'):
                        # Update the left line to the list
                        self.SubPlotList[pos[0]][pos[1]].graph.IntRegionLines[i].set_xdata(
                        [max(self.SubPlotList[spos[0]][spos[1]].graph.i_left_loc, self.SubPlotList[pos[0]][pos[1]].graph.xmin+1),
                        max(self.SubPlotList[spos[0]][spos[1]].graph.i_left_loc, self.SubPlotList[pos[0]][pos[1]].graph.xmin+1)])
                        i+=1
                        # Append the right line of the list
                        self.SubPlotList[pos[0]][pos[1]].graph.IntRegionLines[i].set_xdata(
                        [min(self.SubPlotList[spos[0]][spos[1]].graph.i_right_loc, self.SubPlotList[pos[0]][pos[1]].graph.xmax-1),
                        min(self.SubPlotList[spos[0]][spos[1]].graph.i_right_loc, self.SubPlotList[pos[0]][pos[1]].graph.xmax-1)])
                        i+=1
            else:
                for spos in self.spectral_plot_list:
                    if self.SubPlotList[spos[0]][spos[1]].GetPlotParam('show_electrons'):
                        # Update the left line to the list
                        self.SubPlotList[pos[0]][pos[1]].graph.IntRegionLines[i].set_xdata(
                        [max(self.SubPlotList[spos[0]][spos[1]].graph.e_left_loc, self.SubPlotList[pos[0]][pos[1]].graph.xmin+1),
                        max(self.SubPlotList[spos[0]][spos[1]].graph.e_left_loc, self.SubPlotList[pos[0]][pos[1]].graph.xmin-1)])
                        i+=1
                        # Append the right line of the list
                        self.SubPlotList[pos[0]][pos[1]].graph.IntRegionLines[i].set_xdata(
                        [min(self.SubPlotList[spos[0]][spos[1]].graph.e_right_loc, self.SubPlotList[pos[0]][pos[1]].graph.xmax+1),
                        min(self.SubPlotList[spos[0]][spos[1]].graph.e_right_loc, self.SubPlotList[pos[0]][pos[1]].graph.xmax-1)])
                        i+=1
        self.canvas.draw()
        self.canvas.get_tk_widget().update_idletasks()

        if self.MainParamDict['Recording']:
            self.PrintFig()

    def PrintFig(self, MakingMovie = False, Flag = True):
        if self.movie_dir == '':
            self.movie_dir = os.path.abspath(os.path.join(self.dirname, '..', 'Movie'))
        if Flag:
            if self.MainParamDict['Recording'] or MakingMovie:
                try:
                    os.makedirs(self.movie_dir)

                except (OSError, IOError):
                    if not os.path.isdir(self.movie_dir):
                        self.PrintFig(MakingMovie = MakingMovie, Flag = self.recordProblemsPrompt())

            if MakingMovie and os.path.isdir(self.movie_dir):
                try:
                    os.makedirs(os.path.join(self.movie_dir, '../tmp_erase'))

                except (OSError, IOError):
                    if not os.path.isdir(os.path.join(self.movie_dir, '../tmp_erase')):
                        self.PrintFig(MakingMovie = MakingMovie, Flag = self.recordProblemsPrompt())

            fname = 'iseult_img_'+ str(self.TimeStep.value).zfill(3)+'.png'
            if self.MainParamDict['Recording'] :
                try:
                    self.f.savefig(os.path.join(self.movie_dir, fname))#, dpi = 300)#, facecolor=self.f.get_facecolor())#, edgecolor='none')
                except (OSError, IOError):
                    self.PrintFig(MakingMovie = MakingMovie, Flag = self.recordProblemsPrompt())
            if MakingMovie:
                try:
                    self.f.savefig(os.path.join(self.movie_dir, '../tmp_erase', fname))#, dpi = 300)#, facecolor=self.f.get_facecolor())#, edgecolor='none')
                except (OSError, IOError):
                    self.PrintFig(MakingMovie = MakingMovie, Flag = self.recordProblemsPrompt())


    def recordProblemsPrompt(self):
        if messagebox.askyesno("Recording Problems", "You do not have write access to " +self.movie_dir + ". Would you like record frames to a different directory?"):
            mvdir_opt = {}
            mvdir_opt['initialdir'] = self.dirname
            mvdir_opt['mustexist'] = True
            mvdir_opt['parent'] = self
            self.movie_dir = filedialog.askdirectory(title = 'Please choose a different directory where you have write access to save images.', **self.dir_opt)
            return True
        else:
            return False
            self.MainParamDict['Recording'] = False
            self.playbackbar.RecVar.set(False)

    def MakeAMovie(self, fname, start, stop, step, FPS):
        '''Record a movie'''
        # Where-ever you are create a hidden file and then delete that directory:
        self.PrintFig(MakingMovie= True)
        # Delete all the images in that subdirectory
        if os.path.isdir(os.path.join(self.movie_dir, '../tmp_erase')):
            for name in os.listdir(os.path.join(self.movie_dir, '../tmp_erase')):
                os.remove(os.path.join(self.movie_dir, '../tmp_erase', name))

            # First find the last frame is stop is -1:

            if stop == -1:
                stop = len(self.PathDict['Param'])

            # Now build all the frames we have to visit
            frame_arr = np.arange(start, stop, step)
            if frame_arr[-1] != stop:
                frame_arr = np.append(frame_arr, stop)

            # If total energy plot is showing, we have to loop through everything twice.

            if self.showing_total_energy_plt:
                for k in frame_arr:
                    self.TimeStep.set(k)

            for i in frame_arr:
                self.TimeStep.set(i)
                self.PrintFig(MakingMovie  = True)

            # The ffmpeg command we want to call.
            ## ffmpeg -framerate [FPS] -i [NAME_***].png -c:v prores -pix_fmt yuv444p10le [OUTPUTNAME].mov
            cmdstring = ['xterm', '-e','ffmpeg',
                        '-framerate', str(int(FPS)), # Set framerate to the the user selected option
                        '-pattern_type', 'glob',
                        '-i', os.path.join(self.movie_dir, '../tmp_erase','*.png'),
                        '-c:v',
                        'prores',
                        '-pix_fmt',
                        'yuv444p10le',
                        os.path.join(os.path.join(self.movie_dir,'..'),fname)]#, '&']#, # output name,
                        #'<dev/null', '>dev/null', '2>/var/log/ffmpeg.log', '&'] # run in background
            try:
                subprocess.call(cmdstring)
            except OSError:
                try:
                    subprocess.call(cmdstring[2:])
                except OSError:
                    messagebox.showwarning(
                        "Problems saving a movie",
                        "Please make sure that ffmpeg is installedgg on your machine."
                        )

            for name in os.listdir(os.path.join(self.movie_dir, '../tmp_erase')):
                os.remove(os.path.join(self.movie_dir, '../tmp_erase', name))
            os.rmdir(os.path.join(self.movie_dir, '../tmp_erase'))
            '''
            #THIS METHOD TRIES TO USE SUBPROCCESS AND PIPING TO OUTPUT... LEAVING THIS HERE
            #FOR LATER....
            # Draw frames
            im_list = []
            for i in frame_arr:
                self.TimeStep.set(i)


            # Save the image png as a cString
            ram = cStringIO.StringIO()
            self.f.savefig(ram, format='png', dpi=self.f.dpi, facecolor=self.f.get_facecolor())
            ram.seek(0)
            # write to pipe
            im_list.append(ram.read())
            ram.close()


            # Let's try using CStrings and piping to ffmpeg. Here's what we got from the OIC people
            # ffmpeg -y -f image2 -framerate 8 -pattern_type glob -i '*.png' -codec copy out.mov
            # ffmpeg -y -f image2 -framerate 8 -pattern_type glob -i '*.png' -vcodec libx264 -pix_fmt yuv420p out.mp4
            # OLD DEPRECATED METHOD
            #FFMpegWriter = manimation.writers['ffmpeg']
            #writer = FFMpegWriter(fps=FPS, bitrate = 10000)



            # New Method.... First, let's translate the above command into a string



            cmdstring = ('ffmpeg',
                        '-y', '-f', 'image2', # overwrite, image2 is a colorspace thing..
                        'vcodec', 'png',
                        '-framerate', str(int(FPS)), # Set framerate to the the user selected option
                        '-pattern_type', 'glob', '-i', '-', # Not sure what this does... I am going to get rid of it
                        '-vcodec', 'libx264', '-pix_fmt', 'yuv420p', fname+'.mp4') # output encoding
            p = subprocess.Popen(cmdstring, stdin=subprocess.PIPE)

            # Write the images to the pipe
            for i in range(len(im_list)):
                print i
                p.stdin.write(im_list.pop(0))
            # Finish up
            p.communicate()
            #p.stdout.close()
            '''
    def OpenSaveDialog(self):
        SaveDialog(self)
    def OpenMovieDialog(self):
        MovieDialog(self)

    def onclick(self, event):
        '''After being clicked, we should use the x and y of the cursor to
        determine what subplot was clicked'''

        # Since the location of the cursor is returned in pixels and gs0 is
        # given as a relative value, we must first convert the value into a
        # relative x and y
        if not event.inaxes:
            pass
        if event.button == 1:
            pass
        else:
            fig_size = self.f.get_size_inches()*self.f.dpi # Fig size in px

            x_loc = event.x/fig_size[0] # The relative x position of the mouse in the figure
            y_loc = event.y/fig_size[1] # The relative y position of the mouse in the figure

            sub_plots = self.gs0.get_grid_positions(self.f)
            row_array = np.sort(np.append(sub_plots[0], sub_plots[1]))
            col_array = np.sort(np.append(sub_plots[2], sub_plots[3]))
            i = int((len(row_array)-row_array.searchsorted(y_loc))/2)
            j = int(col_array.searchsorted(x_loc)//2)

            self.SubPlotList[i][j].OpenSubplotSettings()

    def shock_finder(self):
        '''The main idea of the shock finder, is we go to the last timestep
        in the simulation and find where the density is half it's max value.
        We then calculate the speed of the of the shock assuming it is
        traveling at constant velocity. We also calculate the initial B & E fields.'''

        # First load the first field file to find the initial size of the
        # box in the x direction, and find the initial field values
        with h5py.File(os.path.join(self.dirname,self.PathDict['Param'][0]), 'r') as f:
            # Find out what sigma is
            try:
                ''' Obviously the most correct way to do this to to calculate b0 from sigma.
                This is proving more difficult that I thought it would be, so I am calculating it
                as Jaehong did.

                sigma = f['sigma'][0]
                gamma0 = f['gamma0'][0]
                c = f['c'][0]
                btheta = f['btheta'][0]
                bphi = f['bphi'][0]

                ppc0 = f['ppc0'][0]
                mi = f['mi'][0]
                me = f['me'][0]
                print mi, c, ppc0
                if gamma0 <1:
                    beta0 = gamma0
                    gamma0 = 1/np.sqrt(1-gamma0**2)
                else:
                    beta0=np.sqrt(1-gamma0**(-2))


                if sigma <= 1E-10:
                    self.b0 = 1.0
                    self.e0 = 1.0
                else:
                    # b0 in the upstream frame
                    self.b0 = np.sqrt(gamma0*ppc0*.5*c**2*(mi+me)*sigma)
                    # Translate to the downstream frame
                    b_x = self.b0*np.cos(btheta)*np.cos(bphi)
                    b_y = self.b0*np.sin(btheta)*np.cos(bphi)
                    b_z = self.b0*np.sin(btheta)*np.cos(bphi)
                    print 'sigma b0', self.b0
                    '''
                    # Normalize by b0
                if np.abs(f['sigma'][0])==0:
                    self.btheta = np.NaN
                else:
                    self.btheta = f['btheta'][0]
            except KeyError:
                self.btheta = np.NaN


        with h5py.File(os.path.join(self.dirname,self.PathDict['Flds'][0]), 'r') as f:
            by = f['by'][:]
            nxf0 = by.shape[1]
            if np.isnan(self.btheta):
                self.b0 = 1.0
                self.e0 = 1.0
            else:
                # Normalize by b0
                self.bx0 = f['bx'][0,-1,-10]
                self.by0 = by[0,-1,-10]
                self.bz0 = f['bz'][0,-1,-10]
                self.b0 = np.sqrt(self.bx0**2+self.by0**2+self.bz0**2)
                self.ex0 = f['ex'][0,-1,-2]
                self.ey0 = f['ey'][0,-1,-2]
                self.ez0 = f['ez'][0,-1,-2]
                self.e0 = np.sqrt(self.ex0**2+self.ey0**2+self.ez0**2)

        # Load the final time step to find the shock's location at the end.
        with h5py.File(os.path.join(self.dirname,self.PathDict['Flds'][-1]), 'r') as f:
            dens_arr =np.copy(f['dens'][0,:,:])

        with h5py.File(os.path.join(self.dirname,self.PathDict['Param'][-1]), 'r') as f:
            # I use this file to get the final time, the istep, interval, and c_omp
            final_time = f['time'][0]
            istep = f['istep'][0]
            interval = f['interval'][0]
            c_omp = f['c_omp'][0]

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
        self.shock_speed = xshock_final/final_time
        self.prev_shock_loc = np.NaN

    def setKnob(self, value):
        # If the time parameter changes update the plots
        """
        if self.playbackbar.playPressed and not self.MainParamDict['Recording']:
            already_saved = False
            if self.TimeStep.value in self.SavedHashes.keys(): # we have already saved an image for this TimeStep
                # is the current state of Iseult equal to the state when we saved said image?
                already_saved = self.SavedHashes[self.TimeStep.value] ==  self.StateHash

            if not already_saved:
                self.RenewCanvas()

            im = Image.frombuffer('RGBA', self.SavedImgSize[self.TimeStep.value], self.SavedImgStr[self.TimeStep.value], 'raw', 'RGBA', 0, 1)
            self.MovieIm.set_data(im)
            self.MovieCanvas.draw()
#            self.MovieCanvas.get_tk_widget().update_idletasks()
            self.playbackbar.tstep.set(str(value))
            #set the slider
            self.playbackbar.slider.set(value)


        else:
        """
        self.RenewCanvas()

        self.playbackbar.tstep.set(str(value))
        #set the slider
        self.playbackbar.slider.set(value)

    def OpenSettings(self, *args):
        if self.settings_window is None:
            self.settings_window = SettingsFrame(self)
        else:
            self.settings_window.destroy()
            self.settings_window = SettingsFrame(self)

    def TxtEnter(self, e):
        self.playbackbar.TextCallback()

def runMe(cmd_args):
    app = MainApp('Iseult', cmd_args)
    app.mainloop()
"""
    parser = argparse.ArgumentParser(description='Plotting program for Tristan-MP files.')
    #        parser.add_argument('integers', metavar='N', type=int, nargs='+',
    #                        help='The maximum file number to consider')
    #        parser.add_argument('--foo', nargs='?', const='c', default='d')
    #        parser.add_argument('bar', nargs='?', default='d')
    parser.add_argument('-n', nargs = '?',# dest='accumulate', action='store_const',
                        const=-1, default=-1,
                        help='Maximum file # to consider')

    parser.add_argument('-O', nargs = '?',# dest='accumulate', action='store_const',
                        const='', default='',
                        help='Directory Iseult will open. Default is output')

    parser.add_argument('-p', nargs = '?',# dest='accumulate', action='store_const',
                        const='Default', default='Default',
                        help='''Open Iseult with the given saved view.
                              If the name of view contains whitespace,
                              either it must be enclosed in quotation marks or given
                              with whitespace removed. Name is case sensitive.''')
    parser.add_argument("-b", help="Run Iseult from bash script. Makes a movie.",
                        action="store_true")

    #parser.add_argument("--wait", help="Wait until current simulation is finished before making movie.",
    #                    action="store_true")

    cmd_args = parser.parse_args()
    #import sys
    #if cmd_args.wait:
    #    " try to parse stdin"
    #    slurm_num = sys.stdin.read().split[-1]
    #    print(slurm_num)
    #    num = 0
    #    done = False
    #    while num < 2000 and not done:
    #        slurm_queue = subprocess.check_output(["squeue"]    )
    #        if slurm_queue.find(slurm_num) != -1:
    #            num += 1
    #            time.sleep(3E5)
    app = MainApp('Iseult', cmd_args)


    app.mainloop()
"""
