#! /usr/bin/env python

import re # regular expressions
import os, sys # Used to make the code portable
import h5py # Allows us the read the data files
import time,string
import matplotlib
import new_cmaps
import numpy as np
from collections import deque
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from phase_plots import PhasePanel
from fields_plots import FieldsPanel
from density_plots import DensPanel
from spectra import SpectralPanel

import Tkinter as Tk
import ttk as ttk
import tkFileDialog

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

def destroy(e):
    sys.exit()

class MyCustomToolbar(NavigationToolbar2TkAgg):
    def __init__(self, plotCanvas, parent):
        # create the default toolbar
        # plotCanvas is the tk Canvas we want to link to the toolbar,
        # parent is the iseult main app
        NavigationToolbar2TkAgg.__init__(self, plotCanvas, parent)


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
        self.parent = parent
        self.chartType = 'PhasePlot'
        # A dictionary that contains all of the plot types.
        self.PlotTypeDict = {'PhasePlot': PhasePanel,
                             'FieldsPlot': FieldsPanel,
                             'DensityPlot': DensPanel,
                             'SpectraPlot': SpectralPanel}
        # A dictionary that will store where everything is in Hdf5 Files
        self.GenParamDict()
        self.figure = figure
        self.subplot_spec = subplot_spec
        self.pos = pos
        self.pos_in_views = None
        self.graph = graph
        self.Changedto1D = False
        self.Changedto2D = False
        #
    def GetKeys(self):
        return self.graph.set_plot_keys()

    def LoadKey(self, h5key):
        return self.parent.DataDict[h5key]

    def ChangeGraph(self, str_arg):
        # Change the graph type
        self.chartType = str_arg
        # put a list of the previous chart types in iseult

        self.graph = self.PlotTypeDict[self.chartType](self.parent, self)
        self.parent.RenewCanvas(ForceRedraw = True)

    def GenParamDict(self):
        # Generate a dictionary that will store all of the params at dict['ctype']['param_name']
        self.PlotParamsDict = {plot_type: '' for plot_type in self.PlotTypeDict.keys()}
        for elm in self.PlotTypeDict.keys():
            self.PlotParamsDict[elm] = {key: self.PlotTypeDict[elm].plot_param_dict[key] for key in self.PlotTypeDict[elm].plot_param_dict.keys()}

    def SetPlotParam(self, pname, val, ctype = None, update_plot = True):
        NeedsRedraw = False
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
        if ctype is None:
            ctype = self.chartType
        return self.PlotParamsDict[ctype][pname]

    def SetGraph(self, ctype = None):
        if ctype:
            self.chartType = ctype
        self.graph = self.PlotTypeDict[self.chartType](self.parent, self)

    def DrawGraph(self):
        self.graph.draw()

    def RefreshGraph(self):
        self.graph.refresh()

    def OpenSubplotSettings(self):
        self.graph.OpenSettings()
class Knob:
    """
    Knob - simple class with a "setKnob" method.
    A Knob instance is attached to a Param instance, e.g., param.attach(knob)
    Base class is for documentation purposes.
    """
    def setKnob(self, value):
        pass

class Param:

    """
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

    """A Class that will handle the time-stepping in Iseult, and has the
    following, a step left button, a play/pause button, a step right button, a
    playbar, and a settings button."""

    def __init__(self, parent, param, canvas = None):
        Tk.Frame.__init__(self)
        self.parent = parent

        self.skipSize = 5
        self.waitTime = .01
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
        # allow us to select a time
        self.slider = ttk.Scale(self, from_=self.param.minimum, to=self.param.maximum, command = self.ScaleHandler)
        self.slider.set(self.param.value)
        self.slider.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)


        # a measurement button that should lauch a window to take measurements.
        self.MeasuresB= ttk.Button(self, text='Measure', command=self.OpenMeasures)
        self.MeasuresB.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=0)


        # a settings button that should lauch some global settings.
        self.SettingsB= ttk.Button(self, text='Settings', command=self.OpenSettings)
        self.SettingsB.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=0)

        self.RecVar = Tk.IntVar()
        self.RecVar.set(self.parent.recording)
        self.RecVar.trace('w', self.RecChanged)
        self.RecordFrames = ttk.Checkbutton(self, text = 'Record',
                                            variable = self.RecVar)
        self.RecordFrames.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=0)
        #attach the parameter to the Playbackbar
        self.param.attach(self)

    def RecChanged(self, *args):
        if self.RecVar.get() == self.parent.recording:
            pass
        else:
            self.parent.recording = self.RecVar.get()
            if self.parent.recording == 1:
                self.parent.PrintFig()


    def SkipLeft(self, e = None):
        self.param.set(self.param.value - self.skipSize)

    def SkipRight(self, e = None):
        self.param.set(self.param.value + self.skipSize)

    def PlayHandler(self, e = None):
        if not self.playPressed:
            # Set the value of play pressed to true, change the button name to
            # pause, turn off clear_fig, and start the play loop.
            self.playPressed = True
            self.parent.clear_fig = False
            self.playB.config(text='Pause')
            self.after(int(self.waitTime*1E3), self.blink)
        else:
            # pause the play loop, turn clear fig back on, and set the button name back to play
            self.playPressed = False
            self.parent.clear_fig = True
            self.playB.config(text='Play')

    def OpenSettings(self):
        if self.parent.settings_window is None:
            self.parent.settings_window = SettingsFrame(self.parent)
        else:
            self.parent.settings_window.destroy()
            self.parent.settings_window = SettingsFrame(self.parent)

    def OpenMeasures(self):
        if self.parent.measure_window is None:
            self.parent.measure_window = MeasureFrame(self.parent)
        else:
            self.parent.measure_window.destroy()
            self.parent.measure_window = MeasureFrame(self.parent)


    def blink(self):
        if self.playPressed:
            # First check to see if the timestep can get larger
            if self.param.value == self.param.maximum:
                # push pause button
                self.PlayHandler()

            # otherwise skip right by size skip size
            else:
                self.param.set(self.param.value + self.skipSize)

            # start loopin'
            self.after(int(self.waitTime*1E3), self.blink)


    def TextCallback(self):
        try:
            #make sure the user types in a int
            if int(self.tstep.get()) != self.param.value:
                self.param.set(int(self.tstep.get()))
        except ValueError:
            #if they type in random stuff, just set it ot the param value
            self.tstep.set(str(self.param.value))

    def ScaleHandler(self, e):
        # if changing the scale will change the value of the parameter, do so
        if self.param.value != int(self.slider.get()):
            self.param.set(int(self.slider.get()))

    def setKnob(self, value):
        #set the text entry value
        self.tstep.set(str(value))
        #set the slider
        self.slider.set(value)

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
        self.skipSize.set(self.parent.playbackbar.skipSize) # default value
        self.skipSize.trace('w', self.SkipSizeChanged)
        ttk.Label(frm, text="Skip Size:").grid(row=0)
        self.skipEnter = ttk.Entry(frm, textvariable=self.skipSize, width = 6)
        self.skipEnter.grid(row =0, column = 1, sticky = Tk.W + Tk.E)

        # Make an button to change the wait time
        self.waitTime = Tk.StringVar(self)
        self.waitTime.set(self.parent.playbackbar.waitTime) # default value
        self.waitTime.trace('w', self.WaitTimeChanged)
        ttk.Label(frm, text="Playback Wait Time:").grid(row=1)
        self.waitEnter = ttk.Entry(frm, textvariable=self.waitTime, width = 6)
        self.waitEnter.grid(row =1, column = 1, sticky = Tk.W + Tk.E)

        # Have a list of the color maps
        self.cmapList = new_cmaps.cmaps.keys()
        self.cmapvar = Tk.StringVar(self)
        self.cmapvar.set(self.parent.cmap) # default value
        self.cmapvar.trace('w', self.CmapChanged)

        ttk.Label(frm, text="Color map:").grid(row=2)
        cmapChooser = apply(ttk.OptionMenu, (frm, self.cmapvar, self.parent.cmap) + tuple(self.cmapList))
        cmapChooser.grid(row =2, column = 1, sticky = Tk.W + Tk.E)

        # Make an entry to change the number of columns
        self.columnNum = Tk.StringVar(self)
        self.columnNum.set(self.parent.numOfColumns.get()) # default value
        self.columnNum.trace('w', self.ColumnNumChanged)
        ttk.Label(frm, text="# of columns:").grid(row=3)
        self.ColumnSpin = Spinbox(frm,  from_=1, to=self.parent.maxCols, textvariable=self.columnNum, width = 6)
        self.ColumnSpin.grid(row =3, column = 1, sticky = Tk.W + Tk.E)

        # Make an entry to change the number of columns
        self.rowNum = Tk.StringVar(self)
        self.rowNum.set(self.parent.numOfRows.get()) # default value
        self.rowNum.trace('w', self.RowNumChanged)
        ttk.Label(frm, text="# of rows:").grid(row=4)
        self.RowSpin = Spinbox(frm, from_=1, to=self.parent.maxRows, textvariable=self.rowNum, width = 6)
        self.RowSpin.grid(row =4, column = 1, sticky = Tk.W + Tk.E)

        # Control whether or not Title is shown
        self.TitleVar = Tk.IntVar()
        self.TitleVar.set(self.parent.show_title)
        self.TitleVar.trace('w', self.TitleChanged)

        self.LimVar = Tk.IntVar()
        self.LimVar.set(self.parent.xlim[0])
        self.LimVar.trace('w', self.LimChanged)



        self.xleft = Tk.StringVar()
        self.xleft.set(str(self.parent.xlim[1]))
        self.xright = Tk.StringVar()
        self.xright.set(str(self.parent.xlim[2]))


        ttk.Label(frm, text = 'min').grid(row= 6, column = 1, sticky = Tk.N)
        ttk.Label(frm, text = 'max').grid(row= 6, column = 2, sticky = Tk.N)
        cb = ttk.Checkbutton(frm, text ='Set xlim',
                        variable = self.LimVar)
        cb.grid(row = 7, sticky = Tk.N)
        self.eLEnter = ttk.Entry(frm, textvariable=self.xleft, width = 8)
        self.eLEnter.grid(row = 7, column =1, sticky = Tk.N)
        self.eREnter = ttk.Entry(frm, textvariable=self.xright, width = 8)
        self.eREnter.grid(row = 7, column =2, sticky = Tk.N)


        self.yLimVar = Tk.IntVar()
        self.yLimVar.set(self.parent.ylim[0])
        self.yLimVar.trace('w', self.yLimChanged)



        self.yleft = Tk.StringVar()
        self.yleft.set(str(self.parent.ylim[1]))
        self.yright = Tk.StringVar()
        self.yright.set(str(self.parent.ylim[2]))


        cb = ttk.Checkbutton(frm, text ='Set ylim',
                        variable = self.yLimVar)
        cb.grid(row = 8, sticky = Tk.N)
        self.eLEnter = ttk.Entry(frm, textvariable=self.yleft, width = 8 )
        self.eLEnter.grid(row = 8, column =1, sticky = Tk.N)
        self.eREnter = ttk.Entry(frm, textvariable=self.yright, width =8 )
        self.eREnter.grid(row = 8, column =2, sticky = Tk.N)

        cb = ttk.Checkbutton(frm, text = "Show Title",
                        variable = self.TitleVar)
        cb.grid(row = 10, sticky = Tk.N)

        # Control whether or not axes are shared with a radio box:
        self.toLinkList = ['None', 'All spatial', 'All non p-x', 'All 2-D spatial']
        self.LinkedVar = Tk.IntVar()
        self.LinkedVar.set(self.parent.LinkSpatial)

        ttk.Label(frm, text='Share spatial axes:').grid(row = 0, column = 2, sticky = Tk.W)

        for i in range(len(self.toLinkList)):
            ttk.Radiobutton(frm,
                    text=self.toLinkList[i],
                    variable=self.LinkedVar,
                    command = self.RadioLinked,
                    value=i).grid(row = 1+i, column = 2, sticky =Tk.N)

        self.AspectVar = Tk.IntVar()
        self.AspectVar.set(self.parent.plot_aspect)
        self.AspectVar.trace('w', self.AspectVarChanged)

        cb = ttk.Checkbutton(frm, text = "Aspect = 1",
                                variable = self.AspectVar)
        cb.grid(row = 10, column = 1, sticky = Tk.N)

    def AspectVarChanged(self, *args):
        if self.AspectVar.get() == self.parent.plot_aspect:
            pass

        else:
            self.parent.plot_aspect = self.AspectVar.get()
            self.parent.RenewCanvas(ForceRedraw = True)


    def TitleChanged(self, *args):
        if self.TitleVar.get()==self.parent.show_title:
            pass
        else:
            self.parent.show_title = self.TitleVar.get()
            self.parent.RenewCanvas()

    def RadioLinked(self, *args):
        # If the shared axes are changed, the whole plot must be redrawn
        if self.LinkedVar.get() == self.parent.LinkSpatial:
            pass
        else:
            self.parent.LinkSpatial = self.LinkedVar.get()
            self.parent.RenewCanvas(ForceRedraw = True)


    def CmapChanged(self, *args):
    # Note here that Tkinter passes an event object to onselect()
        if self.cmapvar.get() == self.parent.cmap:
            pass
        else:
            self.parent.cmap = self.cmapvar.get()
            if self.parent.cmap in self.parent.cmaps_with_green:
                self.parent.ion_color =  new_cmaps.cmaps['plasma'](0.55)
                self.parent.electron_color = new_cmaps.cmaps['plasma'](0.8)
                self.parent.ion_fit_color = 'r'
                self.parent.electron_fit_color = 'yellow'

            else:
                self.parent.ion_color = new_cmaps.cmaps['viridis'](0.45)
                self.parent.electron_color = new_cmaps.cmaps['viridis'](0.75)
                self.parent.ion_fit_color = 'mediumturquoise'
                self.parent.electron_fit_color = 'lime'


            self.parent.RenewCanvas(ForceRedraw = True)


    def SkipSizeChanged(self, *args):
    # Note here that Tkinter passes an event object to onselect()
        try:
            if self.skipSize.get() == '':
                pass
            else:
                self.parent.playbackbar.skipSize = int(self.skipSize.get())
        except ValueError:
            self.skipSize.set(self.parent.playbackbar.skipSize)

    def RowNumChanged(self, *args):
    # Note here that Tkinter passes an event object to onselect()
        try:
            if self.rowNum.get() == '':
                pass
            elif int(self.rowNum.get())<1:
                self.rowNum.set(1)
            elif int(self.rowNum.get())>self.parent.maxRows:
                self.rowNum.set(self.parent.maxRows)
            else:
                self.parent.numOfRows.set(int(self.rowNum.get()))
        except ValueError:
            self.rowNum.set(self.parent.numOfRows.get())

    def ColumnNumChanged(self, *args):
    # Note here that Tkinter passes an event object to onselect()
        try:
            if self.columnNum.get() == '':
                pass
            elif int(self.columnNum.get())<1:
                self.columnNum.set(1)
            elif int(self.columnNum.get())>self.parent.maxCols:
                self.columnNum.set(self.parent.maxCols)

            else:
                self.parent.numOfColumns.set(int(self.columnNum.get()))
        except ValueError:
            self.columnNum.set(self.parent.numOfColumns.get())

    def WaitTimeChanged(self, *args):
    # Note here that Tkinter passes an event object to onselect()
        try:
            if self.waitTime.get() == '':
                pass
            else:
                self.parent.playbackbar.waitTime = float(self.waitTime.get())
        except ValueError:
            self.waitTime.set(self.parent.playbackbar.waitTime)

    def CheckIfXLimChanged(self):
        to_reload = False
        tmplist = [self.xleft, self.xright]
        for j in range(2):

            try:
            #make sure the user types in a a number and that it has changed.
                if np.abs(float(tmplist[j].get()) - self.parent.xlim[j+1]) > 1E-4:
                    self.parent.xlim[j+1] = float(tmplist[j].get())
                    to_reload += True

            except ValueError:
                #if they type in random stuff, just set it ot the param value
                tmplist[j].set(str(self.parent.xlim[j+1]))
        return to_reload*self.parent.xlim[0]

    def CheckIfYLimChanged(self):
        to_reload = False
        tmplist = [self.yleft, self.yright]
        for j in range(2):

            try:
            #make sure the user types in a int
                if np.abs(float(tmplist[j].get()) - self.parent.ylim[j+1]) > 1E-4:
                    self.parent.ylim[j+1] = float(tmplist[j].get())
                    to_reload += True

            except ValueError:
                #if they type in random stuff, just set it ot the param value
                tmplist[j].set(str(self.parent.ylim[j+1]))
        return to_reload*self.parent.ylim[0]

    def LimChanged(self, *args):
        if self.LimVar.get()==self.parent.xlim[0]:
            pass
        else:
            self.parent.xlim[0] = self.LimVar.get()
            self.parent.RenewCanvas()

    def yLimChanged(self, *args):
        if self.yLimVar.get()==self.parent.ylim[0]:
            pass
        else:
            self.parent.ylim[0] = self.yLimVar.get()
            self.parent.RenewCanvas()

    def SettingsCallback(self, e):
        to_reload = self.CheckIfXLimChanged()
        to_reload += self.CheckIfYLimChanged()

        if to_reload:
            self.parent.RenewCanvas()


    def OnReload(self, event=None):
        self.parent.findDir()

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

        # Make an entry to change the integration region
        # A StringVar for a box to type in a value for the left ion region
        self.ileft = Tk.StringVar()
        # set it to the left value
        self.ileft.set(str(self.parent.i_L.get()))

        # A StringVar for a box to type in a value for the right ion region
        self.iright = Tk.StringVar()
        # set it to the right value
        self.iright.set(str(self.parent.i_R.get()))

        # Now the electrons
        self.eleft = Tk.StringVar()
        self.eleft.set(str(self.parent.e_L.get()))
        self.eright = Tk.StringVar()
        self.eright.set(str(self.parent.e_R.get()))

        ttk.Label(frm, text='Energy Int region:').grid(row = 0, sticky = Tk.W)
        ttk.Label(frm, text='left').grid(row = 0, column = 1, sticky = Tk.N)
        ttk.Label(frm, text='right').grid(row = 0, column = 2, sticky = Tk.N)

        # the ion row
        ttk.Label(frm, text='ions').grid(row= 1, sticky = Tk.N)
        # Make an button to change the wait time

        self.iLEnter = ttk.Entry(frm, textvariable=self.ileft, width=7)
        self.iLEnter.grid(row =1, column = 1, sticky = Tk.W + Tk.E)

        self.iREnter = ttk.Entry(frm, textvariable=self.iright, width=7)
        self.iREnter.grid(row = 1, column =2, sticky = Tk.W + Tk.E)

        ttk.Label(frm, text='electrons').grid(row= 2, sticky = Tk.N)
        self.eLEnter = ttk.Entry(frm, textvariable=self.eleft, width=7)
        self.eLEnter.grid(row = 2, column =1, sticky = Tk.W + Tk.E)
        self.eREnter = ttk.Entry(frm, textvariable=self.eright, width=7)
        self.eREnter.grid(row = 2, column =2, sticky = Tk.W + Tk.E)

        self.RelVar = Tk.IntVar()
        self.RelVar.set(self.parent.e_relative)
        self.RelVar.trace('w', self.RelChanged)
        cb = ttk.Checkbutton(frm, text = "Energy Region relative to shock?",
                        variable = self.RelVar)
        cb.grid(row = 3, columnspan = 3, sticky = Tk.W)

        self.SetTeVar = Tk.IntVar()
        self.SetTeVar.set(self.parent.set_Te)
        self.SetTeVar.trace('w', self.SetTeChanged)
        cb = ttk.Checkbutton(frm, text='Show T_e', variable =  self.SetTeVar)
        cb.grid(row = 5, sticky = Tk.W)

        ttk.Label(frm, text=u'\u0394'+u'\u0263' + ' =').grid(row= 5, column =1, sticky = Tk.N)

        self.SetTpVar = Tk.IntVar()
        self.SetTpVar.set(self.parent.set_Tp)
        self.SetTpVar.trace('w', self.SetTpChanged)

        cb = ttk.Checkbutton(frm, text='Show T_i', variable =  self.SetTpVar)
        cb.grid(row = 6, sticky = Tk.W)
        ttk.Label(frm, text=u'\u0394'+u'\u0263' + ' =').grid(row= 6, column =1, sticky = Tk.N)

        self.delgameVar = Tk.StringVar()
        self.delgameVar.set(str(self.parent.delgam_e))
        self.delgampVar = Tk.StringVar()
        self.delgampVar.set(str(self.parent.delgam_p))


        ttk.Entry(frm, textvariable=self.delgameVar, width = 7).grid(row = 5, column = 2, sticky = Tk.N)
        ttk.Entry(frm, textvariable=self.delgampVar, width = 7).grid(row = 6, column =2, sticky = Tk.N)

        ttk.Label(frm, text='Powerlaw fits:').grid(row = 8, sticky = Tk.W)
        ttk.Label(frm, text='E_min [mc^2]').grid(row = 8, column = 1, sticky = Tk.N)
        ttk.Label(frm, text='E_max [mc^2]').grid(row = 8, column = 2, sticky = Tk.N)

        self.PLFitEVar = Tk.IntVar()
        self.PLFitEVar.set(self.parent.PowerLawFitElectron[0])
        self.PLFitEVar.trace('w', self.PLFitEChanged)
        ttk.Checkbutton(frm, text='Electrons', variable =  self.PLFitEVar).grid(row = 9, sticky = Tk.W)

        self.E1Var = Tk.StringVar()
        self.E1Var.set(str(self.parent.PowerLawFitElectron[1]))
        self.E2Var = Tk.StringVar()
        self.E2Var.set(str(self.parent.PowerLawFitElectron[2]))


        ttk.Entry(frm, textvariable=self.E1Var, width = 7).grid(row = 9, column = 1, sticky = Tk.N)
        ttk.Entry(frm, textvariable=self.E2Var, width = 7).grid(row = 9, column =2, sticky = Tk.N)


        self.PLFitPVar = Tk.IntVar()
        self.PLFitPVar.set(self.parent.PowerLawFitIon[0])
        self.PLFitPVar.trace('w', self.PLFitPChanged)
        ttk.Checkbutton(frm, text='Ions', variable =  self.PLFitPVar).grid(row = 10, sticky = Tk.W)

        self.P1Var = Tk.StringVar()
        self.P1Var.set(str(self.parent.PowerLawFitIon[1]))
        self.P2Var = Tk.StringVar()
        self.P2Var.set(str(self.parent.PowerLawFitIon[2]))

        ttk.Entry(frm, textvariable=self.P1Var, width = 7).grid(row = 10, column = 1, sticky = Tk.N)
        ttk.Entry(frm, textvariable=self.P2Var, width = 7).grid(row = 10, column =2, sticky = Tk.N)


    def PLFitEChanged(self, *args):
        if self.PLFitEVar.get() == self.parent.PowerLawFitElectron[0]:
            pass
        else:
            self.parent.PowerLawFitElectron[0] = self.PLFitEVar.get()
            self.parent.RenewCanvas()

    def PLFitPChanged(self, *args):
        if self.PLFitPVar.get() == self.parent.PowerLawFitIon[0]:
            pass
        else:
            self.parent.PowerLawFitIon[0] = self.PLFitPVar.get()
            self.parent.RenewCanvas()

    def CheckIfTeChanged(self):
        to_reload = False
        try:
        #make sure the user types in a int
            if np.abs(float(self.delgameVar.get()) - self.parent.delgam_e) > 1E-4:
                self.parent.delgam_e = float(self.delgameVar.get())
                to_reload += True*self.parent.set_Te

        except ValueError:
            #if they type in random stuff, just set it ot the param value
            self.delgameVar.set(str(self.parent.delgam_e))
        return to_reload

    def CheckIfTpChanged(self):
        to_reload = False
        try:
        #make sure the user types in a flof
            if np.abs(float(self.delgampVar.get()) - self.parent.delgam_p) > 1E-4:
                    self.parent.delgam_p = float(self.delgampVar.get())
                    to_reload += True*self.parent.set_Tp

        except ValueError:
            #if they type in random stuff, just set it ot the param value
            self.delgampVar.set(str(self.parent.delgam_p))
        return to_reload


    def CheckIfPLChanged(self):
        to_reload = False

        PLList = [self.parent.PowerLawFitElectron, self.parent.PowerLawFitIon]
        VarList = [[self.E1Var, self.E2Var], [self.P1Var, self.P2Var]]

        for j in range(2):
            for k in range(2):
                try:
                    #make sure the user types in a int
                    if np.abs(float(VarList[j][k].get())- PLList[j][k+1])>1E-6:

                        PLList[j][k+1] = float(VarList[j][k].get())
                        to_reload += True

                except ValueError:
                    #if they type in random stuff, just set it to the value
                    VarList[j][k].set(str(PLList[j][k+1]))
        return to_reload

    def CheckIfIntChanged(self, tkVar, valVar):
        to_reload = False
        try:
            #make sure the user types in a int
            if int(tkVar.get()) != valVar.get():
                valVar.set(int(tkVar.get()))
                to_reload = True
            return to_reload
        except ValueError:
            #if they type in random stuff, just set it ot the param value
            tkVar.set(str(valVar.get()))
            return to_reload

    def SetTeChanged(self, *args):
        if self.SetTeVar.get()==self.parent.set_Te:
            pass
        else:
            self.parent.set_Te = self.SetTeVar.get()
            self.parent.RenewCanvas()

    def SetTpChanged(self, *args):
        if self.SetTpVar.get()==self.parent.set_Tp:
            pass
        else:
            self.parent.set_Tp = self.SetTpVar.get()
            self.parent.RenewCanvas()

    def TxtEnter(self, e):
        self.MeasuresCallback()

    def RelChanged(self, *args):
        if self.RelVar.get()==self.parent.e_relative:
            pass
        else:
            self.parent.e_relative = self.RelVar.get()
            self.parent.RenewCanvas()



    def MeasuresCallback(self):
        tkvarIntList = [self.ileft, self.iright, self.eleft, self.eright]
        IntValList = [self.parent.i_L, self.parent.i_R, self.parent.e_L, self.parent.e_R]

        to_reload = False



        for j in range(len(tkvarIntList)):
            to_reload += self.CheckIfIntChanged(tkvarIntList[j], IntValList[j])

        to_reload += self.CheckIfTeChanged()
        to_reload += self.CheckIfTpChanged()
        to_reload += self.CheckIfPLChanged()
        if to_reload:
            self.parent.RenewCanvas()

    def OnClosing(self):
        self.parent.settings_window = None
        self.destroy()


class MainApp(Tk.Tk):
    """ We simply derive a new class of Frame as the man frame of our app"""
    def __init__(self, name):

        Tk.Tk.__init__(self)
        self.update_idletasks()
        menubar = Tk.Menu(self)
        self.wm_title(name)
        self.settings_window = None
        self.measure_window = None
        self.clear_fig = True # A parameter that causes the graph to disappear as soon as something is pressed. Here is the dictionary of the keys:

        self.first_x = None
        self.first_y = None

        self.num_font_size = 11
        self.recording = False
        #
        # Set the number of rows and columns in the figure
        # (As well as the max rows)
        self.maxRows = 5
        self.maxCols = 3

        # a list of cmaps with orange prtl colors
        self.cmaps_with_green = ['viridis', 'Rainbow + White', 'Blue/Green/Red/Yellow', 'Cube YF', 'Linear_L']
        # A param to define whether to share axes
        # 0 : No axes are shared
        # 1 : All axes are shared
        # 2 : All non p-x plots are shared
        # 3 : All 2-D, non p-x plots are shared
        self.LinkSpatial = 2

        # A param that sets the aspect ratio for the spatial, 2d plots
        self.plot_aspect = 0

        # A param that will tell us if we want to set the E_temp
        self.set_Te = False
        self.delgam_e = 0.03

        # A param that will tell us if we want to set the p_temp
        self.set_Tp = False
        self.delgam_p = 0.06


        self.numOfRows = Tk.IntVar(self)
#        self.numOfRows.set(3)
        self.numOfRows.set(2)
        self.numOfRows.trace('w', self.UpdateGridSpec)
        self.numOfColumns = Tk.IntVar(self)
        self.numOfColumns.set(2)
        self.numOfColumns.trace('w', self.UpdateGridSpec)
        self.SubPlotParams = {'left':0.06, 'right':0.95, 'top':.93, 'bottom':0.06, 'wspace':0.15, 'hspace':0.15}
        matplotlib.rc('figure.subplot', **self.SubPlotParams)

        self.show_title = True
        self.xlabel_pad = 0
        self.ylabel_pad = 0
        fileMenu = Tk.Menu(menubar, tearoff=False)
        menubar.add_cascade(label="File", underline=0, menu=fileMenu)
        fileMenu.add_command(label= 'Open Directory', command = self.OnOpen, accelerator='Command+o')
        fileMenu.add_command(label="Exit", underline=1,
                             command=quit, accelerator="Ctrl+Q")
        self.config(menu=menubar)

        self.bind_all("<Control-q>", self.quit)
        self.bind_all("<Command-o>", self.OnOpen)


        # create a bunch of regular expressions used to search for files
        f_re = re.compile('flds.tot.*')
        prtl_re = re.compile('prtl.tot.*')
        s_re = re.compile('spect.*')
        param_re = re.compile('param.*')
        self.re_list = [f_re, prtl_re, s_re, param_re]


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

        # Set the default color map

        self.cmap = 'viridis'

        # Create the figure
        self.f = Figure(figsize = (2,2), dpi = 100)

        # a tk.DrawingArea
        self.canvas = FigureCanvasTkAgg(self.f, master=self)


        # Make the object hold the timestep info
        self.TimeStep = Param(1, minimum=1, maximum=1000)
        self.playbackbar = PlaybackBar(self, self.TimeStep, canvas = self.canvas)

        # Add the toolbar
        self.toolbar =  MyCustomToolbar(self.canvas, self)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=Tk.RIGHT, fill=Tk.BOTH, expand=1)

        # Look for the tristan output files and load the file paths into
        # previous objects
        self.dirname = os.curdir
        self.findDir()

        # Choose the integration region for the particles
        self.e_relative = True
        self.i_L = Tk.IntVar()
        self.i_L.set(-1E4)
        self.i_R = Tk.IntVar()
        self.i_R.set(0)
        self.e_L = Tk.IntVar()
        self.e_L.set(-1E4)
        self.e_R = Tk.IntVar()
        self.e_R.set(0)

        # Whether or not to set a xlim or ylim
        self.xlim = [False, 0, 100]
        self.ylim = [False, 0, 100]

        self.PowerLawFitElectron = [False, 1.0, 10.0]
        self.PowerLawFitIon = [False, 1.0, 10.0]
        # Set the particle colors
        if self.cmap in self.cmaps_with_green:
            self.shock_color = 'w'
            self.ion_color =  new_cmaps.cmaps['plasma'](0.55)
            self.electron_color = new_cmaps.cmaps['plasma'](0.8)
            self.ion_fit_color = 'r'
            self.electron_fit_color = 'yellow'

        else:
            self.shock_color = 'w'
            self.ion_color = new_cmaps.cmaps['viridis'](0.45)
            self.electron_color = new_cmaps.cmaps['viridis'](0.75)
            self.ion_fit_color = 'mediumturquoise'
            self.electron_fit_color = 'limegreen'



        self.TimeStep.attach(self)
        self.DrawCanvas()


        self.playbackbar.pack(side=Tk.TOP, fill=Tk.BOTH, expand=0)
        self.update()
        # now root.geometry() returns valid size/placement
        self.minsize(self.winfo_width(), self.winfo_height())
        self.geometry("1200x700")
        self.bind('<Return>', self.TxtEnter)
        self.bind('<Left>', self.playbackbar.SkipLeft)
        self.bind('<Right>', self.playbackbar.SkipRight)
        self.bind('<space>', self.playbackbar.PlayHandler)

    def quit(self, event):
        print("quitting...")
        sys.exit(0)

    def GenH5Dict(self):
        for pkey in self.PathDict.keys():
            with h5py.File(os.path.join(self.dirname,self.PathDict[pkey][0]), 'r') as f:
                # Because dens is in both spect* files and flds* files,
                for h5key in f.keys():
                    if h5key == 'dens' and pkey == 'Spect':
                        self.H5KeyDict['spect_dens'] = pkey
                    else:
                        self.H5KeyDict[h5key] = pkey

        print self.H5KeyDict

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
        self.PathDict['Flds']= filter(f_re.match, os.listdir(self.dirname))
        self.PathDict['Flds'].sort()

        self.PathDict['Prtl']= filter(prtl_re.match, os.listdir(self.dirname))
        self.PathDict['Prtl'].sort()

        self.PathDict['Spect']= filter(s_re.match, os.listdir(self.dirname))
        self.PathDict['Spect'].sort()

        self.PathDict['Param']= filter(param_re.match, os.listdir(self.dirname))
        self.PathDict['Param'].sort()

        for key in self.PathDict.keys():
            is_okay &= len(self.PathDict[key]) > 0

        if is_okay:
            self.NewDirectory = True
            self.TimeStep.setMax(max(len(self.PathDict['Flds']),1))
            self.playbackbar.slider.config(to =(len(self.PathDict['Flds'])))
            self.shock_finder()

        return is_okay


    def OnOpen(self, e = None):
        """open a file"""


        tmpdir = tkFileDialog.askdirectory(title = 'Choose the directory of the output files', **self.dir_opt)
        if tmpdir == '':
            self.findDir()

        else:
            self.dirname = tmpdir
        if not self.pathOK():
#            p = MyDalog(self, 'Directory must contain either the output directory or all of the following: \n flds.tot.*, ptrl.tot.*, params.*, spect.*', title = 'Cannot find output files')
#            self.wait_window(p.top)
            self.findDir()

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
            tmpdir = tkFileDialog.askdirectory(title = dlgstr, **self.dir_opt)
            if tmpdir != '':
                self.dirname = tmpdir
            if not self.pathOK():
#                p = MyDialog(self, 'Directory must contain either the output directory or all of the following: \n flds.tot.*, ptrl.tot.*, params.*, spect.*', title = 'Cannot find output files')
#                self.wait_window(p.top)
                self.findDir()


    def DrawCanvas(self):
        '''Initializes the figure, and then packs it into the main window.
        Should only be called once.'''

        # figsize (w,h tuple in inches) dpi (dots per inch)
        #f = Figure(figsize=(5,4), dpi=100)

        # Generate all of the subplot wrappers. They are stored in a 2D list
        # where the index [i][j] corresponds to the ith row, jth column

        # divy up the figure into a bunch of subplots using GridSpec.
        self.gs0 = gridspec.GridSpec(self.numOfRows.get(),self.numOfColumns.get())


        # Create the list of all of subplot wrappers
        self.SubPlotList = []
        for i in range(self.maxRows):
            tmplist = [SubPlotWrapper(self, figure = self.f, pos=(i,j)) for j in range(self.maxCols)]
            self.SubPlotList.append(tmplist)
        for i in range(self.maxRows):
            for j in range(self.maxCols):
                self.SubPlotList[i][j].SetGraph('PhasePlot')

        self.SubPlotList[0][1].PlotParamsDict['PhasePlot']['prtl_type'] = 1

        self.SubPlotList[1][1].SetGraph('FieldsPlot')
        self.SubPlotList[2][1].SetGraph('FieldsPlot')
        self.SubPlotList[2][1].PlotParamsDict['FieldsPlot']['field_type'] = 1

        self.SubPlotList[1][0].SetGraph('DensityPlot')
        self.SubPlotList[2][0].SetGraph('SpectraPlot')

        # Make a list that will hold the previous ctype
        self.MakePrevCtypeList()
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        self.ReDrawCanvas()


        self.f.canvas.mpl_connect('button_press_event', self.onclick)
#        self.canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        #self.label = Label(self.top, text = 'Text',bg='orange')
        #self.label.grid()
        # initialize (time index t)

    def UpdateGridSpec(self, *args):
        '''A function that handles updates the gridspec that divides up of the
        plot into X x Y subplots'''
        self.gs0 = gridspec.GridSpec(self.numOfRows.get(),self.numOfColumns.get())
        self.RenewCanvas(keep_view = False, ForceRedraw = True)

    def LoadAllKeys(self):
        ''' A function that will find out will arrays need to be loaded for
        to draw the graphs. If the time hasn't changed, it will only load new keys.'''
        # Make a dictionary that stores all of the keys we will need to load
        # to draw the graphs.
        self.ToLoad = {'Flds': [], 'Prtl': [], 'Param': [], 'Spect': []}
        for i in range(self.numOfRows.get()):
            for j in range(self.numOfColumns.get()):
                # for each subplot, see what keys are needed
                tmpList = self.SubPlotList[i][j].GetKeys()
                # we always load time because it is needed to calculate the shock location
                self.ToLoad[self.H5KeyDict['time']].append('time')
                # We always load enough to calculate xmin, xmax, ymin, ymax:
                self.ToLoad[self.H5KeyDict['c_omp']].append('c_omp')
                self.ToLoad[self.H5KeyDict['istep']].append('istep')
                self.ToLoad[self.H5KeyDict['dens']].append('dens')

                for elm in tmpList:
                    # find out what type of file the key is stored in
                    ftype = self.H5KeyDict[elm]
                    # add the key to the list of that file type
                    self.ToLoad[ftype].append(elm)
        # See if we are in a new Directory
        if self.NewDirectory:
            # Make a list of timesteps we have already loaded.
            self.timestep_visited = []

            # Timestep queue that ensures that we delete the oldest viewed
            # timestep if memory gets too large
            self.timestep_queue = deque()

            # For each timestep we visit, we will load a dictionary and place it in a list
            self.ListOfDataDict = []

            self.NewDirectory = False


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
                                else:
                                    raise

            self.timestep_queue.append(self.TimeStep.value)

        else:
            # The time has changed so we have to reload everything
            self.DataDict = {}
            for pkey in self.ToLoad.keys():
                tmplist = list(set(self.ToLoad[pkey])) # get rid of duplicate keys
                # Load the file
                if len(tmplist) > 0: #check we actually have something to load
                    with h5py.File(os.path.join(self.dirname,self.PathDict[pkey][self.TimeStep.value-1]), 'r') as f:
                        for elm in tmplist:
                        # Load all the keys
                            try:
                                if elm == 'spect_dens':
                                    self.DataDict[elm] = f['dens'][:]
                                else:
                                    self.DataDict[elm] = f[elm][:]
                            except KeyError:
                                if elm == 'sizex':
                                    self.DataDict[elm] = 1
                                else:
                                    raise
            # don't keep more than 15 time steps in memory because of RAM issues
            if len(self.timestep_visited)>30:
                oldest_time = self.timestep_queue.popleft()
                oldest_ind = self.timestep_visited.index(oldest_time)
                self.timestep_visited.remove(oldest_time)
                self.ListOfDataDict.pop(oldest_ind)
            self.timestep_visited.append(self.TimeStep.value)
            self.ListOfDataDict.append(self.DataDict)
            self.timestep_queue.append(self.TimeStep.value)


    def MakePrevCtypeList(self):
        self.prev_ctype_list = []
        for i in range(self.numOfRows.get()):
            tmp_ctype_l = []
            for j in range(self.numOfColumns.get()):
                tmp_ctype_l.append(str(self.SubPlotList[i][j].chartType))
            self.prev_ctype_list.append(tmp_ctype_l)

    def SaveView(self):
        # A function that will make sure our view will stay the same as the
        # plot updates.
        cur_view =  list(self.toolbar._views.__call__())
        # Go to the home view
        self.toolbar._views.home()
        self.toolbar._positions.home()
        home_view =  list(self.toolbar._views.__call__())

        # Filter out the colorbar axes
        cbar_loc = []
        # find the axes they are tuples that are like (0, 1, 0, 1) (but floats)
        for k in range(len(cur_view)):
            first_zero = np.abs(cur_view[k][0])<1E-10
            second_one = np.abs(cur_view[k][1]-1)<1E-10
            third_zero = np.abs(cur_view[k][2])<1E-10
            forth_one = np.abs(cur_view[k][3]-1)<1E-10
            if first_zero and second_one and third_zero and forth_one:
                cbar_loc.append(k)

        # remove all of them
        for j in range(len(cbar_loc))[::-1]:
            cur_view.pop(cbar_loc[j])
            home_view.pop(cbar_loc[j])

        self.is_changed_list = []
        self.old_views = []
        if cur_view is not None:
            for i in range(len(cur_view)):
                is_changed =[]
                for j in range(4):
                    is_changed.append(np.abs(home_view[i][j]-cur_view[i][j])>1E-5)
                self.is_changed_list.append(is_changed)
                self.old_views.append(cur_view[i])

    def LoadView(self):

        self.toolbar._views.clear()
        self.toolbar.push_current()
        next_view = list(self.toolbar._views.__call__())


        # put the parts that have changed from the old view
        # into the proper place in the next view
        m = 0 # a counter that allows us to go from labeling the plots in [i][j] to 1d
        k = 0 # A counter that skips over the cbar axes in next_view
        for i in range(self.numOfRows.get()):
            for j in range(self.numOfColumns.get()):
                tmp_old_view = list(self.old_views.pop(0))
                tmp_new_view = list(next_view[k])
                if self.prev_ctype_list[i][j] == self.SubPlotList[i][j].chartType:
                    is_changed = self.is_changed_list[m]
                    if self.SubPlotList[i][j].Changedto2D or self.SubPlotList[i][j].Changedto1D:
                        # only keep the x values if they have changed
                        for n in range(2):
                            if is_changed[n]:
                                tmp_new_view[n] = tmp_old_view[n]
                    else:
                        # Keep any y or x that is changed
                        for n in range(4):
                            if is_changed[n]:
                                tmp_new_view[n] = tmp_old_view[n]
                next_view[k] = tmp_new_view
                # Handle the counting of the 'views' array in matplotlib
                if self.SubPlotList[i][j].GetPlotParam('twoD') == 1:
                    if self.SubPlotList[i][j].GetPlotParam('show_cbar') == 1:
                        k += 2
                    else:
                        k += 1
                else:
                    k += 1
                m += 1
                self.SubPlotList[i][j].Changedto1D = False
                self.SubPlotList[i][j].Changedto2D = False

        self.toolbar._views.push(next_view)
        self.toolbar.set_history_buttons()
        self.toolbar._update_view()

        # We must now figure out the pos of the charts that changed in the earlier view.


    def RenewCanvas(self, keep_view = True, ForceRedraw = False):

        '''We have two way of updated the graphs: one) by refreshing them using
        self.RefreshCanvas, we don't recreate all of the artists that matplotlib
        needs to make the plot work. self.RefreshCanvas should be fast. Two we
        can ReDraw the canvas using self.ReDrawCanvas. This recreates all the
        artists and will be slow. Sometimes the graph must be redrawn however,
        if the GridSpec changed, more plots are added, the chartype changed, if
        the plot went from 2d to 1D, etc.. If any change occurs that requires a
        redraw, renewcanvas must be called with ForceRedraw = True. '''

        if ForceRedraw:
            self.ReDrawCanvas(keep_view = keep_view)
        else:
            self.RefreshCanvas(keep_view = keep_view)

    def ReDrawCanvas(self, keep_view = True):
        #  We need to see if the user has moved around the zoom level in python.
        # First we see if there are any views in the toolbar
        cur_view =  self.toolbar._views.__call__()

        if cur_view is None:
            keep_view = False
        if self.NewDirectory:
            keep_view = False
        if keep_view:
            self.SaveView()
        self.f.clf()
        #
        if self.clear_fig:
            self.canvas.show()

        self.LoadAllKeys()

        # Calculate the new shock location
        self.shock_loc = self.DataDict['time'][0]*self.shock_speed

        # Calculate the new xmin, and xmax

        # Find the first position with a physical x and y direction:
        self.first_x = None
        self.first_y = None
        k = 0
        for i in range(self.numOfRows.get()):
            for j in range(self.numOfColumns.get()):

                # First handle the axes sharing
                if self.SubPlotList[i][j].chartType == 'SpectraPlot':
                    # The plot type is a spectral plot, which has no spatial dim
                    pass
                elif self.LinkSpatial != 1 and self.SubPlotList[i][j].chartType == 'PhasePlot':
                    # If this is the case we don't care about the phase plots
                    # as we don't want to share the axes
                    pass
                elif self.LinkSpatial == 3 and self.SubPlotList[i][j].GetPlotParam('twoD'):
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

                # Now find the position in the views_list
                if self.SubPlotList[i][j].chartType == 'SpectraPlot':
                    self.SubPlotList[i][j].pos_in_views = k
                    k += 1
                elif self.SubPlotList[i][j].chartType == 'PhasePlot':
                    self.SubPlotList[i][j].pos_in_views = k
                    k += 2
                elif  self.SubPlotList[i][j].GetPlotParam('twoD'):
                    self.SubPlotList[i][j].pos_in_views = k
                    k += 2
                else:
                    self.SubPlotList[i][j].pos_in_views = k
                    k += 1


                # Now... We can draw the graph.
                self.SubPlotList[i][j].DrawGraph()
        if self.show_title:
            self.f.suptitle(os.path.abspath(self.dirname)+ ' at time t = %d $\omega_{pe}$'  % round(self.DataDict['time'][0]), size = 15)

        if keep_view:
            self.LoadView()

        self.MakePrevCtypeList()
        self.canvas.show()
        self.canvas.get_tk_widget().update_idletasks()

        if self.recording:
            self.PrintFig()


    def RefreshCanvas(self, keep_view = True):
        tic = time.time()
        #  We need to see if the user has moved around the zoom level in python.
        # First we see if there are any views in the toolbar
        cur_view =  self.toolbar._views.__call__()

        if cur_view is None:
            keep_view = False
        if self.NewDirectory:
            keep_view = False
        if keep_view:
            self.SaveView()

        self.toolbar._views.clear()

        self.LoadAllKeys()

        # Calculate the new shock location
        self.shock_loc = self.DataDict['time'][0]*self.shock_speed

        # By design, the first_x and first_y cannot change if the graph is
        # being refreshed. Any call that would require this needs a redraw
        # Find the first position with a physical x and y direction:

        k = 0
        for i in range(self.numOfRows.get()):
            for j in range(self.numOfColumns.get()):
                # Now find the position in the views_list
                if self.SubPlotList[i][j].chartType == 'SpectraPlot':
                    self.SubPlotList[i][j].pos_in_views = k
                    k += 1
                elif self.SubPlotList[i][j].chartType == 'PhasePlot':
                    self.SubPlotList[i][j].pos_in_views = k
                    k += 2
                elif  self.SubPlotList[i][j].GetPlotParam('twoD'):
                    self.SubPlotList[i][j].pos_in_views = k
                    k += 2
                else:
                    self.SubPlotList[i][j].pos_in_views = k
                    k += 1


                # Now... We can refresh the graph.
                self.SubPlotList[i][j].RefreshGraph()
        if self.show_title:
            self.f.suptitle(os.path.abspath(self.dirname)+ ' at time t = %d $\omega_{pe}$'  % round(self.DataDict['time'][0]), size = 15)

        if keep_view:
            self.LoadView()

        self.MakePrevCtypeList()
        self.canvas.draw()
#        self.canvas.blit(bbox = self.f.bbox)
#        self.canvas.blit()
        self.canvas.get_tk_widget().update_idletasks()
        toc = time.time()
        print tic-toc
        if self.recording:
            self.PrintFig()

    def PrintFig(self):
        movie_dir = os.path.abspath(os.path.join(self.dirname, '..', 'Movie'))
        try:
            os.makedirs(movie_dir)
        except OSError:
            if not os.path.isdir(movie_dir):
                raise

        fname = 'iseult_img_'+ str(self.TimeStep.value).zfill(3)+'.png'
        self.f.savefig(os.path.join(movie_dir, fname))

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
            i = (len(row_array)-row_array.searchsorted(y_loc))/2
            j = col_array.searchsorted(x_loc)/2

            self.SubPlotList[i][j].OpenSubplotSettings()

    def shock_finder(self):
        '''The main idea of the shock finder, is we go to the last timestep
        in the simulation and find where the density is half it's max value.
        We then calculate the speed of the of the shock assuming it is
        traveling at constant velocity. The function returns an array that contains
        the shock position at every time step.'''

        # First load the first field file to find the initial size of the
        # box in the x direction

#        print os.path.join(self.dirname,self.PathDict['Flds'][0])
        with h5py.File(os.path.join(self.dirname,self.PathDict['Flds'][0]), 'r') as f:
            nxf0 = f['by'][:].shape[1]

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
        jstart = min(10*c_omp/istep, nxf0)
        # build the final x_axis of the plot

        xaxis_final = np.arange(dens_arr.shape[1])/c_omp*istep
        # Find the shock by seeing where the density is 1/2 of it's
        # max value.

        dens_half_max = max(dens_arr[dens_arr.shape[0]/2,jstart:])*.5

        # Find the farthest location where the average density is greater
        # than half max
        ishock_final = np.where(dens_arr[dens_arr.shape[0]/2,jstart:]>=dens_half_max)[0][-1]
        xshock_final = xaxis_final[ishock_final]
        self.shock_speed = xshock_final/final_time


    def setKnob(self, value):
        # If the time parameter changes update the plots
        self.RenewCanvas()

    def TxtEnter(self, e):
        self.playbackbar.TextCallback()

if __name__ == "__main__":
    app = MainApp('Iseult')
    app.mainloop()
