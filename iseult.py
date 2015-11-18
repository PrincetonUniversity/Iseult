#! /usr/bin/env pythonw

import re # regular expressions
import os, sys # Used to make the code portable
import h5py # Allows us the read the data files
from threading import Thread
import time,string
import matplotlib
import new_cmaps
import numpy as np
import  wx.lib.intctrl
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec

matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

import Tkinter as Tk
import ttk as ttk
import tkFileDialog

def destroy(e):
    sys.exit()

colors=[None]*10
for i in range(len(colors)):
    colors[i]=np.random.rand(5,5)


class SubPlotWrapper:
    """A simple class that will eventually hold all of the information
    about each sub_plot in the Figure"""

    def __init__(self, parent, figure=None, sub_plot_spec = None, ctype=None):
        self.parent = parent
        self.chartType = ctype
        # A dictionary that contains all of the plot types.
        self.PlotTypeDict = {'PhasePlot': PhasePanel, 'TestPlot': TestPanel}
        # A dictionary that will store where everything is in Hdf5 Files
        self.GenParamDict()
        self.graph = graph


    def LoadKey(self, h5key):
        pkey = self.parent.H5KeyDict[h5key]
        with h5py.File(os.path.join(self.parent.dirname,self.parent.PathDict[pkey][self.parent.timeStep.value-1]), 'r') as f:
            return f[h5key][:]

    def ChangeGraph(self, str_arg):
        self.chartType = str_arg
        self.parent.ChangeGraph()

    def GenParamDict(self):
        # Generate a dictionary that will store all of the params at dict['ctype']['param_name']
        self.PlotParamsDict = {elm: '' for elm in self.PlotTypeDict.keys() }
        for elm in self.PlotTypeDict.keys():
            self.PlotParamsDict[elm] = {x : '' for x in self.PlotTypeDict[elm].plot_param_list}

    def SetPlotParam(self, pname, val, ctype = None):
        if ctype is None:
            ctype = self.chartType
        self.PlotParamsDict[ctype][pname] = val

    def GetPlotParam(self, pname, ctype = None):
        if ctype is None:
            ctype = self.chartType

        return self.PlotParamsDict[ctype][pname]

    def SetGraph(self, parent, FigWrap, ctype = None, overwrite = True):
        if ctype:
            self.chartType = ctype
        self.graph = self.PlotTypeDict[self.chartType](parent, self, overwrite)

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

    def __init__(self, parent, param):
        Tk.Frame.__init__(self)
        self.parent = parent
#        self.sliderText = wx.TextCtrl(parent, -1, style=wx.TE_PROCESS_ENTER)

        self.skipSize = 1
        self.waitTime = .2
        self.playPressed = False

        self.toolbar_text = ['Play','Pause','Stop']
        self.toolbar_length = len(self.toolbar_text)
        self.toolbar_buttons = [None] * self.toolbar_length

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
        self.v = Tk.StringVar()
        # set it to the param value
        self.v.set(str(self.param.value))

        # the entry box
        self.txtEnter = ttk.Entry(self, textvariable=self.v, width=6)

        # A slider that will show the progress in the simulation as well as
        # allow us to select a time

        self.slider = ttk.Scale(self, from_=self.param.minimum, to=self.param.maximum, command = self.ScaleHandler)
        self.slider.set(self.param.value)

        self.txtEnter.pack(side=Tk.LEFT, fill = Tk.BOTH, expand = 0)
        self.slider.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)
        self.SettingsB= ttk.Button(self, text='Settings')
        self.SettingsB.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=0)

        #attach the parameter to the Playbackbar
        self.param.attach(self)



        # while in start, check if stop is clicked, if not, call blink recursively

    def SkipLeft(self):
        self.param.set(self.param.value - self.skipSize)

    def SkipRight(self):
        self.param.set(self.param.value + self.skipSize)


    def PlayHandler(self):
        if not self.playPressed:
            self.playPressed = True
            self.playB.config(text='Pause')
            self.after(int(self.waitTime*1E3), self.blink)
        else:
            self.playPressed = False
            self.playB.config(text='Play')

    def blink(self):
        if self.playPressed:
            print 'looping',self.param.value

            if self.param.value == self.param.maximum: # push pause button
                self.PlayHandler()
            else:
                self.param.set(self.param.value + self.skipSize)

            self.after(int(self.waitTime*1E3), self.blink)


    def TextCallback(self):
        try:
            if int(self.v.get()) != self.param.value:
                self.param.set(int(self.v.get()))
        except ValueError:
            self.v.set(str(self.param.value))

    def ScaleHandler(self, e):
        if self.param.value != int(self.slider.get()):
            self.param.set(int(self.slider.get()))

    def setKnob(self, value):
        self.v.set(str(value))
        self.slider.set(value)
        self.slider.config(to = self.param.maximum)


class MainApp(Tk.Tk):
    """ We simply derive a new class of Frame as the man frame of our app"""
    def __init__(self, name):

        Tk.Tk.__init__(self)
        self.update_idletasks()
        menubar = Tk.Menu(self)
        self.wm_title(name)




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

        self.PathDict = {'Flds': [], 'Prtl': [], 'Param': [], 'Spect': []}
        self.H5KeyDict ={}
        # Set the default color map

        self.cmap = 'inferno'

        # Make the object hold the time info

        # Look for the tristan output files and load the file paths into
        # previous objects
        self.TimeStep = Param(1, minimum=1, maximum=1000)

        self.dirname = os.curdir
        self.findDir()

        self.TimeStep.attach(self)
        self.DrawCanvas()

        self.playbackbar = PlaybackBar(self, self.TimeStep)
        self.playbackbar.pack(side=Tk.TOP, fill=Tk.BOTH, expand=0)
        self.update()
        # now root.geometry() returns valid size/placement
        self.minsize(self.winfo_width(), self.winfo_height())
        self.geometry("1200x700")
        self.bind('<Return>', self.TxtEnter)

    def quit(self, event):
        print("quitting...")
        sys.exit(0)

    def GenH5Dict(self):
        for pkey in self.PathDict.keys():
            with h5py.File(os.path.join(self.dirname,self.PathDict[pkey][0]), 'r') as f:
                for h5key in f.keys():
                    self.H5KeyDict[h5key] = pkey

    def pathOK(self):
        """ Test to see if the current path contains tristan files
        using regular expressions, then generate the lists of files
        to iterate over"""
        dirlist = os.listdir(self.dirname)
        if 'output' in dirlist:
            self.dirname = os.path.join(self.dirname, 'output')
        is_okay = True
        i = 0
        for key in self.PathDict.keys():
            self.PathDict[key]= (filter(self.re_list[i].match, os.listdir(self.dirname)))
            self.PathDict[key].sort()
            is_okay &= len(self.PathDict[key]) > 0
            i += 1
        self.TimeStep.setMax(len(self.PathDict['Flds']))
        if len(self.H5KeyDict) == 0:
            self.GenH5Dict()
        return is_okay


    def OnOpen(self, e = None):
        """open a file"""
        tmpdir = tkFileDialog.askdirectory(title = 'Choose the directory of the output files', **self.dir_opt)
        if tmpdir != '':
            self.dirname = tmpdir
        if not self.pathOK():
#            p = MyDialog(self, 'Directory must contain either the output directory or all of the following: \n flds.tot.*, ptrl.tot.*, params.*, spect.*', title = 'Cannot find output files')
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
                self.FindDir()

    def DrawCanvas(self):
        # figsize (w,h tuple in inches) dpi (dots per inch)
        #f = Figure(figsize=(5,4), dpi=100)
        self.f = Figure(figsize = (2,2), dpi = 100)
        self.a = self.f.add_subplot(111)
        self.a.pcolor(np.random.rand(5,5))
        # a tk.DrawingArea
        self.canvas = FigureCanvasTkAgg(self.f, master=self)
#        self.canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        #self.label = Label(self.top, text = 'Text',bg='orange')
        #self.label.grid()
        # initialize (time index t)

    def setKnob(self, value):
        self.a.pcolor(colors[self.TimeStep.value-1])
        self.canvas.show()
        self.canvas.get_tk_widget().update_idletasks()

    def TxtEnter(self, e):
        self.playbackbar.TextCallback()


#


        #   refresh the graph

if __name__ == "__main__":
    app = MainApp('Iseult')
    app.mainloop()
