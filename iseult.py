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

class TimeStepper:
    """A Class that will hold the time step info"""
    def __init__(self, initial = 1, minimum = 1, maximum =1):
        self.min_time = min(minimum, maximum)
        self.max_time = max(minimum, maximum)
        self.cur_time = 1
        self.SetTime(initial)

    def GetTime(self):
        return self.cur_time

    def GetMin(self):
        return self.min_time

    def GetMax(self):
        return self.max_time

    def SetMax(self, val):
        if self.max_time < self.min_time:
            self.max_time = self.min_time

        else:
            self.max_time = val

        self.SetTime(self.cur_time)

    def SetMin(self, val):
        if self.min_time > self.max_time:
            self.min_time = self.max_time

        else:
            self.min_time = val

        self.SetTime(self.cur_time)

    def SetTime(self, val):
        if val < self.min_time:
            self.cur_time = self.min_time
        elif val > self.max_time:
            self.cur_time = self.max_time
        else:
            self.cur_time = val

    def Step(self, val):
        self.SetTime(self.cur_time+val)

class MainApp(Tk.Tk):
    """ We simply derive a new class of Frame as the man frame of our app"""
    def __init__(self, name):

        Tk.Tk.__init__(self)
        self.update_idletasks()
        menubar = Tk.Menu(self)
        self.wm_title(name)

        self.skipSize = 1
        self.waitTime = .2


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
        self.TimeStep = TimeStepper(initial = 1, minimum=1, maximum=1000)

        # Look for the tristan output files and load the file paths into
        # previous objects
        self.dirname = os.curdir
        self.findDir()

        self.DrawCanvas()
        self.makeToolbar()

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
        self.TimeStep.SetMax(len(self.PathDict['Flds']))
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
        self.f = Figure()
        self.a = self.f.add_subplot(111)
        self.a.pcolor(np.random.rand(5,5))
        # a tk.DrawingArea
        self.canvas = FigureCanvasTkAgg(self.f, master=self)
#        self.canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        self.canvas.get_tk_widget().grid(row=3,column=0,columnspan=3,sticky= Tk.N+Tk.S+Tk.E+Tk.W)
        self.canvas.show()
        self.bClose = ttk.Button(self, text='Close',command=self.destroy)
        self.bClose.grid()
        #self.label = Label(self.top, text = 'Text',bg='orange')
        #self.label.grid()
        # initialize (time index t)

    def makeToolbar(self):
        self.toolbar_text = ['Play','Pause','Stop']
        self.toolbar_length = len(self.toolbar_text)
        self.toolbar_buttons = [None] * self.toolbar_length

        for toolbar_index in range(self.toolbar_length):
            text = self.toolbar_text[toolbar_index]
            button_id = ttk.Button(self ,text=text)
            button_id.grid(row=4, column=toolbar_index)
            self.toolbar_buttons[toolbar_index] = button_id

            def toolbar_button_handler(event, self=self, button=toolbar_index):
                return self.service_toolbar(button)

            button_id.bind("<Button-1>", toolbar_button_handler)

        # call blink() if start and set stop when stop
    def service_toolbar(self, toolbar_index):
            if toolbar_index == 0:
                self.stop = False
                print self.stop
                self.blink()
            elif toolbar_index == 1:
                self.stop = True
                print self.stop
            elif toolbar_index == 2:
                self.stop = True
                print self.stop
                self.TimeStep.SetTime(self.TimeStep.GetMin())

        # while in start, check if stop is clicked, if not, call blink recursivly
    def blink(self):
        if not self.stop:
            print 'looping',self.stop
            self.a.pcolor(colors[self.TimeStep.GetTime()])

            self.TimeStep.Step(self.skipSize)
            if self.TimeStep.GetTime() == self.TimeStep.GetMax(): # push stop button
                self.service_toolbar(2)

            self.canvas.show()
            self.canvas.get_tk_widget().update_idletasks()
            #self.label.update_idletasks()
            self.after(int(self.waitTime*1E3), self.blink)

if __name__ == "__main__":
    app = MainApp('Iseult')
    app.mainloop()
