#!/usr/bin/env pythonw
import wx # allows us to make the GUI
import re # regular expressions
import os # Used to make the code portable
import h5py # Allows us the read the data files
from thread import start_new_thread
import time
import matplotlib
import new_cmaps

import wx.lib.buttons as buttons
#import icon as icn
import matplotlib.colors as mcolors
from numpy import arange, sin, pi

matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar

class FilesWrapper(object):
    """A simple class wrapper to allow us to make lists of the seperate paths
    of the files and store one of the loaded hdf5 files"""
    def __init__(self, plist=[],f_hdf5=''):
         self.paths = plist
         self.file = f_hdf5

class FigWrapper(object):
    """A simple class wrapper that will eventually hold all of the information
    about each figure in the plot"""
    def __init__(self, ctype='', graph=''):
         self.chartType = ctype
         self.graph = graph


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


class FloatSliderGroup(Knob):
    def __init__(self, parent, label, param):
        self.sliderLabel = wx.StaticText(parent, label=label)
        self.sliderText = wx.TextCtrl(parent, -1, style=wx.TE_PROCESS_ENTER)
        self.slider = wx.Slider(parent, -1)
        self.slider.SetMax(param.maximum*1000)
        self.slider.SetMax(param.minimum*1000)
        self.setKnob(param.value)

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self.sliderLabel, 0, wx.EXPAND | wx.ALIGN_CENTER | wx.ALL, border=2)
        sizer.Add(self.sliderText, 0, wx.EXPAND | wx.ALIGN_CENTER | wx.ALL, border=2)
        sizer.Add(self.slider, 1, wx.EXPAND)
        self.sizer = sizer

        self.slider.Bind(wx.EVT_SLIDER, self.sliderHandler)
        self.sliderText.Bind(wx.EVT_TEXT_ENTER, self.sliderTextHandler)

        self.param = param
        self.param.attach(self)

    def sliderHandler(self, evt):
        value = evt.GetInt() / 1000.
        self.param.set(value)

    def sliderTextHandler(self, evt):
        value = float(self.sliderText.GetValue())
        self.param.set(value)

    def setKnob(self, value):

        self.sliderText.SetValue('%g'%value)
        self.slider.SetValue(value*1000)

def scale_bitmap(bitmap, width, height, flip=False):
    image = wx.ImageFromBitmap(bitmap)
    image = image.Scale(width, height, wx.IMAGE_QUALITY_HIGH)
    if flip:
        image = image.Mirror()
    result = wx.BitmapFromImage(image)
    return result

class PlaybackGroup(Knob):
    def __init__(self, parent, label, param):
        self.Parent = parent
        self.sliderLabel = wx.StaticText(parent, label=label)
        self.sliderText = wx.TextCtrl(parent, -1, style=wx.TE_PROCESS_ENTER)
        self.slider = wx.Slider(parent, -1)
        self.slider.SetMax(param.maximum)
        self.slider.SetMin(param.minimum)
        self.skipSize = 1
        self.waitTime = .2 # 24 fps (may be much slower because of the plotting time)
        self.param = param

        self.setKnob(self.param.value)

        sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.playIcon = scale_bitmap(wx.Bitmap(os.path.join(os.path.dirname(os.path.abspath(__file__)),'icons', 'play.png')), 35,35)
        self.pauseIcon = scale_bitmap(wx.Bitmap(os.path.join(os.path.dirname(os.path.abspath(__file__)),'icons', 'pause.png')), 35,35)
        self.skip_r_Icon = scale_bitmap(wx.Bitmap(os.path.join(os.path.dirname(os.path.abspath(__file__)),'icons', 'skip.png')), 35,35)
        self.skip_l_Icon = scale_bitmap(wx.Bitmap(os.path.join(os.path.dirname(os.path.abspath(__file__)),'icons', 'skip.png')), 35,35, True)
        self.prefIcon = scale_bitmap(wx.Bitmap(os.path.join(os.path.dirname(os.path.abspath(__file__)),'icons', 'params.png')), 35,35)

        self.skiplButton = wx.BitmapButton(parent, -1, self.skip_l_Icon, size = (35,35), style  = wx.NO_BORDER)
        self.playButton = buttons.GenBitmapToggleButton(parent, -1, self.playIcon, size = (35,35),
                            style = wx.NO_BORDER)
        self.playButton.SetBitmapSelected(self.pauseIcon)
        self.skiprButton = wx.BitmapButton(parent, -1, self.skip_r_Icon, size = (35,35), style  = wx.NO_BORDER)
        self.prefButton = wx.BitmapButton(parent, -1, self.prefIcon, size = (35,35), style  = wx.NO_BORDER)

        sizer.Add(self.skiplButton, 0, wx.EXPAND)
        sizer.Add(self.playButton, 0, wx.EXPAND)
        sizer.Add(self.skiprButton, 0, wx.EXPAND)
        sizer.Add(self.sliderLabel, 0, wx.EXPAND | wx.ALIGN_CENTER | wx.ALL, border=2)
        sizer.Add(self.sliderText, 0, wx.EXPAND | wx.ALIGN_CENTER | wx.ALL, border=2)
        sizer.Add(self.slider, 1, wx.EXPAND)
        sizer.Add(self.prefButton, 0, wx.EXPAND)

        self.sizer = sizer

        self.slider.Bind(wx.EVT_SCROLL_THUMBRELEASE, self.sliderHandler)
        self.sliderText.Bind(wx.EVT_TEXT_ENTER, self.sliderTextHandler)

        self.skiplButton.Bind(wx.EVT_BUTTON, self.skipLeft)
        self.skiprButton.Bind(wx.EVT_BUTTON, self.skipRight)
        self.playButton.Bind(wx.EVT_BUTTON, self.onPlay)
        self.prefButton.Bind(wx.EVT_BUTTON, self.openPrefs)

        self.param.attach(self)


    def openPrefs(self, evt):
        win = SettingsFrame(self.Parent, -1, "This is a wx.Frame", size=(350, 200),
                  style = wx.DEFAULT_FRAME_STYLE)
        win.Show(True)

    def skipLeft(self, evt):
        self.param.set(self.param.value - self.skipSize)

    def skipRight(self, evt):
        self.param.set(self.param.value + self.skipSize)

    def onPlay(self, evt):
        start_new_thread(self.playLoop, (self,))

    def playLoop(self, evt):
        start = time.time()
        while self.playButton.GetValue() and self.param.value <self.param.maximum:
            time.sleep(0.02)
            if time.time()-start > self.waitTime:
                self.param.set(self.param.value + self.skipSize)
                start = time.time()
        self.playButton.SetValue(0)

    def sliderHandler(self, evt):
        value = evt.GetInt()
        self.param.set(value)

    def sliderTextHandler(self, evt):
        value = float(self.sliderText.GetValue())
        self.param.set(value)

    def setKnob(self, value):
        self.sliderText.SetValue('%g'%value)
        self.slider.SetValue(value)
        self.slider.SetMax(self.param.maximum)

class CanvasPanel(wx.Window):
    def __init__(self, parent):
        wx.Window.__init__(self, parent)
        self.figure = Figure(figsize=(5, 2), dpi=100)
        self.canvas = FigCanvas(self, -1, self.figure)
        self.Bind(wx.EVT_SIZE, self.sizeHandler)
        self.draw()

    def sizeHandler(self, *args, **kwargs):
        '''Make it so the plot scales with resizing of the window'''
        self.canvas.SetSize(self.GetSize())

    def draw(self):
        self.axes = self.figure.add_subplot(111)
        self.axes.hist2d(self.Parent.prtl.file['xi'][:],self.Parent.prtl.file['ui'][:], bins = [200,200],cmap = new_cmaps.magma, norm = mcolors.PowerNorm(0.4))
        self.canvas.draw()

    def setKnob(self, value):
        self.draw()

class SettingsFrame(wx.Frame):
    def __init__(
            self, parent, ID, title, pos=wx.DefaultPosition,
            size=wx.DefaultSize, style=wx.DEFAULT_FRAME_STYLE
            ):

        wx.Frame.__init__(self, parent, ID, title, pos, size, style)
        panel = wx.Panel(self, -1)
        self.parent = parent
        #Create some sizers
        mainSizer=  wx.BoxSizer(wx.VERTICAL)
        grid =  wx.GridBagSizer(hgap = 5, vgap = 5)
        hSizer=  wx.BoxSizer(wx.HORIZONTAL)
        self.quote = wx.StaticText(self, label='Your Quote: ', pos=(20,30))
        grid.Add(self.quote, pos=(0,0))



        # A multiline TextCtrl - This is here to show how the events
        # work in this program, don't pay too much attention to it
        self.logger = wx.TextCtrl(self, size=(200,300), style= wx.TE_MULTILINE| wx.TE_READONLY)

        # A button
        self.button = wx.Button(self, label='Reload Directory')
        self.Bind(wx.EVT_BUTTON, self.OnReload, self.button)

        # the edit control - one line version.
        self.lblname = wx.StaticText(self, label = 'Your name: ')
        grid.Add(self.lblname, pos=(1,0))
        self.editname = wx.TextCtrl(self, value= 'Enter your name', size = (140,-1))
        grid.Add(self.editname, pos=(1,1))
        self.Bind(wx.EVT_TEXT, self.EvtText, self.editname)
        self.Bind(wx.EVT_CHAR, self.EvtChar, self.editname)
        self.Show()

        # the combobox Control
        self.sampleList = ['Friends', 'advertising', 'web search', 'yellow pages']
        self.lblhear = wx.StaticText(self, label='How did you hear from us?')
        grid.Add(self.lblhear, pos=(3,0))
        self.edithear = wx.ComboBox(self, size = (95,-1), choices = self.sampleList, style = wx.CB_DROPDOWN)
        grid.Add(self.edithear, pos=(3,1))
        self.Bind(wx.EVT_COMBOBOX, self.EvtComboBox, self.edithear)
        self.Bind(wx.EVT_TEXT, self.EvtText, self.edithear)

        # add a spacer to the sizer
        grid.Add((10, 40), pos=(2,0))

        # checkbox
        self.insure = wx.CheckBox(self, label = 'Do you want Insured Shipment?')
        grid.Add(self.insure, pos = (4,0), span=(1,2), flag=wx.BOTTOM, border=5)
        self.Bind(wx.EVT_CHECKBOX, self.EvtCheckBox, self.insure)

        # Radio Boxes
        radioList =['Blue', 'Red', 'Yellow', 'Orange', 'Green', 'Purple', 'Navy Blue', 'Black', 'Gray']
        rb = wx.RadioBox(self, label='What color would you like?', choices=radioList, majorDimension =3, style=wx.RA_SPECIFY_COLS)
        grid.Add(rb, pos=(5,0), span=(1,2))
        self.Bind(wx.EVT_RADIOBOX, self.EvtRadioBox, rb)

        hSizer.Add(grid, 0, wx.ALL, 5)
        hSizer.Add(self.logger)
        mainSizer.Add(hSizer, 0, wx.ALL, 5)
        mainSizer.Add(self.button, 0, wx.LEFT)
        self.SetSizerAndFit(mainSizer)
    # Define functions for the events
    def EvtRadioBox(self, event):
        self.logger.AppendText('EvtRadioBox: %d\n' % event.GetInt())
    def EvtComboBox(self, event):
        self.logger.AppendText('EvtComboBox: %s\n' % event.GetString())
    def OnReload(self, event):
        self.parent.findDir()
    def EvtText(self, event):
        self.logger.AppendText('EvtText: %s\n' % event.GetString())
    def EvtChar(self, event):
        self.logger.AppendText('EvtChar: %d\n' % event.GetKeyCode())
        event.Skip()
    def EvtCheckBox(self, event):
        self.logger.AppendText('EvtCheckBox: %d\n' % event.Checked())

    def OnCloseMe(self, event):
        self.Close(True)

    def OnCloseWindow(self, event):
        self.Destroy()

class MainWindow(wx.Frame):
    """ We simply derive a new class of Frame """
    def __init__(self, parent, title):
        wx.Frame.__init__(self, parent, title =title)

        self.CreateStatusBar() # A statusbar in the bottom of the window
        # intialize the working directory
        # Setting up the menu.
        filemenu = wx.Menu()

        # wx.ID_ABOUT and wx.ID_EXIT are standard IDs provided by wxWidgets.
        menuAbout = filemenu.Append(wx.ID_ABOUT, '&About', ' Information about this program')
        menuExit = filemenu.Append(wx.ID_EXIT,'E&xit', 'Terminate the program')
        menuOpen = filemenu.Append(wx.ID_OPEN, '&Open Directory\tCtrl+o', ' Open the Directory')

        # create the menubar
        menuBar = wx.MenuBar()
        menuBar.Append(filemenu, '&File') # Adding the 'filemenu; to the MenuBar
        self.SetMenuBar(menuBar) # Add the menubar to the frame

        # create a bunch of regular expressions used to search for files
        f_re = re.compile('flds.tot.*')
        p_re = re.compile('prtl.tot.*')
        s_re = re.compile('spect.*')
        par_re = re.compile('param.*')
        self.re_list = [f_re, p_re, s_re, par_re]

        # make a bunch of objects that will store all our file names & values
        self.flds = FilesWrapper()
        self.prtl = FilesWrapper()
        self.param = FilesWrapper()
        self.spect = FilesWrapper()
        self.file_list = [self.flds, self.prtl,  self.spect, self.param]


        # Make the know that will hold the timestep info
        self.timeStep = Param(1, minimum=1, maximum=1000)

        # Look for the tristan output files and load the file paths into
        # previous objects
        self.dirname = os.curdir
        self.findDir()

        # Load the first time step of all the files using h5py
        for elm in self.file_list:
            elm.file = h5py.File(os.path.join(self.dirname,elm.paths[0]), 'r')

        self.timeStep.attach(self)

        # Make some sizers:
        mainsizer = wx.BoxSizer(wx.VERTICAL)
        grid =  wx.GridBagSizer(hgap = 0.5, vgap = 0.5)
        # Make the playback controsl that will control the time slice of the
        # simulation

        self.timeSliderGroup = PlaybackGroup(self, label=' n:', \
            param=self.timeStep)

        # Make the Figures
        self.Fig1 = FigWrapper()
        self.Fig2 = FigWrapper()
        self.Fig3 = FigWrapper()
        self.Fig4 = FigWrapper()
        self.Fig5 = FigWrapper()
        self.Fig6 = FigWrapper()
        self.FigList = [self.Fig1, self.Fig2, self.Fig3, self.Fig4, self.Fig5, self.Fig6]
        col_counter = 0
        for elm in self.FigList:
            elm.graph = CanvasPanel(self)
            self.timeStep.attach(elm.graph)
            grid.Add(elm.graph, pos=(col_counter/2,col_counter%2), flag = wx.EXPAND)
            col_counter += 1

        for i in range(2):
            grid.AddGrowableCol(i)
        for i in range(3):
            grid.AddGrowableRow(i)



        mainsizer.Add(grid,1, wx.EXPAND)
        mainsizer.Add(self.timeSliderGroup.sizer,0, wx.EXPAND | wx.ALIGN_CENTER | wx.ALL, border=5)
        self.SetSizerAndFit(mainsizer)

        # Set events.
        self.Bind(wx.EVT_MENU, self.OnAbout, menuAbout)
        self.Bind(wx.EVT_MENU, self.OnExit, menuExit)
        self.Bind(wx.EVT_MENU, self.OnOpen, menuOpen)

        self.Show(True)

    def setKnob(self, value):
        for elm in self.file_list:
            # close the previous HDF5 file
            elm.file.close()
            # Open the new file
            elm.file = h5py.File(os.path.join(self.dirname,elm.paths[value-1]), 'r')




    # Define the Main Window functions
    def OnAbout(self,e):
        # A message dialog box with an OK buttion. wx.OK is a standardID in wxWidgets.
        dlg = wx.MessageDialog(self, 'A plotting program for Tristan-MP output files', 'About Iseult', wx.OK)
        dlg.ShowModal() # show it
        dlg.Destroy() # destroy it when finished

    def OnExit(self, e):
        self.Close(True)

    def pathOK(self):
        """ Test to see if the current path contains tristan files
        using regular expressions, then generate the lists of files
        to iterate over"""

        is_okay = True
        for i in range(len(self.re_list)):
            self.file_list[i].paths = (filter(self.re_list[i].match, os.listdir(self.dirname)))
            self.file_list[i].paths.sort()
            is_okay &= len(self.file_list[i].paths) > 0
        self.timeStep.setMax(len(self.file_list[i].paths))
        return is_okay

    def OnOpen(self,e):
        """open a file"""
        dlg = wx.DirDialog(self, 'Choose the directory of the output files.', style = wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST)
        if dlg.ShowModal() == wx.ID_OK:
            self.dirname = dlg.GetPath()
        dlg.Destroy()
        if not self.pathOK():
            self.findDir('Directory must contain either the output directory or all of the following: flds.tot.*, ptrl.tot.*, params.*, spect.*')


    def findDir(self, dlgstr = 'Choose the directory of the output files.'):
        """Look for /ouput folder, where the simulation results are
        stored. If output files are already in the path, they are
        automatically loaded"""
        dirlist = os.listdir(self.dirname)
        if 'output' in dirlist:
            self.dirname = os.path.join(self.dirname, 'output')
        if not self.pathOK():
            dlg = wx.DirDialog(self,
                               dlgstr,
                               style = wx.DD_DEFAULT_STYLE
                               | wx.DD_DIR_MUST_EXIST)
            if dlg.ShowModal() == wx.ID_OK:
                self.dirname = dlg.GetPath()
            dlg.Destroy()
            if not self.pathOK() :
                self.findDir('Directory must contain either the output directory or all of the following: flds.tot.*, ptrl.tot.*, params.*, spect.*')

app = wx.App(False)
frame = MainWindow(None, 'Iseult')
app.MainLoop()
