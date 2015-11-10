#!/usr/bin/env pythonw
import wx # allows us to make the GUI
import matplotlib
from validator import MyValidator
import numpy as np
import new_cmaps
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar

class PhasePanel(wx.Window):
    def __init__(self, parent, figwrapper):
        wx.Window.__init__(self, parent)
        self.FigWrap = figwrapper
        self.norm = mcolors.PowerNorm(0.4)
        self.mom_dim = 0
        self.prtl_type = 0 # 1 == electron, 0 == ion
        self.norm_type = "PowerNorm"
        self.pow_num = 0.4

        self.show_cbar = True
        self.weighted = False

        self.figure = Figure(figsize=(3, 1), dpi=100)
        self.canvas = FigCanvas(self, -1, self.figure)
        self.Bind(wx.EVT_SIZE, self.sizeHandler)


        self.draw()
        self.Bind(wx.EVT_ENTER_WINDOW, self.onEnter)
        self.Bind(wx.EVT_LEAVE_WINDOW, self.onLeave)
        self.Bind(wx.EVT_LEFT_UP, self.openGraphPrefs)


    def sizeHandler(self, *args, **kwargs):
        '''Make it so the plot scales with resizing of the window'''
        self.canvas.SetSize(self.GetSize())

    def draw(self, kwargs=''):
        # Choose the normalization
        if self.norm_type =="Linear":
            self.norm = mcolors.Normalize()
        elif self.norm_type =="LogNorm":
            self.norm = mcolors.LogNorm()
        else:
            self.norm = mcolors.PowerNorm(self.pow_num)

        # Choose the particle type and px, py, or pz
        if self.prtl_type == 0:
            self.x_values = self.FigWrap.LoadKey('xi')
            if self.mom_dim == 0:
                self.y_values = self.FigWrap.LoadKey('ui')
            if self.mom_dim == 1:
                self.y_values = self.FigWrap.LoadKey('vi')
            if self.mom_dim == 2:
                self.y_values = self.FigWrap.LoadKey('wi')

        if self.prtl_type == 1:
            self.x_values = self.FigWrap.LoadKey('xe')
            if self.mom_dim == 0:
                self.y_values = self.FigWrap.LoadKey('ue')
            if self.mom_dim == 1:
                self.y_values = self.FigWrap.LoadKey('ve')
            if self.mom_dim == 2:
                self.y_values = self.FigWrap.LoadKey('we')

        self.figure.clf()
        gs = GridSpec(100,100,bottom=0.1,left=0.1,right=0.95, top = 0.95)

        if self.show_cbar:
            self.axes = self.figure.add_subplot(gs[20:,:])
            self.axC = self.figure.add_subplot(gs[:5,:])
            self.cax = self.axes.hist2d(self.x_values,self.y_values, bins = [200,200],cmap = new_cmaps.cmaps[self.Parent.cmap], norm = self.norm)
            self.figure.colorbar(self.cax[3], ax = self.axes, cax = self.axC, orientation = 'horizontal')
        else:
            self.axes = self.figure.add_subplot(gs[5:,:])
            self.cax = self.axes.hist2d(self.x_values,self.y_values, bins = [200,200],cmap = new_cmaps.cmaps[self.Parent.cmap], norm = self.norm)

        self.canvas.draw()

    def onEnter(self, evt):
        if not self.HasCapture():
            self.CaptureMouse()
    def onLeave(self, evt):
        if self.HasCapture():
            self.ReleaseMouse()
    def openGraphPrefs(self, evt):
        win = PhaseSettings(self, -1, "Chart Settings",
                          style = wx.DEFAULT_FRAME_STYLE)
        win.Show(True)

class PhaseSettings(wx.Frame):
    def __init__(
            self, parent, ID, title, pos=wx.DefaultPosition,
            size=wx.DefaultSize, style=wx.DEFAULT_FRAME_STYLE
            ):

        wx.Frame.__init__(self, parent, ID, title, pos, size, style)
        panel = wx.Panel(self, -1)
        self.parent = parent
        #Create some sizers
        self.mainsizer = wx.BoxSizer(wx.VERTICAL)
        grid =  wx.GridBagSizer(hgap = 10, vgap = 10)
        # the Radiobox Control
        self.prtlList = ['ion', 'electron']
        self.rbPrtl = wx.RadioBox(
                self, -1,'Particle', wx.DefaultPosition, wx.DefaultSize,
                self.prtlList,  1, wx.RA_SPECIFY_COLS
                )
        self.rbPrtl.SetSelection(self.parent.prtl_type)
        self.Bind(wx.EVT_RADIOBOX, self.EvtRadioPrtl, self.rbPrtl)

        grid.Add(self.rbPrtl, pos = (1,0))
        self.dimList = ['x-px', 'x-py', 'x-pz']
        self.rbDim = wx.RadioBox(
                self, -1,'Dimension', wx.DefaultPosition, wx.DefaultSize,
                self.dimList,  1, wx.RA_SPECIFY_COLS
                )
        self.rbDim.SetSelection(self.parent.mom_dim)
        grid.Add(self.rbDim, pos = (1,1))

        self.Bind(wx.EVT_RADIOBOX, self.EvtRadioDim, self.rbDim)

        grid1 = wx.GridBagSizer()

        # Group of controls for the cmap normalization:
        self.rbGroupLabel = wx.StaticText(self, label = 'Choose Cmap Norm:') # Title
        self.rbLinear = wx.RadioButton(self, -1, "Linear", style = wx.RB_GROUP)
        self.rbLog = wx.RadioButton(self, -1, "LogNorm" )
        self.rbPow = wx.RadioButton(self, -1, "PowerNorm" )
        self.rbPowNum = wx.TextCtrl(self , -1, str(self.parent.pow_num),
        validator = MyValidator('DIGIT_ONLY') )

        # Set the rb that should be selected from the parent
        if self.parent.norm_type == "Linear":
            self.rbLinear.SetValue(True)
        elif self.parent.norm_type == "LogNorm":
            self.rbLog.SetValue(True)

        else:
            self.rbPow.SetValue(True)

        # Add to Grid
        grid1.Add(self.rbGroupLabel, pos = (0,0), span = (1,2))
        grid1.Add( self.rbLinear, pos = (1,0), flag = wx.ALIGN_LEFT|wx.LEFT|wx.RIGHT|wx.TOP, border = 5)
        grid1.Add( self.rbLog, pos = (2,0), flag = wx.ALIGN_LEFT|wx.LEFT|wx.RIGHT|wx.TOP, border = 5)
        grid1.Add( self.rbPow, pos = (3,0), flag = wx.ALIGN_LEFT|wx.LEFT|wx.RIGHT|wx.TOP, border = 5)
        grid1.Add( self.rbPowNum, pos = (3,1), flag = wx.ALIGN_LEFT|wx.LEFT|wx.RIGHT|wx.TOP, border = 5)

        self.Bind(wx.EVT_RADIOBUTTON, self.EvtRadioNorm, self.rbLinear)
        self.Bind(wx.EVT_RADIOBUTTON, self.EvtRadioNorm, self.rbLog)
        self.Bind(wx.EVT_RADIOBUTTON, self.EvtRadioNorm, self.rbPow)
        self.Bind(wx.EVT_TEXT, self.EvtNormNumSet, self.rbPowNum)

        grid.Add(grid1, pos=(2,0), span=(1,2))

        self.cbColorbar = wx.CheckBox(self, -1, "Show Cbar")
        self.cbColorbar.SetValue(self.parent.show_cbar)
        self.Bind(wx.EVT_CHECKBOX, self.EvtCheckCbar, self.cbColorbar)
        grid.Add(self.cbColorbar, pos=(3,0))

        #Adding WEIGHT *TO DO*
#        self.cbWeight = wx.CheckBox(self, -1, "Weight")
#        self.cbWeight.SetValue(self.parent.weighted)
#        self.Bind(wx.EVT_CHECKBOX, self.EvtCheckWeight, self.cbWeight)
#        grid.Add(self.cbWeight, pos = (3,1))


        self.mainsizer.Add(grid,0, border=15)
#        self.mainsizer.Add(grid1, 0, border = 15)

        self.SetSizerAndFit(self.mainsizer)
    # Define functions for the events

    def EvtCheckCbar(self, evt):
        self.parent.show_cbar = evt.IsChecked()
        self.parent.draw()

    def EvtCheckWeight(self, evt):
        self.parent.weighted = evt.IsChecked()
        self.parent.draw()

    def EvtRadioDim(self, evt):
        self.parent.mom_dim = evt.GetInt()
        self.parent.draw()

    def EvtRadioNorm(self, evt):
        self.parent.norm_type = evt.GetEventObject().GetLabel()
        self.parent.draw()

    def EvtNormNumSet(self, evt):
        tmp_num = evt.GetString()
        if not tmp_num:
            tmp_num = 1E-2
        elif float(tmp_num)<=1E-2:
            tmp_num = 1E-2
        else:
            tmp_num = float(tmp_num)
        self.parent.pow_num = tmp_num
        self.parent.draw()

    def EvtRadioPrtl(self, evt):
        self.parent.prtl_type = evt.GetInt()
        self.parent.draw()

    def OnCloseMe(self, event):
        self.Close(True)

    def OnCloseWindow(self, event):
        self.Destroy()
