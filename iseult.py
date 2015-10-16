#!/usr/bin/env pythonw
import wx # allows us to make the GUI
import re # regular expressions
import os # Used to make the code portable
import h5py # Allows us the read the data files
import matplotlib #plotter

matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar

class MainWindow(wx.Frame):
    """ We simply derive a new class of Frame """
    def __init__(self, parent, title):
        wx.Frame.__init__(self, parent, title =title, size = (200,100))
        self.control = wx.TextCtrl(self, style = wx.TE_MULTILINE)
        self.CreateStatusBar() # A statusbar in the bottome of the window
        # intialize the working directory
        # Setting up the menu.
        filemenu = wx.Menu()
        
        # wx.ID_ABOUT and wx.ID_EXIT are standard IDs provided by wxWidgets.
        menuAbout = filemenu.Append(wx.ID_ABOUT, '&About', ' Information about this program')
        menuExit = filemenu.Append(wx.ID_EXIT,'E&xit', 'Terminate the program')
        menuOpen = filemenu.Append(wx.ID_OPEN, '&Open', ' Open the file')

        # create the menubar
        menuBar = wx.MenuBar()
        menuBar.Append(filemenu, '&File') # Adding the 'filemenu; to the MenuBar
        self.SetMenuBar(menuBar) # Add the menubar to the frame

        # Look for the tristan output files
        self.dirname = os.curdir
        self.findDir()
        print self.dirname
        # Set events.
        self.Bind(wx.EVT_MENU, self.OnAbout, menuAbout)
        self.Bind(wx.EVT_MENU, self.OnExit, menuExit)
        self.Bind(wx.EVT_MENU, self.OnOpen, menuOpen)

        self.Show(True)
    def OnAbout(self,e):
        # A message dialo box with an OK buttion. wx.OK is a standardID in wxWidgets.
        dlg = wx.MessageDialog(self, 'A small text editor', 'About Simple Editor', wx.OK)
        dlg.ShowModal() # show it
        dlg.Destroy() # destroy it when finished

    def OnExit(self, e):
        self.Close(True)
        print self.dirname
        
    def pathOK(self):
        """ Test to see if the current path contains tristan files """
        f_re = re.compile('flds.tot.*')
        p_re = re.compile('prtl.tot.*')
        s_re = re.compile('spect.*')
        par_re = re.compile('param.*')
        re_list = [f_re, p_re, s_re,par_re]
        is_okay = True
        for r in re_list:
            is_okay &= len(filter(r.match, os.listdir(self.dirname))) > 0
        return is_okay
    
    def OnOpen(self,e):
        """open a file"""
        self.dirname = ''

        dlg = wx.DirDialog(self, 'Choose the directory of the output files.', style = wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST)
        if dlg.ShowModal() == wx.ID_OK:
            self.dirname = dlg.GetPath()
        dlg.Destroy()
        if not self.pathOK():
            self.findDir('Directory must contain either the output all of the following: flds.tot.*, ptrl.tot.*, params.*, spect.*') 
            
        
    def findDir(self, dlgstr = 'Choose the directory of the output files.'):
        """Look for /ouput folder, where the simulation results are 
        stored. If output files are already in the path, they are 
        automatically loaded"""
        dirlist = os.listdir(self.dirname)
        if 'output' in dirlist:
            self.dirname = os.path.join(self.dirname, 'output')
        elif not self.pathOK():
            dlg = wx.DirDialog(self,
                               dlgstr,
                               style = wx.DD_DEFAULT_STYLE
                               | wx.DD_DIR_MUST_EXIST)
            if dlg.ShowModal() == wx.ID_OK:
                self.dirname = dlg.GetPath()
            dlg.Destroy()
            if not self.pathOK() :
                self.findDir('Directory must contain all of the following: flds.tot.*, ptrl.tot.*,params.*, spect.*')
    
app = wx.App(False)
frame = MainWindow(None, 'Iseult')
app.MainLoop()

