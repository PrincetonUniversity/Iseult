#!/usr/bin/env pythonw
import Tkinter as Tk
import ttk as ttk
import matplotlib
import numpy as np
import numpy.ma as ma
import new_cmaps
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as PathEffects
from scipy.special import kn # Modified Bessel function

class SpectralPanel:
    # A diction of all of the parameters for this plot with the default parameters

    plot_param_dict = {'spectral_type': 0, #0 dn/dp, 1 = dn/dE
                        'twoD': 0,
                       'show_ions': 1,
                       'show_electrons': 1,
                       'rest_frame': False,
                       'set_ylim': True,
                       'set_xlim': True,
                       'x_min': None,
                       'x_max': None,
                       'spatial_x': False,
                       'spatial_y': False,
                       'show_legend': True,
                       'y_min': -6,
                       'y_max': 0}

    def __init__(self, parent, figwrapper):
        self.settings_window = None
        self.FigWrap = figwrapper
        self.parent = parent
        self.ChartTypes = self.FigWrap.PlotTypeDict.keys()
        self.chartType = self.FigWrap.chartType
        self.figure = self.FigWrap.figure
        if self.GetPlotParam('spectral_type') == 0:
            self.SetPlotParam('x_min', 0.05, update_plot = False)
            self.SetPlotParam('x_max', 200, update_plot = False)
        if self.GetPlotParam('spectral_type') == 1:
            self.SetPlotParam('x_min', 0.0005, update_plot = False)
            self.SetPlotParam('x_max', 100, update_plot = False)


    def ChangePlotType(self, str_arg):
        self.FigWrap.ChangeGraph(str_arg)

    def set_plot_keys(self):
        '''A helper function that will insure that each hdf5 file will only be
        opened once per time step'''
        # First make sure that omega_plasma & istep is loaded so we can fix the
        # integrate over the specified regions

        self.arrs_needed = ['c_omp', 'istep', 'gamma', 'xsl', 'mi', 'me']
        if self.GetPlotParam('rest_frame'):
            # Set the loading of the rest frame spectra
            self.arrs_needed.append('specerest')
            self.arrs_needed.append('specprest')
        else:
            # Load the normal spectra
            self.arrs_needed.append('spece')
            self.arrs_needed.append('specp')


        return self.arrs_needed

    def draw(self):
        # get c_omp and istep to convert cells to physical units
        self.c_omp = self.FigWrap.LoadKey('c_omp')[0]
        self.istep = self.FigWrap.LoadKey('istep')[0]
        self.xsl = self.FigWrap.LoadKey('xsl')/self.c_omp
        self.gamma = self.FigWrap.LoadKey('gamma')

        if self.GetPlotParam('rest_frame'):
            self.spece = np.copy(self.FigWrap.LoadKey('specerest'))
            self.specp = np.copy(self.FigWrap.LoadKey('specprest'))
        else:
            self.spece = np.copy(self.FigWrap.LoadKey('spece'))
            self.specp = np.copy(self.FigWrap.LoadKey('specp'))


        # In output.F90, spece (specp) is defined by the number of electons (ions)
        # divided by gamma in each logrithmic energy bin. So we multiply by gamma.

        for i in range(len(self.xsl)):
            self.spece[:,i] *= self.gamma
            self.specp[:,i] *= self.gamma
        ###############################
        ###### energy spectra, f=(dN/dE)/N
        ###############################
        self.dgamma = np.empty(len(self.gamma))
        delta=np.log10(self.gamma[-1]/self.gamma[0])/len(self.gamma)
        for i in range(len(self.dgamma)):
            self.dgamma[i]=self.gamma[i]*(10**delta-1.)

        # Select the x-range from which to take the spectra
        e_left_loc = self.parent.e_L.get()
        e_right_loc = self.parent.e_R.get()

        if self.parent.e_relative:
            e_left_loc = self.parent.shock_loc/self.c_omp*self.istep+self.parent.e_L.get()
            e_right_loc = self.parent.shock_loc/self.c_omp*self.istep+self.parent.e_R.get()

        eL = self.xsl.searchsorted(e_left_loc)
        eR = self.xsl.searchsorted(e_right_loc, side='right')

        i_left_loc = self.parent.i_L.get()
        i_right_loc = self.parent.i_R.get()

        if self.parent.e_relative:
            i_left_loc = self.parent.shock_loc+self.parent.i_L.get()
            i_right_loc = self.parent.shock_loc+self.parent.i_R.get()

        iL = self.xsl.searchsorted(i_left_loc)
        iR = self.xsl.searchsorted(i_right_loc, side='right')

        if iL >= iR:
            iL = iR
            iR += 1
        if eL >= eR:
            eL = eR
            eR += 1
        # total particles in each linear x bin
        norme = np.copy(self.xsl)
        normp = np.copy(self.xsl)
        for i in range(len(norme)):
            norme[i]=sum(self.spece[:,i])
            normp[i]=sum(self.specp[:,i])

        # energy distribution, f(E)=(dn/dE)/N
        self.fe=np.empty(len(self.gamma))
        self.fp=np.empty(len(self.gamma))


        for k in range(len(self.fe)):
            self.fe[k]=sum(self.spece[k][eL:eR])/(sum(norme[eL:eR])*self.dgamma[k])
            self.fp[k]=sum(self.specp[k][iL:iR])/(sum(normp[iL:iR])*self.dgamma[k])


        #  NOTE: gamma ---> gamma-1 ***
        self.edist=self.gamma*self.fe
        self.pdist=self.gamma*self.fp

        self.momentum=np.sqrt((self.gamma+1)**2-1.)
        self.femom=self.fe/(4*np.pi*self.momentum)/(self.gamma+1)
        self.momedist=self.femom*self.momentum**4
        self.fpmom=self.fp/(4*np.pi*self.momentum)/(self.gamma+1)
        self.mompdist=self.fpmom*self.momentum**4

        # Set the tick color
        tick_color = 'black'

        # Create a gridspec to handle spacing better
        self.gs = gridspec.GridSpecFromSubplotSpec(100,100, subplot_spec = self.parent.gs0[self.FigWrap.pos])#, bottom=0.2,left=0.1,right=0.95, top = 0.95)

        # load the density values
        self.axes = self.figure.add_subplot(self.gs[18:92,:])

        if self.GetPlotParam('spectral_type') == 0: #Show the momentum dist
            if self.GetPlotParam('show_ions'):
                self.axes.plot(self.momentum, self.mompdist, color = self.parent.ion_color)
                # See if we want to plot the electron temperature
                if self.parent.set_Tp:
                    if self.parent.delgam_p >= 0.013:
                        aconst = 1/(self.parent.delgam_p*np.exp(1.0/self.parent.delgam_p)*kn(2, 1.0/self.parent.delgam_p))
                    else:
                        aconst = np.sqrt(2/np.pi)/self.parent.delgam_p**1.5
                        aconst -= 15.0/(4.*np.sqrt(self.parent.delgam_p)*np.sqrt(2*np.pi))
                        aconst += (345*np.sqrt(self.parent.delgam_p))/(64.*np.sqrt(2*np.pi))
                        aconst -= (3285*self.parent.delgam_p**1.5)/(512.*np.sqrt(2*np.pi))
                        aconst += (95355*self.parent.delgam_p**2.5)/(16384.*np.sqrt(2*np.pi))
                        aconst -= (232065*self.parent.delgam_p**3.5)/(131072.*np.sqrt(2*np.pi))

                    fpmommax = self.momentum**4*aconst*(self.gamma+1.0)*np.sqrt((self.gamma+1.0)**2-1)
                    fpmommax *= np.exp(-self.gamma/self.parent.delgam_p)/(4*np.pi*self.momentum)/(self.gamma+1.0)
                    self.T_p, = self.axes.plot(self.momentum, fpmommax, color = self.parent.ion_fit_color, linestyle = '--', linewidth = 1.0)


            if self.GetPlotParam('show_electrons'):

                self.axes.plot(self.momentum, self.momedist, color = self.parent.electron_color)
                # See if we want to plot the electron temperature
                if self.parent.set_Te:
                    self.delgame0=self.parent.delgam_e*self.FigWrap.LoadKey('mi')[0]/self.FigWrap.LoadKey('me')[0]
                    if self.delgame0 >= 0.013:
                        aconst = 1/(self.delgame0*np.exp(1.0/self.delgame0)*kn(2, 1.0/self.delgame0))
                    else:
                        aconst = np.sqrt(2/np.pi)/self.delgame0**1.5
                        aconst -= 15.0/(4.*np.sqrt(self.delgame0)*np.sqrt(2*np.pi))
                        aconst += (345*np.sqrt(self.delgame0))/(64.*np.sqrt(2*np.pi))
                        aconst -= (3285*self.delgame0**1.5)/(512.*np.sqrt(2*np.pi))
                        aconst += (95355*self.delgame0**2.5)/(16384.*np.sqrt(2*np.pi))
                        aconst -= (232065*self.delgame0**3.5)/(131072.*np.sqrt(2*np.pi))

                    femommax = self.momentum**4*aconst*(self.gamma+1.0)*np.sqrt((self.gamma+1.0)**2-1)
                    femommax *= np.exp(-self.gamma/self.delgame0)/(4*np.pi*self.momentum)/(self.gamma+1.0)
                    self.T_e, = self.axes.plot(self.momentum, femommax, color = self.parent.electron_fit_color, linestyle = '--', linewidth = 1.0)



            self.axes.set_xscale("log")
            self.axes.set_yscale("log")
            self.axes.set_axis_bgcolor('lightgrey')
            if self.GetPlotParam('set_xlim'):
                self.axes.set_xlim(self.GetPlotParam('x_min'), self.GetPlotParam('x_max'))
            if self.GetPlotParam('set_ylim'):
                self.axes.set_ylim(10**self.GetPlotParam('y_min'), 10**self.GetPlotParam('y_max'))

            if self.GetPlotParam('show_legend'):
                legend_handles = []
                legend_labels = []
                if self.GetPlotParam('show_electrons') and self.parent.set_Te:
                    legend_handles.append(self.T_e)
                    tmpstr = '%.3f' % self.delgame0
                    legend_labels.append(r'$T_e\ = $' +  ' ' + tmpstr + ' ' + r'$m_e c^2$')
                if self.GetPlotParam('show_ions') and self.parent.set_Tp:
                    legend_handles.append(self.T_p)
                    tmpcon =self.parent.delgam_p*self.FigWrap.LoadKey('mi')[0]/self.FigWrap.LoadKey('me')[0]
                    tmpstr = '%.3f' % tmpcon
                    legend_labels.append(r'$T_p\ = $' +  ' ' + tmpstr + ' ' + r'$m_e c^2$')
                if len(legend_handles)> 0:
                    leg = self.axes.legend(legend_handles, legend_labels, framealpha = .4, fontsize = 11, loc = 'best')
                    leg.get_frame().set_linewidth(0.0)

            self.axes.tick_params(labelsize = self.parent.num_font_size, color=tick_color)

            self.axes.set_xlabel(r'$p\ [mc]$', labelpad = self.parent.xlabel_pad, color = 'black')
            self.axes.set_ylabel(r'$p^4f(p)$', labelpad = self.parent.ylabel_pad, color = 'black')

        if self.GetPlotParam('spectral_type') == 1: #Show the energy dist
            if self.GetPlotParam('show_electrons'):
                self.axes.plot(self.gamma, self.edist, color = self.parent.electron_color)
                if self.parent.set_Te:
                    self.delgame0=self.parent.delgam_e*self.FigWrap.LoadKey('mi')[0]/self.FigWrap.LoadKey('me')[0]
                    if self.delgame0 >= 0.013:
                        aconst = 1/(self.delgame0*np.exp(1.0/self.delgame0)*kn(2, 1.0/self.delgame0))
                    else:
                        aconst = np.sqrt(2/np.pi)/self.delgame0**1.5
                        aconst -= 15.0/(4.*np.sqrt(self.delgame0)*np.sqrt(2*np.pi))
                        aconst += (345*np.sqrt(self.delgame0))/(64.*np.sqrt(2*np.pi))
                        aconst -= (3285*self.delgame0**1.5)/(512.*np.sqrt(2*np.pi))
                        aconst += (95355*self.delgame0**2.5)/(16384.*np.sqrt(2*np.pi))
                        aconst -= (232065*self.delgame0**3.5)/(131072.*np.sqrt(2*np.pi))

                    femax = aconst*self.gamma*(self.gamma+1.0)*np.sqrt((self.gamma+1.0)**2-1)
                    femax *= np.exp(-self.gamma/self.delgame0)
                    self.T_e, = self.axes.plot(self.gamma, femax, color = self.parent.electron_fit_color, linestyle = '--', linewidth = 1.0)

            if self.GetPlotParam('show_ions'):
                self.axes.plot(self.gamma, self.pdist, color = self.parent.ion_color)
                if self.parent.set_Tp:
                    if self.parent.delgam_p >= 0.013:
                        aconst = 1/(self.parent.delgam_p*np.exp(1.0/self.parent.delgam_p)*kn(2, 1.0/self.parent.delgam_p))
                    else:
                        aconst = np.sqrt(2/np.pi)/self.parent.delgam_p**1.5
                        aconst -= 15.0/(4.*np.sqrt(self.parent.delgam_p)*np.sqrt(2*np.pi))
                        aconst += (345*np.sqrt(self.parent.delgam_p))/(64.*np.sqrt(2*np.pi))
                        aconst -= (3285*self.parent.delgam_p**1.5)/(512.*np.sqrt(2*np.pi))
                        aconst += (95355*self.parent.delgam_p**2.5)/(16384.*np.sqrt(2*np.pi))
                        aconst -= (232065*self.parent.delgam_p**3.5)/(131072.*np.sqrt(2*np.pi))

                    fpmax = aconst*self.gamma*(self.gamma+1.0)*np.sqrt((self.gamma+1.0)**2-1)
                    fpmax *= np.exp(-self.gamma/self.parent.delgam_p)
                    self.T_p, = self.axes.plot(self.gamma, fpmax, color = self.parent.ion_fit_color, linestyle = '--', linewidth = 1.0)


            self.axes.set_xscale("log")
            self.axes.set_yscale("log")
            self.axes.set_axis_bgcolor('lightgrey')
            if self.GetPlotParam('set_xlim'):
                self.axes.set_xlim(self.GetPlotParam('x_min'), self.GetPlotParam('x_max'))
            if self.GetPlotParam('set_ylim'):
                self.axes.set_ylim(10**self.GetPlotParam('y_min'), 10**self.GetPlotParam('y_max'))

            if self.GetPlotParam('show_legend'):
                legend_handles = []
                legend_labels = []
                if self.GetPlotParam('show_electrons') and self.parent.set_Te:
                    legend_handles.append(self.T_e)
                    tmpstr = '%.3f' % self.delgame0
                    legend_labels.append(r'$T_e\ = $' +  ' ' + tmpstr + ' ' + r'$m_e c^2$')
                if self.GetPlotParam('show_ions') and self.parent.set_Tp:
                    legend_handles.append(self.T_p)
                    tmpcon =self.parent.delgam_p*self.FigWrap.LoadKey('mi')[0]/self.FigWrap.LoadKey('me')[0]
                    tmpstr = '%.3f' % tmpcon
                    legend_labels.append(r'$T_p\ = $' +  ' ' + tmpstr + ' ' + r'$m_e c^2$')
                if len(legend_handles)> 0:
                    leg = self.axes.legend(legend_handles, legend_labels, framealpha = .4, fontsize = 11, loc = 'best')
                    leg.get_frame().set_linewidth(0.0)

            self.axes.tick_params(labelsize = self.parent.num_font_size, color=tick_color)

            self.axes.set_xlabel(r'$E\ [mc^2]$', labelpad = -2, color = 'black')
            self.axes.set_ylabel(r'$E(dn/dE)/n$', labelpad = 0, color = 'black')

    def GetPlotParam(self, keyname):
        return self.FigWrap.GetPlotParam(keyname)

    def SetPlotParam(self, keyname, value, update_plot = True):
        self.FigWrap.SetPlotParam(keyname, value, update_plot = update_plot)

    def OpenSettings(self):
        if self.settings_window is None:
            self.settings_window = SpectraSettings(self)
        else:
            self.settings_window.destroy()
            self.settings_window = SpectraSettings(self)


class SpectraSettings(Tk.Toplevel):
    def __init__(self, parent):
        self.parent = parent
        Tk.Toplevel.__init__(self)

        self.wm_title('Spectrum (%d,%d) Settings' % self.parent.FigWrap.pos)
        self.parent = parent
        frm = ttk.Frame(self)
        frm.pack(fill=Tk.BOTH, expand=True)
        self.protocol('WM_DELETE_WINDOW', self.OnClosing)
        #Create some sizers

        self.bind('<Return>', self.TxtEnter)

        # Create the OptionMenu to chooses the Chart Type:
        self.ctypevar = Tk.StringVar(self)
        self.ctypevar.set(self.parent.chartType) # default value
        self.ctypevar.trace('w', self.ctypeChanged)

        ttk.Label(frm, text="Choose Chart Type:").grid(row=0, column = 0)
        cmapChooser = apply(ttk.OptionMenu, (frm, self.ctypevar, self.parent.chartType) + tuple(self.parent.ChartTypes))
        cmapChooser.grid(row =0, column = 1, sticky = Tk.W + Tk.E)


        # the Radiobox Control to choose the Field Type
        self.SpectList = ['Momentum', 'Energy']
        self.SpectTypeVar  = Tk.IntVar()
        self.SpectTypeVar.set(self.parent.GetPlotParam('spectral_type'))

        ttk.Label(frm, text='Choose Spectrum Type:').grid(row = 2, sticky = Tk.W)

        for i in range(len(self.SpectList)):
            ttk.Radiobutton(frm,
                text=self.SpectList[i],
                variable=self.SpectTypeVar,
                command = self.RadioSpect,
                value=i).grid(row = 3+i, sticky =Tk.W)


        # show ions
        self.IonVar = Tk.IntVar()
        self.IonVar.set(self.parent.GetPlotParam('show_ions'))
        cb = ttk.Checkbutton(frm, text = "Show ions",
                        variable = self.IonVar,
                        command = lambda:
                        self.parent.SetPlotParam('show_ions', self.IonVar.get()))
        cb.grid(row = 6, column = 0, sticky = Tk.W)

        # show electrons
        self.eVar = Tk.IntVar()
        self.eVar.set(self.parent.GetPlotParam('show_electrons'))
        cb = ttk.Checkbutton(frm, text = "Show electrons",
                        variable = self.eVar,
                        command = lambda:
                        self.parent.SetPlotParam('show_electrons', self.eVar.get()))
        cb.grid(row = 6, column = 1, sticky = Tk.W)
        # show in rest frame
        self.RestVar = Tk.IntVar()
        self.RestVar.set(self.parent.GetPlotParam('rest_frame'))
        cb = ttk.Checkbutton(frm, text = "Show in rest frame",
                        variable = self.RestVar,
                        command = lambda:
                        self.parent.SetPlotParam('rest_frame', self.RestVar.get()))
        cb.grid(row = 7, column = 0, sticky = Tk.W)

        self.xLimVar = Tk.IntVar()
        self.xLimVar.set(self.parent.GetPlotParam('set_xlim'))
        self.xLimVar.trace('w', self.xLimChanged)



        self.xmin = Tk.StringVar()
        self.xmin.set(str(self.parent.GetPlotParam('x_min')))
        self.xmax = Tk.StringVar()
        self.xmax.set(str(self.parent.GetPlotParam('x_max')))


        cb = ttk.Checkbutton(frm, text ='Set xlim',
                        variable = self.xLimVar)
        cb.grid(row = 4, column =3,sticky = Tk.W)
        self.eLEnter = ttk.Entry(frm, textvariable=self.xmin, width=7)
        self.eLEnter.grid(row = 4, column =4)
        self.eREnter = ttk.Entry(frm, textvariable=self.xmax, width=7)
        self.eREnter.grid(row = 4, column =5)


        self.yLimVar = Tk.IntVar()
        self.yLimVar.set(self.parent.GetPlotParam('set_ylim'))
        self.yLimVar.trace('w', self.yLimChanged)



        self.ymin = Tk.StringVar()
        self.ymin.set(str(self.parent.GetPlotParam('y_min')))
        self.ymax = Tk.StringVar()
        self.ymax.set(str(self.parent.GetPlotParam('y_max')))


        cb = ttk.Checkbutton(frm, text ='Set log(y) lim',
                        variable = self.yLimVar)
        cb.grid(row = 5,  column =3, sticky = Tk.W)
        self.eLEnter = ttk.Entry(frm, textvariable=self.ymin, width=7)
        self.eLEnter.grid(row = 5, column =4)
        self.eREnter = ttk.Entry(frm, textvariable=self.ymax, width=7)
        self.eREnter.grid(row = 5, column =5)

    def ctypeChanged(self, *args):
        if self.ctypevar.get() == self.parent.chartType:
            pass
        else:
            self.parent.ChangePlotType(self.ctypevar.get())
            self.destroy()

    def RadioSpect(self):
        if self.SpectTypeVar.get() == self.parent.GetPlotParam('spectral_type'):
            pass
        else:
            self.parent.SetPlotParam('spectral_type', self.SpectTypeVar.get())

    def TxtEnter(self, e):
        self.FieldsCallback()

    def FieldsCallback(self):
        tkvarLimList = [self.xmin, self.xmax, self.ymin, self.ymax]
        plot_param_List = ['x_min', 'x_max', 'y_min', 'y_max']
        tkvarSetList = [self.xLimVar, self.xLimVar, self.yLimVar, self.yLimVar]
        to_reload = False
        for j in range(len(tkvarLimList)):
            try:
            #make sure the user types in a int
                if np.abs(float(tkvarLimList[j].get()) - self.parent.GetPlotParam(plot_param_List[j])) > 1E-4:
                    self.parent.SetPlotParam(plot_param_List[j], float(tkvarLimList[j].get()), update_plot = False)
                    to_reload += True*tkvarSetList[j].get()

            except ValueError:
                #if they type in random stuff, just set it ot the param value
                tkvarLimList[j].set(str(self.parent.GetPlotParam(plot_param_List[j])))
        if to_reload:
            self.parent.SetPlotParam('x_min', self.parent.GetPlotParam('x_min'))
    def xLimChanged(self, *args):
        if self.xLimVar.get() == self.parent.GetPlotParam('set_xlim'):
            pass
        else:
            self.parent.SetPlotParam('set_xlim', self.xLimVar.get())

    def yLimChanged(self, *args):
        if self.yLimVar.get() == self.parent.GetPlotParam('set_ylim'):
            pass
        else:
            self.parent.SetPlotParam('set_ylim', self.yLimVar.get())

    def OnClosing(self):
        self.parent.settings_window = None
        self.destroy()
