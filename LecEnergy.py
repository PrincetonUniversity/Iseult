import matplotlib as mpl
mpl.use('Tkagg')
#mpl.use('agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import numpy.ma as ma
from new_cnorms import PowerNormWithNeg, PowerNormFunc
import matplotlib.colors as mcolors
from scipy.integrate import simps, cumtrapz
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
import h5py
import os, sys
import re
import new_cmaps
import matplotlib.animation as manimation
from scipy.stats import linregress
import matplotlib.patheffects as PathEffects
from matplotlib import transforms
from matplotlib.offsetbox import TextArea, VPacker, AnnotationBbox
import numpy as np
from numba import jit, guvectorize, float64, int32
from math import sqrt

#### A vectorized 2D hist

@guvectorize([(float64[:], # x
               float64[:], # y
               float64[:], # u
               float64[:], # v
               float64[:], # w
               float64[:], # bin_widths
               float64[:], # xmin
               float64[:], # ybin_width
               float64[:], # ymin
               float64[:], # gamma (unbinned)
               float64[:,:], # vx binned
               float64[:,:] # counts
               )],
               '(n),(n),(n),(n),(n),(),(), (),(),(n), (m,p), (m,p)', nopython = True, cache = True, target='parallel')
def Calc2DVxEHists(x,y,  u, v, w, xbin_width, xmin, ybin_width, ymin,  g, vx, counts):
    ybn = ybin_width[0]
    miny = ymin[0]
    maxl = counts.shape[0]
    xbn = xbin_width[0]
    minx = xmin[0]
    maxk = counts.shape[1]

    for i in xrange(len(x)):
        l = int((y[i]-miny)//ybn)
        k = int((x[i]-minx)//xbn)
        if l>=0 and l < maxl and k>=0 and k< maxk:
            g[i] = u[i]*u[i]+v[i]*v[i]+ w[i]*w[i] + 1
            g[i] = sqrt(g[i])
            vx[l,k] += u[i]*g[i]**-1
            counts[l,k] += 1
    for l in xrange(vx.shape[0]):
        for k in xrange(vx.shape[1]):
            if counts[l,k] != 0:
                vx[l,k] *= counts[l,k]**-1
@guvectorize([(float64[:,:], float64[:,:], float64[:,:], float64[:,:],float64[:,:],float64[:,:])], '(m,p),(m,p),(m,p),(m,p),(m,p),(m,p)', cache = True,target='parallel')
def RestFrameBoost(vx_e, ecounts, vx_i, icounts, vx_avg, boost_g):
    for i in xrange(vx_e.shape[0]):
        for j in xrange(vx_e.shape[1]):
            if ecounts[i,j] != 0 or icounts[i,j] != 0:
                vx_avg[i,j] = (vx_e[i,j]*ecounts[i,j]+vx_i[i,j]*icounts[i,j])/(icounts[i,j]+ecounts[i,j])
                boost_g[i,j] = 1/sqrt(1-vx_avg[i,j]*vx_avg[i,j])

@guvectorize([(float64[:], # x
               float64[:], # y
               float64[:], # u
               float64[:], # v
               float64[:], # w
               float64[:], # g
               float64[:,:], # vx_avg
               float64[:,:], # boost_g,
               float64[:], # bin_width
               float64[:], # xmin
               float64[:], # ybin_width
               float64[:], # ymin
               float64[:], # Emax
               float64[:,:], # Counts
               float64[:,:] # T
               )],
               '(n),(n),(n),(n),(n),(n),(m,p),(m,p),(),(),(),(),(),(m,p),(m,p)', nopython = True, cache = True, target='parallel')
def CalcDelGamHists(x, y,u, v, w, g, vx_avg, boost_g, xbin_width, xmin, ybin_width, ymin, Emin, counts, T):
    ybn = ybin_width[0]
    lmax = vx_avg.shape[0]
    miny = ymin[0]
    xbn = xbin_width[0]
    kmax = vx_avg.shape[1]
    minx = xmin[0]
    for i in xrange(len(x)):
        l = int((y[i]-miny)//ybn)
        k = int((x[i]-minx)//xbn)
        if 0 <= l and l <lmax and k>0 and k<kmax: 
            c2 = boost_g[l,k]*(u[i]-g[i]*vx_avg[l,k]) # boosted
            #counts[l,k] += sqrt(c2*c2+v[i]*v[i]+w[i]*w[i]+1)-1

            counts[l,k] +=1 
            T[l,k] += sqrt(c2*c2+v[i]*v[i]+w[i]*w[i]+1)-1
            

    for l in xrange(counts.shape[0]):
        c = 0
        for k in xrange(counts.shape[1]):
            c += counts[l,k]
        if c != 0:
            T[l,:] *= counts.shape[1]*c**-1



@guvectorize([(float64[:], # x
               float64[:], # y
               float64[:], # u
               float64[:], # v
               float64[:], # w
               float64[:], # x_bin_width
               float64[:], # xmin
               float64[:], # y_bin_width
               float64[:], # ymin
               float64[:], # Emin
               float64[:,:], # 2D E_avg binned
               float64[:,:] # 2D counts
               )],
               '(n),(n),(n),(n),(n),(), (),(),(),(),(m,p), (m,p)', nopython = True, cache = True, target='parallel')
def CalcE2DHists(x, y, u, v, w, xbin_width, xmin, ybin_width, ymin, Emin, twoDhist, counts):
    # We follow the convention of numpy.histogram2D where first dimension is y direction and second is x

    ybn = ybin_width[0]
    miny = ymin[0]
    maxl = twoDhist.shape[0]
    xbn = xbin_width[0]
    minx = xmin[0]
    maxk = twoDhist.shape[1]

    for i in xrange(len(x)):
        l = int((y[i]-miny)//ybn)
        k = int((x[i]-minx)//xbn)

        if l>=0 and l < maxl and k>=0 and k< maxk:
            E = sqrt(u[i]*u[i]+v[i]*v[i]+ w[i]*w[i] + 1)-1
            twoDhist[l,k] += E 
            counts[l,k] += 1
    for l in xrange(twoDhist.shape[0]):
        for k in xrange(twoDhist.shape[1]):
            if counts[l,k] != 0:
                twoDhist[l,k] *= counts[l,k]**-1
#### A vectorized 2D hist

@guvectorize([(float64[:], # x
               float64[:], # y
               float64[:], # u
               float64[:], # x_bin_width
               float64[:], # xmin
               float64[:], # y_bin_width
               float64[:], # ymin
               float64[:,:], # 2D E_avg binned
               float64[:,:] # 2D counts
               )],
               '(n),(n),(n),(), (),(),(),(m,p), (m,p)', nopython = True, cache = True, target='parallel')
def Calcpx2DHists(x, y, u,  xbin_width, xmin, ybin_width, ymin, twoDhist, counts):
    # We follow the convention of numpy.histogram2D where first dimension is y direction and second is x

    ybn = ybin_width[0]
    miny = ymin[0]
    maxl = twoDhist.shape[0]
    xbn = xbin_width[0]
    minx = xmin[0]
    maxk = twoDhist.shape[1]

    for i in xrange(len(x)):
        l = int((y[i]-miny)//ybn)
        k = int((x[i]-minx)//xbn)

        if l>=0 and l < maxl and k>=0 and k< maxk:
            twoDhist[l,k] += u[i] 
            counts[l,k] += 1
    #for l in xrange(twoDhist.shape[0]):
    #    for k in xrange(twoDhist.shape[1]):
    #        if counts[l,k] != 0:
    #            twoDhist[l,k] *= counts[l,k]**-1



#plt.style.use('ggplot')
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'


#My own color maps
colors =  ['#D90000','#04756F','#FF8C00', '#988ED5', '#2E0927', '#0971B2']
e_color = '#0971B2'
i_color = '#D90000'
gradient =  np.linspace(0, 1, 256)# A way to make the colorbar display better
gradient = np.vstack((gradient, gradient))


ax_label_size = 30
ticklabel_size = 30

mpl.rcParams['xtick.labelsize'] = ticklabel_size
mpl.rcParams['ytick.labelsize'] = ticklabel_size


# Be sure to call only this with the output directory in the path
# create a bunch of regular expressions used to search for files
f_re = re.compile('flds.tot.*')
prtl_re = re.compile('prtl.tot.*')
s_re = re.compile('spect.*')
param_re = re.compile('param.*')
re_list = [f_re, prtl_re, s_re, param_re]
#root = Tkinter.Tk()



class simulation:
    # a placeholder to hold all the simulation data
    np.nan

dir_names = ['../Gam1.5LowCompSmallerdelY/output']
sim_names = [r'$\theta_B\approx 10^\circ$']
ylims=  [(None,None)]#
xlims = [(None,None)]
sim_time = [-1]


SimList = []
for i in range(len(dir_names)):
    sim = simulation()
    SimList.append(sim)
    sim.name = sim_names[i]
    sim.label = sim_names[i]
    sim.dr = os.path.join(os.curdir,dir_names[i])
    sim.ylims = ylims[i]
    sim.xlims = xlims[i]
    sim.time_slice = sim_time[i]

for sim in SimList:
    sim.PathDict = {}
    #fill all the paths
    sim.PathDict['Flds']= filter(f_re.match, os.listdir(sim.dr))
    sim.PathDict['Flds'].sort()

    sim.PathDict['Prtl']= filter(prtl_re.match, os.listdir(sim.dr))
    sim.PathDict['Prtl'].sort()

    sim.PathDict['Spect']= filter(s_re.match, os.listdir(sim.dr))
    sim.PathDict['Spect'].sort()

    sim.PathDict['Param']= filter(param_re.match, os.listdir(sim.dr))
    sim.PathDict['Param'].sort()


# A dictionary that allows use to see in what HDF5 file each key is stored.
# i.e. {'ui': 'Prtl', 'ue': 'Flds', etc...},

H5KeyDict = {}
for pkey in SimList[0].PathDict.keys():
    with h5py.File(os.path.join(SimList[0].dr,SimList[0].PathDict[pkey][0]), 'r') as f:
    # Because dens is in both spect* files and flds* files,
        for h5key in f.keys():
            if h5key == 'dens' and pkey == 'Spect':
                H5KeyDict['spect_dens'] = pkey
            else:
                H5KeyDict[h5key] = pkey
                H5KeyDict['time'] = 'Param'
# Get the shock speeds:
# First load the first field file to find the initial size of the
# box in the x direction

for sim in SimList:
    with h5py.File(os.path.join(sim.dr,sim.PathDict['Flds'][0]), 'r') as f:
        sim.nxf0 = f['by'][:].shape[1]
    with h5py.File(os.path.join(sim.dr,sim.PathDict['Param'][0]), 'r') as f:
        sim.initial_time = f['time'][0]

    # Load the final time step to find the shock's location at the end.
    with h5py.File(os.path.join(sim.dr,sim.PathDict['Flds'][-1]), 'r') as f:
        dens_arr =np.copy(f['dens'][0,:,:])

    with h5py.File(os.path.join(sim.dr,sim.PathDict['Param'][-1]), 'r') as f:
        # I use this file to get the final time, the istep, interval, and c_omp
        final_time = f['time'][0]
        istep = f['istep'][0]
        c_omp = f['c_omp'][0]

    # Since we're in the piston frame, we have to account for the wall on the left
    istart = 0
    while dens_arr[dens_arr.shape[0]/2,istart]<1E-8:
        istart += 1

    # build the final x_axis of the plot


    jstart = int(min(10*c_omp/istep, sim.nxf0))
    xaxis_final = np.arange(dens_arr.shape[1]-istart)/c_omp*istep
    # Find the shock by seeing where the density is 1/2 of it's
    # max value. First average the density in the y_direction


    dens_half_max = max(dens_arr[dens_arr.shape[0]/2,jstart:])*.5
    ishock_final = np.where(dens_arr[dens_arr.shape[0]/2,jstart:]>=dens_half_max)[0][-1]-istart
    xshock_final = xaxis_final[ishock_final]
    sim.shock_speed = xshock_final/final_time


####
# Get the spectrum and data for t1 and t2
####
for sim in SimList:



    ToLoad = {'Flds': [], 'Prtl': [], 'Param': [], 'Spect': []}
    tmpList = ['c_omp', 'istep', 'dens','ppc0','gamma', 'ppc0','xsl','spece','specp', 'time', 'mi', 'me', 'gamma0', 'bz', 'by','bx','ey', 'ez','c', 'sigma', 'btheta',  'xi', 'yi', 'ui', 'vi', 'wi','xe', 'ye', 'ue', 've', 'we']

    for elm in tmpList:
        # find out what type of file the key is stored in
        ftype = H5KeyDict[elm]
        # add the key to the list of that file type
        ToLoad[ftype].append(elm)

    # Now iterate over each path key and create a datadictionary
    DataDict = {}
    for pkey in ToLoad.keys():
        tmplist = list(set(ToLoad[pkey])) # get rid of duplicate keys
        # Load the file
        with h5py.File(os.path.join(sim.dr,sim.PathDict[pkey][sim.time_slice]), 'r') as f:
            for elm in tmplist:
                # Load all the keys
                if elm == 'spect_dens':
                    DataDict[elm] = f['dens'][:]
                elif elm == 'spece':
                    DataDict[elm] = f['specerest'][:]
                elif elm == 'specp':
                    DataDict[elm] = f['specprest'][:]
                else:
                    DataDict[elm] = f[elm][:]
    dens_arr = DataDict['dens'][0,:,:]
    # Find out where the left wall is:
    sim.mass_ratio = DataDict['mi']/DataDict['me']
    istart = 0
    while dens_arr[dens_arr.shape[0]//2,istart]<1E-8:
        istart += 1
    iend = dens_arr.shape[1]
    while dens_arr[dens_arr.shape[0]//2,iend-1]<1E-8:
        iend += -1
    sim.dens = dens_arr[:,istart:iend]/DataDict['ppc0']
    sim.c_omp = DataDict['c_omp'][0]
    sim.istep = DataDict['istep'][0]
    sim.time = DataDict['time'][0]/np.sqrt(sim.mass_ratio)
    #Let's calculate our 2DHist First let's do xbins =600,
    xbins = 300
    ybins = 100
    xwidth = (sim.dens.shape[1]*sim.istep)//xbins
    ywidth = (sim.dens.shape[0]*sim.istep)//ybins

    MyHist = np.zeros((ybins,xbins))
    MyCounts = np.zeros((ybins,xbins))
    #CalcE2DHists(DataDict['xi'], DataDict['yi'], DataDict['ui'], DataDict['vi'], DataDict['wi'], xwidth, np.array([0]), ywidth, np.array([0]), 2., MyHist, MyCounts)
    Calcpx2DHists(DataDict['xi'], DataDict['yi'], DataDict['ui'], xwidth, np.array([0]), ywidth, np.array([0]), MyHist, MyCounts)
    MyHist *=np.sqrt(DataDict['gamma0'][0]**2-1)**-1

    icounts = np.zeros((ybins,xbins))
    vix = np.zeros((ybins,xbins))
    gi = np.zeros(len(DataDict['xi']))
    Calc2DVxEHists(DataDict['xi'],DataDict['yi'],  DataDict['ui'], DataDict['vi'], DataDict['wi'], xwidth, 0, ywidth, 0,  gi, vix, icounts)

    ecounts = np.zeros((ybins,xbins))
    vex = np.zeros((ybins,xbins))
    ge = np.zeros(len(DataDict['xe']))
    Calc2DVxEHists(DataDict['xe'],DataDict['ye'],  DataDict['ue'], DataDict['ve'], DataDict['we'], xwidth, 0, ywidth, 0,  ge, vex, ecounts)

    vx =np.zeros(ecounts.shape)
    boost_g =np.zeros(ecounts.shape)
    RestFrameBoost(vex, ecounts, vix, icounts, vx, boost_g)
    #MyHist =np.copy(vx)
    MyHist2 = np.zeros((ybins,xbins))
    MyCounts2 = np.zeros((ybins,xbins))
    CalcDelGamHists(DataDict['xi'], DataDict['yi'], DataDict['ui'], DataDict['vi'], DataDict['wi'], gi, vx, boost_g, xwidth, np.array([0]), ywidth, np.array([0]), 1.0, MyCounts2,MyHist2)
    #CalcDelGamHists(DataDict['xi'], DataDict['yi'], DataDict['ui'], DataDict['vi'], DataDict['wi'], gi, vx, boost_g, xwidth, np.array([0]), ywidth, np.array([0]), 4.0, MyCounts2,MyHist2)
    MyHist2 = MyHist2
    MyeHist2 = np.zeros((ybins,xbins))
    MyeCounts2 = np.zeros((ybins,xbins))
    CalcDelGamHists(DataDict['xe'], DataDict['ye'], DataDict['ue'], DataDict['ve'], DataDict['we'], ge, vx, boost_g, xwidth, np.array([0]), ywidth, np.array([0]), 32.0, MyeCounts2,MyeHist2)
    #CalcDelGamHists(DataDict['xe'], DataDict['ye'], DataDict['ue'], DataDict['ve'], DataDict['we'], ge, vx, boost_g, xwidth, np.array([0]), ywidth, np.array([0]), 128.0, MyeCounts2,MyeHist2)
    #MyeHist2 = np.log10(MyeHist2)
    MyeHist2 *= sim.mass_ratio**-1
### MAKE THE PLOTS

gsArgs = {'left':0.1, 'right':0.9, 'top':.98, 'bottom':0.1, 'wspace':0.2, 'hspace':0.1}
MainGS_starting = 45
MainGS_spacing = 40
MainGS_ending = 935

fig = plt.figure(figsize = (12,9))

MainGs = gridspec.GridSpec(1000,1)
graph_extent = int((MainGS_ending - MainGS_starting - 1*MainGS_spacing)//4.)-10
curpos = MainGS_starting
MainGs.update(**gsArgs)
#fig.suptitle('Superluminal and Subluminal shocks at time' + r'$\omega_{pi}t \approx 2600$', size = 30)

ax = fig.add_subplot(MainGs[curpos:curpos+graph_extent,:])
ax.grid(False)
curpos += graph_extent + MainGS_spacing

im = ax.imshow(sim.dens, origin = 'lower')

ax.set_ylabel(r'$y \ [c/\omega_{pi}]$', size = ax_label_size)
im.set_cmap(new_cmaps.cmaps['temperature'])
im.set_extent([0, sim.dens.shape[1]*sim.istep/sim.c_omp/np.sqrt(sim.mass_ratio), 0, sim.dens.shape[0]*sim.istep/sim.c_omp/np.sqrt(sim.mass_ratio)])
ax.set_xlim(sim.xlims)
annotate_kwargs = {'horizontalalignment': 'center',
                   'verticalalignment': 'top',
                   'size' : 30,
                       #'path_effects' : [PathEffects.withStroke(linewidth=1.5,foreground="k")],
                    'xycoords': 'axes fraction',
                    'color' : 'k'}
bbox_props = dict(fc="k", alpha=0.0)

im.norm.vmin = 0
im.norm.vmax = 5
#ax.locator_params(axis='y',nbins=4)
ax.set_ylim(sim.ylims)

ax.tick_params(labelsize = ticklabel_size)#, color=tick_color)
ax.tick_params(axis ='y', which='both', right='off')
ax.tick_params(axis ='x', which='both', top='off')
if i != 1:
    ax.tick_params(axis ='x', which='both', bottom='off')
    #ax.set_xticklabels([])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.1)


# Cbar is vertical
cbar = cax.imshow(np.transpose(gradient)[::-1], aspect='auto',
                                            cmap=new_cmaps.cmaps['temperature'])

cax.tick_params(axis='x',
                which = 'both', # bothe major and minor ticks
                top = 'off', # turn off top ticks
                bottom = 'off', # turn off bottom ticks
                labelbottom = 'off') # turn off bottom ticks
cax.grid(False)

cax.tick_params(axis='y',          # changes apply to the y-axis
                which='both',      # both major and minor ticks are affected
                left='off',      # ticks along the bottom edge are off
                right='on',         # ticks along the top edge are off
                labelleft='off',
                labelright = 'on',
                labelsize = ticklabel_size)


clim = np.copy(im.get_clim())
cbar.set_extent([0,1,clim[0],clim[1]])

# re-create the gradient with the data values
# First make a colorbar in the negative region that is linear in the pow_space
data_range = np.linspace(clim[0],clim[1],512)

cbar.set_data(np.transpose(gradient)[::-1])

cax.set_ylim(clim[0],clim[1])


cax.set_ylabel(r'$n/n_0$', size = ax_label_size, rotation = -90, labelpad = 35)
cax.yaxis.set_label_position("right")
ax.set_xlabel(r'$x\ [c/\omega_{pi}]$', size = ax_label_size)


##### MAKE THE OTHER PLOT
ax = fig.add_subplot(MainGs[curpos:curpos+graph_extent,:], sharex = ax, sharey = ax)
ax.grid(False)
curpos += graph_extent + MainGS_spacing

im = ax.imshow(MyHist,origin = 'lower', aspect = None)#, interpolation = 'Gaussian')
                    #interpolation=self.GetPlotParam('interpolation'))

ax.set_ylabel(r'$y \ [c/\omega_{pi}]$', size = ax_label_size)
im.set_cmap(new_cmaps.cmaps['BuYlRd'])
im.set_extent([0, sim.dens.shape[1]*sim.istep/sim.c_omp/np.sqrt(sim.mass_ratio), 0, sim.dens.shape[0]*sim.istep/sim.c_omp/np.sqrt(sim.mass_ratio)])
ax.set_xlim(sim.xlims)

im.norm.vmin = -500
im.norm.vmax = 500
#ax.locator_params(axis='y',nbins=4)
ax.set_ylim(sim.ylims)

ax.tick_params(labelsize = ticklabel_size)#, color=tick_color)
ax.tick_params(axis ='y', which='both', right='off')
ax.tick_params(axis ='x', which='both', top='off')
if i != 1:
    ax.tick_params(axis ='x', which='both', bottom='off')
    #ax.set_xticklabels([])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.1)


# Cbar is vertical
cbar = cax.imshow(np.transpose(gradient)[::-1], aspect='auto',
                                            cmap=new_cmaps.cmaps['BuYlRd'])

cax.tick_params(axis='x',
                which = 'both', # bothe major and minor ticks
                top = 'off', # turn off top ticks
                bottom = 'off', # turn off bottom ticks
                labelbottom = 'off') # turn off bottom ticks
cax.grid(False)

cax.tick_params(axis='y',          # changes apply to the y-axis
                which='both',      # both major and minor ticks are affected
                left='off',      # ticks along the bottom edge are off
                right='on',         # ticks along the top edge are off
                labelleft='off',
                labelright = 'on',
                labelsize = ticklabel_size)


clim = np.copy(im.get_clim())
cbar.set_extent([0,1,clim[0],clim[1]])

# re-create the gradient with the data values
# First make a colorbar in the negative region that is linear in the pow_space
data_range = np.linspace(clim[0],clim[1],512)

cbar.set_data(np.transpose(gradient)[::-1])


cax.set_ylabel(r'$\langle v_{x}\rangle /c$', size = ax_label_size, rotation = -90, labelpad = 25)
cax.yaxis.set_label_position("right")
ax.set_xlabel(r'$x\ [c/\omega_{pi}]$', size = ax_label_size)

##### MAKE THE OTHER PLOT

ax = fig.add_subplot(MainGs[curpos:curpos+graph_extent,:], sharex = ax, sharey = ax)
ax.grid(False)
curpos += graph_extent + MainGS_spacing

im = ax.imshow(MyHist2,origin = 'lower', aspect = None)#, interpolation = 'Gaussian')
                    #interpolation=self.GetPlotParam('interpolation'))

ax.set_ylabel(r'$y \ [c/\omega_{pi}]$', size = ax_label_size)
#im.set_cmap(new_cmaps.cmaps['Blue/Green/Red/Yellow'])
im.set_cmap(new_cmaps.cmaps['temperature'])
im.set_extent([0, sim.dens.shape[1]*sim.istep/sim.c_omp/np.sqrt(sim.mass_ratio), 0, sim.dens.shape[0]*sim.istep/sim.c_omp/np.sqrt(sim.mass_ratio)])
ax.set_xlim(sim.xlims)

im.norm.vmin = 0
#im.norm.vmax = 1.2
#im.norm.vmax = 1
#ax.locator_params(axis='y',nbins=4)
ax.set_ylim(sim.ylims)

ax.tick_params(labelsize = ticklabel_size)#, color=tick_color)
ax.tick_params(axis ='y', which='both', right='off')
ax.tick_params(axis ='x', which='both', top='off')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.1)


# Cbar is vertical
cbar = cax.imshow(np.transpose(gradient)[::-1], aspect='auto',
                                            #cmap=new_cmaps.cmaps['Blue/Green/Red/Yellow'])
                    cmap=new_cmaps.cmaps['temperature'])

cax.tick_params(axis='x',
                which = 'both', # bothe major and minor ticks
                top = 'off', # turn off top ticks
                bottom = 'off', # turn off bottom ticks
                labelbottom = 'off') # turn off bottom ticks
cax.grid(False)

cax.tick_params(axis='y',          # changes apply to the y-axis
                which='both',      # both major and minor ticks are affected
                left='off',      # ticks along the bottom edge are off
                right='on',         # ticks along the top edge are off
                labelleft='off',
                labelright = 'on',
                labelsize = ticklabel_size)


clim = np.copy(im.get_clim())
cbar.set_extent([0,1,clim[0],clim[1]])

# re-create the gradient with the data values
# First make a colorbar in the negative region that is linear in the pow_space
data_range = np.linspace(clim[0],clim[1],512)

cbar.set_data(np.transpose(gradient)[::-1])


cax.set_ylabel(r'$U^\prime_i$', size = ax_label_size, rotation = -90, labelpad = 10)
#cax.set_yticks([0,.5,1])
cax.yaxis.set_label_position("right")
##### MAKE THE OTHER PLOT

ax = fig.add_subplot(MainGs[curpos:curpos+graph_extent,:], sharex = ax, sharey = ax)
ax.grid(False)
curpos += graph_extent + MainGS_spacing

im = ax.imshow(MyeHist2,origin = 'lower', aspect = None)#, interpolation = 'Gaussian')
                    #interpolation=self.GetPlotParam('interpolation'))

ax.set_ylabel(r'$y \ [c/\omega_{pi}]$', size = ax_label_size)
#im.set_cmap(new_cmaps.cmaps['Blue/Green/Red/Yellow'])
im.set_cmap(new_cmaps.cmaps['temperature'])
im.set_extent([0, sim.dens.shape[1]*sim.istep/sim.c_omp/np.sqrt(sim.mass_ratio), 0, sim.dens.shape[0]*sim.istep/sim.c_omp/np.sqrt(sim.mass_ratio)])
ax.set_xlim(sim.xlims)

#im.norm.vmin = 0
#im.norm.vmax = 1.2
#im.norm.vmax = 250
#ax.locator_params(axis='y',nbins=4)
ax.set_ylim(sim.ylims)

ax.tick_params(labelsize = ticklabel_size)#, color=tick_color)
ax.tick_params(axis ='y', which='both', right='off')
ax.tick_params(axis ='x', which='both', top='off')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.1)


# Cbar is vertical
cbar = cax.imshow(np.transpose(gradient)[::-1], aspect='auto',
                  #cmap=new_cmaps.cmaps['Blue/Green/Red/Yellow'])
                  cmap=new_cmaps.cmaps['temperature'])

cax.tick_params(axis='x',
                which = 'both', # bothe major and minor ticks
                top = 'off', # turn off top ticks
                bottom = 'off', # turn off bottom ticks
                labelbottom = 'off') # turn off bottom ticks
cax.grid(False)

cax.tick_params(axis='y',          # changes apply to the y-axis
                which='both',      # both major and minor ticks are affected
                left='off',      # ticks along the bottom edge are off
                right='on',         # ticks along the top edge are off
                labelleft='off',
                labelright = 'on',
                labelsize = ticklabel_size)


clim = np.copy(im.get_clim())
cbar.set_extent([0,1,clim[0],clim[1]])

# re-create the gradient with the data values
# First make a colorbar in the negative region that is linear in the pow_space
data_range = np.linspace(clim[0],clim[1],512)

cbar.set_data(np.transpose(gradient)[::-1])


#cax.set_ylabel(r'$\langle E^\prime_e \rangle/m_i c^2$', size = ax_label_size, rotation = -90, labelpad = 50)
cax.set_ylabel(r'$U^\prime_e$', size = ax_label_size, rotation = -90, labelpad = 10)
cax.yaxis.set_label_position("right")
#cax.set_yticks([0,.5,1])
ax.set_xlabel(r'$x\ [c/\omega_{pi}]$', size = ax_label_size)
plt.suptitle(r'$\omega_{pi} t = $'+str(int(sim.time)), size = ax_label_size)
plt.show()
