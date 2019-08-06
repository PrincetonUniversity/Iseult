import matplotlib as mpl
mpl.use('Tkagg')
#mpl.use('agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.ticker as ticker
import numpy as np
import numpy.ma as ma
from new_cnorms import PowerNormWithNeg, PowerNormFunc
import matplotlib.colors as mcolors
from scipy.integrate import simps, cumtrapz
from scipy.special import kn
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
import h5py
import os, sys
import re
import new_cmaps
import Tkinter, tkFileDialog
import matplotlib.animation as manimation
from scipy.stats import linregress
import matplotlib.patheffects as PathEffects
from matplotlib import transforms
from matplotlib.offsetbox import TextArea, VPacker, AnnotationBbox

#plt.style.use('ggplot')
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams['xtick.direction']='out'
mpl.rcParams['ytick.direction']='out'
DS = 1 # 1== downstream

#My own color maps
colors =  ['#D90000','#04756F','#FF8C00', '#988ED5', '#2E0927', '#0971B2']
e_color = '#0971B2'
i_color = '#D90000'
gradient =  np.linspace(0, 1, 256)# A way to make the colorbar display better
gradient = np.vstack((gradient, gradient))

Te = 0.003
Tp = 0.00508


ax_label_size = 24
ticklabel_size = 20

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

#dir_names = ['./Gam1.5', './Gam3','./Gam5','./Gam5']#,'../SuperLarge10deg/output']
dir_names = ['../Gam1.5XXXL/output', '../Gam3XXL/output','../Gam5DS/output']#,'../SuperLarge10deg/output']
#dir_names = ['../Gam3XXL/output', '../Gam3XXL/output','../Gam5DS/output']#,'../SuperLarge10deg/output']

sim_names = [r'$\Gamma_0 =1.5$',r'$\Gamma_0 = 3$', r'$\Gamma_0 = 5$', '']
#sim_i_regions = [(-320,-40),(-640,-80), (-1010,-125),(-640,-80)]
sim_i_regions = [(40,80),(80,160), (125,250),(125,250)]
#sim_e_regions = [(-320, -40),(-640,-80),(-1010,-125),(-640,-80)]
sim_e_regions = [(40,80),(80,160), (125,250),(125,250)]
if DS:
    sim_i_regions = [(-1500, -500),(-1500,-500),(-1500,-500),(-1000,-80)]
    sim_e_regions = [(-1500, -500),(-1500,-500),(-1500,-500),(-1000,-80)]

#sim_time = [28, 27, 44]
#ylims=  [(0,79), (210,289),(0,79)]
ylims=  [(0,79), (0,79),(0,79), (0,79)]
#ylims=  [(0,79), (490,569),(0,79)]
sim_time = [-1, -1, -1, -1]

ls = ['-', '-','-',':']
lcolor  = ['#FF8C00',  '#2E0927','#04756F','#04756F']


SimList = []
for i in range(len(dir_names)):
    sim = simulation()
    SimList.append(sim)
    sim.name = sim_names[i]
    sim.label = sim_names[i]
    sim.dr = os.path.join(os.curdir,dir_names[i])
    sim.i_region = sim_i_regions[i]
    sim.e_region = sim_e_regions[i]
    sim.ylims = ylims[i]
    sim.time_slice = sim_time[i]
    sim.ls = ls[i]
    sim.color = lcolor[i]
    sim.dens_norm = 0.75
    sim.dens_midpoint = 1.0
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
    dens_arr = np.average(dens_arr.reshape(-1, dens_arr.shape[-1]), axis = 0)
    while dens_arr[istart]<1E-8:
        istart += 1

    # build the final x_axis of the plot


    jstart = int(min(10*c_omp/istep, sim.nxf0))
    xaxis_final = np.arange(dens_arr.shape[0]-istart)/c_omp*istep
    # Find the shock by seeing where the density is 1/2 of it's
    # max value. First average the density in the y_direction


    dens_half_max = max(dens_arr[jstart:])*.5
    ishock_final = np.where(dens_arr[jstart:]>=dens_half_max)[0][-1]-istart
    xshock_final = xaxis_final[ishock_final]
    sim.shock_speed = xshock_final/final_time


####
# Get the spectrum and data for t1 and t2
####
for sim in SimList:



    ToLoad = {'Flds': [], 'Prtl': [], 'Param': [], 'Spect': []}
    tmpList = ['c_omp', 'istep', 'dens','ppc0','gamma', 'ppc0','xsl','spece','specp', 'time', 'mi', 'me', 'gamma0', 'bz', 'by','bx','ey', 'ez','c', 'sigma', 'btheta']

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
                #elif elm =='gamma0':
                #    DataDict[elm] = 1.5
                #elif elm == 'spece':
                #    DataDict[elm] = f['specerest'][:]
                #elif elm == 'specp':
                #    DataDict[elm] = f['specprest'][:]
                else:
                    DataDict[elm] = f[elm][:]

    dens_arr = DataDict['dens'][0,:,:]
    # Find out where the left wall is:
    sim.mass_ratio = DataDict['mi']/DataDict['me']
    istart = 0
    while dens_arr[dens_arr.shape[0]//2,istart]<1E-8:
        istart += 1
    sim.dens = dens_arr[:,istart:]/DataDict['ppc0']
    sim.p_shock_loc = DataDict['time'][0]*sim.shock_speed
    # All the data for this time step is now stored in the datadict
    # Calculate the current shock location:
    sim.c_omp = DataDict['c_omp'][0]
    sim.istep = DataDict['istep'][0]
    sim.time = r'$\omega_{pi}t=$'+str(int(DataDict['time'][0]//np.sqrt(sim.mass_ratio)))

    shock_loc = DataDict['time'][0]*sim.shock_speed+istart/sim.c_omp*sim.istep
    ### CALCULATE EPS B IN LEFT WALL FRAME
    bx_prime = DataDict['bx'][0,:,istart:]
    by_prime = 1.5*(DataDict['by'][0,:,istart:]+np.sqrt(1-1.5**-2)*DataDict['ez'][0,:,istart:])
    bz_prime = 1.5*(DataDict['bz'][0,:,istart:]-np.sqrt(1-1.5**-2)*DataDict['ey'][0,:,istart:])
    eps_B = (bx_prime**2+by_prime**2+bz_prime**2)/((1.5-1)*DataDict['ppc0']*.5*DataDict['c']**2*(DataDict['mi']+DataDict['me']))
    sim.eps_B = np.mean(eps_B, axis = 0)
    sim.x_array = np.arange(len(sim.eps_B))*sim.istep/sim.c_omp/np.sqrt(sim.mass_ratio)
    sim.eps_B_arr = np.log10(eps_B)
    #### Get all the data:

    c_omp = DataDict['c_omp'][0]
    istep = DataDict['istep'][0]
    xsl = DataDict['xsl']/c_omp
    gamma = DataDict['gamma']
    sim.LF = DataDict['gamma0']
    spece = DataDict['spece']
    specp = DataDict['specp']
    # In output.F90, spece (specp) is defined by the number of electons (ions)
    # divided by gamma in each logrithmic energy bin. So we multiply by gamma.

    for j in range(len(xsl)):
        spece[:,j] *= gamma
        specp[:,j] *= gamma
    dgamma = np.empty(len(gamma))
    delta=np.log10(gamma[-1]/gamma[0])/len(gamma)
    for j in range(len(dgamma)):
        dgamma[j]=gamma[j]*(10**delta-1.)

    # Select the x-range from which to take the spectra
    e_left_loc = shock_loc + sim.e_region[0]
    e_right_loc = shock_loc + sim.e_region[1]

    eL = xsl.searchsorted(e_left_loc)
    eR = xsl.searchsorted(e_right_loc, side='right')
    if eL == eR:
        eR += 1
    i_left_loc = shock_loc + sim.i_region[0]
    i_right_loc = shock_loc + sim.i_region[1]

    iL = xsl.searchsorted(i_left_loc)
    iR = xsl.searchsorted(i_right_loc, side='right')
    #    print iL, eL, iR, eR
    if iL == iR:
        iR += 1

    # total particles in each linear x bin
    norme = np.empty(len(xsl))
    normp = np.empty(len(xsl))

    for k in range(len(norme)):
        norme[k]=sum(spece[:,k])
        normp[k]=sum(specp[:,k])

    # energy distribution, f(E)=(dn/dE)/N
    fe=np.empty(len(gamma))
    fp=np.empty(len(gamma))

    for k in range(len(fe)):
        fe[k]=sum(spece[k][eL:eR])/(sum(norme[eL:eR]))
        fp[k]=sum(specp[k][iL:iR])/(sum(normp[iL:iR]))


    fe[fe<=0]=1E-100
    fp[fp<=0]=1E-100

    # MASK OUT THE ZEROS


    #  NOTE: gamma ---> gamma-1 ***
    edist = np.copy(fe)
    pdist = np.copy(fp)

    masked_edist = np.ma.masked_where(fe < 1E-20, fe)
    masked_pdist = np.ma.masked_where(fp < 1E-20, fp)


    sim.exdata = gamma*(DataDict['me'][0]/DataDict['mi'][0])/(sim.LF-1)
    sim.eydata = masked_edist*sim.exdata
    sim.ixdata = gamma/(sim.LF-1)
    sim.iydata = masked_pdist*sim.ixdata



    delgame0=Te*DataDict['mi'][0]/DataDict['me'][0]
    if delgame0 >= 0.013:
        aconst = 1/(delgame0*np.exp(1.0/delgame0)*kn(2, 1.0/delgame0))
    else:
        aconst = np.sqrt(2/np.pi)/delgame0**1.5
        aconst -= 15.0/(4.*np.sqrt(delgame0)*np.sqrt(2*np.pi))
        aconst += (345*np.sqrt(delgame0))/(64.*np.sqrt(2*np.pi))
        aconst -= (3285*delgame0**1.5)/(512.*np.sqrt(2*np.pi))
        aconst += (95355*delgame0**2.5)/(16384.*np.sqrt(2*np.pi))
        aconst -= (232065*delgame0**3.5)/(131072.*np.sqrt(2*np.pi))



    femax = aconst*gamma*(gamma+1.0)*np.sqrt((gamma+1.0)**2-1)
    femax *= np.exp(-gamma/delgame0)
    maxg=(delgame0)*40+1. #!40 times the temperature is usually enough
    gamma_table = np.empty(200)
    for i in range(len(gamma_table)):
        gamma_table[i]=(maxg-1.)/(len(gamma_table)-1)*(i)#!+1. remove 1. from definition to avoid underflows.
    tmp_arr=(gamma_table+1.)*np.sqrt(gamma_table*(gamma_table+2.))*np.exp(-(gamma_table)/delgame0)
    pdf_table = tmp_arr*gamma_table#len(gamma_table))
    pdf_table=pdf_table/max(pdf_table)#tmp_arr)

    sim.etempy =  pdf_table# femax#*sum(normp[eL:eR])#*gamma**2
    sim.etempx = gamma_table/(sim.mass_ratio)/(sim.LF-1)

    delgami0=Tp
    if delgami0 >= 0.013:
        aconst = 1/(delgami0*np.exp(1.0/delgami0)*kn(2, 1.0/delgami0))
    else:
        aconst = np.sqrt(2/np.pi)/delgami0**1.5
        aconst -= 15.0/(4.*np.sqrt(delgami0)*np.sqrt(2*np.pi))
        aconst += (345*np.sqrt(delgami0))/(64.*np.sqrt(2*np.pi))
        aconst -= (3285*delgami0**1.5)/(512.*np.sqrt(2*np.pi))
        aconst += (95355*delgami0**2.5)/(16384.*np.sqrt(2*np.pi))
        aconst -= (232065*delgami0**3.5)/(131072.*np.sqrt(2*np.pi))

    fimommax = aconst*gamma*(gamma+1.0)*np.sqrt((gamma+1.0)**2-1)
    fimommax *= np.exp(-gamma/delgami0)

    sim.itempy =  fimommax#/np.sqrt(sim.LF**2-1)
    sim.itempx = gamma/(sim.LF-1)



#####
#
# We've now loaded all of the data. Let's make the 2 spectra:
#
#####


# Make the figure
fig = plt.figure(figsize = (7,8))
gsArgs = {'left':0.12, 'right':0.95, 'top':.95, 'bottom':0.1, 'wspace':0.2, 'hspace':0.0}
MainGs = gridspec.GridSpec(1000,1000)
MainGs.update(**gsArgs)
i_ax = fig.add_subplot(MainGs[35:490,50:990]) # axes to plot the downstream spectra
e_ax = fig.add_subplot(MainGs[510:990,50:990]) # axes to plot the upstream spectra


for sim in SimList:
    #Plot the eps_e and eps_p lines
    sim.iplot, = i_ax.loglog(sim.ixdata, sim.iydata, ls = sim.ls, linewidth = 1.5, color = sim.color, label = sim.label)
    sim.eplot, = e_ax.loglog(sim.exdata, sim.eydata, ls = sim.ls, linewidth = 1.5, color = sim.color, label = sim.time)
#    sim.Te, = e_ax.loglog(sim.etempx, sim.etempy, ls = '--', linewidth = 2,color = '#FF2D00' )#, label = sim.label)
#    sim.Tp, = i_ax.loglog(sim.itempx, sim.itempy, ls = '--', linewidth = 2,color = '#FF2D00' )#, label = sim.label)

if DS:
    i_ax.text(8.,.07,r'Ions',size = ax_label_size)
    e_ax.text(4.,.07,r'Electrons',size = ax_label_size)
else:
    i_ax.text(4.,.035,r'Ions',size = ax_label_size)
    e_ax.text(2.,.035,r'Electrons',size = ax_label_size)

ax_list = [i_ax, e_ax]
k = 0
for ax in ax_list:
    # MAKE THE LEGEND OMITING THE 4RD SIM
    #if DS and k == 0:
    #    leg = ax.legend(handlelength=0, handletextpad=0, fontsize = ax_label_size-1, bbox_to_anchor=(0.36,1.05))
    #if k ==0:
        #leg = ax.legend(handlelength=0, handletextpad=0, fontsize = ax_label_size, bbox_to_anchor=(0.36,0.52))
  

    leg = ax.legend(handlelength=0, handletextpad=0, fontsize = ax_label_size, loc='best')        
    leg.get_frame().set_linewidth(0)
    leg.get_frame().set_alpha(0.0)
    for item in leg.legendHandles:
        item.set_visible(False)
    for color,text in zip(lcolor,leg.get_texts()):
        text.set_color(color)

    k +=1

    ax.tick_params(axis='x',
                     which = 'both', # bothe major and minor ticks
                      top = 'off') # turn off top ticks

    ax.tick_params(axis='y',          # changes apply to the y-axis
                   which='major',      # both major and minor ticks are affected
                   left='on',      # ticks along the bottom edge are off
                   right='off')         # ticks along the top edge are off)

    ax.tick_params(axis='y',          # changes apply to the y-axis
                   which='minor',      # both major and minor ticks are affected
                   left='off',      # ticks along the bottom edge are off
                   right='off')         # ticks along the top edge are off


    ax.set_ylabel(r'$E^2 {\rm d} n/{\rm d}E \ [{\rm arb. unit}]$', size = ax_label_size)

    ax.set_yticks([1E-6, 1E-4, 1E-2])
    ax.set_ylim(1E-6,5E-1)
    if DS:
        ax.set_xlim(5E-3,5E1)
        ax.set_ylim(1E-6,5E-1)
    else:
        ax.set_xlim(2E-4,1E2)

#i_ax.add_artist(l1)
#earlyLine = mlines.Line2D([], [], color='K', ls='-', label=r'$\omega_{pi}t= 840$')
#lateLine = mlines.Line2D([], [], color='k', ls=':', label=r'$\omega_{pi}t = 1970$')
#e_ax.legend(handles=[earlyLine, lateLine], fontsize = 20, handletextpad=0, handlelength = .75, loc = 'upper left', fancybox = False, edgecolor = 'w', framealpha = 0)
#i_ax.set_title(r'$\omega_{pi}t/\gamma_0 \approx1330$', size = ax_label_size, loc = 'right')
e_ax.set_xlabel(r'${\rm KE}/(m_i(\gamma_0-1))$', size = ax_label_size)
i_ax.set_xticklabels([])
i_ax.tick_params(axis='x',
                 which = 'both', # bothe major and minor ticks
                 bottom = 'off') # turn off top ticks
fig.suptitle(r'$\theta_B =0,\quad \sigma_0 = 0.01, \quad m_i/m_e = 64$', size = ax_label_size)
plt.savefig('SpectVGam.pdf')
plt.show()
