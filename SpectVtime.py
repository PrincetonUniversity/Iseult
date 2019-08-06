import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import h5py
import new_cmaps
import os, sys
import re

import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FuncFormatter
from scipy.signal import savgol_filter
from scipy.stats import linregress
from scipy.integrate import simps, cumtrapz
from sklearn.isotonic import IsotonicRegression


from scipy.special import kn # Modified Bessel function

#plt.style.use('ggplot')
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams['xtick.direction']= 'out'
mpl.rcParams['ytick.direction']= 'out'
#gradient = np.vstack((gradient, gradient))


#Tp = 0.15
Tp = 0.32
Te = .8*Tp
ax_label_size = 22
ticklabel_size = 22

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
skip_size = 5
max_num = 83
sim = simulation()
##sim.dr = os.path.join(os.curdir,'../StampedeGamBeta0_33_sigma0_1/')
#sim.dr = os.path.join(os.curdir,'../StampedeGamBeta0_33/')
sim.dr = os.path.join(os.curdir,'../StampedeMorePPC/gambeta0_33/')
#sim.dr = os.path.join(os.curdir,'../StampedeGamBeta3/')
#sim.dr = os.path.join(os.curdir,'../StampedeGamBeta1/savedRun/')
#sim.dr = os.path.join(os.curdir,'../Gam3LowComp/output')
#sim.dr = os.path.join(os.curdir,'../Gam1.5LowCompSmallerdelY/output')

relative = True

#sim.i_region = (-100000,100000)
#sim.i_region = (300,700)
sim.i_region = (-1000,-100)
#sim.i_region = (250,500)
#sim.e_region = (-100000,0)
#sim.e_region = (300,700)
sim.e_region = (-1000,-100)
#sim.e_region = (550, 500)
#sim.name = r'$\sigma = 0.01\quad \theta \approx 10^\circ \quad \gamma_0 = 1.5$'
sim.ls = '-'

cmap = new_cmaps.cmaps['deep']

#######
#
# Get all the file paths
#
########

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
for pkey in sim.PathDict.keys():
    with h5py.File(os.path.join(sim.dr,sim.PathDict[pkey][0]), 'r') as f:
    # Because dens is in both spect* files and flds* files,
        for h5key in f.keys():
            if h5key == 'dens' and pkey == 'Spect':
                H5KeyDict['spect_dens'] = pkey
            else:
                H5KeyDict[h5key] = pkey
H5KeyDict['time'] = 'Param'


#which time slices?
if max_num == -1:
    sim.files = (np.arange(len(sim.PathDict['Param'])//skip_size)+1)*skip_size -1# Every 4th timestep, but since it wasn't run until 60 we have to cheat the last time
    sim.files[-1] = -1
else:
    sim.files = (np.arange(max_num//skip_size)+1)*skip_size -1# Every 4th timestep, but since it wasn't run until 60 we have to cheat the last time
    sim.files[-1] = max_num

#sim.files = sim.files[2:]
print sim.files
ToLoad = {'Flds': [], 'Prtl': [], 'Param': [], 'Spect': []}
tmpList = ['c_omp', 'istep', 'dens','ppc0','gamma', 'xsl','spece','specp', 'time', 'mi', 'me', 'gamma0', 'c', 'sigma', 'btheta']
sim.ix = []
sim.iy = []
sim.ex = []
sim.ey = []
sim.t  = []

for i in sim.files:

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
        with h5py.File(os.path.join(sim.dr,sim.PathDict[pkey][i]), 'r') as f:
            for elm in tmplist:
                # Load all the keys
                if elm == 'spect_dens':
                    DataDict[elm] = f['dens'][:]
                elif elm =='gamma0':
                    #DataDict[elm] = 1.41
                    DataDict[elm] = 1.055
                    #DataDict[elm] = 3.31
                elif elm == 'spece':
                    #DataDict[elm] = f['specerest'][:]
                    DataDict[elm] = f['spece'][:]
                elif elm == 'specp':
                    #DataDict[elm] = f['specprest'][:]
                    DataDict[elm] = f['specp'][:]
                else:
                    DataDict[elm] = f[elm][:]


    sim.mass_ratio = DataDict['mi']/DataDict['me']
    istart = 0


    x = np.arange(DataDict['dens'][:,:,:].shape[-1])#[::-1]
    y = np.average(DataDict['dens'][:,:,:].reshape(-1, x.shape[0]), axis = 0)[::-1]



    # #############################################################################
    # Fit IsotonicRegression and LinearRegression models

    ir = IsotonicRegression()

    y_ = ir.fit_transform(x, y)



        
    # All the data for this time step is now stored in the datadict
    # Calculate the current shock location:
    sim.c_omp = DataDict['c_omp'][0]
    sim.istep = DataDict['istep'][0]

    shock_loc = (len(y_)-y_.searchsorted(max(y_)*.5))*sim.istep/sim.c_omp
    print sim.mass_ratio
    sim.t.append(DataDict['time'][0]/np.sqrt(sim.mass_ratio))
    # figure out the indices of xL and x
    xmax = (DataDict['dens'].shape[2]-istart)*sim.istep/sim.c_omp
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
    e_left_loc = sim.e_region[0]
    e_right_loc = sim.e_region[1]
    if relative:
        e_left_loc += shock_loc
        e_right_loc += shock_loc
    eL = xsl.searchsorted(e_left_loc)
    eR = xsl.searchsorted(e_right_loc, side='right')
    if eL == eR:
        eR += 1
    i_left_loc = sim.i_region[0]
    i_right_loc = sim.i_region[1]
    if relative:
        i_left_loc += shock_loc
        i_right_loc += shock_loc

    iL = xsl.searchsorted(i_left_loc)
    iR = xsl.searchsorted(i_right_loc, side='right')
    #    print iL, eL, iR, eR
    if iL == iR:
        iR += 1
    if iR == len(xsl):
        iR = -1
    print xsl[iR]/8.
    # total particles in each linear x bin
    norme = np.copy(xsl)
    normp = np.copy(xsl)

    for k in range(len(norme)):
        norme[k]=sum(spece[:,k])
        normp[k]=sum(specp[:,k])

    # energy distribution, f(E)=(dn/dE)/N
    fe=np.empty(len(gamma))
    fp=np.empty(len(gamma))


    for k in range(len(fe)):
        fe[k]=sum(spece[k][eL:eR])/(sum(norme[eL:eR])*dgamma[k])
        #fe[k]=sum(spece[k][eL:eR])/(dgamma[k])
        
        fp[k]=sum(specp[k][iL:iR])/(sum(normp[iL:iR])*dgamma[k])
        #fp[k]=sum(specp[k][iL:iR])/(dgamma[k])

    fe[fe<=0]=1E-100
    fp[fp<=0]=1E-100

    # MASK OUT THE ZEROS

    #  NOTE: gamma ---> gamma-1 ***
    edist = gamma*fe
    pdist = gamma*fp

    masked_edist = np.ma.masked_where(edist < 1E-20, edist)
    masked_pdist = np.ma.masked_where(pdist < 1E-20, pdist)
        
    momentum=np.sqrt((gamma+1)**2-1.)
    femom=fe/(4*np.pi*momentum)/(gamma+1)
    momedist=femom*momentum**4
    fpmom=fp/(4*np.pi*momentum)/(gamma+1)
    mompdist=fpmom*momentum**4
            
    masked_momedist = np.ma.masked_where(edist < 1E-40, momedist)
    masked_mompdist = np.ma.masked_where(pdist < 1E-40, mompdist)
            
    if shock_loc > i_right_loc:
        sim.ex.append(gamma*DataDict['me'][0]/DataDict['mi'][0]/(DataDict['gamma0']-1))
        sim.ey.append(masked_edist/np.sqrt(sim.LF-1))
        sim.ix.append(gamma/(DataDict['gamma0']-1))
        sim.iy.append(masked_pdist/np.sqrt(sim.LF-1))
    else:
        sim.ex.append([])
        sim.ey.append([])
        sim.ix.append([])
        sim.iy.append([])

        """
       if i == sim.files[-1]:
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
            
        femommax = momentum**4*aconst*(gamma+1.0)*np.sqrt((gamma+1.0)**2-1)
        femommax *= np.exp(-gamma/delgame0)/(4*np.pi*momentum)/(gamma+1.0)
        sim.etempy =  femommax/sim.mass_ratio/np.sqrt(sim.LF**2-1)
        sim.etempx = momentum*DataDict['me'][0]/DataDict['mi'][0]/np.sqrt(DataDict['gamma0']**2-1)
        
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
            
        fimommax = momentum**4*aconst*(gamma+1.0)*np.sqrt((gamma+1.0)**2-1)
        fimommax *= np.exp(-gamma/delgami0)/(4*np.pi*momentum)/(gamma+1.0)
        sim.itempy =  fimommax/np.sqrt(sim.LF**2-1)
        sim.itempx = momentum/np.sqrt(DataDict['gamma0']**2-1)
        """
# Make the figure
fig = plt.figure(figsize = (7,7))
gsArgs = {'left':0.18, 'right':0.95, 'top':.95, 'bottom':0.1, 'wspace':0.2, 'hspace':0.0}
MainGs = gridspec.GridSpec(100,100)
MainGs.update(**gsArgs)

i_ax = fig.add_subplot(MainGs[2:49,:80])
#e_ax = fig.add_subplot(MainGs[51:98,:80], sharex=i_ax)
e_ax = fig.add_subplot(MainGs[51:98,:80])
cax =  fig.add_subplot(MainGs[2:98,82:87])

for i in range(len(sim.ex)):
    ## Calculate the color for this time-step
    t_step_c = cmap((sim.t[i]/sim.t[-1])**1.0)[0]

    #Plot the ion and electron spectra
    if i != len(sim.ex)-1:
        i_ax.loglog(sim.ix[i], sim.iy[i], ls = sim.ls, linewidth = 1.0, color = t_step_c, alpha = 0.75)
        e_ax.loglog(sim.ex[i], sim.ey[i], ls = sim.ls, linewidth = 1.0, color = t_step_c, alpha = 0.75)
    else:
        sim.iplot, = i_ax.loglog(sim.ix[i], sim.iy[i], ls = sim.ls, linewidth = 3.5, color = t_step_c, label = 'Ions')
        sim.eplot, = e_ax.loglog(sim.ex[i], sim.ey[i], ls = sim.ls, linewidth = 3.5, color = t_step_c, label = 'Electrons')
        

#e_ax.loglog(sim.etempx, sim.etempy, ls = '--', linewidth = 2, color = '#FF2D00')
#i_ax.loglog(sim.itempx, sim.itempy, ls = '--', linewidth = 2, color = '#FF2D00')
#i_ax.text(.3,.0003,r'$k T_p\approx0.32 m_i c^2$',size = 0.75*ax_label_size, color = '#FF2D00')
#e_ax.text(.3,.0003,r'$T_e\approx 80\%T_p$',size = 0.75*ax_label_size, color = '#FF2D00')

ax_list = [i_ax, e_ax]

for ax in ax_list:
    # MAKE THE LEGEND OMITING THE 4RD SIM
    leg = ax.legend(handlelength=0, handletextpad=0, fontsize = ax_label_size, loc="lower left")
    leg.get_frame().set_linewidth(0)
    for item in leg.legendHandles:
        item.set_visible(False)
    #for color,text in zip(lcolor,leg.get_texts()):
    #    text.set_color(color)

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
    

    ax.set_ylabel(r'$E {\rm d} n/{\rm d}E \ [{\rm arb. unit}]$', size = ax_label_size)


    #ax.set_yticks([1E-6, 1E-4, 1E-2])

    #ax.set_ylim(1E-5,5)
    ax.set_ylim(1E3,3E8)
    ax.set_xlim(5E-3,40)


### Make the Colorbar
#i_ax.add_artist(l1)

#i_ax.set_title(r'$\omega_{pi}t/\gamma_0 \approx1330$', size = ax_label_size, loc = 'right')

e_ax.set_xlabel(r'${\rm KE}/(m_i(\gamma_0-1))$', size = ax_label_size)

i_ax.set_xticklabels([])
i_ax.tick_params(axis='x',
                 which = 'both', # bothe major and minor ticks
                 bottom = 'off') # turn off top ticks


gradient =  np.linspace(0, 1, 256)# A way to make the colorbar display better
gradient = np.vstack((gradient, gradient))


cb = cax.imshow(np.transpose(gradient)[::-1], cmap=cmap, aspect = 'auto')

cb.set_extent([0,1.0, 0, sim.t[-1]])
for t in sim.t:
    cax.axhline(t, color = 'w', alpha = .5)
cax.set_yticks(np.append(np.append([0],sim.t)[::2],sim.t[-1]))
cax.tick_params(axis='x',
                which = 'both', # bothe major and minor ticks
                top = 'off', # turn off top ticks
                bottom = 'off',
                labelbottom = 'off')
cax.tick_params(axis='y',          # changes apply to the y-axis
                which='both',      # both major and minor ticks are affected
                left='off',      # ticks along the bottom edge are off
                right='on',         # ticks along the top edge are off
                labelleft='off',
                pad = -3,
                labelsize = .8*ax_label_size,
                labelright='on')


labels = cax.get_yticklabels()
plt.setp(labels, rotation=-45)

import types
SHIFT = 60. # Data coordinates
for label in labels:
    label.customShiftValue = SHIFT
    label.set_y = types.MethodType( lambda self, y: mpl.text.Text.set_y(self, y+self.customShiftValue), 
                                    label, mpl.text.Text )
for tick in cax.yaxis.get_majorticklabels():
    tick.set_verticalalignment("top")
cax.set_ylabel(r'$\omega_{pi}t$', labelpad = 25, rotation = -90,size = ax_label_size*1.2)
cax.yaxis.set_label_position("right")
cax.grid(False)
#plt.savefig('gam1vtime.png', dpi=200)
#i_ax.set_title(r'$\gamma_0\beta_0=1$' + ' Downstream sim', size=20)
i_ax.set_title(r'$\gamma_0=1.5$' + ' Upstream sim', size=20)
plt.show()
