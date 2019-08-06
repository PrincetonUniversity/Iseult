import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import h5py
import new_cmaps
import os, sys
import re
import ConfigParser
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FuncFormatter
from scipy.signal import savgol_filter
from scipy.stats import linregress
from scipy.integrate import simps, cumtrapz

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams['lines.linewidth'] = .75
mpl.rcParams['xtick.labelsize']=15
mpl.rcParams['ytick.labelsize']=15
mpl.rcParams['axes.labelsize']=18
mpl.rcParams['figure.subplot.left']= 0.15  # the left side of the subplots of the figure
mpl.rcParams['figure.subplot.bottom']= 0.15    # the bottom of the subplots of the figure
mpl.rcParams['figure.subplot.top'] =  0.9    # the top of the subplots of the figure

#My own color maps
colors =  ['#D90000','#04756F','#FF8C00', '#988ED5', '#2E0927', '#0971B2']

e_color = '#0971B2'
i_color = '#D90000'

gradient =  np.linspace(0, 1, 256)# A way to make the colorbar display better
gradient = np.vstack((gradient, gradient))


ax_label_size = 20
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
skip_size = 3
sim = simulation()
sim.dr = os.path.join(os.curdir,'../StampedeGamBeta0_33/')
#sim.dr = os.path.join(os.curdir,'../StampedeGamBeta3/')
#sim.dr = os.path.join(os.curdir,'../StampedeGamBeta1/savedRun/')
#sim.dr = os.path.join(os.curdir,'../Gam3LowComp/output')
#sim.dr = os.path.join(os.curdir,'../Gam1.5LowCompSmallerdelY/output')

relative = False

sim.e_region = (300,700)

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
sim.files = (np.arange(len(sim.PathDict['Param'])//skip_size)+1)*skip_size -1# Every 4th timestep, but since it wasn't run until 60 we have to cheat the last time
sim.files[-1] = -1
#sim.files = sim.files[2:]
print sim.files
ToLoad = {'Flds': [], 'Prtl': [], 'Param': [], 'Spect': []}
tmpList = ['c_omp', 'istep', 'dens','ppc0','gamma', 'xsl','spece','specp', 'time', 'mi', 'me', 'gamma0', 'c', 'sigma', 'btheta', 'ue', 've', 'we', 'qi', 'ui', 'vi', 'wi', 'stride', 'bx', 'by', 'bz', 'ex', 'ey', 'ez','xe', 'xi']
sim.KE_e = []
sim.KE_i = []
sim.B = []
sim.E = []
sim.t  = []

for i in sim.files:
    print(i)
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


        
    # All the data for this time step is now stored in the datadict
    # Calculate the current shock location:
    sim.c_omp = DataDict['c_omp'][0]
    sim.istep = DataDict['istep'][0]

    sim.t.append(DataDict['time'][0]/np.sqrt(sim.mass_ratio)) 
    GoodLecs = DataDict['xe']<=sim.e_region[1]*DataDict['c_omp'][0]
    GoodLecs *= DataDict['xe']>sim.e_region[0]*DataDict['c_omp'][0]

    TotalElectronKE = DataDict['ue']*DataDict['ue']
    TotalElectronKE += DataDict['ve']*DataDict['ve']
    TotalElectronKE += DataDict['we']*DataDict['we']+1
    TotalElectronKE = np.sum(np.sqrt(TotalElectronKE[GoodLecs])-1)
    TotalElectronKE *= DataDict['stride'][0] # multiply by the stride.
    TotalElectronKE *= np.abs(DataDict['qi'][0])*DataDict['c'][0]**2 # * m_e c^2, mass of particle is its charge, qe/me=1
    print(TotalElectronKE)
    sim.KE_e.append(TotalElectronKE)

    GoodIons = DataDict['xi']<=sim.e_region[1]*DataDict['c_omp'][0]
    GoodIons *= DataDict['xi']>sim.e_region[0]*DataDict['c_omp'][0]
    TotalIonKE = DataDict['ui']*DataDict['ui']
    TotalIonKE += DataDict['vi']*DataDict['vi']
    TotalIonKE += DataDict['wi']*DataDict['wi']+1
    TotalIonKE = np.sqrt(TotalIonKE)-1
    TotalIonKE = np.sum(TotalIonKE[GoodIons])
    TotalIonKE *= DataDict['stride'][0] # multiply by the stride
    TotalIonKE *= DataDict['mi'][0]/DataDict['me'][0]*np.abs(DataDict['qi'][0])*DataDict['c'][0]**2 #mass of particle is its charge, qe/me=1

    sim.KE_i.append(TotalIonKE)
    
    fL = int(sim.e_region[0]*DataDict['c_omp'][0]/DataDict['istep'][0])
    fR = int(sim.e_region[1]*DataDict['c_omp'][0]/DataDict['istep'][0])
    print(fL, fR)
    BzEnergy = DataDict['bz'][:,:,fL:fR]*DataDict['bz'][:,:,fL:fR]
    TotalBEnergy = DataDict['bx'][:,:,fL:fR]*DataDict['bx'][:,:,fL:fR]
    TotalBEnergy += DataDict['by'][:,:,fL:fR]*DataDict['by'][:,:,fL:fR]
    TotalBEnergy += BzEnergy
    # sum over the array and then multiply by the number of points len(x)*len(y)

    BzEnergy = np.sum(BzEnergy)
    BzEnergy *= DataDict['istep'][0]**2*.5
    TotalBEnergy = np.sum(TotalBEnergy) #*DataDict['bx'][0,:,:].shape[1]*DataDict['bx'][0,:,:].shape[0]
    TotalBEnergy *= DataDict['istep'][0]**2*.5

    TotalEEnergy = DataDict['ex'][:,:,fL:fR]*DataDict['ex'][:,:,fL:fR]
    TotalEEnergy += DataDict['ey'][:,:,fL:fR]*DataDict['ey'][:,:,fL:fR]
    TotalEEnergy += DataDict['ez'][:,:,fL:fR]*DataDict['ez'][:,:,fL:fR]

    # sum over the array and then divide by the number of points len(x)*len(y)
    TotalEEnergy = np.sum(TotalEEnergy) #*DataDict['ex'][0,:,:].shape[1]*DataDict['ex'][0,:,:].shape[0]
    TotalEEnergy *= DataDict['istep'][0]**2*.5

    sim.B.append(TotalBEnergy)
    sim.E.append(TotalEEnergy)

# Make the figure
fig = plt.figure()

plt.plot(sim.t, np.array(sim.KE_i)+np.array(sim.KE_e), '-', c=colors[2], label = "prtl total")
plt.plot(sim.t, np.array(sim.KE_i), 'x-', c=i_color, label = "ions")
plt.plot(sim.t, np.array(sim.KE_e), '.-', c=e_color, label = "electrons")
plt.plot(sim.t, np.array(sim.B)+np.array(sim.E), 'd-', c=colors[1], label="fields")
plt.xlabel(r'$\omega_{pi} t$')
plt.ylabel('Energy')
plt.ylim(0,None)
leg = plt.legend(handlelength=0, handletextpad=0, fontsize = 20, fancybox=False, edgecolor = 'w', markerscale=0)
for item in leg.legendHandles:
    item.set_visible(False)
for color,text in zip([colors[2], i_color, e_color, colors[1]],leg.get_texts()):
    #for color,text in zip([colors[2], i_color, e_color],leg.get_texts()):
    text.set_color(color)

plt.title(r'$\gamma_0\beta_0=0.33$' + ' Upstream sim', size=20)
#plt.title(r'$\gamma_0\beta_0=1$' + ' Downstream sim', size=20)
plt.show()
