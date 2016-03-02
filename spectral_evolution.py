import matplotlib as mpl
from matplotlib import pyplot
import numpy as np
import new_cmaps
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import h5py
import os, sys
import re
import Tkinter, tkFileDialog

def pathOK(dirname):
    """ Test to see if the current path contains tristan files
    using regular expressions, then generate the lists of files
    to iterate over"""

    f_re = re.compile('flds.tot.*')
    is_okay = len(filter(f_re.match, os.listdir(dirname)))>0

    return is_okay

def findDir(parent, dlgstr = 'Choose the directory of the output files.'):
    """Look for /ouput folder, where the simulation results are
    stored. If output files are already in the path, they are
    automatically loaded"""
    # defining options for opening a directory


    dir_opt = {}
    dir_opt['initialdir'] = os.curdir
    dir_opt['mustexist'] = True



    tmpdir = tkFileDialog.askdirectory(parent = parent, title = dlgstr, **dir_opt)
    dirlist = os.listdir(tmpdir)
    if 'output' in dirlist:
        tmpdir = os.path.join(tmpdir, 'output')
    print tmpdir
    return tmpdir


# Be sure to call this with the output directory in the path
# create a bunch of regular expressions used to search for files
f_re = re.compile('flds.tot.*')
prtl_re = re.compile('prtl.tot.*')
s_re = re.compile('spect.*')
param_re = re.compile('param.*')
re_list = [f_re, prtl_re, s_re, param_re]
#root = Tkinter.Tk()
dirname = os.path.join(os.curdir,'output')
#root.destroy()

PathDict = {'Flds': [], 'Prtl': [], 'Param': [], 'Spect': []}

#fill all the paths
PathDict['Flds']= filter(f_re.match, os.listdir(dirname))
PathDict['Flds'].sort()

PathDict['Prtl']= filter(prtl_re.match, os.listdir(dirname))
PathDict['Prtl'].sort()

PathDict['Spect']= filter(s_re.match, os.listdir(dirname))
PathDict['Spect'].sort()

PathDict['Param']= filter(param_re.match, os.listdir(dirname))
PathDict['Param'].sort()


# A dictionary that allows use to see in what HDF5 file each key is stored.
# i.e. {'ui': 'Prtl', 'ue': 'Flds', etc...},

H5KeyDict = {}
for pkey in PathDict.keys():
    with h5py.File(os.path.join(dirname,PathDict[pkey][0]), 'r') as f:
    # Because dens is in both spect* files and flds* files,
        for h5key in f.keys():
            if h5key == 'dens' and pkey == 'Spect':
                H5KeyDict['spect_dens'] = pkey
            else:
                H5KeyDict[h5key] = pkey

# Get the shock location:
# First load the first field file to find the initial size of the
# box in the x direction

#        print os.path.join(self.dirname,self.PathDict['Flds'][0])

with h5py.File(os.path.join(dirname,PathDict['Flds'][0]), 'r') as f:
    nxf0 = f['by'][:].shape[1]
with h5py.File(os.path.join(dirname,PathDict['Param'][-1]), 'r') as f:
    initial_time = f['time'][0]

# Load the final time step to find the shock's location at the end.
with h5py.File(os.path.join(dirname,PathDict['Flds'][-1]), 'r') as f:
    dens_arr =np.copy(f['dens'][0,:,:])

with h5py.File(os.path.join(dirname,PathDict['Param'][-1]), 'r') as f:
    # I use this file to get the final time, the istep, interval, and c_omp
    final_time = f['time'][0]
    istep = f['istep'][0]
    interval = f['interval'][0]
    c_omp = f['c_omp'][0]

# Find out where the shock is at the last time step.
jstart = min(10*c_omp/istep, nxf0)
# build the final x_axis of the plot

xaxis_final = np.arange(dens_arr.shape[1])/c_omp*istep
# Find the shock by seeing where the density is 1/2 of it's
# max value. First average the density in the y_direction


dens_half_max = max(dens_arr[dens_arr.shape[0]/2,jstart:])*.5
ishock_final = np.where(dens_arr[dens_arr.shape[0]/2,jstart:]>=dens_half_max)[0][-1]
xshock_final = xaxis_final[ishock_final]

print xshock_final
shock_speed = xshock_final/final_time

ion_e_region = (-1E5,0) # relative to the shock
e_e_region = (-1E5,0)



# Create figure where the lines will be plotted
# First make the gridspec
gs = gridspec.GridSpec(100,100)
gsArgs = {'left':0.15, 'right':0.9, 'top':.95, 'bottom':0.06, 'wspace':0.2, 'hspace':0.2}
gs.update(**gsArgs)

fig = pyplot.figure()
ax1 = fig.add_subplot(gs[:2,:])
ion_axes = fig.add_subplot(gs[10:45,:])
e_axes = fig.add_subplot(gs[60:95,:])

cmap = new_cmaps.cmaps['viridis']
#for i in range(len(PathDict['Spect'])):
for i in range(len(PathDict['Spect'])):
    t_step_color = cmap(i*1.0/(len(PathDict['Spect'])-1))

    ToLoad = {'Flds': [], 'Prtl': [], 'Param': [], 'Spect': []}
    tmpList = ['c_omp', 'istep', 'gamma', 'xsl','spece','specp', 'time']
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
        with h5py.File(os.path.join(dirname,PathDict[pkey][i]), 'r') as f:
            for elm in tmplist:
                # Load all the keys
                if elm == 'spect_dens':
                    DataDict[elm] = f['dens'][:]
                else:
                    DataDict[elm] = f[elm][:]

    # All the data for this time step is now stored in the datadict
    # Calculate the current shock location:
    shock_loc = DataDict['time'][0]*shock_speed

    c_omp = DataDict['c_omp'][0]
    istep = DataDict['istep'][0]
    xsl = DataDict['xsl']/c_omp
    gamma = DataDict['gamma']

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
    e_left_loc = shock_loc + e_e_region[0]
    e_right_loc = shock_loc + e_e_region[1]

    eL = xsl.searchsorted(e_left_loc)
    eR = xsl.searchsorted(e_right_loc, side='right')
    if eL == eR:
        eR += 1
    i_left_loc = shock_loc + ion_e_region[0]
    i_right_loc = shock_loc + ion_e_region[1]

    iL = xsl.searchsorted(i_left_loc)
    iR = xsl.searchsorted(i_right_loc, side='right')
    print iL, eL, iR, eR
    if iL == iR:
        iR += 1

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
        fp[k]=sum(specp[k][iL:iR])/(sum(normp[iL:iR])*dgamma[k])

    #  NOTE: gamma ---> gamma-1 ***
    edist = gamma*fe
    pdist = gamma*fp

    momentum=np.sqrt((gamma+1)**2-1.)
    femom=fe/(4*np.pi*momentum)/(gamma+1)
    momedist=femom*momentum**4
    fpmom=fp/(4*np.pi*momentum)/(gamma+1)
    mompdist=fpmom*momentum**4


    # Create a gridspec to handle spacing better

    e_axes.plot(momentum, momedist, color = t_step_color, linewidth = 0.5)


    ion_axes.plot(momentum, mompdist, color = t_step_color, linewidth = 0.5)


# Make the colorbar up top:
norm = mcolors.Normalize(vmin=initial_time, vmax=final_time)

# ColorbarBase derives from ScalarMappable and puts a colorbar
# in a specified axes, so it has everything needed for a
# standalone colorbar.  There are many more kwargs, but the
# following gives a basic continuous colorbar with ticks
# and labels.
cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                norm=norm,
                                orientation='horizontal')
cb1.ax.tick_params(labelsize = 8)
#cb1.set_label(r'$\omega_p$',labelpad =25, size =18 )

e_axes.set_xscale("log")
e_axes.set_yscale("log")
e_axes.set_xlim(-0.05,500)
#e_axes.set_ylim(1E-6,1)
e_axes.tick_params(labelsize = 10)
e_axes.set_xlabel(r'$p(m_e c)$', color = 'black')
e_axes.set_ylabel(r'$p^4f_e (p)$')

ion_axes.set_xscale("log")
ion_axes.set_yscale("log")
ion_axes.set_xlim(-0.05,500)
#ion_axes.set_ylim(1E-6,1)
ion_axes.tick_params(labelsize = 10)
ion_axes.set_xlabel(r'$p(m_i c)$', color = 'black')
ion_axes.set_ylabel(r'$p^4 f_p (p)$')
pyplot.show()
