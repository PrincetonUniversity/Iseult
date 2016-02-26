# Iseult

A TKinter python 2.7 GUI for visualizing Tristan-MP plots. A work in progress.

Written by:
Patrick Crumley, patrick.crumley@gmail.com

UPDATES:
-------
The code as stands is also ready enough to send out to others.
What's left before I feel comfortable sharing:
1) Spectral plots + settings
2) Density plots + settings (Easy to do)
3) Ability to set the minimum and maximum x and y values
4) Ability to take some measurements, Maxwellian fits and eps_e/b



Dependencies:
-------------

Python packages required: matplotlib 1.4, python 2.7, h5py

This should work on Windows, MacOS & Linux.

To get this working on orbital type the following:
```bash
$ module load anaconda/2.4.0
$ source activate root
```

To run:

```bash
$ cd /path/to/Iseult/
$ chmod +x ./iseult.py
$ export PATH="$PATH:/path/to/Iseult"
$ iseult.py
```

If you are in the path to your tristan output files or the parent directory of
the output directory, Iseult automatically loads your Tristan-MP data.
Otherwise you can open it through the file menu or by pressing crtl + O.

Enjoy!

| Implemented: |
| ------------ |
| Time stepping |
| Movies |
| Basic plotting |
| ability to modularly change plots. |
| plot control panel to edit things about indv. plots (not finished for all plots)|
| shock-finding |


| Left to Implement:|
| ------------------ |
| Density plot type |
| Spectrum plot type |
| Ability to take measurements |
| gifs |
| figure saving |
| ability to save iseult settings as a json object|
| Longer term goals |

Resources:
----------
| Useful links |
| ----------------------- |
| http://python.org |
| http://effbot.org/tkinterbook/ |
| http://matplotlib.org |
| http://h5py.org |
