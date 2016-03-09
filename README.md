# Iseult

A TKinter python 2.7 GUI for visualizing Tristan-MP plots. A work in progress.

Written by:

Patrick Crumley, patrick.crumley@gmail.com, based on jaehong's Tristan analysis
IDL script.

UPDATES:
-------
The code as stands is almost ready enough to send out to others.
What's left before I feel comfortable sharing:
1) Spectral plots (DONE!) + settings (not yet!)
2) Density plots + settings (DONE!)   
3) Ability to set the minimum and maximum x and y values
4) Ability to take some measurements, Maxwellian fits and eps_e/b
5) Save plots (Done!).


Dependencies:
-------------

Python packages required: matplotlib 1.4, python 2.7, h5py

This should work on Windows, MacOS & Linux.

To run iseult on orbital type the following:
```bash
$ module load anaconda/2.4.0
$ source activate root
$ cd /path/to/Iseult/
$ chmod +x ./iseult.py
$ iseult.py
```

When iseult is started, it prompts you to select the directory of where you
Tristan-MP is saved. To edit/change any of the plots, just click on the subplot
directly. You can also change the number of columns, the colormap, and spacing
of the plots by clicking the settings button. The measure button allows you to
take measurements like T_i, T_e, etc.

Enjoy!


| Implemented: |
| ------------ |
| Time stepping |
| Movie (without recording) |
| Basic plotting |
| ability to modularly change plots. |
| plot control panel to edit things about indv. plots |
| shock-finding |
| figure saving |
| zooming |


| Left to Implement:|
| ------------------ |
| Ability to take measurements |
| gifs |
| ability to save iseult settings as a json object|
| Longer term goals (???)|

Resources:
----------
| Useful links |
| ----------------------- |
| http://python.org |
| http://effbot.org/tkinterbook/ |
| http://matplotlib.org |
| http://h5py.org |
