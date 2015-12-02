# Iseult

A python GUI for visualizing Tristan-MP plots. A work in progress.

Written by:
Patrick Crumley, patrick.crumley@gmail.com

UPDATE:
-------
I have decided that I should re-write the program in Tkinter instead of wxpython
for two reasons. One, I need to significantly change the wxpython code anyway so
that there is only one figure with 6 subplots (instead of 6 separate figures).
This means it will be much much easier to change the number of subplots and to
save the figure. Two, the code as written will not run on orbital, and the
Tkinter code is guaranteed to run as it is part of the python stack. Wish me
luck!

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
$ chmod +x ./iseult
$ ./iseult
```

If you are in the path to your tristan output files or the parent directory of
the output directory, Iseult automatically loads your Tristan-MP data.
Otherwise you can open it through the file menu.

Enjoy!

| Implemented: |
| ------------ |
| Time stepping |
| Basic plotting |
| ability to modularly change plots. |
| plot control panel to edit things about indv. plots|


| Left to Implement:|
| ------------------ |
| All the plot types |
| shock-finding |
| gifs |
| figure saving |
| Much much more... |

Resources:
----------
| Useful links |
| ----------------------- |
| http://python.org
| http://www.wxpython.org |
| http://matplotlib.org |
| http://h5py.org |
