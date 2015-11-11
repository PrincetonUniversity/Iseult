# Iseult

A python GUI for visualizing Tristan-MP plots. A work in progress.

Written by:
Patrick Crumley, patrick.crumley@gmail.com

Dependencies:
-------------

Python packages required: matplotlib, wxpython, h5py

This should work on Windows, MacOS & Linux.

To get this working on orbital type the following:
```bash
$ module load python/2.7
$ module load hdf5/gcc/1.8.12
$ module load h5py27/2.2.1/hdf5-1.8.12
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
