# Iseult

[![Princeton RSE Badge](https://img.shields.io/badge/Princeton_RSE-2024-%23F58025.svg?logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIGlkPSJMYXllcl8xIiB2aWV3Qm94PSIwIDAgMzkyLjkgNTAwIj48ZGVmcz48c3R5bGU+LmNscy0xLC5jbHMtM3tzdHJva2Utd2lkdGg6MH0uY2xzLTN7ZmlsbDojZmZmfTwvc3R5bGU+PC9kZWZzPjxwYXRoIGQ9Ik0zODIgMjExYzE3LTgwIDktMTM4IDktMTM4QzI1OCAxMDkgMTk2IDAgMTk2IDBTMTM1IDEwOSAyIDczYzAgMC05IDU4IDkgMTM4aDM3MVoiIGNsYXNzPSJjbHMtMSIvPjxwYXRoIGQ9Ik0xMSAyMTFhNDgxIDQ4MSAwIDAgMCAxODUgMjg5IDU0NSA1NDUgMCAwIDAgMTM0LTE1M2MyOC00OCA0My05NSA1Mi0xMzZIMTFaIiBzdHlsZT0iZmlsbDojZjU4MDI1O3N0cm9rZS13aWR0aDowIi8+PHBhdGggZD0iTTE5NiAyMTggNjMgMzQ3YzEyIDIxIDI2IDQzIDQzIDY0bDkwLTg5IDkxIDg5YzE3LTIxIDMxLTQzIDQzLTY0TDE5NiAyMThaIiBjbGFzcz0iY2xzLTEiLz48cGF0aCBkPSJtMTc2IDg2LTE4IDNjLTkgMi0yMCAzLTMwLTF2ODZjNCAzIDExIDMgMTcgM2wzMC00YzctMSAxNC0xIDIwIDNWOTFjLTUtNi0xMi02LTE5LTVaTTIzNSA4OWwtMTktM2MtNi0xLTEzLTEtMTggNXY4NWM2LTQgMTMtNCAxOS0zbDMxIDRjNSAwIDEzIDAgMTctM1Y4OGMtMTAgNC0yMSAzLTMwIDFaIiBjbGFzcz0iY2xzLTMiLz48cGF0aCBkPSJtMTk2IDE5MSAxMC0xdi02YTM0NyAzNDcgMCAwIDAgNzAgMFY5OGwtOC0xdjc4cy0xIDMtNSA0Yy05IDMtMTUgMi0yMyAxbC0xNy0zYy05LTEtMTktNC0yNSA1di0xaC0zdjFjLTItMy00LTQtNi01LTYtMi0xMyAwLTE5IDBsLTE4IDNjLTcgMS0xMyAyLTIyLTEtNC0xLTUtNC02LTRWOTdsLTggMXY4NmM0IDEgMTcgMyAzMCAzIDIwIDAgMzktMyA0MS0zdjZsOSAxWk0xMDUgMTA5aDh2MTZoLTh6TTEwNSAxNTFoOHYxNWgtOHpNMjg4IDEyNWgtOHYtMTZoOHpNMjg4IDE2NmgtOHYtMTVoOHoiIGNsYXNzPSJjbHMtMyIvPjxwYXRoIGQ9Im0xNTIgMTA0LTMgMS01IDEyLTEgMXYtMWwtMy01LTEtNC0yLTMtMi0xdi0xaDh2MWwtMyAxaDFsMyA5IDMtNyAxLTItMi0xLTEtMWg3djFNMTY3IDExN2gtMTR2LTFoMmwxLTJ2LTlsLTMtMXYtMWgxM3Y0bC0zLTNoLTV2NWgybDItMXYtMWgxdjZoLTFsLTEtM2gtM3Y0bDEgMXYxaDVsMy0zdjRNMTg0IDEwN2gtMWwtMS0yLTItMWgtMnYxMGMwIDIgMSAyIDMgMnYxaC04di0xbDMtMXYtOWwtMS0yaC0ybC0zIDN2LTRoMTR2NE0xNTEgMTMwaC0xbC0xLTItMi0yaC0ydjExbDMgMmgtOGMyIDAgMi0xIDMtMnYtMTFoLTNsLTMgM3YxLTVoMTR2NU0xNjcgMTQwbC03LTFoLTcgMWwyLTF2LTEwYzAtMi0xLTItMi0ydi0xaDEybDEgMXYzaC0xbC0yLTItMS0xaC00djVoMmwxLTEgMS0xaDFsLTEgMyAxIDNoLTFsLTItMi0xLTFoLTF2NmwxIDFoM2wyLTFzMi0xIDItM2gxdjRsLTEgMU0xNzcgMTQwbC0zLTFoLTJ2LTRoMWwxIDIgMyAyIDItMSAxLTJ2LTFsLTItMS0zLTEtMi0yLTEtMmMwLTIgMi00IDQtNGw0IDF2LTFoMXY1aC0xbC0xLTMtMy0xLTIgMiAyIDIgMiAxIDIgMWMyIDEgMiAyIDIgNHMtMiA0LTUgNE0xNTcgMTYyaC0xNHYtMWgybDEtM3YtOGwtMy0yaDEzdjNoLTFsLTItMmgtNXY1aDJsMi0xdi0yaDF2N2wtMS0xLTEtMmgtM3Y1bDEgMWg1bDMtM3Y0TTE3NyAxNDhjLTEgMC0yIDEtMiAzdjExaC0xbC0xMC0xMXY4bDEgMmgydjFoLTd2LTFsMy0xdi05bC0xLTF2LTFsLTItMWg1bDcgOCAyIDJ2LTdsLTEtMi0yLTFoNy0xTTIyNCAxMDRjLTEgMC0yIDAtMiAzdjEwaC0xbC0xMC0xMXY4bDEgMmgydjFoLTd2LTFsMy0xdi04bC0xLTJ2LTFoLTJ2LTFoNWw3IDggMiAydi03bC0xLTJoLTJ2LTFoN3YxTTIzNyAxMDVsLTMtMi0yIDItMiA1IDEgNSAzIDJjMiAwIDMtMSAzLTJsMS01LTEtNW0yIDEwYy0xIDItMyAzLTUgM2wtNS0zLTItNSAyLTVjMi0yIDQtMiA1LTJsNCAxYzIgMSA0IDMgNCA2IDAgMi0xIDQtMyA1TTI2MCAxMDRsLTMgMS01IDEyLTEgMXYtMWwtMy01LTEtNC0yLTMtMi0xdi0xaDh2MWwtMiAxIDMgOSAzLTcgMS0yLTItMS0xLTFoN3YxTTIyMCAxMzBoLTFsLTEtMi0yLTJoLTJ2MTFsMyAyaC04YzIgMCAyLTEgMy0ydi05bC0xLTJoLTJsLTMgM3YxLTVoMTR2NU0yMjggMTI4bC0xIDItMSAzaDRsLTItNW05IDEyaC0xbC0zLTEtMyAxdi0xbDItMS0xLTItMS0yaC00bC0yIDN2MWwyIDFoLTYgMWwyLTEgMS00IDMtNSAxLTJ2LTJoMWwxIDEgMiA2IDIgMyAxIDMgMiAxdjFNMjU5IDE0MGwtNi0xaC0xIDFsMS0ydi0xMGwtMyA2LTEgMy0xIDNoLTFsLTEtMi01LTExdjEwYzAgMiAxIDMgMyAzdjFoLTFsLTMtMWgtM2wyLTEgMS0ydi04bC0xLTJoLTJ2LTFoNmwxIDEgNCA5IDEtMiAzLTYgMS0yaDV2MWwtMiAxdjhsMSA0aDF2MU0yMjAgMTUyaC0xbC0xLTItMy0xaC0xdjEwYzAgMiAxIDIgMyAydjFoLTh2LTFsMi0xdi0xMWgtMmwtMyAzaC0xdi00aDE1djRNMjM3IDE0OGwtMiAyLTYgMTJoLTFsLTItNi0yLTQtMS0yLTItMmg3bC0yIDF2MWw0IDkgMy03di0ybC0yLTJoNk0yNTkgMTYyaC03di0xaDFsMS0ydi0xMGwtMyA3LTEgMi0xIDNoLTFjLTEgMCAwIDAgMCAwbC0xLTItNS0xMHY5YzAgMiAxIDMgMyAzdjFoLTZ2LTFsMi0xdi0xMWwtMi0xaDZsNCAxMCAxLTIgMy03IDEtMWg1bC0yIDF2OWwxIDNoMXYxIiBjbGFzcz0iY2xzLTEiLz48L3N2Zz4=&labelColor=%235A575B)](https://researchcomputing.princeton.edu/services/research-software-engineering)

A TKinter based python GUI for visualizing Tristan-MP plots. A work in progress.

![Iseult Set-up](https://raw.githubusercontent.com/pcrumley/Iseult/gh-pages/images/IseultPanels.png)
An example visualization of a Tristan-MP simulation.

Written by:

- Patrick Crumley, patrick.crumley@gmail.com, based on Jaehong's Tristan analysis IDL script.
- Robert Caddy, rcaddy@princeton.edu

## UPDATES:

July 2024: Ported to Python 3.11, see environment.yml for details on all dependencies.
Added support for Tristan v2 data.

May 8th 2019: Ported code to python 3.7.x & Matplotlib 3.0.x There may be a few bugs
here and there, but I think it is working.

The code is now it's beta phase. See the implemented column for what has
already been implemented.

## Dependencies:


Python packages required: See the environment.yml file

To use the movie saving feature: ffmpeg & xterm.

Iseult should work on Windows, MacOS & Linux.

## Setup & Running

There are two different methods to setup what you need for Iseult to run. The recommended method is to setup a python environment with conda, pyenv, etc and a conda `environment.yml` file is provided. The other method for those with access to the Stellar cluster at Princeton is to simply load the proper Anaconda module.

### Recommended Method

```bash
# Recommended setup
# If on Stellar
$ module load anaconda3/2024.6 # or the latest version
# If not on Stellar then install/load Anaconda
$ cd /path/to/Iseult/
$ conda env create
$ conda activate iseult
$ chmod +x ./iseult.py

# After initial setup
# make sure the proper environment is activated
$ conda activate iseult
# Then run iseult. This method opens a GUI to allow you to search for the data you want to open
$ cd /path/to/Iseult/
$ ./iseult.py
# Or if you want to open data directly
$ cd /directory/that/contains/data # usually the `output` directory
$ /path/to/iseult.py # you can append an ampersand (&) to open iseult in the background
```

### Stellar Only Method

The previous method is still recommended on Stellar, this secondary method is simplier but potentially less reliable if the software stack on Stellar is changed.

```bash
# Setup
$ module load anaconda3/2023.3
$ cd /path/to/Iseult/
$ chmod +x ./iseult.py

# After initial setup
$ module load anaconda3/2023.3
# Then run iseult. This method opens a GUI to allow you to search for the data you want to open
$ cd /path/to/Iseult/
$ ./iseult.py
# Or if you want to open data directly
$ cd /directory/that/contains/data # usually the `output` directory
$ /path/to/iseult.py # you can append an ampersand (&) to open iseult in the background
```

When Iseult is started, it checks to see if Tristan-MP data is located at the
current path, if it isn't Iseult prompts you to select the directory of where
your Tristan-MP data is saved. To edit/change any of the plots, just right click
on the subplot directly. You can also change the number of columns, the
colormap, and other general settings by clicking the settings button. The
measure button allows you to take measurements like T_i, T_e, measure spectra,
take 1-D FFTs, etc. The matplotlib interactive toolbar is beneath the playbar,
it allows you to save the figure, use your mouse to zoom around, etc.

If you get a set-up you like, you can go to file menu and choose Save Iseult
State. It will give you an option to name the 'view.' To replace the default
state of Iseult, name the view Default. The views are saved as .cfg files in
.iseult_configs folder. You must restart Iseult to see the saved config in the
preset views menu.

Enjoy!


| Implemented: |
| ------------ |
| Time stepping |
| Movie (without recording) |
| Basic plotting |
| ability to modularly change plots. |
| plot control panels to edit things about indv. plots |
| shock-finding |
| figure saving |
| Ability to take measurements |
| ability to save Iseult settings in a config file|
| zooming |


| Left to Implement:|
| ------------------ |
| gifs/movies |
| Longer term goals (???)|

Resources:
----------
| Useful links |
| ----------------------- |
| http://python.org |
| http://effbot.org/tkinterbook/ |
| http://matplotlib.org |
| http://matplotlib.org/users/navigation_toolbar.html |
| http://h5py.org |
