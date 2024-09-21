#!/usr/bin/python3
#--------  Code Dependencies   ----------
#\__________General Utilities____________/
import pprint
import numpy as np
#\__________matplotlib functions__________/
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import mlab
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
import matplotlib.dates as mdates
#\_____________________________________/
#\__________Scripts List__________/
"""
Script         Description
pltTrSp        Plot two graphs one on the topo of the other
gplot
"""
#
# -------------- pltTrSp  ----------------------------
"""
Plotting script for 1-D time series

"""
def pltTrSp(x1, y1, x2, y2,
            x1label=False, y1label=False, y1log=False, clr1 = 'k', 
            x2label=False, y2label=False, y2log=False, clr2 = 'r'
            ):
#
    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.plot(x1, y1, clr1)  #, label="Timeseries data"
    #ax1 = plt.gca()
    #ax1.set_xlim([xmin, xmax])
    #ax1.set_ylim([ymin, ymax])
    if x1label: ax1.set_xlabel(x1label)
    if y1label: ax1.set_ylabel(y1label)
    if y1log:   ax1.set_yscale('log')
    ax1.grid(True)
    #ax1.legend()
    #
    ax2.plot(x2, y2, clr2)  #, label="Amplitude spectrum"
    if x2label: ax1.set_xlabel(x2label)
    if y2label: ax1.set_ylabel(y2label)
    if y2log:   ax2.set_yscale('log')
    ax2.grid(True)
    #ax2.legend()
    plt.show()
#
# -------------- End of function   ---------------------
#
# -------------- pltTrSp  ----------------------------
"""
Trace-like plotting script
tr  -> A trace
"""
def gplot(t, tr, trlog=False, xlbl=None, ylbl=None, title=None, clr = 'k'):
#------------ Plot trace
    fig, ax = plt.subplots()
    ax.plot(t, tr.data, clr, alpha=0.7)
    if xlbl: ax.set_xlabel(xlbl)
    if ylbl: ax.set_ylabel(ylbl)
    if ylbl: ax.set_yscale('log')
    if title: ax.suptitle(title)
    ax.grid(True)
    plt.show()
#
#    plt.legend()
#    plt.show(block=False)
#
# -------------- End of function   ---------------------