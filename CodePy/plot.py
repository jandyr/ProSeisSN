#!/usr/bin/python3
#--------  Code Dependencies   ----------
#\__________General Utilities____________/
import os, sys
import pprint
import numpy as np
import math
#\__________ObsPy functions__________/
from obspy import UTCDateTime
from obspy.imaging.cm import obspy_sequential
#\__________matplotlib functions__________/
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import mlab
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
import matplotlib.dates as mdates
from matplotlib.colors import Normalize
matplotlib.use('nbagg')
plt.style.use('ggplot')
#\_____________________________________/
#\__________Scripts List__________/
"""
Script         Description
pgather        Plot gather
pltTrSp        Plot two graphs one on the topo of the other
gplot
"""
#
# -------------- Plot FK O/P  -------------------
def pltbaz(out, stime, etime, title='TTB'):
  """ 
  Plot
  """
# -------------- Warning on plots
  sys.stdout.write('\n')
  print(f'\n>> Warning: pltbaz works with VAGAROSITY.')
  sys.stdout.write('\n')
# -------------- Plot graphs
  labels = ['rel.power', 'abs.power', 'baz', 'slow']
#    xlocator = mdates.AutoDateLocator()
  fig = plt.figure()
  for i, lab in enumerate(labels):
      ax = fig.add_subplot(4, 1, i + 1)
#-- Work with velocity instead of vagarosity
#      val = 1./ out[:, i + 1] if lab == 'vel' else out[:, i + 1]
#-- Change x-coord
      dummy = out[-1, 0] - out[0, 0]
      dummy = np.linspace(stime, etime,
                          num=len(out))
      out[:, 0] = dummy.copy()
      ax.scatter(out[:, 0], out[:, i + 1], c=out[:, 1], alpha=0.6,
                 edgecolors='none', cmap=obspy_sequential)
#-- Axis limits
      ax.set_ylabel(lab)
      if lab == 'baz':
#        baz[baz < 0.0] += 360                 #backazimuth to values between 0 and 360
        ax.set_ylim(-180,180)
      else:
        ax.set_ylim(out[:, i + 1].min(), out[:, i + 1 ].max())
      ax.set_xlim(out[0, 0], out[-1, 0])

#        ax.xaxis.set_major_locator(xlocator)
#        ax.xaxis.set_major_formatter(mdates.AutoDateFormatter(xlocator))
#  fig.suptitle(title+' %s' % (
#      stime.strftime('%Y-%m-%d'), ))
#    fig.autofmt_xdate()
  ax.set_xlabel('s')
  fig.subplots_adjust(left=0.15, top=0.95, right=0.95, bottom=0.2, hspace=0.5)
#    plt.show()
#
# -------------- Polar plot
#-- Sums the relative power in gridded bins, each defined by backazimuth 
#    and slowness of the analyzed signal part. The backazimuth is counted
#    clockwise from north, the slowness limits can be set by hand.
#
  cmap = obspy_sequential
#-- Make output human readable, adjust backazimuth to values between 0 and 360
  t, rel_power, abs_power, baz, slow = out.T
  baz[baz < 0.0] += 360
# choose number of fractions in plot (desirably 360 degree/N is an integer!)
  N = 36
  N2 = 30
  abins = np.arange(N + 1) * 360. / N
  sbins = np.linspace(0, 3, N2 + 1)
#-- Sum rel power in bins given by abins and sbins
  hist, baz_edges, sl_edges = \
    np.histogram2d(baz, slow, bins=[abins, sbins], weights=rel_power)
#-- Transform to radian
  baz_edges = np.radians(baz_edges)
#-- Add polar and colorbar axes
  fig = plt.figure(figsize=(8, 8))
  cax = fig.add_axes([0.85, 0.2, 0.05, 0.5])
  ax = fig.add_axes([0.10, 0.1, 0.70, 0.7], polar=True)
  ax.set_theta_direction(-1)
  ax.set_theta_zero_location("N")
#
  dh = abs(sl_edges[1] - sl_edges[0])
  dw = abs(baz_edges[1] - baz_edges[0])

  # circle through backazimuth
  for i, row in enumerate(hist):
      bars = ax.bar((i * dw) * np.ones(N2),
                    height=dh * np.ones(N2),
                    width=dw, bottom=dh * np.arange(N2),
                    color=cmap(row / hist.max()))
#
  ax.set_xticks(np.linspace(0, 2 * np.pi, 4, endpoint=False))
  ax.set_xticklabels(['N', 'E', 'S', 'W'])

  # set slowness limits
  ax.set_ylim(0, 3)
  [i.set_color('grey') for i in ax.get_yticklabels()]
  ColorbarBase(cax, cmap=cmap,
               norm=Normalize(vmin=hist.min(), vmax=hist.max()))
  plt.show()
#
# -------------- End of function   ---------------------
#
#-------------- Plot gather Program --------------
def pgather(xg, yg, annotg, coord='polar', degrees=True):
  """
  Program to plot Gather.
   zero = True  -> Annotate 0 to first point.
   zero = False -> Annotate 1 to first point.
  """
#-- Transform to cartezians: xg-> radial distance; yg-> azimuth in radians
  if coord == 'polar': xg, yg = util.Pol_Car(xg, yg, degrees = degrees)
#-- Plots
  fig, ax = plt.subplots()
  ax.scatter(xg, yg)
  for xi, yi, text in zip(xg, yg, annotg):
    ax.annotate(text,
      xy=(xi, yi), xycoords='data',
      xytext=(2.5, 2.5), textcoords='offset points')
  plt.xlabel('X(m)')
  plt.ylabel('Y(m)')
  plt.suptitle('Gather')
#
#    plt.legend()
  plt.show()
#  plt.show(block=False)
  return ''


#
# -------------- pltTrSp  ----------------------------
"""
Plotting script for 1-D time series with two plots

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
tr0 and tr1 -> Two 1-D data sets
"""
def gplot(t, tr0, tr1=None, clr0 = 'k', clr1 = 'r', trlog=False, xlbl=None, ylbl=None, title=None):
#------------ Plot trace
    fig, ax = plt.figure()                #subplots()
    ax.plot(t, tr0, clr0, alpha=0.7)
    if tr1 is not None: ax.plot(t, tr1, clr1, alpha=0.7)
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
#
# -------------- Plot spectrogram  -------------------
def Pspect(time, trZ, wl=None):
    """
    tr: Trace object
    wlen(int or float): fft window length(s). None=128 samples
    cmap: Colormap instance
    show=False: Do not call plt.show() at end of routine
    clip=[0.0, 1.0]: Clip colormap to [lower,upper] ends.
    """
#
    fig = plt.figure()
    ax1 = fig.add_axes([0.1, 0.75, 0.7, 0.2]) #[left bottom width height]
    ax2 = fig.add_axes([0.1, 0.1, 0.7, 0.60], sharex=ax1)
    ax3 = fig.add_axes([0.83, 0.1, 0.03, 0.6])

    #plot waveform (top subfigure)    
    ax1.plot(time, trZ.data, 'k')
    cmap = mpl.cm.jet

    fig = trZ.spectrogram(wlen=wl,log=True, dbscale=True, per_lap=0.9,
                    cmap=cmap, show=False,
                    title='Spectrogram(dB)', axes=ax2)

    plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap), cax=ax3, label='dB')

    plt.show()
#
# -------------- End of function   ---------------------
#
# -------------- nearest_pow_2(x)  -------------------
def _nearest_pow_2(x):
    """
    Find power of two nearest to x
    :type x: float
    :param x: Number
    :rtype: Int
    :return: Nearest power of 2 to x
    """
    a = math.pow(2, math.ceil(np.log2(x)))
    b = math.pow(2, math.floor(np.log2(x)))
    if abs(a - x) < abs(b - x):
        return a
    else:
        return b
#
# -------------- End of function   ---------------------
