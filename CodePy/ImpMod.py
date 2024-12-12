#------ Import External Libraries
import os, sys, glob
import numpy as np
import scipy
from scipy.fftpack import fft,ifft,next_fast_len
import pandas as pd
import matplotlib

matplotlib.use('nbagg')
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import Normalize
plt.style.use('ggplot')
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, AutoLocator, FixedLocator
from matplotlib.ticker import ScalarFormatter, FuncFormatter
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator
#------ Import specialized ObsPy packages + sanity.
try:
    import obspy
    print('ObsPy version ==>', obspy.__version__)
except TypeError:
    print('Stopping RUNTIME. ObsPy not found.')
    exit()

from obspy import read
from obspy import Stream
from obspy import UTCDateTime
from obspy import read, Stream
from obspy import Trace
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import FDSNException
from obspy.signal.invsim import cosine_taper
import obspy.signal        
import obspy.signal.filter             # To estimate envelpe
from obspy.signal.array_analysis import array_processing
from obspy.signal.array_analysis import clibsignal, cosine_taper, get_geometry, get_timeshift
from obspy.geodetics.base import gps2dist_azimuth
from obspy.core.util import AttribDict
from obspy.imaging.cm import obspy_sequential
from obspy.signal.invsim import corn_freq_2_paz
from obspy.signal.trigger import plot_trigger
from obspy.signal.trigger import classic_sta_lta
from obspy.geodetics.base import gps2dist_azimuth
from obspy import UTCDateTime, read, read_inventory
from obspy.signal.util import next_pow_2
from obspy.core import AttribDict

import utm
#------ Import local routines.

import plot as p
import utils as u