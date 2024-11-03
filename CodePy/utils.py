#!/usr/bin/python3
#--------  Code Dependencies   ----------
#\__________General Utilities____________/
import pprint
import time
import timeit
#import argparse
import sys
import numpy as np
from scipy import fftpack
from scipy import signal
from scipy.signal import find_peaks, peak_widths
import scipy.fft
#\__________ObsPy Utilities____________/
from obspy import read, Stream, Trace
import obspy.signal                    # To estimate envelpe
from obspy.signal.array_analysis import array_processing
from obspy.core.util import AttribDict
from obspy import UTCDateTime
from obspy.imaging.cm import obspy_sequential
#\_____________________________________/
#
#\__________Scripts List__________/
#
"""
Script         Description
otrstr         Preprocess the stream/trace st. Calls Proc
ObsPyStream    Create a new stream and loop over traces. Add the location to the “header”.
Proc           Preprocess the stream/trace st
BeamFK         Beamforming - FK Analysis
TrGain         apply a gain
nearest_pow_2  Find power of two nearest to x
AuxReset       saves result for the next cell.
lmt_ValInd     Limits 1-D array a1 to a given value and saves the indexes to limit another 1-D array a2.
"""
#\__________Scripts__________/
#
# ---------- Preprocess the stream/trace  ----------
def ObsPyStream(st):
    """ 
    Create a new ObsPy stream 
    <Arguments>
    st                -> A stream
    Returns a new st
    """ 
#------ Initialization
    gather = Stream()
    lon = 0.
    lat = 0.
    i   = 0
#
#------ Loop through traces
    for t in st:
        tr = Trace(data=t.data)
        tr.stats.sampling_rate = t.stats.sampling_rate
        tr.stats.station       = f"gloc[i,0]"                           # Assign station name   
        tr.stats.starttime     = t.stats.starttime
        tr.stats.network       = "TTB22"                                # Assign network code
        tr.stats.channel       = "HHZ"                                  # Assign channel code
        tr.stats.location      = "0"                                    # Assign location code
        tr.stats.distance      = t.stats.seg2.RECEIVER_LOCATION         # Distance along cable
        tr.stats.coordinates = \
                                 AttribDict({'latitude': ttb_gloc[i,1],
                                             'longitude': ttb_gloc[i,2],
                                             'elevation': 10.})
#
        lon += ttb_gloc[i,2]
        lat += ttb_gloc[i,1]
        i += 1
        gather += tr                                                      # gather.append(tr)
#
#------ Baricenter
    lon /= float(i)
    lat /= float(i)
    print(f">> Gather baricenter is at: lat= {lat}, lon= {lon}.")
#
    return st
#
# -------------- End of function   ---------------------
#
def Proc(st):
    """ 
    Preprocess the stream/trace st
    <Arguments>
    st                -> A stream
    """ 
#------ Copy I/P stream
    st0 = st.copy()
#------ Remove the mean and trend
    st.detrend("linear")
    st.detrend("demean") 
# 
#------ Notch 60Hz spectral line.
    print(f">> Notch the data at [bs 59.2, 60.8]")
    st, ftype, flims = u.TrFlt(st, ['bs', 59.2, 60.8])
#
#-------- Bandpass filter the data.
    print("\r", end="")
    print(f">> Bandpass filter the data")
    ent = input(f' Enter dflt [bp 5. 50.], or enter your choice: ')
    ent = None if ent else ['bp', 5., 50.]
    st, ftype, flims = u.TrFlt(st, ent=ent)
#-------- Taper the data: Hanning
    print(f'>> Taper window ends with 0.1')
    st.taper(type = 'hann', max_percentage = 0.1)
#---------- Amplitude loss
    gain = np.round(np.median(env0) / np.median(env1), 0)
    dbLoss = np.round(10.0 * np.log10( 1. / gain ), 0)
    ent = input(f'>> Stream loss in amplitude is {dbLoss}dB. Apply a gain of {gain} (dftl = no, any I/P = yes) ')
    gain = gain if ent else 1.
    st.data *= gain
#
    return st
#
# -------------- End of function   ---------------------
#
# ---------- Beamforming - FK Analysis  ----------
def BeamFK(st, MTparam, phone_geo, **kwargs):
    """ 
    Perform FK analysis with ObsPy. The data is bandpass filtered, prewhitening disabled.
    <Arguments>
    st                -> A stream
    MTparam           -> A list [fmin, fmax]
    phone_geo         -> A list of [[phone#, lat, lon], ...]
    stime, etime      -> Relative time limits (s) on the time window for event,
                          later to be transformed to UTCDateTime format.
    sll_o, sl_s       -> Slowness grid (s/km) and step fraction.
                 slm_y +─────+
                       │     │
                       │     │
                 sll_y +─────+    
                    sll_x   slm_x
    win_len, win_frac -> window length and step fraction (s).
    semb_thres, vel_thres -> infinitesimally small numbers; must not be changed.
    timestamp         -> written in 'mlabday', read directly by plotting routine.
    """ 
#
#------------ Select event times --------------
    print(f'\n>> Select start and end times (s) for beanforming')
    print(f'  stime etime')
    print(f'    └─────│──> Initial event time')
    print(f'          └──> Final event time')
    ent = input(f'   Enter t0 t1:\n')
    if not ent: raise ValueError("t0 t1 mandatory.")
    stime0, etime0 = np.array(ent.rstrip().split(' '), dtype=float)
    dummy = UTCDateTime(st[0].stats.starttime)
    stime = dummy + stime0
    etime = dummy + etime0
#------------ Relevant parameters --------------
    kwargs = dict(
# slowness grid : Xmin , Xmax , Ymin , Ymax , Slow Step
    sll_x =-3.0, slm_x =3.0, sll_y =-3.0, slm_y =3.0, sl_s =0.03,
# Changed for TTB 4/8/23
#    sll_x =-5.0, slm_x =5.0, sll_y =-5.0, slm_y =5.0, sl_s =0.025,
# sliding window properties
    win_len =1.0, win_frac =0.1,
# frequency properties
    frqlow =MTparam[0], frqhigh =MTparam[1], prewhiten =0,
# restrict output
    semb_thres=-1.e9, vel_thres=-1.e9 , timestamp='julsec',
    stime=stime , etime=etime
    )
#
    print(f'\n>>  F-K parameters')
    pprint.pprint(kwargs)
#    pline(['\n>> F-K parameters.'], width=13, align='^')
#    pline(['    1) Slowness grid (km).'], width=13, align='^')
#    pline(['        y=['+str(kwargs.get(sll_y))+', '+str(kwargs.get(slm_y))+']'], width=6, align='^')
#    pline(['        x=['+str(kwargs.get(sll_x))+', '+str(kwargs.get(slm_x))+']'], width=6, align='^')
#    pline(['        Slow Step= '+str(kwargs.get(sl_s))], width=6, align='^')
#    pline(['    2) Sliding window.'], width=13, align='^')
#    pline(['        win_len= '+str(kwargs.get(win_len))+', win_frac'+str(kwargs.get(win_frac))], width=6, align='^')
#-- Geophone numbering was INVERTED as for Geode headers. Sort files by phone now.
    dummy = phone_geo
    dummy = dummy[np.argsort(dummy[:,0])]
    phone_geo = dummy
#        └────> phone# (int), x, y (floats)
#
#------------- FK processing (obspy) --------------------
    for ind, tr in enumerate(st):
#-- Normalize
        tr.data = tr.data/max(tr.data)
#------------ Geophone positions in geographic. TTB altitude: 10m
        tr.stats.coordinates = AttribDict({
                'latitude':  phone_geo[ind,0],
                'elevation': float(10),
                'longitude': phone_geo[ind,1]})
#------------ Execute array_processing
    out = array_processing(st, **kwargs)
#-- return
    return out, stime0, etime0
#
# -------------- End of function   ----------------------------
#
#
# ---------- Stream/Trace Pre-processing  ----------
def otrstr(tr, MTparam, verbose=True):
    """ 
    Process the stream/trace tr. A simpler version of the original otrstr.
    <Arguments>
    tr                -> A stream or trace
    MTparam           -> A list with I/P parameters
     │─────>  dtr line ftype Fmin Fmax taper gain
     │   i =   0    1    2    3    4     5     6
     └─> dtr       -> Remove trends: 0 = no; 1 = yes  
         line      -> Notch 60Hz:    0 = no; 1 = yes
         ftype     -> Filter type:   lp=lowpass, hp=highpass, bp=bandpass, bs=bandstop  
         Fmin Fmax -> Corner frequencies.
         taper     -> Taper ends:    0 = no; 1 = yes
         gain      -> Gain data:     0 = no; 1 = yes
    """
#-- Fix # of corners for filters
    nc = 4
#-- Copy of original trace/stream
    tr0 = tr.copy()
#
#------------ Detrend
    if MTparam[0] == int(1):
        tr.detrend("linear")
        tr.detrend("demean")
        if verbose:
            dummy = np.mean(tr0[0].data, dtype=np.float64) if hasattr(tr0, 'traces') else np.mean(tr0.data, dtype=np.float64)
            print(f">> The mean of 1st original trace is {dummy}")
#------------  Filter 60Hz line with fixed limits for TTB
    if MTparam[1] == int(1):
        dummy = ['bs', 59.2, 60.8]
        tr, ftype, flims = TrFlt(tr, ent=dummy)
        if verbose:
            print(f">> Notched original trace from {dummy[1]} to {dummy[-1]}")
#
#------------  Filter the data.
    ent = [MTparam[2], MTparam[3], MTparam[4]]
    tr, ftype, flims = TrFlt(tr, ent=ent)
    if verbose: print(f">> Useful range due to {ftype} filter: {flims[0]} to {flims[1]}Hz.")
#
#------------  Taper the data with 10% Hanning
    if MTparam[5] == int(1):
        tr.taper(type = 'hann', max_percentage = 0.1)
#
#------------  Gain trace
    if MTparam[6] == int(1):
        tr, dummy = TrGain(tr0, tr)
        if verbose:
            print(f">> Applyied a maximum gain of {dummy[0]}dB to compensate processing losses.")
#
#------------  return processed trace
    return tr
#
# -------------- End of function   ---------------------
#
# ---------- Filter trace  ----------
"""
Find the loss of amplitude between two trace instances and apply a gain
tr....... Original trace (modified in place!!)
"""
def TrFlt(trZ, ent=None):
#------------ Asks for input
    if ent is None:
        print(f"Enter the filter, lower freq., upper freq., filter order, zerophase: ftype, f0, f1, nc, zP")
        print(f" Filter minimum options:")
        print(f"                lp (lowpass)  -> f0         is required")
        print(f"                hp (highpass) -> f0         is required")
        print(f"                bp (bandpass) -> f0 and f1 are required")
        print(f"                bs (bandstop) -> f0 and f1 are required")
        ent = input(f' Enter ftype, f0, f1, nc and zP (dflt: Nc=4, zP=True): ')
        ent = ent.rstrip().split(' ')
#
    if ent[0] == 'lp' or ent[0] == 'hp':
        ent[0] = 'lowpass' if ent[0] == 'lp' else 'highpass'
        ent[1] = float(ent[1])
        ent[len(ent):] = [int(4), True] if len(ent) == 2 else [int(ent[2]), True] if len(ent) == 3 else [int(ent[2]), bool(ent[3])]

    elif ent[0] == 'bp' or ent[0] == 'bs':
        ent[0]  = 'bandpass' if ent[0] == 'bp' else 'bandstop'
        ent[1:3] = [ float(dummy) for dummy in ent[1:3] ]
        ent[len(ent):] = [int(4), True] if len(ent) == 3 else [int(ent[3]), True] if len(ent) == 3 else [int(ent[3]), bool(ent[4])]
#
    """    trace data
    trZ=> +──+──...+──+
          |  |     |  └──> fNy, Nyquist  frequency
          |  |     └─────> f1, Higher cutoff frequency
          |  └───────────> f0, Lower cutoff frequency
          └──────────────> fmin        Minumum frequency
    """
#------- Use aliases to make life easier: unpack ent[:]
    if len(ent) == 5:
        ftype, f0, f1, nc , zP = ent[:]
    else:
        ftype, f0,     nc , zP = ent[:]
        f1 =   f0
#
#------- Acausal filter/zero phase if zP=True, otherwise a causal filter.
    if ftype in ['bandpass', 'bandstop']:
        trZ.filter(ftype, freqmin=f0, freqmax=f1, zerophase=zP, corners=nc)
        dummy = str(f0)+' '+str(f1)
    else:
        trZ.filter(ftype, freq=f0, zerophase=zP, corners=nc)
        dummy = str(f0)
#
    ftype = ftype+' '+dummy
#
#------- Return filtered trace, ftype and [fmin, fmax]
    return trZ, ftype, [f0, f1]
#
# -------------- End of function   ---------------------
#
# ---------- Trace gain  ----------
"""
Find the loss of amplitude between two trace instances and apply a gain.
Use the median of the envelope.
tr0..... Original trace
tr1..... Processed trace
"""
def TrGain(tr0, tr1):
#------------  Gain stream, which has attribute 'traces'
#               => Each trace is independently gained!
    if hasattr(tr0, 'traces'):
        maxg = float('-inf')
        ming = float('inf')
#
        for ind in range(len(tr0)):
            gain = np.median(obspy.signal.filter.envelope(tr0[ind].data))                     # max( )
            gain = np.round(gain/np.median(obspy.signal.filter.envelope(tr1[ind].data)),0)
#
            tr1[ind].data *= gain
            maxg = np.maximum(gain, maxg)
            ming = np.minimum(gain, ming)
#------------  Gain trace, which has attribute 'data'
    elif hasattr(tr0, 'data'):
            gain = np.median(obspy.signal.filter.envelope(tr0.data))
            gain = np.round(gain / np.median(obspy.signal.filter.envelope(tr1.data)),0)
#
            tr1.data *= gain
#
            maxg = gain; ming = gain
#-- Sanity
    else:
        raise ValueError("Bad hasattr")
#
#---------- Amplitude gain in dB
    maxg = np.round(10.0 * np.log10(maxg ), 0)
    ming = np.round(10.0 * np.log10(ming ), 0)
#
#------------  Returns gained trace
    return tr1, [maxg, ming]
#
# -------------- End of function   ----------------------------
#
# -------------- End of function   ---------------------
# ---------- nearest_pow_2  ----------
"""
Find power of two nearest to x
:type x: float
:param x: Number
:rtype: Int
:return: Nearest power of 2 to x
Call as:
nfft = int(_nearest_pow_2(wlen * samp_rate))
"""
def _nearest_pow_2(x):

    a = math.pow(2, math.ceil(np.log2(x)))
    b = math.pow(2, math.floor(np.log2(x)))
    if abs(a - x) < abs(b - x):
        return a
    else:
        return b
#
# -------------- End of function   ---------------------
#
# ---------- AuxReset  ----------
"""
An auxiliary code to aid in reseting cell to a previous cell value and
 saving result for the next cell.

 tr0 = previous cell value
 tr1 = present cell value for future
 lst= an existing list
 app = a *list* object to append to list *lst*
"""
def AuxReset(tr0, tr1, lst=None, app=None):
#------------------- Reset the trace from the last cell ------------------
    ent = input(f' Run this cell again (rtn= No)?: ') or False
    if ent:
        tr1 = tr0.copy()
        print(f' Original trace copied back. ==> RUN CELL AGAIN.')
    else:
        tr0 = tr1.copy()
        print(f' Resuls saved + a safety trace copy was created.')
        if isinstance(lst, list):
            lst[len(lst):] = app
#              +^^^^^^^^^+─> The slicing takes the space after the last item and unpacks "app"
#
    return tr0, tr1, lst
#
# -------------- End of function   ---------------------
#
# ---------- lmt_ValInd  ----------
"""
Limits 1-D array a1 to a given value and saves the indexes to limit another  1-D array a2.

Args:
a1: The numpy array to limit.
a2: The numpy array to limit based on the indexes of a1.
limit_value: The value to limit a1 to.

Returns:
The limited a1 and a2.
"""
def lmt_ValInd(a1, a2, limit):
#
    indexes = np.where(a1 <= limit)
    limited_a1 = a1[indexes]
    limited_a2 = a2[indexes]
    return limited_a1, limited_a2
#
# -------------- End of function   ---------------------
