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



TrGain         apply a gain
nearest_pow_2  Find power of two nearest to x
AuxReset       saves result for the next cell.
lmt_ValInd     Limits 1-D array a1 to a given value and saves the indexes to limit another 1-D array a2.
"""
#\__________Scripts__________/
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
# ---------- Filter trace  ----------
"""
Find the loss of amplitude between two trace instances and apply a gain
tr....... Original trace (modified in place!!)
"""
def TrFlt(trZ):
#
    print(f"Enter the filter, lower freq., upper freq., filter order, zerophase: ftype, f0, f1, nc, zP")
    print(f" Filter minimum options:")
    print(f"                lp (lowpass)  -> f0         is required")
    print(f"                hp (highpass) -> f0         is required")
    print(f"                bp (bandpass) -> f0 and f1 are required")
    print(f"                bs (bandstop) -> f0 and f1 are required")

    ent = input(f' Enter ftype, f0, f1, nc and zP (dflt: Nc=4, zP=True): ')
    ent = ent.rstrip().split(' ')
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
Find the loss of amplitude between two trace instances and apply a gain
env0..... Original envelope
env1..... Processed envelope
gain.... if True apply gain
=False
"""
def TrGain(env0, env1):
#---------- Amplitude loss
    gain = np.round(np.median(env0) / np.median(env1), 0)
    dbLoss = np.round(10.0 * np.log10( 1. / gain ), 0)
    print(f'==> Loss in amplitude of {dbLoss}dB, requiring a gain of {gain}')
#---------- Returns multiplicative gain
    return gain
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
