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
#\_____________________________________/
#
#\__________Scripts List__________/
#
"""
Script         Description
nearest_pow_2  Find power of two nearest to x


"""
#\__________Scripts__________/

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
"""
def AuxReset(tr0, tr1):
#------------------- Reset the trace from the last cell ------------------
    ent = input(f' Do you want to run this cell again?(rtn= Not):\n') or False
    if ent:
        tr1 = tr0.copy()
        print(f' Original trace copied back. Run code again.')
    else:
        tr0 = tr1.copy()
        print(f' A safety copy was created for the next cell.')
#\___A safety copy for an eventual reset in the next cell____/
#
    return tr0, tr1
#
# -------------- End of function   ---------------------