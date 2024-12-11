#!/usr/bin/python3
#--------  Code Dependencies   ----------
#\__________General Utilities____________/
import pprint
import time
import timeit
#import argparse
import math 
import sys
import numpy as np
from numpy.linalg import inv
from scipy import fftpack
from scipy import signal
from scipy.signal import find_peaks, peak_widths
from scipy.fftpack import fft,ifft,next_fast_len
import scipy.fft
import warnings, utm
warnings.filterwarnings("ignore")
#\__________ObsPy Utilities____________/
import utm
from obspy import read, Stream, Trace
import obspy.signal                    # To estimate envelpe
from obspy.signal.array_analysis import array_processing
from obspy.core.util import AttribDict
from obspy import UTCDateTime
from obspy.imaging.cm import obspy_sequential
from obspy.core.util.base import _get_function_from_entry_point
from obspy.signal.util import next_pow_2
#\__________Local___________________/
import plot as p
#
#\__________Scripts List__________/
#
"""
Script         Description

pairs              Pairwise distances between points
divisors           Divisors  with zero reminder
get_coords         Retrieves coordinates
cut_trace           Cuts data into time segments
taper               Applies taper to time segments
whiten              Spectral spectral normalization and whitening
noise_processing    Wrapper for spectral normalization and whitening
moving_ave          Running smooth average
correlate           Does the cross-correlation in freq domain
extract_dispersion  Takes the dispersion image from CWT as input
sldtw_fk            FK analysis
array_lsq           Pairwise cross-correlation and least-squares inversion
sldtw               sliding time-window porcessing array
creastrm            Create a new stream
cat_wfrm            Constructs a waveform catalog
RGloc               Read a cvs file
otrstr              Preprocess the stream/trace st. Calls Proc
ObsPyStream         Create a new stream and loop over traces. Add the location to the “header”.
Proc                Preprocess the stream/trace st
BeamFK              Beamforming - FK Analysis
TrGain              apply a gain
nearest_pow_2       Find power of two nearest to x
AuxReset            saves result for the next cell.
lmt_ValInd          Limits 1-D array a1 to a given value and saves the indexes to limit another 1-D array a2.
"""
#\__________Scripts__________/

#
# ---------- Pairwise distances between points ----------
"""
  <Args>
    matrix -> A 3-column matrix: [[phone, x, y], ...] (x, y = cartesian coordinates)
    ainc   -> An angle increment
  <Returns>
    maxdist -> A 4-column matrix: [[phone1, phone2, dist,  angle], ...]
                                      int     int   float  float
                    angle -> (degrees) polar angle
                    dist  -> distance between phone1 and phone2)
                    maxdist is filtered by the largest distance near a given angle

"""
def pairs(matrix, angle_increment):
#
#------  Convert to numpy array if necessary
    if not isinstance(matrix, np.ndarray):
        matrix = np.array(matrix)
#
#------  Calculate the baricenter
    baricenter_x = np.mean(matrix[:, 1])
    baricenter_y = np.mean(matrix[:, 2])
#
#------ Create an empty list to store results
    results = []
#
#------ Iterate through angles 
    for angle in range(0, 180, angle_increment):
      # Convert angle to radians
        rad = np.radians(angle)
        # Define the diameter
        diameter_slope = np.tan(rad)
#
#------  Find the pair of points with minimum distance from the current diameter
        min_distance = float('inf')
        point1_index = -1
        point2_index = -1
        for i in range(len(matrix)):
            for j in range(i + 1, len(matrix)):
                # Calculate distance from the line
                dist_i = abs((matrix[i, 2] - baricenter_y) - diameter_slope * (matrix[i, 1] - baricenter_x)) / (np.sqrt(1 + diameter_slope**2))
                dist_j = abs((matrix[j, 2] - baricenter_y) - diameter_slope * (matrix[j, 1] - baricenter_x)) / (np.sqrt(1 + diameter_slope**2))
                total_dist = dist_i + dist_j

                if total_dist < min_distance:
                    
                   # print(point1_index,point2_index,total_dist,min_distance)
                    
                    min_distance = total_dist
                    point1_index = int(matrix[i,0])
                    point2_index = int(matrix[j,0])
#
#------ Calculate distance between the selected pair
        distance = np.sqrt((matrix[int(np.where(matrix[:,0] == point1_index)[0]), 1] - matrix[int(np.where(matrix[:,0] == point2_index)[0]), 1])**2 + (matrix[int(np.where(matrix[:,0] == point1_index)[0]), 2] - matrix[int(np.where(matrix[:,0] == point2_index)[0]), 2])**2)
#
#------ Angle relative to vertical axis
        results.append([int(point1_index), int(point2_index), distance, angle])  #90. - angle
    return np.array(results, dtype=object)
#
# -------------- End of function   ---------------------

#
# ---------- Temporal normalization ----------
"""
Temporal normalization the traces. Traces must be already demeaned, detrended and filtered.

    <Methods>
    1) "clipping"
        The signal is clipped to ->  clip_factor * sdev
    2) "clipping_iter"
        The signal is clipped iteratively based on clip_factor. Hint: clip_factor = 6

         values above 'clip_factor * std' 
        are divided by 'clip_weight'. until the whole signal is below 
        'clip_factor * std'
        clip_factor recommended: 6 (times std)

    3) "ramn"
        Running absolute mean normalization: a sliding window runs along the 
        signal. The values within the window are used to calculate a 
        weighting factor, and the center of the window is scaled by this 
        factor. 
            weight factor: w = np.mean(np.abs(tr.data[win]))/(2. * norm_win + 1) 
        finally, the signal is tapered with a tukey window (alpha = 0.2).

        norm_win: running window length, in seconds.
          recommended: half the longest period

    4) "1bit"
        Rtains only the sign of the signal.

    <Args>
    tr          -> A trace.
    clip_factor -> A multiplicative factor to the standard deviation (sdev). Hint: 6
                                    = clip_factor * sdev
    clip_weight -> A multiplicative factor to the standard deviation (sdev). Hint: 10
                                    = clip_factor * sdev
    norm_win    -> Running window length (s). Hint: half the longest period
  <Returns>
    tr: Normalized trace
"""
#
def normalize(tr, clip_factor=6, clip_weight=10, norm_win=None, norm_method="1bit"): 
#
    if norm_method == 'clipping':
        lim = clip_factor * np.std(tr.data)
        tr.data[tr.data > lim] = lim
        tr.data[tr.data < -lim] = -lim

    elif norm_method == "clipping_iter":
        lim = clip_factor * np.std(np.abs(tr.data))
        
        # as long as still values left above the waterlevel, clip_weight
        while tr.data[np.abs(tr.data) > lim] != []:
            tr.data[tr.data > lim] /= clip_weight
            tr.data[tr.data < -lim] /= clip_weight

    elif norm_method == 'ramn':
        lwin = tr.stats.sampling_rate * norm_win
        st = 0                                               # starting point
        N = lwin                                             # ending point

        while N < tr.stats.npts:
            win = tr.data[st:N]

            w = np.mean(np.abs(win)) / (2. * lwin + 1)
            
            # weight center of window
            tr.data[st + lwin / 2] /= w

            # shift window
            st += 1
            N += 1

        # taper edges
        taper = get_window(tr.stats.npts)
        tr.data *= taper

    elif norm_method == "1bit":
        tr.data = np.sign(tr.data)
        tr.data = np.float32(tr.data)

    return tr
#
# -------------- End of function   ---------------------
#
# ---------- get_coords ----------
"""
  Retrieves the coordinates (X, Y) of a point with a given index from a 3-column matrix.
  <Args>
    matrix: A NumPy array representing the 3-column matrix (index, X, Y).
    index: The integer index of the point to find.
  <Returns>
    A tuple (X, Y) representing the coordinates of the point.
  <Raises>
    ValueError: If the index is not found in the matrix.
"""
#
def get_coords(matrix, index):
#
  # Find the row corresponding to the given index
  row_index = np.where(matrix[:, 0] == index)

  if row_index[0].size == 0 :
    raise ValueError(f"Index {index} not found in the matrix.")
  
  #--- Extract X and Y coordinates
  x = matrix[row_index[0][0], 1]
  y = matrix[row_index[0][0], 2]

  return (x, y)
#
# -------------- End of function   ---------------------
#
# ---------- divisors ----------
"""
  Finds all divisors of a given number with zero reminder
  <Args>
    num:   The number for which to find divisors
    thres: A cap on the divisors
  <Returns>
    A list of integers
"""
#
def divisors(num, thres):
  div = []
  for i in range(1, num + 1):
    if num % i == 0:
      div.append(i)
  div =  div[1:-2] if len(div) >= 3 else None
#--- Limits list values
  if div is not None:
      div =  [div for div in div if div >= thres]
  return div
#
# -------------- End of function   ---------------------
#
# ---------- Cuts data into time segments ----------
"""
Cuts continous noise data into user-defined segments, estimate statistics for each segment and keep timestamps for later use.
  <Args>
    tr        -> A trace, an 1-dimensional timeseries array
    step      -> % overlapping (0, 1) between two sliding windows
    cc_len    -> Segment length (sec) to cut trace data
  <Returns>
    trace_stdS: standard deviation of the noise amplitude of each segment
    trace_madS: mad of the noise amplitude of each segment
    dataS_t:    timestamps of each segment
    dataS:      2D matrix of the segmented data
"""
def cut_trace(tr, cc_len, step):
#
#------ Trace stats and initialize return variables
#--- sampling_rate
    sps       = int(tr.stats.sampling_rate)
#--- trace window length (sec)
    tw        = round(tr.stats.endtime - tr.stats.starttime, 0)
#--- Date and time of the first data sample given in UTC (default value is “1970-01-01T00:00:00.0Z”)
    starttime = tr.stats.starttime - obspy.UTCDateTime(1970,1,1)
#--- Overlap in sec
    step = cc_len * step
#--- Number of segments and of points in each segment
    nseg      = int(np.floor((tw - cc_len)/step))    #tw/24*86400-cc_len
    npts = int(cc_len * sps)
#--- Copy data into array and check if data is shorter than the trim
    data      = tr.data
    if data.size < sps * nseg: raise ValueError("cc_len must be << trace.")
#--- Overlapping (sec) between two sliding windows
    step      = np.ceil(cc_len * step).astype(int)
#
#------ initialize variables
    dataS_t = []; dataS = []
#--- Relevant arrays
    dataS      = np.zeros(shape=(nseg, npts),dtype=np.float32)
    trace_madS = np.zeros(nseg,dtype=np.float32)
    trace_stdS = np.zeros(nseg,dtype=np.float32)
    trace_madS = np.zeros(nseg,dtype=np.float32)
    dataS_t    = np.zeros(nseg,dtype=np.float32)
#
#------ Statistic to detect segments that may be associated with a given event
#--- median absolute deviation over all trace
    all_madS = np.median(np.absolute(data - np.median(data)))
#--- standard deviation over all noise window
    all_stdS = np.std(data)
    if all_madS==0 or all_stdS==0 or np.isnan(all_madS) or np.isnan(all_stdS):
        print("Warn: madS or stdS == 0 for %s" % tr)
        return None, dataS_t, dataS
#

#    tw, cc_len, step 60.0 3.0 1.5
#    data.size, sps, nseg, npts 3000 50 38 150
    print('tw, cc_len, step', tw, cc_len, step)
    
    indx1 = 0
    for iseg in range(nseg):
        indx2 = indx1+npts

        
        print('data, iseg, indx1,indx2', data[indx1:indx2].size, iseg, indx1, indx2)

        
        dataS[iseg, :] = data[indx1:indx2]
        trace_madS[iseg] = (np.max(np.abs(dataS[iseg]))/all_madS)
        trace_stdS[iseg] = (np.max(np.abs(dataS[iseg]))/all_stdS)
        dataS_t[iseg]    = starttime+step*iseg
        indx1 = indx1+step*sps
#
#------ Data conditioning. It is assumed data is already demeaned and detrended
    dataS = taper(dataS)
#
#------  Returns
    return trace_stdS, trace_madS, dataS_t, dataS
#
# -------------- End of function   ---------------------
#
# ---------- Applies taper to time segments ----------
"""
Applies a cosine taper using obspy functions
  <Args>
    data -> An input data matrix
  <Returns>
    data -> A tapered data matrix
"""
def taper(data):
#ndata = np.zeros(shape=data.shape,dtype=data.dtype)
    if data.ndim == 1:
        npts = data.shape[0]
        # window length
        if npts*0.05>20:wlen = 20
        else:wlen = npts*0.05
        # taper values
        func = _get_function_from_entry_point('taper', 'hann')
        if 2*wlen == npts:
            taper_sides = func(2*wlen)
        else:
            taper_sides = func(2*wlen+1)
        # taper window
        win  = np.hstack((taper_sides[:wlen], np.ones(npts-2*wlen),taper_sides[len(taper_sides) - wlen:]))
        data *= win
    elif data.ndim == 2:
        npts = data.shape[1]
        # window length
        if npts*0.05>20:wlen = 20
        else:wlen = npts*0.05
        # taper values
        func = _get_function_from_entry_point('taper', 'hann')
        if 2*wlen == npts:
            taper_sides = func(2*wlen)
        else:
            taper_sides = func(2*wlen + 1)
        # taper window
        win  = np.hstack((taper_sides[:wlen], np.ones(npts-2*wlen),taper_sides[len(taper_sides) - wlen:]))
        for ii in range(data.shape[0]):
            data[ii] *= win
#
#------  Returns
    return data
#
# -------------- End of function   ---------------------
#
# ---------- Spectral spectral normalization and whitening ----------
"""
Transforms to frequency domain, whitens the amplitude spectrum in the frequency domain between *freqmin* and *freqmax*,
 and returns the whitened fft.
  <Args>
    tr: A trace: 1-dimensional timeseries array
    dt: The sampling space of the `data`
    freqmin: The lower frequency bound
    freqmax: The upper frequency bound
    smooth_N: integer, it defines the half window length to smooth
    freq_norm: whitening method between 'one-bit' and 'RMA'
  <Returns>
    FFTRawSign: The FFT (numpy.ndarray) of the whitened input trace in [freqmin, freqmax]
"""
#
def whiten(data, delta, freqmin, freqmax, smooth_N, freq_norm):
#
#------ Speed up FFT by padding to optimal size for FFTPACK
    if data.ndim == 1:
        axis = 0
    elif data.ndim == 2:
        axis = 1
#
    Nfft = int(next_fast_len(int(data.shape[axis])))
#
#------ Apodization number to the left and to the right
    Napod = 100
#
    Nfft = int(Nfft)
    freqVec = scipy.fftpack.fftfreq(Nfft, d=delta)[:Nfft // 2]
    J = np.where((freqVec >= freqmin) & (freqVec <= freqmax))[0]
    low = J[0] - Napod
    if low <= 0:
        low = 1
#
    left = J[0]
    right = J[-1]
    high = J[-1] + Napod
    if high > Nfft/2:
        high = int(Nfft//2)
#
    FFTRawSign = scipy.fftpack.fft(data, Nfft,axis=axis)
#
#------  Left tapering:
    if axis == 1:
        FFTRawSign[:,0:low] *= 0
        FFTRawSign[:,low:left] = np.cos(
            np.linspace(np.pi / 2., np.pi, left - low)) ** 2 * np.exp(
            1j * np.angle(FFTRawSign[:,low:left]))
#
#------ Pass band:
        if freq_norm == 'phase_only':
            FFTRawSign[:,left:right] = np.exp(1j * np.angle(FFTRawSign[:,left:right]))
        elif freq_norm == 'rma':
            for ii in range(data.shape[0]):
                tave = moving_ave(np.abs(FFTRawSign[ii,left:right]),smooth_N)
                FFTRawSign[ii,left:right] = FFTRawSign[ii,left:right]/tave
#
#------  Right tapering:
        FFTRawSign[:,right:high] = np.cos(
            np.linspace(0., np.pi / 2., high - right)) ** 2 * np.exp(
            1j * np.angle(FFTRawSign[:,right:high]))
        FFTRawSign[:,high:Nfft//2] *= 0
#
#------  Hermitian symmetry -> input is real
        FFTRawSign[:,-(Nfft//2)+1:] = np.flip(np.conj(FFTRawSign[:,1:(Nfft//2)]),axis=axis)
    else:
        FFTRawSign[0:low] *= 0
        FFTRawSign[low:left] = np.cos(
            np.linspace(np.pi / 2., np.pi, left - low)) ** 2 * np.exp(
            1j * np.angle(FFTRawSign[low:left]))
#
#------  Pass band:
        if freq_norm == 'phase_only':
            FFTRawSign[left:right] = np.exp(1j * np.angle(FFTRawSign[left:right]))
        elif freq_norm == 'rma':
            tave = moving_ave(np.abs(FFTRawSign[left:right]),smooth_N)
            FFTRawSign[left:right] = FFTRawSign[left:right]/tave
#
#------  Right tapering:
        FFTRawSign[right:high] = np.cos(
            np.linspace(0., np.pi / 2., high - right)) ** 2 * np.exp(
            1j * np.angle(FFTRawSign[right:high]))
        FFTRawSign[high:Nfft//2] *= 0
#
#------  Hermitian symmetry  -> input is real
        FFTRawSign[-(Nfft//2)+1:] = FFTRawSign[1:(Nfft//2)].conjugate()[::-1]
#
#------  Returns
    return FFTRawSign
#
# -------------- End of function   ---------------------
#
# ---------- Wrapper for spectral normalization and whitening ----------
"""
Transforms to frequency domain, whitens the amplitude spectrum in the frequency domain between *freqmin* and *freqmax*,
 and returns the whitened fft.
  <Args>
    data: A 2D matrix of all segmented noise data
    dt: The sampling distance in seconds
    freqmin:   The lower frequency bound
    freqmax:   The upper frequency bound
    smooth_N:  Integer, it defines the half window length to smooth
    time_norm: Time-domain normalization -> 'one_bit' or 'rma'
    freq_norm: Whitening method -> 'one-bit' and 'RMA'
  <Returns>
    FFTRawSign: The FFT (numpy.ndarray) of the whitened input trace in [freqmin, freqmax]
"""
def noise_processing(data, delta, freqmin, freqmax, 
                     smooth_N = 100, time_norm = False, freq_norm = False):
#
    N = data.shape[0]
#
#------  Time normalization
    if time_norm:
#--- Sign normalization
        if time_norm == 'one_bit':
            white = np.sign(data)
#--- Normalization over smoothed absolute average; running mean
        elif time_norm == 'rma':
            white = np.zeros(shape=data.shape,dtype=data.dtype)
            for kkk in range(N):
                white[kkk,:] = data[kkk,:]/moving_ave(np.abs(data[kkk,:]),smooth_N)
#--- Don't normalize
    else:
        white = data
#
#------ whiten
    if freq_norm:
#--- Whiten and return FFT
        White = whiten(white, delta, freqmin, freqmax, smooth_N, freq_norm)
    else:
#--- Return FFT
        Nfft = int(next_fast_len(int(dataS.shape[1])))
        White = scipy.fftpack.fft(white, Nfft, axis=1)
#
#------  Returns
    return White
#
# -------------- End of function   ---------------------
#
# ---------- Running smooth average ----------
"""
This function does running smooth average for an array (use numba for performance)
  <Args>
    A: 1-D array of data to be smoothed
    N: Defines the half window length to smooth (integer)
  <Returns>
    B: 1-D array with smoothed data
"""
def moving_ave(A,N):
#
    A = np.concatenate((A[:N],A,A[-N:]),axis=0)
    B = np.zeros(A.shape,A.dtype)

    tmp=0.
    for pos in range(N,A.size-N):
        # do summing only once
        if pos==N:
            for i in range(-N,N+1):
                tmp+=A[pos+i]
        else:
            tmp=tmp-A[pos-N-1]+A[pos+N]
        B[pos]=tmp/(2*N+1)
        if B[pos]==0:
            B[pos]=1
    return B[N:-N]
#
# -------------- End of function   ---------------------
#
# ---------- Cross-correlation ----------
"""
Does the cross-correlation in freq domain. Keeps the sub-stacks of the cross-correlation if needed, taking advantage of
 the linear relationship of ifft. Stacking is performed in spectrum domain, reducing the total number of ifft.
 
  <Args>
    fft1:          Source station power (smoothed) spectral density
    fft2:          Receiver station raw FFT spectrum
    dt:            sampling rate (in s)
    freqmin:       minimum frequency (Hz)
    freqmax:       maximum frequency (Hz)
    Nfft:          number of frequency points for ifft
    maxlag:        maximum lags to keep in the cross correlation
    dataS_t:       matrix of datetime object
    method:        'xcorr' or "coherency"
    substack:      sub-stack cross-correlationS or not
    substack_len:  multiples of cc_len to stack ove
    smoothspect_N: number of points to be smoothed for running-mean average (freq-domain, method=="coherency")
  <Returns>
    s_corr: 1D or 2D matrix of the averaged or sub-stacks of cross-correlation functions in time domain
    t_corr: timestamp for each sub-stack or averaged function
    n_corr: number of included segments for each sub-stack or averaged function
"""
#
def correlate(fft1, fft2, dt, dataS_t,
              freqmin, freqmax, maxlag,
              method = 'xcorr', substack = False, substack_len  = None, smoothspect_N = None):
#
#------ check on fft1/2 dimensions: put a cap on them, assuning same diensions in both
    if len(fft1) != len(fft2): raise ValueError("fft1 and fft2 must have the same length")
#--- nwin = number of segments in the 2D fft(i) matrix
#
#------ fft1 & fft2 are 2-D
    if fft1.ndim == 2:
        nwin  = fft1.shape[0]
        Nfft = fft1.shape[1]
        Nfft2 = Nfft//2
#--- Cap second dimension to Nfft2
        fft1 = fft1[:, :Nfft2]
        fft2 = fft2[:, :Nfft2]
#--- convert 2D array into 1D
        corr = np.zeros(nwin*Nfft2,dtype=np.complex64)
#        corr = fft1.reshape(fft1.size,)*fft2.reshape(fft2.size,)
        corr = fft1.reshape(-1) * fft2.reshape(-1)
#
#------ fft1 & fft2 are 1-D
    elif fft1.ndim == 1:
        nwin  = 1
        Nfft = fft1.shape[0]
        Nfft2 = Nfft//2
#--- Cap to Nfft2
        fft1 = fft1[:Nfft2]
        fft2 = fft2[:Nfft2]
        corr = np.zeros(nwin*Nfft2,dtype=np.complex64)
        corr = fft1 * fft2
    else:
        raise ValueError("fft must be either 1 or 2-dimensions")
#
#------ method == coherency
    if method == "coherency":
        temp = moving_ave(np.abs(fft2.reshape(fft2.size,)),smoothspect_N)
        corr /= temp
    corr  = corr.reshape(nwin,Nfft2)            #it was Nfft2
#
#------ method != coherency
    if substack:
        if substack_len == cc_len:
            # choose to keep all fft data for a day
            s_corr = np.zeros(shape=(nwin,Nfft),dtype=np.float32)   # stacked correlation
            ampmax = np.zeros(nwin,dtype=np.float32)
            n_corr = np.zeros(nwin,dtype=np.int16)                  # number of correlations for each substack
            t_corr = dataS_t                                        # timestamp
            crap   = np.zeros(Nfft,dtype=np.complex64)
            for i in range(nwin):
                n_corr[i]= 1
                crap[:Nfft2] = corr[i,:]            #it was Nfft2
                crap[:Nfft2] = crap[:Nfft2]-np.mean(crap[:Nfft2])   # remove the mean in freq domain (spike at t=0). It was Nfft2
                crap[-(Nfft2)+1:] = np.flip(np.conj(crap[1:(Nfft2)]),axis=0)            #it was Nfft2
                crap[0]=complex(0,0)
                s_corr[i,:] = np.real(np.fft.ifftshift(scipy.fftpack.ifft(crap, Nfft, axis=0)))

            # remove abnormal data
            ampmax = np.max(s_corr,axis=1)
            tindx  = np.where( (ampmax<20*np.median(ampmax)) & (ampmax>0))[0]
            s_corr = s_corr[tindx,:]
            t_corr = t_corr[tindx]
            n_corr = n_corr[tindx]

        else:
            # get time information
            Ttotal = dataS_t[-1]-dataS_t[0]             # total duration of what we have now
            tstart = dataS_t[0]

            nstack = int(np.round(Ttotal/substack_len))
            ampmax = np.zeros(nstack,dtype=np.float32)
            s_corr = np.zeros(shape=(nstack,Nfft),dtype=np.float32)
            n_corr = np.zeros(nstack,dtype=np.int)
            t_corr = np.zeros(nstack,dtype=np.float)
            crap   = np.zeros(Nfft,dtype=np.complex64)

            for istack in range(nstack):
                # find the indexes of all of the windows that start or end within
                itime = np.where( (dataS_t >= tstart) & (dataS_t < tstart+substack_len) )[0]
                if len(itime)==0:tstart+=substack_len;continue

                crap[:Nfft2] = np.mean(corr[itime,:],axis=0)   # linear average of the correlation. It was Nfft2
                crap[:Nfft2] = crap[:Nfft2]-np.mean(crap[:Nfft2])   # remove the mean in freq domain (spike at t=0). It was Nfft2
                crap[-(Nfft2)+1:]=np.flip(np.conj(crap[1:(Nfft2)]),axis=0)            #it was Nfft2
                crap[0]=complex(0,0)
                s_corr[istack,:] = np.real(np.fft.ifftshift(scipy.fftpack.ifft(crap, Nfft, axis=0)))
                n_corr[istack] = len(itime)               # number of windows stacks
                t_corr[istack] = tstart                   # save the time stamps
                tstart += substack_len
                #print('correlation done and stacked at time %s' % str(t_corr[istack]))

            # remove abnormal data
            ampmax = np.max(s_corr,axis=1)
            tindx  = np.where( (ampmax<20*np.median(ampmax)) & (ampmax>0))[0]
            s_corr = s_corr[tindx,:]
            t_corr = t_corr[tindx]
            n_corr = n_corr[tindx]
#
#------ Not substack
    else:
        # average daily cross correlation functions
        ampmax = np.max(corr,axis=1)
        tindx  = np.where( (ampmax<20*np.median(ampmax)) & (ampmax>0))[0]
        n_corr = nwin
        s_corr = np.zeros(Nfft,dtype=np.float32)
        t_corr = dataS_t[0]
        crap   = np.zeros(Nfft,dtype=np.complex64)
        crap[:Nfft2] = np.mean(corr[tindx],axis=0)            #it was Nfft2
        crap[:Nfft2] = crap[:Nfft2]-np.mean(crap[:Nfft2],axis=0)            #it was Nfft2
        crap[-(Nfft2)+1:]=np.flip(np.conj(crap[1:(Nfft2)]),axis=0)            #it was Nfft2
        s_corr = np.real(np.fft.ifftshift(scipy.fftpack.ifft(crap, Nfft, axis=0)))

    # trim the CCFs in [-maxlag maxlag]
    t = np.arange(-Nfft2+1, Nfft2)*dt            #it was Nfft2
    ind = np.where(np.abs(t) <= maxlag)[0]
    if s_corr.ndim==1:
        s_corr = s_corr[ind]
    elif s_corr.ndim==2:
        s_corr = s_corr[:,ind]
#
#------  Returns
    return s_corr, t_corr, n_corr
#
# -------------- End of function   ---------------------


#
# ---------- extract the dispersion from the image ----------
"""
Takes the dispersion image from CWT as input, tracks the global maxinum on
    the wavelet spectrum amplitude and extract the sections with continous and high quality data
  <Args>
    amp:   2D amplitude matrix of the wavelet spectrum
    phase: 2D phase matrix of the wavelet spectrum
    per:   period vector for the 2D matrix
    vel:   vel vector of the 2D matrix
  <Returns>
    per:  central frequency of each wavelet scale with good data
    gv:   group velocity vector at each frequency
"""
def extract_dispersion(amp,per,vel):
#
    maxgap = 5
    nper = amp.shape[0]
    gv   = np.zeros(nper,dtype=np.float32)
    dvel = vel[1]-vel[0]
#
#------  find global maximum
    for ii in range(nper):
        maxvalue = np.max(amp[ii],axis=0)
        indx = list(amp[ii]).index(maxvalue)
        gv[ii] = vel[indx]
#
#------  check the continuous of the dispersion
    for ii in range(1,nper-15):
        # 15 is the minumum length needed for output
        for jj in range(15):
            if np.abs(gv[ii+jj]-gv[ii+1+jj])>maxgap*dvel:
                gv[ii] = 0
                break
#
#------  remove the bad ones
    indx = np.where(gv>0)[0]
#
#------  Returns
    return per[indx],gv[indx]
#
# -------------- End of function   ---------------------
#
# ---------- FK analysis ----------
""" 
    FK analysis with ObsPy. The data is bandpass filtered, prewhitening disabled.
    <Arguments>
    st                -> A stream
    MTparam           -> 
    tstart, tend      -> Relative time limits (s) on the time window for event,
                          later to be transformed to UTCDateTime format.
    sll_o, sl_s       -> Slowness grid (s/km) and step fraction.
                 slm_y +─────+
                       │     │
                       │     │
                 sll_y +─────+    
                    sll_x   slm_x
    win_len, win_frac -> window length and the overlap fraction (s).
    semb_thres, vel_thres -> infinitesimally small numbers; must not be changed.
    timestamp         -> written in 'mlabday', read directly by plotting routine.
""" 
def sldtw_fk(st, tstart, tend, MTparam, **kwargs):
#
#------ Coordinates
    for st_i in st:
        st_i.stats.coordinates = AttribDict({
            'latitude': st_i.stats.coordinates['latitude'],
            'elevation': st_i.stats.coordinates['elevation'],
            'longitude': st_i.stats.coordinates['longitude']})
#
#------------ Relevant parameters --------------
    kwargs = dict(
#
#------ slowness grid : Xmin , Xmax , Ymin , Ymax , Slow Step
    sll_x =-3.0, slm_x =3.0, sll_y =-3.0, slm_y =3.0, sl_s =0.03,
# Changed for TTB 4/8/23
#    sll_x =-5.0, slm_x =5.0, sll_y =-5.0, slm_y =5.0, sl_s =0.025,
#
#------ sliding window properties
    win_len =1.0, win_frac =0.1,
#
#------ Frequency properties from Filter parameters
    frqlow = MTparam[3], frqhigh = MTparam[4], prewhiten =0,
#
#------ restrict output
    semb_thres=-1.e9, vel_thres=-1.e9 , timestamp='julsec',
#
#------ Time 
    stime=st[0].stats.starttime+tstart, etime=st[0].stats.starttime+tend, verbose=False
    )
#
#------ Calculations
    out = array_processing(st, **kwargs)
#
#------ returns 
    return out
#
# -------------- End of function   ---------------------
#
# ---------- array coordinates ----------
"""
    Returns the array coordinates for an array, in km with respect to the reference array provided

    <Arguments>
    st          -> An ObsPy Stream object: the array data
    ref_station -> A String containing the name of the reference station or point

    <Returns>
    X           -> [Nx2] NumPy array of array coordinates in km
    stnm        -> [Nx1] list of names
"""
#
def acoord(st, ref_station, km=True):
#
    X = np.zeros((len(st), 2))
    stnm = []
#
#------ Transforms gather st coordinates to utm
    for i in range(0, len(st)):
        E, N, _, _ = utm.from_latlon(st[i].stats.coordinates['latitude'], st[i].stats.coordinates['longitude'])
        X[i,0] = E; X[i,1] = N
        stnm.append(st[i].stats.station)
#
#------ Transforms ref_station coordinates to utm
        Xref = np.zeros((2))
        Xref[0], Xref[1], _, _ = utm.from_latlon(ref_station[0].stats.coordinates['latitude'], ref_station[0].stats.coordinates['longitude'])
#
#------ Adjusting to the reference station, and converting to m or km:
#    ref_station_ix = np.where(np.array(stnm) == ref_station)[0][0]    # index of reference station
    km = 1000. if km else 1.
    X[:,0] = (X[:,0] - Xref[0])
    X[:,1] = (X[:,1] - Xref[1])
    X = X / km
#
#------ returns 
    return X, stnm
#
# -------------- End of function   ---------------------
#
# ---------- Adds beam channel to a Stream ----------
"""
    Adds the beam channel to an ObsPy Stream using the time of the reference station
    <Arguments>
    st     -> An ObsPy Stream object: the array data
    beam   -> NumPy array containing beam data
    refst  -> Str with the name of the reference station
    reftr  -> Trace of refst
    <Returns>
    st  -> ObsPy Stream object including beam channel with station name = 'BEAM'
    WARNING: st is modifyied in place!
"""
#
def add_beam(st, beam, refst, reftr=0):
#
#------ Obtain trace for reference station
    tr = st.select(station=refst)[reftr].copy()
    tr.data = beam[0:len(tr.data)]
    tr.stats.station = 'BEAM'
    st.append(tr)
#------ returns station name = 'BEAM'
    return st
#
# -------------- End of function   ---------------------
#
# ---------- sliding time-window array ----------
"""
    Performs sliding time-window array processing using the least-squares array processing
    method in array_lsq (by Stephen Arrowsmith)
    <Call>
    array_lsq
    <Arguments>
    st     -> An ObsPy Stream object: the array data
    X      -> Array coordinates
    └───────> [gloc[i,1], gloc[i,2]]
    tstart -> Start time,in seconds after Stream start-time
    tend   -> End time, in seconds after Stream start-time
    twin   -> Time window (s) for array processing 
    overlap - Overlap (s) for array processing
    <Returns>
plt.plot(X[:,0], X[:,1], 'o')
    
"""
#
def sldtw(st, X, tstart, tend, twin, overlap):
#
    time_start = st[0].stats.starttime + tstart
    time_end = time_start+tend

    time_start_i = time_start
    time_end_i = time_start_i+twin

    t = tstart; T = []; V = []; B = []
    while time_end_i < time_end:
        st_win = st.slice(time_start_i, time_end_i)
        vel, baz = array_lsq(st_win, X)
        T.append(t + twin/2); V.append(vel); B.append(baz)
        t = t + overlap
        time_start_i = time_start_i + overlap
        time_end_i = time_end_i + overlap
    T = np.array(T); V = np.array(V); B = np.array(B)
    
    return T, V, B
#
# -------------- End of function   ---------------------
#
# ---------- Pairwise cross-correlation ----------
"""
    Performs pairwise cross-correlation on each trace in st, and least-squares inversion
    for the slowness vector corresponding to the best-fitting plane wave (by Stephen Arrowsmith)
    
    <Arguments>
    st     -> An ObsPy Stream object: the array data in a time window suitable for array processing
    X      -> Array coordinates: [Nx2] NumPy array of array coordinates (in km relative to a reference point)
    └───────> [phone, lat, lon]
    <Returns>
    baz     -> Backazimuth (in degrees from North)
    vel     -> Phase velocity (in km/s)    
"""
#
def array_lsq(st, X):
#
#------ Initializing empty arrays for array distances and delay times
    N = len(st)           # Number of elements
    M = int(N*(N-1)/2)    # Number of pairs of elements
    R = np.zeros((M,2))   # Array to hold relative coordinates between elements
    tau = np.zeros((M,1)) # Array to hold delay times

    k = 0
    for i in range(0,N):
        for j in range(i+1,N):

            tr1 = st[i]; tr2 = st[j]
            C = np.correlate(tr1.data, tr2.data, mode='full')
            lags = np.arange(-np.floor(len(C)/2), np.floor(len(C)/2)+1, 1)*tr1.stats.delta

            # Computing lag corresponding to maximum correlation:
            ix = np.argmax(C); tau[k] = lags[ix]

            # Computing vector of distances between array coordinates:
            R[k,:] = X[i,:] - X[j,:]

            k = k + 1
#
#------ Performing least squares inversion:
    R = np.matrix(R); tau = np.matrix(tau)
    u = (inv(np.transpose(R)*R)*np.transpose(R))*tau
    v = 1/np.sqrt(u[0]**2 + u[1]**2)
    azimut = 180 * math.atan2(u[0], u[1]) / math.pi
    baz = (azimut % -360 + 180) % 360
    
    return float(v), float(baz)
#
# -------------- End of function   ---------------------
#
# ---------- Create a stream ----------
"""
    Create a new stream.  Loop over tracesand add the distance information
    in the “header” and then add that trace to the stream.
    <Arguments>
    st     -> An original (SEG2) stream
    gloc   -> A numpy array with phone number and coordinates in degrees
      └─────> [phone, lat, lon]
    hght   -> A fixed height (m)
    <Returns>
    ganther -> a gather
    bcenter -> Gather baricenter (field: lat= -1.201925, lon= -48.506498; Datum  WGS84)
"""
def creastrm(st, gloc, hght = 5.):
#------ Initialization
    gather = Stream()
    lon = 0.; lat = 0.; i   = 0
#
#------ Loop through traces
    for t in st:
        tr = Trace(data=t.data)
        tr.stats.sampling_rate = t.stats.sampling_rate
        tr.stats.station       = f"{gloc[i,0]}"                         # Assign station name   
        tr.stats.starttime     = t.stats.starttime
        tr.stats.network       = "TTB22"                                # Assign network code
        tr.stats.channel       = "HHZ"                                  # Assign channel code
        tr.stats.location      = "0"                                    # Assign location code
        tr.stats.distance      = t.stats.seg2.RECEIVER_LOCATION         # Distance along cable
        tr.stats.coordinates = \
                                 AttribDict({'latitude':  gloc[i,1],
                                             'longitude': gloc[i,2],
                                             'elevation': hght})
#
        lon += gloc[i,2]
        lat += gloc[i,1]
        i += 1
        gather += tr                                                      # gather.append(tr)
#
#------ Gather baricenter
    lon /= float(i)
    lat /= float(i)
    bcenter = Stream()
    trace = Trace(data=np.zeros(len(st[0].data)))
    trace.stats.sampling_rate = st[0].stats.sampling_rate 
    trace.stats.station       = f"bcenter"  
    trace.stats.starttime     = st[0].stats.starttime 
    trace.stats.network       = "TTB22"
    trace.stats.channel       = "HHZ"  
    trace.stats.location      = "0"    
    trace.stats.distance      = -1.
    trace.stats.coordinates = \
                             AttribDict({'latitude':  lat,
                                         'longitude': lon,
                                         'elevation': hght})
    bcenter += trace
#
    print(f">> A new gather with baricenter {bcenter[0].stats.station} is created at {bcenter[0].stats.coordinates}")
#------ Return gather
    return gather, bcenter
#
# -------------- End of function   ---------------------
#
# ---------- Construct a waveform catalog ----------
"""
Construct catalog of template waveforms to generate characteristic functions
 to be used for event detection.
    <Arguments>
    trZ    -> A trace
    cat    -> A catalog

    Returns an extended catalog
"""
def cat_wfrm(trZ, cat, tname, phone):
    while True:
#------ Read waveform limits, otherwise return
        trTplt = trZ.copy()
        ent = input(f'<< Enter pck0 and pck1 in s (rtn=exit): ')
        if not ent: return cat, tname
#
        ent = ent.rstrip().split(' ')
        pck0 = trTplt.stats.starttime + float(ent[0])
        pck1 = trTplt.stats.starttime + float(ent[1])
#------ Trim the template
        trTplt.slice(pck0, pck1)
        print(f">> Template is at [{pck0}, {pck1}] with {len(trTplt)} data points.")
#------ Relative time: nummpy array
        dummy = trTplt.times(type="relative")        
#------ Check with wavefrom and fft. IF OK add to catalog
        fNy = 0.8 * trTplt.stats.sampling_rate / 2.
        FtrTplt = np.fft.rfft(trTplt.data)
        freq = np.linspace(0, fNy, len(FtrTplt))
        p.pltTrSp( dummy, trTplt.data, freq, abs(FtrTplt),
                  x1label='s', y1label='Ampl.', y1log=False, clr1 = 'k', 
                  x2label='Hz', y2label='Spec.', y2log=True, clr2 = 'r' )        

        ent = input(f'<< Accept waveform (rtn=yes): ')
        if not ent:
            cat.append(trTplt)
            tname.append(f"p{phone:02}t{len(cat):03}")
#                            └──────+...+────│─────+──> wavefrom catalog
#                                            └────────> waveform name: p..t...
#                                                        p..=phone number, t... template number
#
# -------------- End of function   ---------------------
#
# ---------- A snippet to read a CSV file  ----------
"""
Read a cvs file and stores the information in an object numpy array
"""
def RGloc(filename):
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
            data = []
            for line in lines[1:]:               # Skip header row
                parts = line.strip().split(',')
                data.append([int(parts[0]), float(parts[1]), float(parts[2])])
        return np.array(data, dtype=object)
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return None                                 
#
#
# -------------- End of function   ---------------------
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
    st, ftype, flims = TrFlt(st, ['bs', 59.2, 60.8])
#
#-------- Bandpass filter the data.
    print("\r", end="")
    print(f">> Bandpass filter the data")
    ent = input(f' Enter dflt [bp 5. 50.], or enter your choice: ')
    ent = None if ent else ['bp', 5., 50.]
    st, ftype, flims = TrFlt(st, ent=ent)
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
def otrstr(tr0, MTparam, line = ['bs', 59.2, 60.8], verbose=True):
    """ 
    Process the stream/trace tr. A simpler version of the original otrstr.
    <Arguments>
    tr                -> A stream or trace
    MTparam           -> A list with I/P parameters
     │─────>  dtr line ftype Fmin Fmax taper gain
     │   i =   0    1    2    3    4     5     6
     └─> dtr       -> Remove trends: 0 = no; 1 = yes  
         line      -> Notch 60Hz:    0 = no; 1 = yes
         ftype     -> Filter type:   lp=lowpass, hp=highpass, bp=bandpass, bs=bandstop, no= no filter
         Fmin Fmax -> Corner frequencies.
         taper     -> Taper ends:    0 = no; 1 = yes
         gain      -> Gain data:     0 = no; 1 = yes
    """
#-- Fix # of corners for filters
    nc = 4
#-- Copy of original trace/stream
    tr = tr0.copy()
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

#        dummy = ['bs', 59.2, 60.8]
        tr, ftype, flims = TrFlt(tr, ent=line)
        if verbose:
            print(f">> Notched original trace from {line[1]} to {line[2]}")
#
#------------  Filter the data.
    ent = [MTparam[2], MTparam[3], MTparam[4]]
    if MTparam[2] != 'no':
        tr, ftype, flims = TrFlt(tr, ent=ent)
        if verbose: print(f">> Useful range due to {ftype} filter: {flims[0]} to {flims[1]}Hz.")
#
#------------  Taper the data with 10% Hanning
    if MTparam[5] == int(1):
        tr.taper(type = 'hann', max_percentage = 0.1)
#
#------------  Gain trace
    if MTparam[6] == int(1):
        tr, dummy = TrGain(tr, tr0)
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
