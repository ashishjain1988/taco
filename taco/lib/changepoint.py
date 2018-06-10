'''
TACO: Multi-sample transcriptome assembly from RNA-Seq
'''
from collections import namedtuple

import numpy as np

from taco.lib.stats import mannwhitneyu
from taco.lib.cchangepoint import mse as mse_cython


__author__ = "Matthew Iyer, Yashar Niknafs, and Balaji Pandian"
__copyright__ = "Copyright 2012-2018"
__credits__ = ["Matthew Iyer", "Yashar Niknafs", "Balaji Pandian"]
__license__ = "MIT"
__version__ = "0.7.3"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"

ChangePoint = namedtuple('ChangePoint', ['pos', 'start', 'end',
                                         'pvalue', 'foldchange', 'sign'])


def mse(x):
    mse_min = None
    mse_i = -1
    for i in xrange(1, len(x)):
        x1 = x[:i]
        x2 = x[i:]
        mse1 = np.power(x1 - x1.mean(), 2).sum()
        mse2 = np.power(x2 - x2.mean(), 2).sum()
        mse = mse1 + mse2
        if (mse_i == -1) or (mse < mse_min):
            mse_min = mse
            mse_i = i
    return mse_min, mse_i


def mwu(a1, a2):
    if len(a1) == 0 or len(a2) == 0:
        return (None, 1)
    elif (len(np.unique(a1)) == 1 and
          np.array_equal(np.unique(a1), np.unique(a2))):
        return (None, 1)
    else:
        U, p = mannwhitneyu(a1, a2)
        return (U, p)


def slope_extract(slope_a, i):
    slope_a = np.sign(slope_a)
    sign = slope_a[i]
    if slope_a[i] == 0:
        j = 0
        k = 0
    else:
        j = 0
        while (i - j) >= 0 and slope_a[i - j] == slope_a[i]:
            j += 1
        k = 0
        while (i + k) < len(slope_a) and slope_a[i + k] == slope_a[i]:
            k += 1
    return (j, k, sign)


def bin_seg_slope(a, s_a, pval=0.05, fc_cutoff=0.80, size_cutoff=20,
                  cp_func=mse_cython, cps=None, offset=0):
    '''
    a: numpy input array
    s_a: numpy array containing slope of changes in array 'a'
    pval: mann-whitney-u p-value threshold
    fc_cutoff: fold change threshold across change point
    size_cutoff: stop searching for change points when vector
                 length < size_cutoff
    cp_func: function to choose change point
    cps: (recursion) list of change points
    offset: (recursion) offset into vectors
    '''
    if cps is None:
        cps = []
    if a.shape[0] < size_cutoff:
        return cps
    # condense expression array to indices where expression values differ
    diff_indexes = np.ediff1d(a).nonzero()[0]
    diff_arr = a[diff_indexes]
    # choose candidate change point
    stat, diff_i = cp_func(diff_arr)
    if diff_i < 0:
        return cps
    i = diff_indexes[diff_i]
    if i <= 1:
        return cps
    # mann-whitney-u statistic to assess significance
    U, p = mwu(diff_arr[:diff_i], diff_arr[diff_i:])
    if p >= pval:
        return cps
    # fold change is ratio of smaller mean to larger mean
    m1 = a[:i].mean()
    m2 = a[i:].mean()
    if m1 > m2:
        m1, m2 = m2, m1
    fc = m1 / m2
    if fc >= fc_cutoff:
        return cps
    # significant change point found - compute interval of slope change
    j, k, sign = slope_extract(s_a, offset + i)
    # ensure slope interval length is nonzero
    if j == 0 or k == 0:
        return cps
    # save changepoint
    cps.append(ChangePoint(pos=offset+i, start=offset+i-j, end=offset+i+k,
                           pvalue=p, sign=sign, foldchange=fc))
    # test left segment
    if (offset+i-j) > offset:
        b1 = a[:i-j]
        cps = bin_seg_slope(b1, s_a, pval, fc_cutoff, size_cutoff,
                            cp_func, cps=cps, offset=offset)
    # test right segment
    if (offset+i+k) < offset+len(a):
        b2 = a[i+k:]
        cps = bin_seg_slope(b2, s_a, pval, fc_cutoff, size_cutoff,
                            cp_func, cps=cps, offset=(offset + i + k))
    return cps


def smooth(x, window_len=11, window='hanning'):
    """
    http://scipy-cookbook.readthedocs.org/items/SignalSmooth.html

    smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd
                    integer
        window: the type of window from 'flat', 'hanning', 'hamming',
                'bartlett', 'blackman', flat window will produce a moving
                average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    TODO: the window parameter could be the window itself if an array instead
          of a string
    NOTE: length(output) != length(input), to correct this:
          return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len < 3:
        return x
    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is not one of 'flat', 'hanning', 'hamming', "
                         "'bartlett', 'blackman'")

    s = np.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]
    # moving average
    if window == 'flat':
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.'+window+'(window_len)')

    y = np.convolve(w/w.sum(), s, mode='valid')
    return y[(window_len/2-1):-(window_len/2)]


def run_changepoint(a, pval=0.05, fc_cutoff=0.80, size_cutoff=20,
                    cp_func=mse_cython, smooth_window="hanning",
                    smooth_window_len=11):
    '''
    Detects change points in 1D array 'a' by minimizing mean-squared error
    (MSE). Selected change points must have mannwhitneyu pvalue < 'pval' and
    have a relative fold change of at least 'fc_cutoff'.

    a: array with expression data
    pval: mann-whitney-u test critical value for selecting significant change
          points
    fc_cutoff: fold-change cutoff
    size_cutoff: minimum length of interval to allow recursive change point
                 searches
    cp_func: distance function to compute change point location
    smooth_window: numpy smoothing window type (must be 'flat', 'hanning',
                   'hamming', 'bartlett', 'blackman')
    smooth_window_len: size of smoothing window

    returns list of ChangePoint namedtuple objects
    '''
    if a.shape[0] < smooth_window_len:
        return []
    # get smoothing window
    s_a = np.gradient(smooth(a, window_len=smooth_window_len,
                             window=smooth_window))
    cps = bin_seg_slope(a, s_a, pval=pval, fc_cutoff=fc_cutoff,
                        size_cutoff=size_cutoff,
                        cp_func=cp_func)
    return cps


def find_threshold_points(a, start=0, threshold=0):
    '''
    a - numpy array of expression data (all values >= 0)
    start - genomic start of 'a'
    threshold - expression threshold to detect changes
    '''
    if a.shape[0] == 0:
        return []
    if a.shape[0] == 1:
        return []
    change_points = []
    for i in xrange(1, a.shape[0]):
        if (a[i-1] > threshold) and (a[i] <= threshold):
            change_points.append(start + i)
        elif (a[i-1] <= threshold) and (a[i] > threshold):
            change_points.append(start + i)
    return change_points


def permute(a, nperms=10):
    d = np.ediff1d(a)
    positions = np.arange(a.shape[0])
    start_values = d[d > 0]
    stop_values = d[d < 0]

    e = np.zeros((nperms, a.shape[0]))
    for i in xrange(nperms):
        print i
        e[i, :] += a[0]
        start_pos = np.random.choice(positions, size=len(start_values), replace=False)
        for j in xrange(len(start_values)):
            e[i, :start_pos[j]] += start_values[j]
        stop_pos = np.random.choice(positions, size=len(stop_values), replace=False)
        for j in xrange(len(stop_values)):
            e[i, stop_pos[j]:] += stop_values[j]

    return e
