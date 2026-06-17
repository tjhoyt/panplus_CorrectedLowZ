#!/usr/bin/env python

import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.stats import bootstrap, gaussian_kde
# TODO: masked=False arg in sigma_clip instead of current two-liner.
try:
    from astropy.stats import sigma_clip
except:
    print('astopy.stats.sigma_clip() not found. Cannot initialize calc_median,mean,dispersion()')
try:
    import statsmodels.api as sm
except:
    print('No statsmodels module installed. Continuing regardless.')
    pass
def __nan_avg_old(x, xerr):
    '''
    Calculate error on weighted avg of x with uncertainties xerr.
    
    This uses the formula sigma_avg = 1/sqrt(sum(xerr))

    Returns
    -------
    wavg

    avgerr
    '''
    # trim NaN; notify user
    N1 = len(x)
    x    = x[~np.isnan(x)]
    xerr = xerr[~np.isnan(xerr)]
    N2 = len(x)
    print('Calculating weighted average and error for {0:d} non-NaN points '
          'out of {1:d} total points.'.format(N2, N1))
    
    # do calculation
    invVars = 1/xerr**2
    w = invVars/np.sum(invVars)
    wavg   = np.sum(x*w)
    avgerr = 1/np.sqrt(np.sum(invVars))/np.sqrt(N2)
    return wavg, avgerr


### Convenience wrapper for NMAD ###
from scipy.stats import median_abs_deviation
def nmad(x):
    return median_abs_deviation(x, scale='normal')

def nan_avg(x, xerr, std_err_method='sum_inverse_square_errors', verbose=False):
    '''
    Calculate error on weighted avg of x with uncertainties xerr.
 
    Std. Err. Methods:
    ------------------
    sum_inverse_variance   : Sum of inverse variances method. Straightforward.
    
    observed_variance_rootn: Simply the weighted variance of the input data array divided
                             by RootN, but with some actual statistical rigor presented in:
                             https://seismo.berkeley.edu/~kirchner/Toolkits/Toolkit_12.pdf.
                             This should agree exactly with sigma/sqrt(n) in the limit of equal variances.

    Returns
    -------
    wavg

    avgerr
    '''
    # trim NaN; notify user
    N1 = len(x)
    nanmask = np.array(np.isnan(x) | np.isnan(xerr))
    x    = x[~nanmask]
    xerr = xerr[~nanmask]
    N2 = len(x)
    if verbose:
        print('Calculating weighted average and error for {0:d} non-NaN points '
              'out of {1:d} total points.'.format(N2, N1))

    # Return nada if no valid data in the array after dropping NaN values
    if N2 == 0:
        if verbose:
            print('No valid numbers in array. Returning (NaN,NaN).')
        return (np.nan, np.nan)
    
    # Calculate inverse variances and normalized weights
    n = N2
    invVars = 1/xerr**2
    w = invVars/np.sum(invVars)
    
    # Compute weighted average
    xAvg = np.sum(w*x)/np.sum(w)
    
    # Compute std_err for sum of inverse variance
    if std_err_method == 'sum_inverse_square_errors':
        xStderr = np.sqrt(1/np.sum(invVars))
    elif std_err_method == 'observed_variance_rootn':
        sample_var = (np.sum(w*(x-xAvg)**2)/np.sum(w)) * (n/(n-1))
        if xVar < 0:
            print('negative variance on Weighted average. Probable bug in code')
        xStderr = np.sqrt(xVar/n)
    else:
        print('Not a supported std_err_method, defaulting to sum of inverse square errors')
        xStderr = np.sqrt(1/np.sum(invVars))
    return xAvg, xStderr
    
def std_err(x):
    '''Calculate standard error as simply sigma/sqrt(N) where N is
    the length of x'''
    n = len(x)
    sigma = np.std(x, ddof=1)
    stderr = sigma/np.sqrt(n)
    return stderr

def quad(arr, axis=0):
    return np.sqrt(np.sum(np.array(arr)**2,axis=axis))

def sample_weighted_avg(mus, sigs, ns):
    weights = ns/np.sum(ns) * 1/sigs**2
    return weights*mus

def combine_normal_samples(mu1, mu2, sig1, sig2, n1, n2):
    '''
    FULL CREDIT: Mothas on SE question #94685
    This can frequently return NaN values. pooled_variance() is preferred.
    '''
    var = (n1*sig1**2/(n1+n2) + n2*sig2**2/(n1+n2)) +  (n1*mu1**2/(n1+n2) + n2*mu2**2/(n1+n2) - np.mean([mu1,mu2])**2)
    return np.sqrt(var)

def pooled_variance(mu1, mu2, sig1, sig2, n1, n2):
    '''
    FULL CREDIT: Ben on SE question #43159
    '''
    var = (n1+n2-1)**-1 *( (n1-1)*sig1**2 + (n2-1)*sig2**2 + (n1*n2)*(mu1-mu2)**2/(n1+n2) )
    return np.sqrt(var)

def simple_regplot(
    x, y, n_std=2, n_pts=100, ax=None, scatter_kws=None, line_kws=None, ci_kws=None
):
    '''Draw a regression line with error interval.
    *FULL CREDIT*: Legendre17 (SO question #22852244)
    '''
    ax = plt.gca() if ax is None else ax

    # calculate best-fit line and interval
    x_fit = sm.add_constant(x)
    fit_results = sm.OLS(y, x_fit).fit()

    eval_x = sm.add_constant(np.linspace(np.min(x), np.max(x), n_pts))
    pred = fit_results.get_prediction(eval_x)

    # draw the fit line and error interval
    ci_kws = {} if ci_kws is None else ci_kws
    ax.fill_between(
        eval_x[:, 1],
        pred.predicted_mean - n_std * pred.se_mean,
        pred.predicted_mean + n_std * pred.se_mean,
        **ci_kws,
    )
    line_kws = {} if line_kws is None else line_kws
    h = ax.plot(eval_x[:, 1], pred.predicted_mean, **line_kws)

    # draw the scatterplot
    scatter_kws = {} if scatter_kws is None else scatter_kws
    ax.scatter(x, y, c=h[0].get_color(), **scatter_kws)

    return fit_results


def hwhm_to_sigma(hw):
    return hw/1.1775

def fit_slope_leastsq(x, y, yerr=None):
    '''
    use the statsmodels api to do an unweighted leastsq fit

    automatically drop nan values

    if y errors are supplied, do a weighted fit
    '''
    X = sm.add_constant(x)
    if yerr is None:
        model = sm.OLS(y, X, missing='drop')
    else:
        model = sm.WLS(y, X, weights=1/yerr**2, missing='drop')
    result = model.fit()
    return result


def calc_median(a, clip=True, sigclip=5, vectorized=True, nresample=500, scipy_method='percentile',
               confidence_level=0.95):
    '''
    compute the median of input array and bootstrapped error on that median.

    sigma clip at the sigclip threshold if clip is True.

    Returns the stdev of the (clipped) residuals, the bootstrap-computed uncertainty on that,
    as well as a conveniently formatted string to print the results, if desired.
    '''
    if clip:
        a = sigma_clip(a, sigma=sigclip)
        # drop masked values instead of masking
        a = a.data[~a.mask]
    else:
        pass

    N = len(a)
    medi = np.median(a)
    if N < 2:
        bs = None
        stderr = np.nan
    else:
        bs = bootstrap((a,), np.median, method=scipy_method, vectorized=vectorized,
                   n_resamples=nresample, confidence_level=confidence_level)
        stderr = bs.standard_error
    statstr = '${0:6.3f} \pm {1:0.3f}$\n$N={2:d}$'.format(medi, stderr, N)
    return medi, bs, N, statstr

def calc_mean(a, clip=True, sigclip=5, vectorized=True, nresample=500, scipy_method='percentile',
               confidence_level=0.95):

    '''
    compute the mean of input array and bootstrapped error on that mean.

    sigma clip at the sigclip threshold if clip is True.

    Returns the stdev of the (clipped) residuals, the bootstrap-computed uncertainty on that,
    as well as a conveniently formatted string to print the results, if desired.
    '''
    if clip:
        a = sigma_clip(a, sigma=sigclip)
        # drop masked values instead of masking
        a = a.data[~a.mask]
    else:
        pass

    N = len(a)
    mean = np.mean(a)
    if N < 2:
        bs = None
        stderr = np.nan
    else:
        bs = bootstrap((a,), np.mean, method=scipy_method, vectorized=vectorized,
                   n_resamples=nresample, confidence_level=confidence_level)
        stderr = bs.standard_error
    statstr = '${0:6.3f} \pm {1:0.3f}$\n$N={2:d}$'.format(mean, stderr, N)
    return mean, bs, N, statstr

def calc_stdev(a, weights=None, clip=True, sigclip=5, vectorized=True, nresample=500, scipy_method='percentile',
               confidence_level=0.95):
    '''
    compute the dispersion of input array and error on that dispersion.

    sigma clip at the sigclip threshold if clip is True.

    Returns the stdev of the (clipped) residuals, the bootstrap-computed uncertainty on that,
    as well as a conveniently formatted string to print the results, if desired.
    '''
    if clip:
        a = sigma_clip(a, sigma=sigclip)
        # drop masked values instead of masking
        a = a.data[~a.mask]
    else:
        pass

    N = len(a)
    if N < 2:
        std = np.nan
        bs = None
        stderr = np.nan
    else:
        std = np.std(a)
        bs = bootstrap((a,), np.std, method=scipy_method, vectorized=vectorized,
                   n_resamples=nresample, confidence_level=confidence_level)
        stderr = bs.standard_error
    statstr = '${0:6.3f} \pm {1:0.3f}$\n$N={2:d}$'.format(std, stderr, N)
    return std, bs, N, statstr

def robust_wavg(a, weights, vectorized=True, nresample=5000, scipy_method='percentile', confidence_level=0.95):
    '''
    The incorporation of a custom stats function into scipy's bootstrap from J-J-J and slightly modified
    https://stackoverflow.com/questions/75921938/
    '''
    N = len(a)
    # Assign unique index to each ith ordered pair of the a and weights
    dict_x_w =  {ind:(ai,wi) for ind, (ai, wi) in enumerate(zip(a,weights))}
    indices = np.array(list(dict_x_w.keys()))
    def weighted_mean(z, axis=0):
        #creating an array of weights, by mapping z to dict_x_w    
        data = np.vectorize(dict_x_w.get)(z)
        return np.average(data[0], weights=data[1],axis=axis)
    bs = bootstrap( (indices,), weighted_mean, method=scipy_method, vectorized=vectorized, 
                   n_resamples=nresample, confidence_level=confidence_level)
    wavg = weighted_mean(indices)
    stderr = bs.standard_error
    statstr = f'${wavg:6.3f} \pm {stderr:0.3f}$\n$N={N:d}$'
    return wavg, bs, N, statstr

def shuffle_array(x, nbin=None):
    if nbin is None:
        nbin = len(x)
    else:
        nbin = nbin

    initial_inds = np.linspace(0, len(x)-1, len(x), dtype=int)
    initial_inds.sort()
    # compute the shuffled indices
    shuffled_inds = np.random.choice(initial_inds, size=nbin, replace=False)
    # permute redshift and scale mu according to the prediction
    x_shuf = x[shuffled_inds]
    return x_shuf

def flip_sign(x):
    sign = np.random.uniform(-1, 1, size=len(x))
    sign_flips = np.where(sign<0, -1, 1)
    return x * sign_flips

def bootstrap_array(x):
    # Define original array indices based on length
    initial_inds = np.linspace(0, len(x)-1, len(x), dtype=int)
    # Resample them with replacement (i.e. bootstrap)
    resampled_inds = np.random.choice(initial_inds, size=len(x), replace=True)
    return x[resampled_inds], resampled_inds

def subsample_inds(a, frac, Nsub=None, aindex=None):
    '''
    generate list of indices that form a subset of a with length f * len(a).

    if aindex is None, then a is assumed to be indexed from 0 to len(a)-1.

    Nsub overrides kwarg frac if specified
    '''
    if aindex is None:
        aindex = range(len(a))

    if Nsub is not None:
        num_rows_to_select = Nsub
    else:
        # Calculate the number of rows to select
        num_rows_to_select = int(len(a) * frac)

    # Randomly select 'num_rows_to_select' rows
    random_indices = random.sample(aindex, num_rows_to_select)
    
    return random_indices

def mad_to_sig(mad):
    return 1.4826 * mad

def gaussian_model(sig):
    '''
    return gaussian function evaluatable about x
    '''
    def gaussian(x):
        return (
            1.0 / (np.sqrt(2.0 * np.pi) * sig) * np.exp(-np.power(x / sig, 2.0) / 2)
        )
    return gaussian

def gaussian_kde_stats(x):
    '''
    Use scipy's gaussian_kde to estimate statistics of dataset
    '''

    gk = gaussian_kde(x)
    std = np.sqrt(gk.covariance/gk.factor)

    return gk, std

def compute_CI(x, interv_size, return_latex_str=False):
    half_size = interv_size/2
    qvals = [50-half_size, 50, 50+half_size]
    #interv = [np.percentile(x, q=q) for q in qvals]
    interv = np.percentile(x, q=qvals)
    if return_latex_str:
        left, right = interv[1] - interv[0], interv[2] - interv[1]
        s = f'{interv[1]:0.3f}_{{-{left:0.3f}}}^{{+{right:0.3f}}}'
        return interv, s
    return interv 

def weighted_avg_and_std(values, weights):
    """
    From Eric O. Lebigot: https://stackoverflow.com/questions/2413522/weighted-standard-deviation-in-numpy
    Return the weighted average and standard deviation.

    They weights are in effect first normalized so that they
    sum to 1 (and so they must not all be 0).

    values, weights -- NumPy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))

def clipped_percentiles(x, pctiles=[16,50,84], sig2clip=5):
    '''
    Compute percentiles after sigma-clipping the input data array.

    Uses astropy's sigma_clip() function.

    Defaults to the 16-50-84 percentiles (68% + median)
    '''
    clipped_array = sigma_clip(x, masked=False, sigma=sig2clip)
    return list(np.percentile(clipped_array, q=pctiles)) + [len(clipped_array)]

def compute_all_BS_stats(x, sig2clip=5, CI=0.68):
    '''
    Convenient wrapper of the three bootstrapped statistics defined above. 

    Returns as a 3x3 array the (central value, low BS bound, high BS bound) for each of (mean, median, stdev)
    '''
    stats = []
    for stat_func in [calc_mean, calc_median, calc_stdev]:
        stat_bs_results = stat_func(x, sigclip=sig2clip, confidence_level=CI)
        val = stat_bs_results[0]
        lo = stat_bs_results[1].confidence_interval.low
        hi = stat_bs_results[1].confidence_interval.high
        N = stat_bs_results[2]
        stats.append([val, lo, hi, N])
    return np.array(stats)

