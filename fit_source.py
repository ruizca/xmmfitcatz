#!/usr/bin/env python
"""
Created on Wed Nov 22 12:15:29 2017

@author: ruizca
"""

import logging
import os
import time
import json
import argparse
import warnings

from itertools import product, ifilter, izip, repeat

#import matplotlib
#matplotlib.use('AGG') 
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import numpy as np
import scipy.stats
import pymultinest

from astropy.io import fits
from astropy.cosmology import Planck15 as cosmo
from astropy import units as u
from astroML.density_estimation import bayesian_blocks
from scipy.interpolate import interp1d
from sklearn.model_selection import train_test_split

from sherpa.astro.ui import *
import bxa.sherpa as bxa


from models import *
from priors import *

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes.formatter', limits=(-3,3))
plt.style.use('bmh')
plt.ioff()

def mode(inputData, axis=None, dtype=None):
    """
    Robust estimator of the mode of a data set using the half-sample mode.
    Based on functions from Henry Freudenriech (Hughes STX) statistics 
    library (called ROBLIB) from AstroIDL User's Library.

     https://fornax.phys.unm.edu/lwa/subversion/trunk/lsl/lsl/statistics/robust.py
    .. versionadded: 1.0.3
    """
	
    if axis is not None:
        fnc = lambda x: mode(x, dtype=dtype)
        dataMode = np.apply_along_axis(fnc, axis, inputData)
    else:
        # Create the function that we can use for the half-sample mode
        def _hsm2(data):
            if len(data) < 3:
                return data.mean()

            # Round up to include the middle value, in the case of an odd-length array
            half_idx = int((len(data) + 1)/2) 

            # Calculate all interesting ranges that span half of all data points
            ranges = data[-half_idx:] - data[:half_idx]
            smallest_range_idx = np.argmin(ranges)

            # Now repeat the procedure on the half that spans the smallest range
            data_subset = data[smallest_range_idx : (smallest_range_idx + half_idx)]

            return _hsm2(data_subset)
				
        data = inputData.ravel()
        if type(data).__name__ == "MaskedArray":
            data = data.compressed()
        if dtype is not None:
            data = data.astype(dtype)
			
        # Sort and remove repeated values 
        #(the algorithm doesn't work with unsorted or repeated values)
        data = np.unique(data)
		
        # Find the mode
        dataMode = _hsm2(data)

    return dataMode

def bootstrap(data, num_samples, alpha):
    """
    Returns the results from num_samples bootstrap samples for an input test 
    statistic and its 100*(1-alpha)% confidence level interval.

    From https://github.com/ilanfri/StatsML/blob/master/Bayesian_bootstrap.ipynb
    """
    # Generate the indices for the required number of permutations/(resamplings 
    # with replacement) required
    idx = np.random.randint(0, len(data), (num_samples, len(data)))

    # Generate the multiple resampled data set from the original one
    samples = data[idx]

    # Estimate the mode for each of the data sets produced by 
    # the resampling and order the resulting statistic by decreasing size.
    stats = np.sort(mode(samples, axis=1))
    stat = stats.mean()

    p95 = np.percentile(samples, 95, axis=1)
    p95 = p95.mean()

    p05 = np.percentile(samples, 5, axis=1)
    p05 = p05.mean()

    # Estimate the median for the resamplings (no sorting because the 
    # CI is estimated for the mode)
    stats2 = np.median(samples, axis=1)
    stat2 = stats.mean()

    # Return the value of the computed statistic at the upper and lower 
    # percentiles specified by the alpha parameter given. These are, by 
    # definition, the boundaries of the Confidence Interval for that value of 
    # alpha. E.g. alpha=0.05 ==> CI 95%
    low_ci = stats[int((alpha / 2.0) * num_samples)]
    high_ci = stats[int((1 - alpha / 2.0) * num_samples)]

    return stat, stat2, p05, p95
 
def update_rmf(specfile, rmf_folder):
    """
    Read specfile (XMM spectra in fits format) and returns the corresponding rmf file
    defined in the spectrum header
    """
    #epoch_limits
    d1 = time.strptime("2007-01-01", "%Y-%m-%d")
    d2 = time.strptime("2014-01-01", "%Y-%m-%d")

    hdulist = fits.open(specfile)
    obs_date = hdulist[0].header['DATE-OBS']
    rmf_old = hdulist[1].header['RESPFILE']
    camera = hdulist[1].header['INSTRUME']
    deltaE = hdulist[1].header['SPECDELT']
    
    hdulist.close()
    
    if camera == 'EPN':
        obs_date = time.strptime(obs_date, "%Y-%m-%dT%H:%M:%S")
        
        if obs_date < d1:
            epoch = "e1"
            
        elif obs_date >= d1 and obs_date < d2:
            epoch = "e2"
                    
        else:
            epoch = "e3"

        rmf_split = rmf_old.split('_',1)
        if (rmf_split[1].startswith("e1") or 
            rmf_split[1].startswith("e2")  or 
            rmf_split[1].startswith("e3")):

            rmf_new = "{0}{1}_{2}".format(rmf_folder, rmf_split[0], 
                                          rmf_split[1].replace(".rmf", "_v16.0.rmf"))

        else:
            rmf_new = "{0}{1}_{2}_{3}".format(rmf_folder, rmf_split[0], epoch,
                                              rmf_split[1].replace(".rmf", "_v16.0.rmf"))            
    else:
        if deltaE == 5 :
            rmf_new = "{0}5ev/{1}".format(rmf_folder, rmf_old)
        else :
            rmf_new = "{0}{1}".format(rmf_folder, rmf_old)

    return rmf_new, camera


def calc_bkg_model(sid) :
    """
    Estimate the background model for the background spectra of the data set sid,
    as used in the wstat statistics .
    [see https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixStatistics.html]
    """
    
    notice_id(sid)
    S = get_data(sid).counts
    B = get_bkg(sid).counts
    m = get_model_plot(sid).y * (get_model_plot(sid).xhi - get_model_plot(sid).xlo)
    texps = get_data(sid).exposure * get_areascal(sid)
    texpb = get_bkg(sid).exposure * get_areascal(sid, 1) / get_bkg_scale(sid)

    T = texps + texpb
    d = np.sqrt((T*m - S - B)**2 + 4*T*B*m)
    f = (S + B - T*m + d)/(2*T)
    
    return f


def calc_srcbkg_model(sid) :
    
    f = calc_bkg_model(sid)
    m = get_model_plot(sid).y
    de = get_model_plot(sid).xhi - get_model_plot(sid).xlo
    
    return m + f/de


def ks_goodness(ids, emin=0, emax=np.inf) :
    """
    KS test between data and model. If the null hypothesis is not rejected with
    high probability, the fit is good.
    """

    model = []
    data = []
    for sid in ids :
        notice_id(sid)
        texps = get_data(sid).exposure * get_areascal(sid)
        e = (get_model_plot(sid).xhi + get_model_plot(sid).xlo)/2
        de = get_model_plot(sid).xhi - get_model_plot(sid).xlo
        mask = np.logical_and(e >= emin, e <= emax)
      
        sb_ctr = calc_srcbkg_model(sid)
        model = np.concatenate((model, texps*sb_ctr[mask]*de[mask]))
        data = np.concatenate((data, get_data(sid).counts[mask]))
        
        ignore_id(sid, ":{:f},{:f}:".format(emin, emax))
        
    ks = ks_2samp(data, model)
    en = np.sqrt(1.0*len(data)*len(model)/(len(data) + len(model))) 
    prob = scipy.stats.distributions.kstwobign.sf((en + 0.12 + 0.11/en) * ks)

    return ks, prob


def ks_2samp(sample1, sample2):
    """
    Kolmogorov-Smirnov test for two samples.
    """
    cdf1 = sample1.cumsum()/sample1.sum()
    cdf2 = sample2.cumsum()/sample2.sum()

    return np.abs(cdf1 - cdf2).max()


def ad_2samp(sample1, sample2):
    """
    Implementation of the Anderson-Darling test for two samples.

    See Scholz, F. W and Stephens, M. A. (1987), K-Sample Anderson-Darling Tests, 
    Journal of the American Statistical Association, Vol. 82, pp. 918-924    
    """
    cdf1 = sample1.cumsum()/sample1.sum()
    cdf2 = sample2.cumsum()/sample2.sum()

    m = len(sample1)
    n = len(sample2)
    N = m + n
    HN = (m*cdf1 + n*cdf2)/N

    msk = np.logical_and(HN>0, HN<1)    
    I = (cdf1[msk] - cdf2[msk])**2/(HN[msk]*(1 - HN[msk]))
    dHN = np.diff(HN[msk])
    
    return np.sum(I[:-1]*dHN)*(m*n/N)


def goodness(ids, emin=0, emax=np.inf, niter=1000):
    """
    KS and Anderson-Darling test between data and model. If the null hypothesis 
    is not rejected with high probability, the fit is good. 

    Note that in this kind of tests the parameters of the model should be
    independent of the data. If they are estimated by fitting the data,
    the probabilities of the KS test (and other non-paramatric test based on 
    comparison of the empirical cumulative distributions) are wrong!!!
    Fortunately, we can estimate the p-value using non-parametric bootstrap 
    resampling. 

    To do this we used a permutation test, spliting the data+model
    sample in two equal size sample and estimating the statistic.
    Doing this N times, we can estimate the p-value as how often the
    statistic of the permutation is larger than the observed value in the
    original data/model test. 

    See https://asaip.psu.edu/Articles/beware-the-kolmogorov-smirnov-test
    and Babu, G. J.  & Rao, C. R. (2004) and
    https://stats.stackexchange.com/questions/59774/test-whether-variables-follow-the-same-distribution/59875#59875
    """

    # Obtain data and model counts
    model = []
    data = []
    for sid in ids:
        notice_id(sid)
        texps = get_data(sid).exposure * get_areascal(sid)
        e = (get_model_plot(sid).xhi + get_model_plot(sid).xlo)/2
        de = get_model_plot(sid).xhi - get_model_plot(sid).xlo
        mask = np.logical_and(e >= emin, e <= emax)
      
        sb_ctr = calc_srcbkg_model(sid)
        model = np.concatenate((model, texps*sb_ctr[mask]*de[mask]))
        data = np.concatenate((data, get_data(sid).counts[mask]))
        
        ignore_id(sid, ":{:f},{:f}:".format(emin, emax))

    # Estimate KS and AD stats
    ks = ks_2samp(data, model)
    ad = ad_2samp(data, model)
    
    # Estimate the p-values by bootstraping
    count = np.array([0, 0])
    bs_data = np.concatenate((data,model))

    for _ in repeat(None, niter):
        perm1, perm2 = train_test_split(bs_data, test_size=len(data))
        bs_ks = ks_2samp(perm1, perm2)
        bs_ad = ad_2samp(perm1, perm2)

        count[0] += (bs_ks >= ks)
        count[1] += (bs_ad >= ad)
        
    pvals = 1.*count/niter

    return ks, ad, pvals[0], pvals[1]

def goodness2(ids, emin=0, emax=np.inf, niter=1000):
    """
    KS and Anderson-Darling test between data and model. If the null hypothesis 
    is not rejected with high probability, the fit is good. 

    Note that in this kind of tests the parameters of the model should be
    independent of the data. If they are estimated by fitting the data,
    the probabilities of the KS test (and other non-paramatric test based on 
    comparison of the empirical cumulative distributions) are wrong!!!
    Fortunately, we can estimate the p-value using non-parametric bootstrap 
    resampling. 

    To do this we used a permutation test, spliting the data+model
    sample in two equal size sample and estimating the statistic.
    Doing this N times, we can estimate the p-value as how often the
    statistic of the permutation is larger than the observed value in the
    original data/model test. 

    See https://asaip.psu.edu/Articles/beware-the-kolmogorov-smirnov-test
    and Babu, G. J.  & Rao, C. R. (2004) and
    https://stats.stackexchange.com/questions/59774/test-whether-variables-follow-the-same-distribution/59875#59875
    """

    
    ks_all = []
    ad_all = []
    counts_all = []
    pvals_total = np.array([1,1])
    
    for sid in ids:
        # Obtain data and model counts
        notice_id(sid)
        texps = get_data(sid).exposure * get_areascal(sid)
        e = (get_model_plot(sid).xhi + get_model_plot(sid).xlo)/2
        de = get_model_plot(sid).xhi - get_model_plot(sid).xlo
        mask = np.logical_and(e >= emin, e <= emax)
      
        sb_ctr = calc_srcbkg_model(sid)
        model = texps*sb_ctr[mask]*de[mask]
        data = get_data(sid).counts[mask]
        
        ignore_id(sid, ":{:f},{:f}:".format(emin, emax))

        # Estimate KS and AD stats
        ks = ks_2samp(data, model)
        ad = ad_2samp(data, model)
        
        ks_all = np.concatenate((ks_all, [ks]))
        ad_all = np.concatenate((ad_all, [ad]))
        counts_all = np.concatenate((counts_all, [data.sum()]))
    
        # Estimate the p-values by bootstraping
        count = np.array([0, 0])
        bs_data = np.concatenate((data,model))

        for _ in repeat(None, niter):
            perm1, perm2 = train_test_split(bs_data, test_size=len(data))
            bs_ks = ks_2samp(perm1, perm2)
            bs_ad = ad_2samp(perm1, perm2)
    
            if bs_ks > ks:
                c = 1.0
            elif bs_ks == ks:
                c = 0.5
            else :
                c = 0.0

            count[0] += c #(bs_ks >= ks)
            count[1] += (bs_ad >= ad)
            
        pvals = 1.*count/niter
        pvals_total = pvals_total*pvals
        print pvals, pvals_total

    ks_final = np.max(ks_all)
    ad_final = np.average(ad_all, weights=counts_all)

    print ks_final, ad_final, pvals_total[0], pvals_total[1]
    
    return ks_final, ad_final, pvals_total[0], pvals_total[1]



def ticks_format(value, index):
    """
    get the value and returns the value as:
       integer: [0,99]
       1 digit float: [0.1, 0.99]
       n*10^m: otherwise
    To have all the number of the same size they are all returned as latex strings
    """
    exp = np.floor(np.log10(value))
    base = value/10**exp
    if exp == 0 or exp == 1:   
        return '${0:d}$'.format(int(value))
    if exp == -1 :
        return '${0:.1f}$'.format(value)
    else:
        if base == 1 :
            return '$10^{{{0:d}}}$'.format(int(exp))
        else :
            return '${0:d}\\times10^{{{1:d}}}$'.format(int(base), int(exp))
        
        
def my_plot_bkg_spec(ids, plot_folder, emin=0, emax=np.inf, save=True, label_model="") :

    for sid in ids :
        notice_id(sid)
        e = (get_model_plot(sid).xhi + get_model_plot(sid).xlo)/2
        texpb = get_bkg(sid, 1).exposure * get_areascal(sid, 1) / get_bkg_scale(sid)
        bkg_data = get_bkg(sid, 1).counts
        bkg_model = texpb*calc_bkg_model(sid)
        mask = np.logical_and(e >= emin, e <= emax)

        fig = plt.figure()
        plt.scatter(e[mask], bkg_data[mask], s=5)
        plt.loglog(e[mask], bkg_model[mask], c='red', alpha=0.5)
        plt.xlabel("Energy / keV")
        plt.ylabel("counts")
        
        ax = plt.gca()        
        ax.xaxis.set_major_locator(matplotlib.ticker.LogLocator(subs=(1.0, 2.0, 5.0)))
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticks_format))
        
        plt.tight_layout()
        
        if save :
            label = get_bkg(sid, 1).name.split('/')[-1]
            fig.savefig("{0}{1}_bkgspec_{2}.png".format(
                        plot_folder, label_model, label[11:17]))
        else :
            fig.show()

        plt.close(fig)


def my_plot_fit(ids, plot_folder, emin=0, emax=np.inf, save=True, label_model="", z=None) :

    all_model = []
    all_emodel = []
    all_data = []
    all_dataxerr = []    
    all_datayerr = []
    all_edata = []
    all_ratio = []
    all_ratioerr = []

    # One plot per spectrum        
    for sid in ids :
        group_snr(sid, 3)
        d = get_data_plot(sid)
        
        e = (get_model_plot(sid).xhi + get_model_plot(sid).xlo)/2
        model = calc_srcbkg_model(sid)
        model_de = model * (get_model_plot(sid).xhi - get_model_plot(sid).xlo)

        bins = np.concatenate((d.x - d.xerr/2, [d.x[-1]+d.xerr[-1]]))
        model_binned, foo1, foo2 = scipy.stats.binned_statistic(e, model_de, bins=bins, statistic='sum')
        model_binned = model_binned/d.xerr

        #delchi = resid/d.yerr
        ratio = d.y/model_binned
        
        mask_data = np.logical_and(d.x+d.xerr/2 >= emin, d.x-d.xerr/2 <= emax)
        mask_model = np.logical_and(e >= emin, e <= emax)

        ymin = np.min(d.y[mask_data])
        ylim = 10**np.floor(np.log10(ymin))        
                
        fig = plt.figure()
        gs = matplotlib.gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 

        # Data and model
        ax0 = plt.subplot(gs[0])                
        ax0.errorbar(d.x[mask_data], d.y[mask_data], 
                     xerr=d.xerr[mask_data]/2, yerr=d.yerr[mask_data], 
                     fmt="o", ms=2, elinewidth=0.8, ls="None", zorder=1000)
        ax0.loglog(e[mask_model], model[mask_model], c='red', alpha=0.5)
        ax0.loglog(e[mask_model], get_model_plot(sid).y[mask_model])
                                    
                                    
        ax0.set_ylabel("count rate / $\mathrm{s}^{-1}\:\mathrm{keV}^{-1}$")
        ax0.set_ylim(bottom=ylim)

        # Ratio
        ax1 = plt.subplot(gs[1], sharex = ax0)
        ax1.axhline(1, ls="--", c="gray")
        ax1.errorbar(d.x[mask_data], ratio[mask_data],
                     xerr=d.xerr[mask_data]/2, 
                     yerr=d.yerr[mask_data]/model_binned[mask_data], 
                     elinewidth=0.8, ls="None", zorder=1000)
                     
        plt.setp(ax0.get_xticklabels(), visible=False)

        ax1.set_yscale("log")
        ax1.set_xlabel("Energy / keV")
        ax1.set_ylabel("ratio")
        
        ax1.xaxis.set_major_locator(matplotlib.ticker.LogLocator(subs=(1.0, 2.0, 5.0)))
        ax1.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticks_format))
        ax1.yaxis.set_major_locator(matplotlib.ticker.LogLocator(subs=(1.0, 3.0, 5.0)))        
        ax1.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticks_format))
        ax1.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        
        plt.xlim(0.5,10)        
        plt.tight_layout()

        # Inset for the FeKa line region (only if enough data available)
        if z is not None :
            resid = d.y - model_binned
            msk_inset = np.where(np.logical_and(d.x > 6/(1+z), d.x<7/(1+z)))
            
            if len(msk_inset[0]) > 10 :
                inset_ymax = np.max(resid[msk_inset] + 2*d.yerr[msk_inset])
                inset_ymin = np.min(resid[msk_inset] - 2*d.yerr[msk_inset])

                axins = inset_axes(ax0, width="44%", height="35%", loc="lower left")
                axins.errorbar(d.x[mask_data], resid[mask_data],
                               xerr=d.xerr[mask_data]/2, yerr=d.yerr[mask_data], 
                               elinewidth=0.8, ls="None", zorder=1000)

                axins.axhline(0, ls="--", c="gray")
                axins.axvline(6.4/(1+z), ls="--", c="gray")

                axins.set_xlim(6/(1+z), 7/(1+z))
                axins.set_ylim(inset_ymin, inset_ymax)
                axins.xaxis.tick_top()
                axins.yaxis.tick_right()
                #axins.set_ylabel("resid.")
                #axins.yaxis.set_label_position("right")
                axins.tick_params(axis='both', which='major', labelsize=8)


        if save :
            label = get_data(sid).name.split('/')[-1]
            fig.savefig("{0}{1}_srcspec_{2}.png".format(
                        plot_folder, label_model, label[11:17]))
        else :
            fig.show()

        plt.close(fig)


        if len(ids) > 1 :
            all_model.append(model[mask_model])
            all_emodel.append(e[mask_model])
            all_data.append(d.y[mask_data])
            all_dataxerr.append(d.xerr[mask_data])
            all_datayerr.append(d.yerr[mask_data])
            all_edata.append(d.x[mask_data])
            all_ratio.append(ratio[mask_data])
            all_ratioerr.append(d.yerr[mask_data]/model_binned[mask_data])
        
        group_counts(sid, 1)

    # Show all spectra in one plot
    if len(ids) > 1 :        
        fig = plt.figure()
        gs = matplotlib.gridspec.GridSpec(2, 1, height_ratios=[3, 1])

        ax0 = plt.subplot(gs[0])
        ylim_total=1
        for i in range(len(ids)) :
            ax0.errorbar(all_edata[i], all_data[i], 
                         xerr=all_dataxerr[i]/2, yerr=all_datayerr[i], 
                         fmt="o", ms=2, elinewidth=0.8, ls="None", zorder=1000)
            ax0.loglog(all_emodel[i], all_model[i], c='red', alpha=0.5)
            ymin = np.min(all_data[i])
            ylim = 10**np.floor(np.log10(ymin))        
            if ylim < ylim_total :
                ylim_total = ylim

        ax0.set_ylim(bottom=ylim_total)
        ax0.set_ylabel("count rate / $\mathrm{s}^{-1}\:\mathrm{keV}^{-1}$")
        ax0.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticks_format))
        
        ax1 = plt.subplot(gs[1], sharex = ax0)
        for i in range(len(ids)) :
            ax1.errorbar(all_edata[i], all_ratio[i],
                         xerr=all_dataxerr[i]/2, 
                         yerr=all_ratioerr[i], 
                         elinewidth=0.8, ls="None", zorder=1000)
        
        ax1.axhline(1, ls="--", c="gray")

        plt.setp(ax0.get_xticklabels(), visible=False)

        ax1.set_yscale("log")
        ax1.set_xlabel("Energy / keV")
        ax1.set_ylabel("ratio")
        
        ax1.xaxis.set_major_locator(matplotlib.ticker.LogLocator(subs=(1.0, 2.0, 5.0)))
        ax1.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticks_format))
        ax1.yaxis.set_major_locator(matplotlib.ticker.LogLocator(subs=(1.0, 3.0, 5.0)))        
        ax1.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(ticks_format))
        ax1.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        
        plt.xlim(0.5,10)        
        plt.tight_layout()

        if save :
            fig.savefig("{0}{1}_srcspec_ALL.png".format(plot_folder, label_model))
        else :
            fig.show()
            
        plt.close(fig)


def qq_plot(ids, plot_folder, emin=0, emax=np.inf, save=True, label_model="") :

    top = 0
    fig = plt.figure()
    for sid in ids :
        notice(sid)
        texps = get_data(sid).exposure * get_areascal(sid)
        e = (get_model_plot(sid).xhi + get_model_plot(sid).xlo)/2
        de = get_model_plot(sid).xhi - get_model_plot(sid).xlo
        mask = np.logical_and(e >= emin, e <= emax)

        sb_model = calc_srcbkg_model(sid)
        model = (texps*sb_model[mask]*de[mask])
        data = get_data(sid).counts[mask]

        top = np.max([top, np.max(data.cumsum()), np.max(model.cumsum())])
        plt.plot(data.cumsum(), model.cumsum())        

    plt.plot([0, top], [0, top], ":")
    plt.xlim(0, top)
    plt.ylim(0, top)    

    plt.xlabel("Data (source + background) / counts")
    plt.ylabel("Model (source + background) / counts")
    plt.tight_layout()

    if save :
        fig.savefig("{0}{1}_qqplot.png".format(plot_folder, label_model))        
    else :
        fig.show()

    plt.close(fig)


def plot_marginals(analyzer, parameters, stats, plot_folder, label_model="", 
    save=True, bbahist=False) :
    p = pymultinest.PlotMarginal(analyzer)
    values = analyzer.get_equal_weighted_posterior()
    warnings.filterwarnings("ignore")

    n_params = len(parameters)
    assert n_params == len(stats['marginals'])
    modes = stats['modes']


    fig = plt.figure(figsize=(5*n_params, 5*n_params))
    for i in range(n_params):
        plt.subplot(n_params, n_params, i + 1)
        plt.xlabel(parameters[i].fullname)

        m = stats['marginals'][i]
        plt.xlim(m['5sigma'])

        oldax = plt.gca()

        if bbahist:
            bb_edges = bayesian_blocks(values[:,i], fitness='events', p0=0.01)
            # Find and remove spikes from bayesian blocks (they apear sometimes)
            counts, edges = np.histogram(values[:,i], bins=bb_edges, density=True)
            mcounts = np.median(counts)
            scounts = np.median(np.abs(counts - mcounts))
            bad_bins = 1 + np.where(counts > mcounts + 5*scounts)[0]
            bb_edges = np.delete(edges, bad_bins)

            oldax.hist(values[:,i], bins=bb_edges, histtype='step', 
                       color='red', normed=True, linewidth=1.5)
        
        x, w, patches = oldax.hist(values[:,i], bins='auto', 
                                 edgecolor='black', color='black', 
                                 histtype='stepfilled', alpha=0.2, normed=True)
        
        newax = plt.gcf().add_axes(oldax.get_position(), sharex=oldax, frameon=False)
        p.plot_marginal(i, ls='-', color='blue', linewidth=3)
        newax.set_ylim(0, 1)

        ylim = newax.get_ylim()
        y = ylim[0] + 0.05*(ylim[1] - ylim[0])
        center = m['median']
        low1, high1 = m['1sigma']
        newax.errorbar(x=center, y=y,
                       xerr=np.transpose([[center - low1, high1 - center]]), 
                       color='blue', linewidth=2, marker='s')
        oldax.set_yticks([])
        newax.set_ylabel("Probability")
        ylim = oldax.get_ylim()
        newax.set_xlim(m['5sigma'])
        oldax.set_xlim(m['5sigma'])

        for j in range(i):
            plt.subplot(n_params, n_params, n_params * (j + 1) + i + 1)
            p.plot_conditional(i, j, bins=20, cmap = plt.cm.gray_r)
            for n in modes:
                plt.errorbar(x=n['mean'][i], y=n['mean'][j], 
                             xerr=n['sigma'][i], yerr=n['sigma'][j])
            plt.xlabel(parameters[i].fullname)
            plt.ylabel(parameters[j].fullname)
            #plt.xlim(m['5sigma'])

    fig.tight_layout()
    
    if save :
        fig.savefig("{}{}_param_marginals.png".format(plot_folder, label_model))
    else :
        plt.show()

    plt.close(fig)
    
        
def define_model(model, ids, nhgal, z=None) :

    if model == "po" :
        define_model_po(ids, z, nhgal)
    
    elif model == "wapo" :
        define_model_wapo(ids, z, nhgal)

    elif model == "mekal" :
        define_model_mekal(ids, z, nhgal)
        
    elif model == "wamekal" :
        define_model_wamekal(ids, z, nhgal)

    elif model == "bb" :
        define_model_bb(ids, z, nhgal)
        
    elif model == "wabb" :
        define_model_wabb(ids, z, nhgal)

    elif model == "wamekalpo" :
        define_model_wamekalpo(ids, z, nhgal)

    elif model == "wapopo" :
        define_model_wapopo(ids, z, nhgal)

    elif model == "wabbpo" :
        define_model_wabbpo(ids, z, nhgal)
        
    else :
        raise ValueError("unknown model!!!")


def define_priors(model, ids, zcdf=None) :

    if model == "po" :
        params, priorfun = define_priors_po(ids, zcdf)
    
    elif model == "wapo" :
        params, priorfun = define_priors_wapo(ids, zcdf)

    elif model == "mekal" :
        params, priorfun = define_priors_mekal(ids, zcdf)
        
    elif model == "wamekal" :
        params, priorfun = define_priors_wamekal(ids, zcdf)

    elif model == "bb" :
        params, priorfun = define_priors_bb(ids, zcdf)
        
    elif model == "wabb" :
        params, priorfun = define_priors_wabb(ids, zcdf)

    elif model == "wamekalpo" :
        params, priorfun = define_priors_wamekalpo(ids, zcdf)

    elif model == "wapopo" :
        params, priorfun = define_priors_wapopo(ids, zcdf)

    elif model == "wabbpo" :
        params, priorfun = define_priors_wabbpo(ids, zcdf)
        
    else :
        raise ValueError("unknown model!!!")
        
    return params, priorfun
        

def calc_flux_pars(params, parvals, ids, ebands=[(0.5,10.0)], 
                      z=0, thawedz=False) :

    for p, v in zip(params, parvals):
        p.val = v

        if thawedz and p.fullname.endswith("redshift") :
            z = p.val

    ebands = np.array(ebands)/(1+z)
    
    f = [[calc_energy_flux(id=sid, lo=band[0], hi=band[1]) for sid in ids]
         for band in ebands]

    return np.mean(f, axis=1)

        
def get_distribution_with_fluxes(chains, parameters, ids, ebands=[(0.5,10.0)]) :

    nrows, ncols_chains = np.shape(chains)
    ncols = ncols_chains + len(ebands)

    r = np.full((nrows, ncols), np.nan)

    first_colflux = ncols_chains
    r[:, 0:first_colflux] = chains
    r[:, first_colflux:] = [calc_flux_pars(parameters, row, ids, ebands=ebands)
                            for row in chains]
    
    return r


def get_distribution_with_intfluxes(chains, parameters, ids, ebands=[(0.5,10.0)], 
                                    z=0, abscorr=True, thawedz=False) :
                                        
    if abscorr :
        for p in parameters:
            if p.fullname.startswith("intabs") :
                set_par("{}.nH".format(p.fullname.split(".")[0]), val=0)

    nrows, ncols_chains = np.shape(chains)
    ncols = ncols_chains + len(ebands)
    
    r = np.full((nrows, ncols), np.nan)

    first_colflux = ncols_chains
    r[:, 0:first_colflux] = chains
    r[:, first_colflux:] = [calc_flux_pars(parameters, row, ids, z=z, 
                                           ebands=ebands, thawedz=thawedz)
                            for row in chains]

    return r

def flux_stats(wF, ebands) :

    fmedian = [None]*len(ebands)
    fmode = [None]*len(ebands)
    flimits = np.full((len(ebands), 2), np.nan)

    for i in range(len(ebands)) :
        # group by flux and add all points in the bins within the group
        edges = np.linspace(np.min(wF[:,i+1]), np.max(wF[:,i+1]), num=11)
        wF_binned = np.zeros(len(edges)-1)
        for k in range(len(edges[:-1])):
            msk = np.logical_and(wF[:,i+1] > edges[k], wF[:,i+1] <= edges[k+1])
            wF_binned[k] = np.sum(wF[msk,0])
    
        # estimate cumulative distribution of fluxes
        mid_edges = (edges[:-1]+edges[1:])/2            
        wF_cum = np.cumsum(wF_binned)/np.sum(wF[:,0])
        wF_interp = interp1d(wF_cum, mid_edges, bounds_error=False, fill_value='extrapolate')

        fmedian[i] = np.asscalar(wF_interp(0.5))
        fmode[i] = mid_edges[np.argmax(wF_binned)]
        flimits[:,i] = wF_interp(0.05), wF_interp(0.95)

    return fmode, fmedian, flimits

def calc_fluxes_fast(chain_params, params, ids, ebands=[(0.5,10.0)], 
                     abscorr=True, z=0, thawedz=False):

    # The idea is to estimate the N-dimensional histogram of the parameter 
    # posterior and calculate the flux only for the non-empty bins. Then use 
    # the number of sources in each bin to estimate the flux probability 
    # distribution. It works reasonably ok and reduces significantly the number 
    # of fluxes we need to estimate, compared to using the whole posterior 
    # distribution of the parameters
    #
    # Note: it doesn't work for high dimensions or large counts, the histogram
    # doesn't fit in memory :-(

    fluxmedian = [None]*len(ebands)*3
    fluxmode = [None]*len(ebands)*3
    fluxlimits = np.full((2, len(ebands)*3), np.nan)

    H, edges = np.histogramdd(chain_params, bins=10)
    mid = np.array([(edge[:-1]+edge[1:])/2 for edge in edges])
    idx_nz = np.nonzero(H)    
    pars = np.transpose([row[idx_nz[i]] for i, row in enumerate(mid)]) 

    ## Observed fluxes
    # Estimate flux for the combination of parameters in non-empty bins
    wF = np.array([np.append(h, calc_flux_pars(params, p, ids, ebands=ebands)) 
                   for h, p in izip(H[idx_nz], pars)])

    (fluxmode[:len(ebands)], 
     fluxmedian[:len(ebands)], 
     fluxlimits[:,:len(ebands)]) = flux_stats(wF, ebands)

    ## Intrinsic fluxes
    if abscorr :
        for p in params:
            if p.fullname.startswith("intabs") :
                set_par("{}.nH".format(p.fullname.split(".")[0]), val=0)

    # Estimate flux for the combination of parameters in non-empty bins
    wFint = np.array([np.append(h, calc_flux_pars(params, p, ids, ebands=ebands, 
                                                  z=z, thawedz=thawedz)) 
                      for h, p in izip(H[idx_nz], pars)])

    (fluxmode[len(ebands):2*len(ebands)], 
     fluxmedian[len(ebands):2*len(ebands)], 
     fluxlimits[:,len(ebands):2*len(ebands)]) = flux_stats(wFint, ebands)

    # Luminosities
    if thawedz:
        for j,p in enumerate(params) :
                if p.fullname.endswith("redshift") :
                    z_idx=j
                    break
        dL = cosmo.luminosity_distance(pars[:,z_idx]).to(u.cm).value

    else:
        dL = np.array([cosmo.luminosity_distance(z).to(u.cm).value]*len(wFint[:,0]))

    wL = np.c_[wFint[:,0], 4*np.pi*dL[:,None]**2 * wFint[:,1:]]

    (fluxmode[2*len(ebands):], 
     fluxmedian[2*len(ebands):], 
     fluxlimits[:,2*len(ebands):]) = flux_stats(wL, ebands)

    return fluxmedian, fluxmode, fluxlimits

def save_bestfit_pars(parameters, parmedian, parmode, parlimits, filename) :

    params_dict = {}
    for i,p in enumerate(parameters) :
        parval_dict = {"median" : parmedian[i],
                       "mode" : parmode[i],
                       "parmin" : parlimits[0,i],
                       "parmax" : parlimits[1,i],
                       "errmedianmin" : parmedian[i] - parlimits[0,i],
                       "errmedianmax" : parlimits[1,i] - parmedian[i],
                       "errmodemin" : parmode[i] - parlimits[0,i],
                       "errmodenmax" : parlimits[1,i] - parmode[i]}
                       
        params_dict[p.fullname] = parval_dict

    json.dump(params_dict, file(filename, 'w'), indent=2)


def save_bestfit_fluxes(ebands, parmedian, parmode, parlimits, firstcol, 
                        filename, bs_fluxes) :

    if bs_fluxes:
       method = 'sampling'
    else:
       method = 'full'

    fluxes_dict = {}
    for i,e in enumerate(ebands) :
        flux_dict = {"median" : parmedian[firstcol+i],
                     "mode" : parmode[firstcol+i],
                     "fluxmin" : parlimits[0,firstcol+i],
                     "fluxmax" : parlimits[1,firstcol+i],
                     "errmedianmin" : parmedian[firstcol+i] - parlimits[0,firstcol+i],
                     "errmedianmax" : parlimits[1,firstcol+i] - parmedian[firstcol+i],
                     "errmodemin" : parmode[firstcol+i] - parlimits[0,firstcol+i],
                     "errmodemax" : parlimits[1,firstcol+i] - parmode[firstcol+i],
                     "emin" : e[0],
                     "emax" : e[1],
                     "method" : method}

        fluxes_dict["flux" + str(i)] = flux_dict
    
    json.dump(fluxes_dict, file(filename, 'w'), indent=2)


def save_goodness(ks, p, ad, adp, npar, prefix) :

    stats = get_stat_info()[-1]
    wstat = stats.statval
    dof = stats.numpoints - npar

    with open('{}stats.json'.format(prefix)) as data_file:
        multinest_stats = json.load(data_file)
        
    dict_goodness = {"wstat" : wstat,
                     "dof" : dof,
                     "npar" : npar,
                     "ks" : ks,
                     "pvalue" : p,
                     "ad" : ad,
                     "ad_pvalue" : adp,
                     "lnZ" : multinest_stats["global evidence"],                         
                     "lnZ_err" : multinest_stats["global evidence error"],
                     "lnev" : multinest_stats["nested sampling global log-evidence"],
                     "lnev_err" : multinest_stats["nested sampling global log-evidence error"]}

    json.dump(dict_goodness, file("{}bestfit_goodness.json".format(prefix), 'w'), indent=2)


def save_counts(tcounts, bcounts, ebands, filename) :

    bcounts = np.sum(bcounts, axis=0)
    
    counts_dict = {}
    for i,e in enumerate(ebands) :
        band_dict = {"counts" : bcounts[i],
                     "emin" : e[0],
                     "emax" : e[1]}
                     
        counts_dict["band" + str(i+1)] = band_dict

    counts_dict["total"] = {"counts" : np.sum(tcounts),
                            "emin" : ebands[0][0],
                            "emax" : ebands[-1][1]}
                            
    json.dump(counts_dict, file(filename, 'w'), indent=2)


def counts_spectra(spec_files, ebands, rmf, inst) :
    counts_full = np.full((len(spec_files),), np.nan)
    counts_ebands = np.full((len(spec_files), len(ebands)), np.nan)
    spec_info = {}

    clean()    
    # Get counts per band and total counts
    for j, spec in enumerate(spec_files) :
        # Load spectrum and auxiliary files (ARF, RMF, BACKG)
        load_pha(spec)        
        load_rmf(rmf[j])
        subtract()

        ungroup(id=1)
        for k,band in enumerate(ebands) :
            notice()
            ignore(":{:f},{:f}:".format(band[0], band[1]))            
            counts_ebands[j,k] = sum(get_counts(filter=True))
        
        notice()
        ignore(":{:f},{:f}:".format(ebands[0][0], ebands[-1][-1]))
        
        counts_full[j] = sum(get_counts(filter=True))
        clean()
        
        exp_dict = {}
        exp_dict["exposure"] = spec.split("/")[-1][13:17]
        exp_dict["detector"] = inst[j]
        exp_dict["counts"] = counts_full[j]
        if counts_full[j] > 50 :
            exp_dict["good"] = 1
        else :
            exp_dict["good"] = 0
        
        spec_info["exp{:d}".format(j+1)] = exp_dict

    return counts_full, counts_ebands, spec_info


def src_fitall(ids, galnH, z, zcdf, ebands, models_to_fit, 
               results_folder=None, plots_folder=None,
               use_zpdf=True, make_plots=True, calc_fluxes=True, 
               bs_fluxes=False) :


    for model_name in models_to_fit :
        print model_name
        ignore(":{:f}, {:f}:".format(ebands[0][0], ebands[-1][-1]))
                
        if use_zpdf :
            # Define models
            define_model(model_name, ids, galnH)

            # Define priors
            parameters, priorfunction = define_priors(model_name, ids, zcdf=zcdf)

            prefix = "{}chains_{}_withzpdf_".format(results_folder, model_name)
            
        else :
            # Define models
            define_model(model_name, ids, galnH, z=z)

            # Define priors
            parameters, priorfunction = define_priors(model_name, ids)

            prefix = "{}chains_{}_".format(results_folder, model_name)

        
        # Run analysis
        set_stat("wstat")
        bxa.nested_run(otherids=ids,
                       prior=priorfunction, 
                       parameters=parameters,
                       resume=True, verbose=False,
                       sampling_efficiency=0.3,  # 0.3
                       n_live_points=400,         # 400
                       evidence_tolerance=0.5,     # 0.1
                       outputfiles_basename=prefix)


        # Estimate posterior distributions and stats
        # (this generates the file "post_equal_weights.dat" and "stats.json")
        a = pymultinest.analyse.Analyzer(n_params=len(parameters), 
                                         outputfiles_basename=prefix)
        stats = a.get_stats()
        json.dump(stats, open(prefix + 'stats.json', 'w'), indent=2)
        
        # Set best fit parameters
        params = a.get_best_fit()['parameters']
        for p,v in zip(parameters, params):
            p.val = v        

        # Estimate and save goodness of fit
        #ks, p = ks_goodness(ids, emin=ebands[0][0], emax=ebands[-1][-1])

        ks, ad, pks, pad = goodness(ids, emin=ebands[0][0], emax=ebands[-1][-1])
        save_goodness(ks, pks, ad, pad, len(parameters), prefix)

        print ks, ad, pks, pad
        goodness2(ids, emin=ebands[0][0], emax=ebands[-1][-1])

        # Make plots
        if make_plots :
            if use_zpdf :
                label_model = "{}_withzpdf".format(model_name)
                
                # Finde the redshift of the best fit
                for p in parameters :
                    if p.fullname.endswith("redshift") :
                        zplot = p.val
                        break
            else :
                label_model = "{}".format(model_name)
                zplot = z

            my_plot_bkg_spec(ids, plots_folder, 
                             emin=ebands[0][0], emax=ebands[-1][-1],
                             label_model=label_model)

            my_plot_fit(ids, plots_folder, 
                        emin=ebands[0][0], emax=ebands[-1][-1],
                        label_model=label_model, z=zplot)

            qq_plot(ids, plots_folder, 
                    emin=ebands[0][0], emax=ebands[-1][-1],
                    label_model=label_model)

            plot_marginals(a, parameters, stats, plots_folder, 
                           label_model=label_model)

        # Load posterior dist.
        chain = a.get_equal_weighted_posterior()
        
        # Estimate best fit parameters and 90% errors from the posterior
        chain_params = chain[:,0:len(parameters)]

        parmode = mode(chain_params, axis=0)
        parmedian = np.median(chain_params, axis=0)
        parlimits = np.percentile(chain_params, [5,95], axis=0)

        # Save best fit parameters and errors
        save_bestfit_pars(parameters, parmedian, parmode, parlimits, 
                          '{}bestfit_pars.json'.format(prefix))

        galabs.nH = 0
        if calc_fluxes:
            if bs_fluxes:
                # Estimate fluxes drawing a subsample from the full posterior
                # distribution of parameters
                idx = np.random.randint(chain.shape[0], size=300)
                subchain = chain[idx,:]

            else :
                subchain = chain

            # Estimate flux distributions (observed and intrinsic)
            chain = get_distribution_with_fluxes(subchain, parameters, 
                                                 ids, ebands=ebands)

            abscorr = model_name.startswith("wa")
            chain = get_distribution_with_intfluxes(chain, parameters,
                                                    ids, ebands=ebands, z=z,
                                                    abscorr=abscorr, 
                                                    thawedz=use_zpdf)
            
            # Estimate luminosity distributions
            if use_zpdf :
                for j,p in enumerate(parameters) :
                        if p.fullname.endswith("redshift") :
                            z_idx=j
                            break
    
                dL = cosmo.luminosity_distance(chain[:,z_idx]).to(u.cm).value
    
            else :
                dL = cosmo.luminosity_distance(z).to(u.cm).value
        
            for j in range(len(ebands)) :
                flux = chain[:, len(parameters) + len(ebands) + 1 + j]
                lumin = 4 * np.pi * dL**2 * flux
                chain = np.append(chain, lumin[:,np.newaxis], axis=1)

            if not bs_fluxes:
                np.savetxt(prefix + "post_equal_weights_fluxeslumins.csv", 
                           chain, delimiter=",")
            
            # Estimate fluxes/luminosities and 90% errors
            chain_flux = chain[:,len(parameters)+1:]
            ncols = 3*len(ebands)
            fluxmode = np.full((ncols), np.nan)
            fluxmedian = np.full((ncols), np.nan)
            fluxlimits = np.full((2,ncols), np.nan)
    
            for col in range(ncols) :
                mask = chain_flux[:,col] > 0
                logflux = np.log10(chain_flux[mask,col])

                #print logflux
                fluxmode[col] = 10**mode(logflux)
                fluxmedian[col] = 10**np.median(logflux)
                fluxlimits[:,col] = 10**np.percentile(logflux, [5,95])
                
            del chain_flux
            del logflux

            # Save best fit fluxes and 90% errors
            first_fluxcol = 0
            save_bestfit_fluxes(ebands, fluxmedian, fluxmode, fluxlimits, first_fluxcol, 
                                '{}bestfit_obsflux.json'.format(prefix), bs_fluxes)
            
            first_fluxcol = len(ebands)
            save_bestfit_fluxes(ebands, fluxmedian, fluxmode, fluxlimits, first_fluxcol, 
                                '{}bestfit_intflux.json'.format(prefix), bs_fluxes)

            # Save best fit luminosities and 90% errors
            first_fluxcol = 2*len(ebands)
            save_bestfit_fluxes(ebands, fluxmedian, fluxmode, fluxlimits, first_fluxcol, 
                                '{}bestfit_lumin.json'.format(prefix), bs_fluxes)

            del a
            del chain


def main(args) :
    logger = logging.getLogger("sherpa")
    logger.setLevel(logging.ERROR)
    set_conf_opt('numcores', 2)

    energy_limits = np.array([0.5, 2.0, 10.0])
    ebands = [(x,y) for x,y in zip(energy_limits[0:-1], energy_limits[1:])]
        
    # Get rmfs
    rmf = []
    inst = []
    for spec in args.specfiles :
        # Get response matrix        
        rmft, instt = update_rmf(spec, rmf_folder=args.rmfs_folder)
        rmf.append(rmft)
        inst.append(instt)

    sorted_inst = np.array(inst)
    sorted_specfiles = np.array(args.specfiles)
    sorted_rmf = np.array(rmf)
    
    idx_instsorted = np.flip(np.argsort(inst), 0)
    sorted_inst = sorted_inst[idx_instsorted]
    sorted_specfiles = sorted_specfiles[idx_instsorted]
    sorted_rmf = sorted_rmf[idx_instsorted]
        
    # Get counts per each available spectrum
    counts_full, counts_ebands, spec_info = counts_spectra(sorted_specfiles, ebands, 
                                                           sorted_rmf, sorted_inst)

    # Identify good spectra (counts in 0.5-10 keV band > 50)
    msk = counts_full > 50
    good_spec_files = np.array(sorted_specfiles)[msk]
    good_rmf = np.array(sorted_rmf)[msk]

    if len(good_spec_files) == 0 :
        print "No good spectra available for this detection!!!"
        return
    
    # Define and make plots folder
    plots_folder = "./plots/{}/".format(args.detid)
    if not os.path.exists(plots_folder) :
        os.makedirs(plots_folder)        

    # Define and make results folder
    results_folder = "./bxa_results/{}/".format(args.detid)
    if not os.path.exists(results_folder) :
        os.makedirs(results_folder)

    # Save counts results 
    save_counts(counts_full, counts_ebands, ebands,
                "{}/counts.json".format(results_folder))

    # Save spec_info
    json.dump(spec_info, 
              file("{}/spec_info.json".format(results_folder), 'w'), 
              indent=2)

    # Load zcdf if needed
    if args.use_zpdf :
        zcdf_object = np.load("temp_zcdf_{}.npz".format(args.detid))
        zcdf = zcdf_object["zcdf"]
    else :
        zcdf = None

    # Load good spectra for analysis
    clean()
    ignore_filter = ":{:f},{:f}:".format(energy_limits[0], energy_limits[-1])

    for j, spec in enumerate(good_spec_files, 1) :
        load_pha(j, str(spec))
        load_rmf(j, str(good_rmf[j-1]))
        load_bkg_rmf(j, str(good_rmf[j-1]))

        # Ignore data outside the 0.5-10 keV range, and group for 1 count per bin
        ungroup(id=j)
        ignore_id(j, ignore_filter)
        group_counts(j, 1)
    
    # Select models for fitting
    if np.sum(counts_full) > 500 :
        models_to_fit = args.simple_models + args.complex_models
    else :
       models_to_fit = args.simple_models

    # Bayesian spectral analysis (+ plots and errors estimation)
    src_fitall(list_data_ids(), args.galnH, args.z, zcdf, ebands, models_to_fit, 
               results_folder=results_folder, 
               plots_folder=plots_folder,
               use_zpdf=args.use_zpdf, 
               make_plots=args.make_plots, 
               calc_fluxes=args.calc_fluxes,
               bs_fluxes=args.bs_fluxes)

    #print counts_full


if __name__ == '__main__' :
    # Parser for shell parameters
    parser = argparse.ArgumentParser(description='X-ray spectral fitting for xmmfitcatz using BXA.')
                                             
    parser.add_argument('--specfiles', dest='specfiles', action='store',
                        default=None, nargs='*',
                        help='List with full routes to the X-ray spectra.')
    
    parser.add_argument('--rmfs', dest='rmfs_folder', action='store',
                        default='./data/rmfcanned/',
                        help='Full route to the folder with the RMFs.')
    
    parser.add_argument('--models_simple', dest='simple_models', action='store',
                        default=['wapo'], nargs='*',
                        help='List of simple models (po, wapo, mekal, wamekal, bb, wabb).')
    
    parser.add_argument('--models_complex', dest='complex_models', action='store',
                        default=['wapopo'], nargs='*',
                        help='List of complex models (wapopo, wamekalpo, wabbpo).')
    
    parser.add_argument('--use_zpdf', dest='use_zpdf', action='store_true',
                        default=False,
                        help='Use the full pdf of photo-z instead of single values.')
    
    parser.add_argument('--make_plots', dest='make_plots', action='store_true',
                        default=False,
                        help='Make all plots for the results (spectra, qq plots, etc).')
    
    parser.add_argument('--calc_fluxes', dest='calc_fluxes', action='store_true',
                        default=False,
                        help='Estimate fluxes from posterior distribution.')

    parser.add_argument('--bs_fluxes', dest='bs_fluxes', action='store_true',
                        default=False,
                        help='Estimate fluxes using a subsample.')

    parser.add_argument('-detid', dest='detid', action='store',
                        default=None, help='Detection ID of the source')

    parser.add_argument('-z', dest='z', action='store', type=float,
                        default=None, help='Redshift of the source')

    parser.add_argument('-nh', dest='galnH', action='store', type=float,
                        default=0.01, help='Galactic nH')
    
    args = parser.parse_args()

    main(args)
