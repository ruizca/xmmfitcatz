#!/usr/bin/env python
"""
Created on Wed Nov 22 12:15:29 2017

@author: ruizca
"""

import glob
import logging
import os
import time
import json
import argparse
import itertools
import sys
import gc
import joblib

import numpy as np
import scipy.stats
import pymultinest

from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import Planck15 as cosmo
from astropy import units as u

from sherpa.astro.ui import *
import bxa.sherpa as bxa

import matplotlib
matplotlib.use('AGG') 
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes.formatter', limits=(-3,3))
plt.style.use('bmh')

from models import *
from priors import *


def get_rmf(specfile, rmf_folder) :
    """
    Read specfile (XMM spectra in fits format) and returns the corresponding rmf file
    defined in the spectrum header
    """
    #epoch_limits
    d1 = time.strptime("2007-01-01", "%Y-%m-%d")
    d2 = time.strptime("2014-01-01", "%Y-%m-%d")

    hdulist = fits.open(specfile)
    obs_date = hdulist[0].header['DATE-OBS']
    rmf_old  = hdulist[1].header['RESPFILE']
    camera   = hdulist[1].header['INSTRUME']
    deltaE   = hdulist[1].header['SPECDELT']
    
    hdulist.close()

    
    if camera == 'EPN' :
        obs_date = time.strptime(obs_date, "%Y-%m-%dT%H:%M:%S")
        
        if obs_date < d1 :
            epoch = "e1"
            
        elif obs_date >= d1 and obs_date < d2 :
            epoch = "e2"
                    
        else :
            epoch = "e3"

        rmf_split = rmf_old.split('_',1)
        if rmf_split[1].startswith("e1") or rmf_split[1].startswith("e3")  or rmf_split[1].startswith("e3") :
            rmf_new = rmf_folder + rmf_split[0] + "_" + rmf_split[1].replace(".rmf", "_v16.0.rmf")                
        else :
            rmf_new = rmf_folder + rmf_split[0] + "_" + epoch + "_" + rmf_split[1].replace(".rmf", "_v16.0.rmf")

    else :
        if deltaE == 5 :
            rmf_new = rmf_folder + "5ev/" + rmf_old
        else :
            rmf_new = rmf_folder + rmf_old

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
        
        ignore_filter = ":" + str(emin) + "," + str(emax) + ":"
        ignore_id(sid, ignore_filter)

    modelc = model.cumsum()/model.sum()
    datac = data.cumsum()/data.sum()
        
    ks = np.abs(modelc - datac).max()
    en = np.sqrt(1.0*len(data)*len(model)/(len(data) + len(model))) 
    prob = scipy.stats.distributions.kstwobign.sf((en + 0.12 + 0.11/en) * ks)

    return ks, prob


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
            label = label[11:17]
            filename = plot_folder + label_model + "_bkgspec_" + label + ".png"
            fig.savefig(filename)
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
                     fmt="o", ms=2, elinewidth=0.8, ls="None")
        ax0.loglog(e[mask_model], model[mask_model], c='red', alpha=0.5)
                                     
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
            label = label[11:17]
            filename = plot_folder + label_model + "_srcspec_" + label + ".png"
            
            fig.savefig(filename)
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
                         fmt="o", ms=2, elinewidth=0.8, ls="None")
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
            filename = plot_folder + label_model + "_srcspec_ALL.png"
            fig.savefig(filename)
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
        filename = plot_folder + label_model + "_qqplot.png"
        fig.savefig(filename)
    else :
        fig.show()

    plt.close(fig)


def plot_marginals(analyzer, parameters, stats, plot_folder, label_model="", save=True) :
    p = pymultinest.PlotMarginal(analyzer)
    values = analyzer.get_equal_weighted_posterior()
    
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
        x,w,patches = oldax.hist(values[:,i], bins='auto', 
                                 edgecolor='black', color='black', 
                                 histtype='stepfilled', alpha=0.2)

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
            for m in modes:
                plt.errorbar(x=m['mean'][i], y=m['mean'][j], xerr=m['sigma'][i], yerr=m['sigma'][j])
            plt.xlabel(parameters[i].fullname)
            plt.ylabel(parameters[j].fullname)

    fig.tight_layout()
    
    if save :
        filename = plot_folder + label_model + "_param_marginals.png"
        fig.savefig(filename, dpi=50)
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
                modelname = p.fullname.split(".")[0]
                par = modelname + ".nH"
                set_par(par, val=0)

    nrows, ncols_chains = np.shape(chains)
    ncols = ncols_chains + len(ebands)
    
    r = np.full((nrows, ncols), np.nan)

    first_colflux = ncols_chains
    r[:, 0:first_colflux] = chains
    r[:, first_colflux:] = [calc_flux_pars(parameters, row, ids,
                                           ebands=ebands, thawedz=thawedz)
                            for row in chains]

    return r


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
    
    return params_dict


def save_bestfit_fluxes(ebands, parmedian, parmode, parlimits, firstcol, filename) :

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
                     "emax" : e[1]}

        fluxes_dict["flux" + str(i)] = flux_dict
    
    json.dump(fluxes_dict, file(filename, 'w'), indent=2)
    
    return fluxes_dict


def save_goodness(ks, p, npar, prefix) :

    stats = get_stat_info()[-1]
    wstat = stats.statval
    dof = stats.numpoints - npar

    with open(prefix + 'stats.json') as data_file:
        multinest_stats = json.load(data_file)
        
    dict_goodness = {"wstat" : wstat,
                     "dof" : dof,
                     "npar" : npar,
                     "ks" : ks,
                     "pvalue" : p,
                     "lnZ" : multinest_stats["global evidence"],                         
                     "lnZ_err" : multinest_stats["global evidence error"],
                     "lnev" : multinest_stats["nested sampling global log-evidence"],
                     "lnev_err" : multinest_stats["nested sampling global log-evidence error"]}

    json.dump(dict_goodness, file(prefix + "bestfit_goodness.json", 'w'), indent=2)
    
    return dict_goodness


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
    
    return counts_dict


def make_fit(i, detid, obsid, srcnum, galnH, z, use_zpdf=True, make_plots=True, calc_fluxes=True) :

    #fp = open("last_source.dat", "w")
    #fp.write(str(i+first_source))
    #fp.close()

    clean()
    print "\n\nDet. " + str(detid) + " (" + str(i + first_source + 1) + " out of " + str(nsrcs) + ")"
    spec_folder = args.spectra_folder + str(obsid) + "/product/"

    # find spectra of interest for this observation and this source
    spec_files = glob.glob(spec_folder + '*SRSPEC' 
                           + format(srcnum, '04x').upper() + '.FTZ')

    if len(spec_files) == 0 :
        print "No spectra available for this detection!!!"
        return

    
    counts_full = np.full((len(spec_files),), np.nan)
    counts_ebands = np.full((len(spec_files), len(ebands)), np.nan)
    rmf = []
    spec_info = {}
    
    # Identify good spectra (counts in 0.5-10 keV band > 50)
    for j, spec in enumerate(spec_files) :
        # Get response matrix        
        rmft, inst = get_rmf(spec, rmf_folder=args.rmfs_folder)
        rmf.append(rmft)
                
        # Load spectrum and auxiliary files (ARF, RMF, BACKG)
        load_pha(spec)
        load_rmf(rmf[j])
        subtract()

        ungroup(id=1)
        for k,band in enumerate(ebands) :
            notice()
            ignore(":" + str(band[0]) + ", " + str(band[1]) + ":")
            counts_ebands[j,k] = sum(get_counts(filter=True))
        
        notice()
        ignore(ignore_filter)
        
        counts_full[j] = sum(get_counts(filter=True))
        clean()
        
        exp_dict = {}
        exp_dict["exposure"] = spec.split("/")[-1][13:17]
        exp_dict["detector"] = inst
        exp_dict["counts"] = counts_full[j]
        if counts_full[j] > 50 :
            exp_dict["good"] = 1
        else :
            exp_dict["good"] = 0
        
        exp_name = "exp" + str(j+1)
        spec_info[exp_name] = exp_dict

    msk = counts_full > 50
    good_spec_files = np.array(spec_files)[msk]
    good_rmf = np.array(rmf)[msk]

    if len(good_spec_files) == 0 :
        print "No good spectra available for this detection!!!"
        return

    
    # Load good spectra for analysis
    for j, spec in enumerate(good_spec_files, 1) :
        load_pha(j, str(spec))
        load_rmf(j, str(good_rmf[j-1]))
        load_bkg_rmf(j, str(good_rmf[j-1]))

        # Ignore data outside the 0.5-10 keV range, and group for 1 count per bin
        ungroup(id=j)
        ignore_id(j, ignore_filter)
        group_counts(j, 1)

    
    # Define and make plots folder
    plots_folder = "./plots/" + str(xmm_detid[i]) + "/"
    if not os.path.exists(plots_folder) :
        os.makedirs(plots_folder)        

    # Define and make results folder
    results_folder = "./bxa_results/" + str(xmm_detid[i])
    if not os.path.exists(results_folder) :
        os.makedirs(results_folder)


    # Save counts results 
    counts_json = results_folder + "/counts.json"
    save_counts(counts_full, counts_ebands, ebands, counts_json)

    # Save spec_info
    spec_json = results_folder + "/spec_info.json"
    json.dump(spec_info, file(spec_json, 'w'), indent=2)
    

    # Select models for fitting
    if np.sum(counts_full) > 500 :
        models_to_fit = args.simple_models + args.complex_models
    else :
       models_to_fit = args.simple_models


    # Estimate cumulative distribution of z if needed
    if use_zpdf :
        if not np.isnan(z) :
            use_zpdf_src = False

        else :
            use_zpdf_src = True
            mask_z = ~np.isnan(zrange[i])
            cdf = zpdf[i][mask_z].cumsum()
            zcdf = np.array([zrange[i][mask_z], cdf/np.max(cdf)])
    else :
        use_zpdf_src = False
    

    ### Bayesian spectral analysis
    for model_name in models_to_fit :
        print model_name
        ignore(ignore_filter)
                
        if use_zpdf_src :
            # Define models
            define_model(model_name, list_data_ids(), nHgal[i])

            # Define priors
            parameters, priorfunction = define_priors(model_name, list_data_ids(), zcdf=zcdf)

            prefix = results_folder + "/chains_" + model_name + "_withzpdf_"
            
        else :
            # Define models
            define_model(model_name, list_data_ids(), nHgal[i], z=z)

            # Define priors
            parameters, priorfunction = define_priors(model_name, list_data_ids())

            prefix = results_folder + "/chains_" + model_name + "_"

        
        # Run analysis
        set_stat("wstat")
        bxa.nested_run(otherids = list_data_ids(), 
                       prior = priorfunction, 
                       parameters = parameters,
                       resume = True, verbose = False,
                       sampling_efficiency = 0.3,  # 0.3
                       n_live_points = 400,         # 400
                       evidence_tolerance = 0.5,     # 0.1
                       outputfiles_basename = prefix)


        # Estimate posterior distributions and stats
        # (this generates the file "post_equal_weights.dat" and "stats.json")
        a = pymultinest.analyse.Analyzer(n_params = len(parameters), 
                                         outputfiles_basename = prefix)
        stats = a.get_stats()
        json.dump(stats, open(prefix + 'stats.json', 'w'), indent=4)
        
        # Set best fit parameters
        params = a.get_best_fit()['parameters']
        for p,v in zip(parameters, params):
            p.val = v        

        # Estimate and save goodness of fit
        ks, p = ks_goodness(list_data_ids(), 
                            emin=energy_limits[0], 
                            emax=energy_limits[-1])

        save_goodness(ks, p, len(parameters), prefix)


        # Make plots
        if make_plots :
            if use_zpdf_src :
                label_model = model_name + "_withzpdf"
                
                # Finde the redshift of the best fit
                for p in parameters :
                    if p.fullname.endswith("redshift") :
                        zplot = p.val
                        break
            else :
                label_model = model_name
                zplot = z

            my_plot_bkg_spec(list_data_ids(), plots_folder, 
                             emin=energy_limits[0], emax=energy_limits[-1],
                             label_model=label_model)

            my_plot_fit(list_data_ids(), plots_folder, 
                        emin=energy_limits[0], emax=energy_limits[-1],
                        label_model=label_model, z=zplot)

            qq_plot(list_data_ids(), plots_folder, 
                    emin=energy_limits[0], emax=energy_limits[-1],
                    label_model=label_model)

            plot_marginals(a, parameters, stats, plots_folder, 
                           label_model=label_model)


        # Load posterior dist.
        chain = a.get_equal_weighted_posterior()
        
        # Estimate best fit parameters and 90% errors from the posterior
        chain_params = chain[:,0:len(parameters)]
        parmode = np.full((len(parameters)), np.nan)
        for col in range(len(parameters)) :
            hist, xedges = np.histogram(chain[:,col], bins="auto")
            x = (xedges[:-1] + xedges[1:])/2            
            maxp = np.argsort(hist)[-1]  
            parmode[col] = x[maxp]
            
        parmedian = np.median(chain_params, axis=0)
        parlimits = np.percentile(chain_params, [5,95], axis=0)


        if calc_fluxes :
            # Estimate flux distributions (observed and intrinsic)
            galabs.nH = 0            
            chain = get_distribution_with_fluxes(chain, parameters, 
                                                 list_data_ids(), ebands=ebands)
    
            abscorr = model_name.startswith("wa")
            chain = get_distribution_with_intfluxes(chain, parameters,
                                                    list_data_ids(),
                                                    ebands=ebands, z=z,
                                                    abscorr=abscorr, thawedz=use_zpdf_src)
            
            # Save posterior distributions with fluxes
            np.savetxt(prefix + "post_equal_weights_fluxes.csv", chain, delimiter=",")


        # Estimate luminosity distributions
        if not calc_fluxes :
            chain = np.loadtxt(prefix + "post_equal_weights_fluxes.csv", delimiter=",")

        if use_zpdf_src :
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

        np.savetxt(prefix + "post_equal_weights_fluxeslumins.csv", chain, delimiter=",")

        
        # Estimate fluxes/luminosities and 90% errors from the posterior
        chain_flux = chain[:,len(parameters)+1:]
        ncols = 3*len(ebands)
        fluxmode = np.full((ncols), np.nan)
        fluxmedian = np.full((ncols), np.nan)
        fluxlimits = np.full((2,ncols), np.nan)

        for col in range(ncols) :
            mask = chain_flux[:,col] > 0
            logflux = np.log10(chain_flux[mask,col])
            hist, xedges = np.histogram(logflux, bins="auto")
            x = (xedges[:-1] + xedges[1:])/2
            maxp = np.argsort(hist)[-1]  

            fluxmode[col] = 10**x[maxp]
            fluxmedian[col] = 10**np.median(logflux)
            fluxlimits[:,col] = 10**np.percentile(logflux, [5,95])

        del chain
        del chain_flux
        del logflux


        # Save best fit parameters and errors
        fitpars_json = prefix + 'bestfit_pars.json'
        save_bestfit_pars(parameters, parmedian, parmode, parlimits, fitpars_json)

        # Save best fit fluxes and 90% errors
        first_fluxcol = 0
        fitflux_json = prefix + 'bestfit_obsflux.json'
        save_bestfit_fluxes(ebands, fluxmedian, fluxmode, fluxlimits, first_fluxcol, fitflux_json)
        
        first_fluxcol = len(ebands)
        fitflux_json = prefix + 'bestfit_intflux.json'
        save_bestfit_fluxes(ebands, fluxmedian, fluxmode, fluxlimits, first_fluxcol, fitflux_json)

        # Save best fit luminosities and 90% errors
        first_fluxcol = 2*len(ebands)
        fitlumin_json = prefix + 'bestfit_lumin.json'
        save_bestfit_fluxes(ebands, fluxmedian, fluxmode, fluxlimits, first_fluxcol, fitlumin_json)
        
        gc.collect()
        del gc.garbage[:]


# Parser for shell parameters
parser = argparse.ArgumentParser(description='X-ray spectral fitting for xmmfitcatz using BXA.')
                 
parser.add_argument('--catalogue', dest='sources_table', action='store',
                    default='./data/3xmmdr6_detwithphotoz_nhgal.fits',
                    help='Full route to the detections catalogue.')

parser.add_argument('--catalogue_zpdf', dest='zpdfs_table', action='store',
                    default='./data/3xmmdr6_detwithphotoz_zpdfs.fits',
                    help='Full route to the catalogue of photo-z pdfs.')
# Both tables (catalogue and catalogue_zpdf) 
# must be sorted exactly the same way!!!!
                    
parser.add_argument('--spectra', dest='spectra_folder', action='store',
                    default='./data/3xmmspectra/',
                    help='Full route to the folder with the X-ray spectra.')

parser.add_argument('--rmfs', dest='rmfs_folder', action='store',
                    default='./data/rmfcanned/',
                    help='Full route to the folder with the RMFs.')

parser.add_argument('--models_simple', dest='simple_models', action='store',
                    default=['wapo'], nargs='*',
                    help='List of simple models (po, wapo, mekal, wamekal, bb, wabb).')

parser.add_argument('--models_complex', dest='complex_models', action='store',
                    default=['wapopo'], nargs='*',
                    help='List of simple models (po, wapo, mekal, wamekal, bb, wabb).')

parser.add_argument('--use_zpdf', dest='use_zpdf', action='store_true',
                    default=False,
                    help='Use the full pdf of photo-z instead of single values.')

parser.add_argument('--make_plots', dest='make_plots', action='store_true',
                    default=False,
                    help='Make all plots for the results (spectra, qq plots, etc).')

parser.add_argument('--calc_fluxes', dest='calc_fluxes', action='store_true',
                    default=False,
                    help='Estimate fluxes from posterior distribution.')

args = parser.parse_args()


logger = logging.getLogger("sherpa")
logger.setLevel(logging.ERROR)
set_conf_opt('numcores', 2)


energy_limits = np.array([0.5, 2.0, 10.0])
ebands = [(x,y) for x,y in zip(energy_limits[0:-1], energy_limits[1:])]        
ignore_filter = ":" + str(energy_limits[0]) + "," + str(energy_limits[-1]) + ":"    

# Find last fitted source
try :
    fp = open("last_source.dat", "r")
    first_source = int(fp.readline())
    fp.close()
except :
    first_source = 0


## Load detections table
xmmdet_cat = Table.read(args.sources_table, format="fits")
xmm_detid = xmmdet_cat.columns["DETID"][first_source:]
xmm_obsid = xmmdet_cat.columns["OBS_ID"][first_source:]
xmm_srcnum = xmmdet_cat.columns["SRC_NUM"][first_source:]
nHgal = xmmdet_cat.columns["NHGAL"][first_source:]
photoz = xmmdet_cat.columns["PHOT_Z"][first_source:]
specz = xmmdet_cat.columns["SPEC_Z"][first_source:]

nsrcs = len(xmmdet_cat.columns["DETID"])

del xmmdet_cat


# Load (if needed) zpdfs table
if args.use_zpdf :
    zpdfs_cat = Table.read(args.zpdfs_table, format="fits")
    zrange = zpdfs_cat.columns["PDF_z"][first_source:]
    zpdf = zpdfs_cat.columns["PDF values"][first_source:]
    del zpdfs_cat
    
    # Change photometric redshifts for spectroscopic redshifts?
    #zvals = specz
    zvals = np.full((len(specz)), np.nan)

else :
    # Change photometric redshifts for spectroscopic redshifts?
    zvals = photoz
    mask_spec = ~np.isnan(specz)
    zvals[mask_spec] = specz[mask_spec]



joblib.Parallel(n_jobs=4)(joblib.delayed(make_fit)(i, detid, obsid, srcnum, galnH, z, 
                                                   use_zpdf=args.use_zpdf, make_plots=args.make_plots, 
                                                   calc_fluxes=args.calc_fluxes)
                          for i, detid, obsid, srcnum, galnH, z in itertools.izip(
                                   itertools.count(), xmm_detid, xmm_obsid, xmm_srcnum, nHgal, zvals)
                         )



