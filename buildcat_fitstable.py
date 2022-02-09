#!/usr/bin/env python
"""
Created on Tue Dec  5 12:46:05 2017

@author: ruizca
"""

import json
import numpy as np

from astropy.table import Table, hstack
from astropy import units as u

sources_table = "./data/3xmmdr6_detwithphotoz_test_nhgal.fits"
final_table = "./data/3xmmdr6_spectralfits_test.fits"


xmmdet_cat = Table.read(sources_table, format="fits")
xmm_iauname = xmmdet_cat.columns["IAUNAME"]
xmm_srcra = xmmdet_cat.columns["SC_RA"]
xmm_srcdec = xmmdet_cat.columns["SC_DEC"]
xmm_srcid = xmmdet_cat.columns["xmm_SRCID"]
xmm_detid = xmmdet_cat.columns["DETID"]
xmm_obsid = xmmdet_cat.columns["OBS_ID"]
xmm_srcnum = xmmdet_cat.columns["SRC_NUM"]
nHgal = xmmdet_cat.columns["NHGAL"]
zph = xmmdet_cat.columns["PHOT_Z"]
zphErr = xmmdet_cat.columns["PHOT_ZERR"]
zspc = xmmdet_cat.columns["SPEC_Z"]

nsrcs = len(xmm_iauname)

models = np.array(["po", "wapo", "mekal", "wamekal", "bb", "wabb", "wabbpo", "wamekalpo", "wapopo"])

# Source to get the model parameters
test_detid = "101010401010001"
prefix_pars = "./bxa_results/" + test_detid + "/chains_"


### Define new columns
afit = np.full((nsrcs), np.nan)
best_model = [None]*nsrcs
good_models = [None]*nsrcs

print len(best_model)
counts_total = np.full((nsrcs), np.nan)
counts_soft = np.full((nsrcs), np.nan)
counts_hard = np.full((nsrcs), np.nan)

models_dict = {}
for model in models :
    # Goodness
    models_dict[model + "_wstat"] = np.full((nsrcs), np.nan)
    models_dict[model + "_dof"] = np.full((nsrcs), np.nan)
    models_dict[model + "_ks"] = np.full((nsrcs), np.nan)
    models_dict[model + "_pvalue"] = np.full((nsrcs), np.nan)
    models_dict[model + "_logZ"] = np.full((nsrcs), np.nan)
    
    # Parameters
    paramfile = prefix_pars + model + '_withzpdf_bestfit_pars.json'
    with open(paramfile) as data_file:
        params = json.load(data_file)
    
    for p in params.keys() :
        parname = p.split(".")[-1]
        if not (parname.startswith("c") or parname.startswith("lognorm")) :
            key = model + "_" + p
            models_dict[key] = np.full((nsrcs), np.nan)
            models_dict[key + "_min"] = np.full((nsrcs), np.nan)
            models_dict[key + "_max"] = np.full((nsrcs), np.nan)
        
    # Fluxes
    models_dict[model + "_flux_soft"] = np.full((nsrcs), np.nan)
    models_dict[model + "_fluxmin_soft"] = np.full((nsrcs), np.nan)
    models_dict[model + "_fluxmax_soft"] = np.full((nsrcs), np.nan)
    models_dict[model + "_flux_hard"] = np.full((nsrcs), np.nan)
    models_dict[model + "_fluxmin_hard"] = np.full((nsrcs), np.nan)
    models_dict[model + "_fluxmax_hard"] = np.full((nsrcs), np.nan)
    
    # Intrinsic fluxes (rest frame band, absorbtion corrected)
    models_dict[model + "_intflux_soft"] = np.full((nsrcs), np.nan)
    models_dict[model + "_intfluxmin_soft"] = np.full((nsrcs), np.nan)
    models_dict[model + "_intfluxmax_soft"] = np.full((nsrcs), np.nan)
    models_dict[model + "_intflux_hard"] = np.full((nsrcs), np.nan)
    models_dict[model + "_intfluxmin_hard"] = np.full((nsrcs), np.nan)
    models_dict[model + "_intfluxmax_hard"] = np.full((nsrcs), np.nan)

    # Luminosities (rest frame band, absorbtion corrected)
    models_dict[model + "_lumin_soft"] = np.full((nsrcs), np.nan)
    models_dict[model + "_luminmin_soft"] = np.full((nsrcs), np.nan)
    models_dict[model + "_luminmax_soft"] = np.full((nsrcs), np.nan)
    models_dict[model + "_lumin_hard"] = np.full((nsrcs), np.nan)
    models_dict[model + "_luminmin_hard"] = np.full((nsrcs), np.nan)
    models_dict[model + "_luminmax_hard"] = np.full((nsrcs), np.nan)



for i,detid in enumerate(xmm_detid) :
    print i,detid

    data_folder = "./bxa_results/" + str(detid) + "/"
    
    # Counts
    countsfile = data_folder + "counts.json"    
    with open(countsfile) as data_file :
        counts = json.load(data_file)

    counts_total[i] = counts["total"]["counts"]
    counts_soft[i] = counts["band1"]["counts"]
    counts_hard[i] = counts["band2"]["counts"]
    
    logZ = np.full((len(models)), np.nan)
    
    afit[i] = False
    
    # Model related columns
    for j,model in enumerate(models) :
        try :
            # Goodnes of fit
            goodnessfile = data_folder + "chains_" \
                           + model + "_withzpdf_bestfit_goodness.json"                           

            with open(goodnessfile) as data_file :
                goodness = json.load(data_file)
    
            models_dict[model + "_wstat"][i] = goodness["wstat"]
            models_dict[model + "_dof"][i] = goodness["dof"]
            models_dict[model + "_ks"][i] = goodness["ks"]
            models_dict[model + "_pvalue"][i] = goodness["pvalue"]
            models_dict[model + "_logZ"][i] = goodness["lnZ"]
            logZ[j] = goodness["lnZ"]
            
            if goodness["pvalue"] > 0.1 :
                afit[i] = True


            # Best fit parameters
            paramsfile = data_folder + "chains_" \
                           + model + "_withzpdf_bestfit_pars.json"

            with open(paramsfile) as data_file :
                params = json.load(data_file)

            for p in params.keys() :
                parname = p.split(".")[-1]
                if not (parname.startswith("c") or parname.startswith("lognorm")) :
                    key = model + "_" + p
                    models_dict[key][i] = params[p]['mode']
                    models_dict[key + "_min"][i] = params[p]['parmin']
                    models_dict[key + "_max"][i] = params[p]['parmax']
                
                if parname == "redshift" :
                    z = params[p]['mode']

            
            # Observed fluxes
            fluxfile = data_folder + "chains_" \
                           + model + "_withzpdf_bestfit_obsflux.json"

            with open(fluxfile) as data_file :
                obsflux = json.load(data_file)
                           
            models_dict[model + "_flux_soft"][i] = obsflux["flux0"]["mode"] 
            models_dict[model + "_fluxmin_soft"][i] = obsflux["flux0"]["fluxmin"]
            models_dict[model + "_fluxmax_soft"][i] = obsflux["flux0"]["fluxmax"]
            models_dict[model + "_flux_hard"][i] = obsflux["flux1"]["mode"]
            models_dict[model + "_fluxmin_hard"][i] = obsflux["flux1"]["fluxmin"]
            models_dict[model + "_fluxmax_hard"][i] = obsflux["flux1"]["fluxmax"]


            # Intrinsic fluxes
            fluxfile = data_folder + "chains_" \
                           + model + "_withzpdf_bestfit_intflux.json"

            with open(fluxfile) as data_file :
                obsflux = json.load(data_file)
                           
            models_dict[model + "_intflux_soft"][i] = obsflux["flux0"]["mode"]
            models_dict[model + "_intfluxmin_soft"][i] = obsflux["flux0"]["fluxmin"]
            models_dict[model + "_intfluxmax_soft"][i] = obsflux["flux0"]["fluxmax"]
            models_dict[model + "_intflux_hard"][i] = obsflux["flux1"]["mode"]
            models_dict[model + "_intfluxmin_hard"][i] = obsflux["flux1"]["fluxmin"]
            models_dict[model + "_intfluxmax_hard"][i] = obsflux["flux1"]["fluxmax"]


            # Luminosities
            luminfile = data_folder + "chains_" \
                           + model + "_withzpdf_bestfit_lumin.json"

            with open(luminfile) as data_file :
                lumin = json.load(data_file)
                           
            models_dict[model + "_lumin_soft"][i] = lumin["flux0"]["mode"]
            models_dict[model + "_luminmin_soft"][i] = lumin["flux0"]["fluxmin"]
            models_dict[model + "_luminmax_soft"][i] = lumin["flux0"]["fluxmax"]
            models_dict[model + "_lumin_hard"][i] = lumin["flux1"]["mode"]
            models_dict[model + "_luminmin_hard"][i] = lumin["flux1"]["fluxmin"]
            models_dict[model + "_luminmax_hard"][i] = lumin["flux1"]["fluxmax"]

            
        except IOError :
            print "Model " + model + " wasn't fitted for this source!"

    # Finde best and good models
    # The most likely model is used as normalization.
    # Uniform model priors are assumed, with a cut of 
    # log10(limit) to rule out models.'
    limit = 30 # for example, Jeffreys scale for the Bayes factor
    
    mask_fitted = ~np.isnan(logZ)
    models_fit = models[mask_fitted]
    logZ = logZ[mask_fitted]    

    best = np.argsort(logZ)[-1]
    Zrel = logZ - logZ[best]
    Ztotal = np.log(np.sum(np.exp(Zrel)))
    
    good_mask = np.logical_and(Zrel >= Ztotal - np.log(limit), Zrel<0)
    best_model[i] = models_fit[best]
    good_models[i] = " ".join(models_fit[good_mask])
    
    

table_srcdata = Table((xmm_iauname, 
                       xmm_srcra, 
                       xmm_srcdec,
                       xmm_srcid,
                       xmm_detid,
                       xmm_obsid,
                       xmm_srcnum,
                       zph,
                       zphErr,
                       zspc,
                       nHgal))

table_counts = Table(data=(counts_total, counts_soft, counts_hard, afit, best_model, good_models),
                     names=("T_COUNTS","S_COUNTS","H_COUNTS", "A_FIT", "P_MODEL", "A_MODELS"))

cols_goodness = ["_wstat",
                "_dof",
                "_ks",
                "_pvalue",
                "_logZ"]

cols_flux = [
    "_flux_soft", "_fluxmin_soft", "_fluxmax_soft",
    "_flux_hard", "_fluxmin_hard", "_fluxmax_hard",
    "_intflux_soft", "_intfluxmin_soft", "_intfluxmax_soft",
    "_intflux_hard", "_intfluxmin_hard", "_intfluxmax_hard",
    "_lumin_soft", "_luminmin_soft", "_luminmax_soft",
    "_lumin_hard", "_luminmin_hard", "_luminmax_hard"]

cols_pars_po = [
    "_po.PhoIndex", "_po.PhoIndex_min", "_po.PhoIndex_max",
    "_po.redshift", "_po.redshift_min", "_po.redshift_max"]

cols_pars_wapo = [
    "_intabs.logNH", "_intabs.logNH_min", "_intabs.logNH_max",
    "_po.PhoIndex", "_po.PhoIndex_min", "_po.PhoIndex_max",
    "_po.redshift", "_po.redshift_min", "_po.redshift_max"]
                  
cols_pars_thermal = [
    "_thermal.kT", "_thermal.kT_min", "_thermal.kT_max",
    "_thermal.redshift", "_thermal.redshift_min", "_thermal.redshift_max"]
                  
cols_pars_wathermal = [
    "_intabs.logNH", "_intabs.logNH_min", "_intabs.logNH_max",
    "_thermal.kT", "_thermal.kT_min", "_thermal.kT_max",
    "_thermal.redshift", "_thermal.redshift_min", "_thermal.redshift_max"]

cols_pars_wabbpo = [
    "_intabs.logNH", "_intabs.logNH_min", "_intabs.logNH_max",
    "_thermal.kT", "_thermal.kT_min", "_thermal.kT_max",
    "_po.PhoIndex", "_po.PhoIndex_min", "_po.PhoIndex_max",
    "_po.redshift", "_po.redshift_min", "_po.redshift_max"]

cols_pars_wamekalpo = [
    "_intabs1.logNH1", "_intabs1.logNH1_min", "_intabs1.logNH1_max",
    "_intabs2.logNH2", "_intabs2.logNH2_min", "_intabs2.logNH2_max",
    "_thermal.kT", "_thermal.kT_min", "_thermal.kT_max",
    "_po.PhoIndex", "_po.PhoIndex_min", "_po.PhoIndex_max",
    "_po.redshift", "_po.redshift_min", "_po.redshift_max"]

cols_pars_wapopo = [
    "_intabs1.logNH1", "_intabs1.logNH1_min", "_intabs1.logNH1_max",
    "_intabs2.logNH2", "_intabs2.logNH2_min", "_intabs2.logNH2_max",
    "_po1.PhoIndex", "_po1.PhoIndex_min", "_po1.PhoIndex_max",
    "_po2.PhoIndex", "_po2.PhoIndex_min", "_po2.PhoIndex_max",    
    "_po2.redshift", "_po2.redshift_min", "_po2.redshift_max"]


columns_data = []
columns_names = []
for model in models :

    if model == "po" :
        cols_pars = cols_pars_po
    
    elif model == "wapo" :
        cols_pars = cols_pars_wapo

    elif model == "mekal" or model == "bb" :
        cols_pars = cols_pars_thermal
        
    elif model == "wamekal" or model == "wabb" :
        cols_pars = cols_pars_wathermal

    elif model == "wamekalpo" :
        cols_pars = cols_pars_wamekalpo

    elif model == "wapopo" :
        cols_pars = cols_pars_wapopo

    elif model == "wabbpo" :
        cols_pars = cols_pars_wabbpo
        
    columns = np.concatenate((cols_goodness, cols_pars, cols_flux))
    for col in columns :
        name = model + col
        columns_data.append(models_dict[name])
        columns_names.append(name)

table_models = Table(data=columns_data, names=columns_names)        

table_all = hstack([table_srcdata, table_counts, table_models], join_type='exact')
table_all.write(final_table, format="fits", overwrite=True)