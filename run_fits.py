#!/usr/bin/env python
"""
Created on Mon Dec 11 17:08:18 2017

@author: ruizca
"""

import glob
import itertools
import argparse
import time
import subprocess
import os

import numpy as np

from astropy.table import Table


def main(args):
        
    # Find last fitted source
    try :
        fp = open(args.file_lastsource, "r")
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
        zvals = specz
        #zvals = np.full((len(specz)), np.nan)
    
    else :
        # Change photometric redshifts for spectroscopic redshifts?
        zvals = photoz
        mask_spec = ~np.isnan(specz)
        zvals[mask_spec] = specz[mask_spec]
        

    # Define fixed arguments to pass to fit_source.py
    if args.make_plots :
        fit_args_plots = "--make_plots "
    else :
        fit_args_plots = ""

    if args.calc_fluxes :
        fit_args_flux = "--calc_fluxes "
    else :
        fit_args_flux = ""

    if args.bs_fluxes :
        fit_args_bsflux = "--bs_fluxes "
    else :
        fit_args_bsflux = ""

        
    fit_args_rmf = "--rmfs {} ".format(args.rmfs_folder)
    fit_args_smodels = "--models_simple {} ".format(" ".join(args.simple_models))
    fit_args_cmodels = "--models_complex {} ".format(" ".join(args.complex_models))


    start_total = time.time()
    for i, detid, obsid, srcnum, galnH, z in itertools.izip(itertools.count(), xmm_detid,
                                                         xmm_obsid, xmm_srcnum, nHgal, zvals):
        fp = open(args.file_lastsource, "w")
        fp.write(str(i+first_source))
        fp.close()

        print "\n\nDet. {:d} ({:d} out of {:d})".format(detid, i + first_source + 1, nsrcs)

        # Find spectra of interest for this observation and this source
        spec_folder = os.path.join(args.spectra_folder, obsid, 'product/')
        spec_files = glob.glob("{}*SRSPEC{:04X}.FTZ".format(spec_folder, srcnum))

        if len(spec_files) == 0 :
            print "No spectra available for this detection!!!"
            continue

        fit_args_specfiles = "--specfiles {} ".format(" ".join(spec_files))
        
        # Estimate cumulative distribution of z if needed
        if args.use_zpdf :
            if (not np.isnan(z)) and (z > 0) :
                fit_args_zpdf = ""
    
            else :
                mask_z = ~np.isnan(zrange[i])
                cdf = zpdf[i][mask_z].cumsum()
                zcdf_file = "temp_zcdf_{}".format(detid)
                np.savez(zcdf_file, zcdf=np.array([zrange[i][mask_z], 
                                                   cdf/np.max(cdf)]))
                fit_args_zpdf = "--use_zpdf "
        else :
            fit_args_zpdf = ""


        # Fit source
        fit_command = "./fit_source.py "
        fit_args = "-z {:f} -nh {} -detid {} ".format(z, galnH, detid)
        fit_args = fit_args + fit_args_zpdf + fit_args_flux + fit_args_bsflux + \
                  fit_args_plots + fit_args_rmf + fit_args_specfiles + \
                  fit_args_smodels + fit_args_cmodels             
        print fit_command + fit_args
        
        start_time = time.time()
        try :
            subprocess.call(fit_command + fit_args, shell=True)
            
        except :
            print "Something went wrong fitting det. {}".format(detid)

                    
        print "Time of fit {:d}: {:f}".format(i+1, time.time() - start_time )
        try :
            os.remove(zcdf_file + ".npz")
            
        except :
            print "No zcdf file!"

        break
        #if i==10 :
        #    break

    print time.time() - start_total
    print "All done!!!"  

                                                             
if __name__ == '__main__':
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

    parser.add_argument('--lsf', dest='file_lastsource', action='store',
                        default='last_source.dat',
                        help='File to store the last fitted source.')
                        
    args = parser.parse_args()

    main(args)
