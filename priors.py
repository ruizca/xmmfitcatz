#!/usr/bin/env python
# Filename: priors.py

import scipy.interpolate

import bxa.sherpa as bxa
from sherpa.astro.ui import *
from sherpa.models.parameter import Parameter


def get_inv_zpdf_func(zcdf) :

    cdf = scipy.interpolate.interp1d(zcdf[1,:], zcdf[0,:], 
                                     bounds_error=False, fill_value=0.001)
    return lambda x : cdf(x)

    
def define_priors_po(ids, zcdf=None) :

    # Define array with the thawed parameters from the defined 
    # model of the first spectrum
    lognorm = Parameter(modelname='po', name='lognorm', val=-5, min=-30, max=0)
    po.norm = 10**lognorm

    parameters = [po.PhoIndex, lognorm]
        
    priors = []
    # photon index (Distribution from n Nandra & Pounds (1994))
    priors += [bxa.create_gaussian_prior_for(parameters[0], 1.9, 0.15)]
    #priors += [bxa.create_uniform_prior_for(parameters[1])]  
    # power-law norm
    priors += [bxa.create_uniform_prior_for(parameters[1])] 

    if zcdf is not None :
        parameters += [po.redshift]
        priors += [get_inv_zpdf_func(zcdf)]
    
    # Add norm. constants to parameters and define priors
    for sid in ids :
        if sid > 1 :
            cpar = get_model(sid).pars[0]
            parameters += [cpar]
            priors += [bxa.create_jeffreys_prior_for(cpar)]
            
    return parameters, bxa.create_prior_function(priors = priors)

       
def define_priors_wapo(ids, zcdf=None) :
    # Define array with the thawed parameters from the defined 
    # model of the first spectrum
    lognH = Parameter(modelname='intabs', name='logNH', val=22, min=20, max=25)
    intabs.nH = 10**(lognH - 22)    
    lognorm = Parameter(modelname='po', name='lognorm', val=-5, min=-30, max=0)
    po.norm = 10**lognorm

    parameters = [lognH, po.PhoIndex, lognorm]
        
    priors = []
    # intrinsic nH
    priors += [bxa.create_uniform_prior_for(parameters[0])]  
    # photon index (Distribution from n Nandra & Pounds (1994))
    priors += [bxa.create_gaussian_prior_for(parameters[1], 1.9, 0.15)]
    #priors += [bxa.create_uniform_prior_for(parameters[1])]  
    # power-law norm
    priors += [bxa.create_uniform_prior_for(parameters[2])] 

    if zcdf is not None :
        parameters += [po.redshift]
        priors += [get_inv_zpdf_func(zcdf)]
    
    # Add norm. constants to parameters and define priors
    for sid in ids :
        if sid > 1 :
            cpar = get_model(sid).pars[0]
            parameters += [cpar]
            priors += [bxa.create_jeffreys_prior_for(cpar)]
            
    return parameters, bxa.create_prior_function(priors = priors)


def define_priors_mekal(ids, zcdf=None) :
    # Define array with the thawed parameters from the defined 
    # model of the first spectrum

    lognorm = Parameter(modelname='thermal', name='lognorm', val=-5, min=-30, max=0)
    thermal.norm = 10**lognorm

    parameters = [thermal.kT, lognorm]
        
    priors = []
    # plasma temperature in keV
    priors += [bxa.create_uniform_prior_for(parameters[0])]
    # power-law norm
    priors += [bxa.create_uniform_prior_for(parameters[1])] 

    if zcdf is not None :
        parameters += [thermal.redshift]
        priors += [get_inv_zpdf_func(zcdf)]
        
    # Add norm. constants to parameters and define priors
    for sid in ids :
        if sid > 1 :
            cpar = get_model(sid).pars[0]
            parameters += [cpar]
            priors += [bxa.create_jeffreys_prior_for(cpar)]
            
    return parameters, bxa.create_prior_function(priors = priors)
    
    
def define_priors_wamekal(ids, zcdf=None) :
    # Define array with the thawed parameters from the defined 
    # model of the first spectrum

    lognH = Parameter(modelname='intabs', name='logNH', val=22, min=20, max=25)
    intabs.nH = 10**(lognH - 22)    
    lognorm = Parameter(modelname='thermal', name='lognorm', val=-5, min=-30, max=0)
    thermal.norm = 10**lognorm

    parameters = [lognH, thermal.kT, lognorm]
        
    priors = []
    # intrinsic nH
    priors += [bxa.create_uniform_prior_for(parameters[0])]  
    # plasma temperature in keV
    priors += [bxa.create_uniform_prior_for(parameters[1])]
    # power-law norm
    priors += [bxa.create_uniform_prior_for(parameters[2])] 

    if zcdf is not None :
        parameters += [thermal.redshift]
        priors += [get_inv_zpdf_func(zcdf)]
    
    # Add norm. constants to parameters and define priors
    for sid in ids :
        if sid > 1 :
            cpar = get_model(sid).pars[0]
            parameters += [cpar]
            priors += [bxa.create_jeffreys_prior_for(cpar)]
            
    return parameters, bxa.create_prior_function(priors = priors)
    

def define_priors_bb(ids, zcdf=None) :
    # Define array with the thawed parameters from the defined 
    # model of the first spectrum

    lognorm = Parameter(modelname='thermal', name='lognorm', val=-5, min=-30, max=0)
    thermal.norm = 10**lognorm

    parameters = [thermal.kt, lognorm]
        
    priors = []
    # plasma temperature in keV
    priors += [bxa.create_uniform_prior_for(parameters[0])]
    # norm
    priors += [bxa.create_uniform_prior_for(parameters[1])] 

    if zcdf is not None :
        parameters += [thermal.redshift]
        priors += [get_inv_zpdf_func(zcdf)]
        
    # Add norm. constants to parameters and define priors
    for sid in ids :
        if sid > 1 :
            cpar = get_model(sid).pars[0]
            parameters += [cpar]
            priors += [bxa.create_jeffreys_prior_for(cpar)]
            
    return parameters, bxa.create_prior_function(priors = priors)


def define_priors_wabb(ids, zcdf=None) :
    # Define array with the thawed parameters from the defined 
    # model of the first spectrum

    lognH = Parameter(modelname='intabs', name='logNH', val=22, min=20, max=25)
    intabs.nH = 10**(lognH - 22)    
    lognorm = Parameter(modelname='thermal', name='lognorm', val=-5, min=-30, max=0)
    thermal.norm = 10**lognorm

    parameters = [lognH, thermal.kt, lognorm]
        
    priors = []
    # intrinsic nH
    priors += [bxa.create_uniform_prior_for(parameters[0])]  
    # plasma temperature in keV
    priors += [bxa.create_uniform_prior_for(parameters[1])]
    # power-law norm
    priors += [bxa.create_uniform_prior_for(parameters[2])] 

    if zcdf is not None :
        parameters += [thermal.redshift]
        priors += [get_inv_zpdf_func(zcdf)]
    
    # Add norm. constants to parameters and define priors
    for sid in ids :
        if sid > 1 :
            cpar = get_model(sid).pars[0]
            parameters += [cpar]
            priors += [bxa.create_jeffreys_prior_for(cpar)]
            
    return parameters, bxa.create_prior_function(priors = priors)
    

def define_priors_wamekalpo(ids, zcdf=None) :
    # Define array with the thawed parameters from the defined 
    # model of the first spectrum

    lognH1 = Parameter(modelname='intabs1', name='logNH1', val=22, min=20, max=25)
    intabs1.nH = 10**(lognH1 - 22)
    lognH2 = Parameter(modelname='intabs2', name='logNH2', val=22, min=20, max=25)
    intabs2.nH = 10**(lognH2 - 22)    
    lognormt = Parameter(modelname='thermal', name='lognormt', val=-5, min=-30, max=0)
    thermal.norm = 10**lognormt
    lognormp = Parameter(modelname='po', name='lognormp', val=-5, min=-30, max=0)
    po.norm = 10**lognormp

    parameters = [lognH1, thermal.kT, lognH2, po.PhoIndex, lognormt, lognormp]
        
    priors = []
    # intrinsic nH 1
    priors += [bxa.create_uniform_prior_for(parameters[0])]  
    # plasma temperature in keV
    priors += [bxa.create_uniform_prior_for(parameters[1])]
    # intrinsic nH 2
    priors += [bxa.create_uniform_prior_for(parameters[2])]  
    # photon index
    priors += [bxa.create_gaussian_prior_for(parameters[3], 1.9, 0.15)]
    # thermal norm
    priors += [bxa.create_uniform_prior_for(parameters[4])]
    # po norm
    priors += [bxa.create_uniform_prior_for(parameters[5])]

    if zcdf is not None :
        parameters += [po.redshift]
        priors += [get_inv_zpdf_func(zcdf)]
    
    # Add norm. constants to parameters and define priors
    for sid in ids :
        if sid > 1 :
            cpar = get_model(sid).pars[0]
            parameters += [cpar]
            priors += [bxa.create_jeffreys_prior_for(cpar)]
            
    return parameters, bxa.create_prior_function(priors = priors)


def define_priors_wapopo(ids, zcdf=None) :
    # Define array with the thawed parameters from the defined 
    # model of the first spectrum

    lognH1 = Parameter(modelname='intabs1', name='logNH1', val=22, min=20, max=25)
    intabs1.nH = 10**(lognH1 - 22)
    lognH2 = Parameter(modelname='intabs2', name='logNH2', val=22, min=20, max=25)
    intabs2.nH = 10**(lognH2 - 22)        
    lognorm1 = Parameter(modelname='po1', name='lognorm1', val=-5, min=-30, max=0)
    po1.norm = 10**lognorm1
    lognorm2 = Parameter(modelname='po2', name='lognorm2', val=-5, min=-30, max=0)
    po2.norm = 10**lognorm2

    parameters = [lognH1, po1.PhoIndex, lognH2, po2.PhoIndex, lognorm1, lognorm2]
        
    priors = []
    # intrinsic nH1
    priors += [bxa.create_uniform_prior_for(parameters[0])]  
    # photon index 1
    priors += [bxa.create_gaussian_prior_for(parameters[1], 1.9, 0.15)]
    # intrinsic nH2
    priors += [bxa.create_uniform_prior_for(parameters[2])] 
    # photon index 2
    priors += [bxa.create_gaussian_prior_for(parameters[3], 1.9, 0.15)]
    # po norm 1
    priors += [bxa.create_uniform_prior_for(parameters[4])] 
    # po norm 2
    priors += [bxa.create_uniform_prior_for(parameters[5])] 

    if zcdf is not None :
        parameters += [po2.redshift]
        priors += [get_inv_zpdf_func(zcdf)]
    
    # Add norm. constants to parameters and define priors
    for sid in ids :
        if sid > 1 :
            cpar = get_model(sid).pars[0]
            parameters += [cpar]
            priors += [bxa.create_jeffreys_prior_for(cpar)]
            
    return parameters, bxa.create_prior_function(priors = priors)
    
    
def define_priors_wabbpo(ids, zcdf=None) :
    # Define array with the thawed parameters from the defined 
    # model of the first spectrum

    lognH = Parameter(modelname='intabs', name='logNH', val=22, min=20, max=25)
    intabs.nH = 10**(lognH - 22)
    lognormt = Parameter(modelname='thermal', name='lognormt', val=-5, min=-30, max=0)
    thermal.norm = 10**lognormt
    lognormp = Parameter(modelname='po', name='lognormp', val=-5, min=-30, max=0)
    po.norm = 10**lognormp

    parameters = [lognH, thermal.kt, po.PhoIndex, lognormt, lognormp]
        
    priors = []
    # intrinsic nH
    priors += [bxa.create_uniform_prior_for(parameters[0])]  
    # plasma temperature in keV
    priors += [bxa.create_uniform_prior_for(parameters[1])]
    # photon index
    priors += [bxa.create_gaussian_prior_for(parameters[2], 1.9, 0.15)]
    # thermal norm
    priors += [bxa.create_uniform_prior_for(parameters[3])]
    # po norm
    priors += [bxa.create_uniform_prior_for(parameters[4])]

    if zcdf is not None :
        parameters += [po.redshift]
        priors += [get_inv_zpdf_func(zcdf)]
    
    # Add norm. constants to parameters and define priors
    for sid in ids :
        if sid > 1 :
            cpar = get_model(sid).pars[0]
            parameters += [cpar]
            priors += [bxa.create_jeffreys_prior_for(cpar)]
            
    return parameters, bxa.create_prior_function(priors = priors)
    

