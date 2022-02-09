#!/usr/bin/env python
# Filename: models.py

from sherpa.astro.ui import *


def define_model_po(ids, z, nhgal) :

    for sid in ids :
        if sid == 1 :
            model_str = "xstbabs.galabs*xszpowerlw.po"
            set_source(sid, model_str)
        else :
            model_str = "const1d.c" + str(sid) + "*galabs*po"
            set_source(sid, model_str)
            par = "c" + str(sid) + ".c0"
            set_par(par, val=1, min=0, max=10000)
    
    # Define freezed parameters and limits (for priors)
    galabs.nH = nhgal / 1e22
    freeze(galabs.nH)
    
    po.PhoIndex = 2
    po.PhoIndex.min = 0
    po.PhoIndex.max = 10
    po.norm.min = 1e-30
    po.norm.max = 1

    if z is None :
        thaw(po.redshift)
    else :
        po.redshift = z


def define_model_wapo(ids, z, nhgal) :

    for sid in ids :
        if sid == 1 :
            model_str = "xstbabs.galabs*xsztbabs.intabs*xszpowerlw.po"
            set_source(sid, model_str)
        else :
            model_str = "const1d.c" + str(sid) + "*galabs*intabs*po"
            set_source(sid, model_str)
            par = "c" + str(sid) + ".c0"
            set_par(par, val=1, min=0, max=10000)
    
    # Define freezed parameters and limits (for priors)
    galabs.nH = nhgal / 1e22
    freeze(galabs.nH)
    
    intabs.nH.min = 0
    intabs.nH.max = 1000

    po.PhoIndex = 2
    po.PhoIndex.min = 0
    po.PhoIndex.max = 10
    po.norm.min = 1e-30
    po.norm.max = 1

    if z is None :
        intabs.redshift = po.redshift
        thaw(po.redshift)
    else :
        intabs.redshift = z
        po.redshift = z


def define_model_wamekal(ids, z, nhgal) :

    for sid in ids :
        if sid == 1 :
            model_str = "xstbabs.galabs*xsztbabs.intabs*xsmekal.thermal"
            set_source(sid, model_str)
        else :
            model_str = "const1d.c" + str(sid) + "*galabs*intabs*thermal"
            set_source(sid, model_str)
            par = "c" + str(sid) + ".c0"
            set_par(par, val=1, min=0, max=10000)
    
    # Define freezed parameters and limits (for priors)
    galabs.nH = nhgal / 1e22
    freeze(galabs.nH)
    
    intabs.nH.min = 0
    intabs.nH.max = 1000

    thermal.kT = 0.5
    thermal.kT.min = 0.08
    thermal.kT.max = 20
    thermal.norm.min = 1e-30
    thermal.norm.max = 1

    if z is None :
        intabs.redshift = thermal.redshift
        thaw(thermal.redshift)
    else :
        intabs.redshift = z
        thermal.redshift = z
        
        
def define_model_mekal(ids, z, nhgal) :

    for sid in ids :
        if sid == 1 :
            model_str = "xstbabs.galabs*xsmekal.thermal"
            set_source(sid, model_str)
        else :
            model_str = "const1d.c" + str(sid) + "*galabs*thermal"
            set_source(sid, model_str)
            par = "c" + str(sid) + ".c0"
            set_par(par, val=1, min=0, max=10000)
    
    # Define freezed parameters and limits (for priors)
    galabs.nH = nhgal / 1e22
    freeze(galabs.nH)

    thermal.kT = 0.5
    thermal.kT.min = 0.08
    thermal.kT.max = 20
    thermal.norm.min = 1e-30
    thermal.norm.max = 1

    if z is None :
        thaw(thermal.redshift)
    else :
        thermal.redshift = z
        

def define_model_wabb(ids, z, nhgal) :

    for sid in ids :
        if sid == 1 :
            model_str = "xstbabs.galabs*xsztbabs.intabs*xszbbody.thermal"
            set_source(sid, model_str)
        else :
            model_str = "const1d.c" + str(sid) + "*galabs*intabs*thermal"
            set_source(sid, model_str)
            par = "c" + str(sid) + ".c0"
            set_par(par, val=1, min=0, max=10000)
    
    # Define freezed parameters and limits (for priors)
    galabs.nH = nhgal / 1e22
    freeze(galabs.nH)
    
    intabs.nH.min = 0
    intabs.nH.max = 1000

    thermal.kt = 0.5
    thermal.kt.min = 0.01
    thermal.kt.max = 10
    thermal.norm.min = 1e-30
    thermal.norm.max = 1

    if z is None :
        intabs.redshift = thermal.redshift
        thaw(thermal.redshift)
    else :
        intabs.redshift = z
        thermal.redshift = z


def define_model_bb(ids, z, nhgal) :

    for sid in ids :
        if sid == 1 :
            model_str = "xstbabs.galabs*xszbbody.thermal"
            set_source(sid, model_str)
        else :
            model_str = "const1d.c" + str(sid) + "*galabs*thermal"
            set_source(sid, model_str)
            par = "c" + str(sid) + ".c0"
            set_par(par, val=1, min=0, max=10000)
    
    # Define freezed parameters and limits (for priors)
    galabs.nH = nhgal / 1e22
    freeze(galabs.nH)
    
    thermal.kt = 0.5
    thermal.kt.min = 0.01
    thermal.kt.max = 10
    thermal.norm.min = 1e-30
    thermal.norm.max = 1

    if z is None :
        thaw(thermal.redshift)
    else :
        thermal.redshift = z
        

def define_model_wamekalpo(ids, z, nhgal) :

    for sid in ids :
        if sid == 1 :
            model_str = "xstbabs.galabs*xsztbabs.intabs1*(xsmekal.thermal + xsztbabs.intabs2*xszpowerlw.po)"
            set_source(sid, model_str)
        else :
            model_str = "const1d.c" + str(sid) + "*galabs*intabs1*(thermal + intabs2*po)"
            set_source(sid, model_str)
            par = "c" + str(sid) + ".c0"
            set_par(par, val=1, min=0, max=10000)
    
    # Define freezed parameters and limits (for priors)
    galabs.nH = nhgal / 1e22
    freeze(galabs.nH)
    
    intabs1.nH.min = 0
    intabs1.nH.max = 1000

    intabs2.nH.min = 0
    intabs2.nH.max = 1000

    thermal.kT = 0.5
    thermal.kT.min = 0.08
    thermal.kT.max = 20
    thermal.norm.min = 1e-30
    thermal.norm.max = 1

    po.PhoIndex = 2
    po.PhoIndex.min = 0
    po.PhoIndex.max = 10
    po.norm.min = 1e-30
    po.norm.max = 1
    
    if z is None :
        intabs1.redshift = po.redshift
        intabs2.redshift = po.redshift
        thermal.redshift = po.redshift
        thaw(po.redshift)
    else :
        intabs1.redshift = z
        intabs2.redshift = z
        thermal.redshift = z
        po.redshift = z


def define_model_wapopo(ids, z, nhgal) :

    for sid in ids :
        if sid == 1 :
            model_str = "xstbabs.galabs*xsztbabs.intabs1*(xszpowerlw.po1 + xsztbabs.intabs2*xszpowerlw.po2)"
            set_source(sid, model_str)
        else :
            model_str = "const1d.c" + str(sid) + "*galabs*intabs1*(po1 + intabs2*po2)"
            set_source(sid, model_str)
            par = "c" + str(sid) + ".c0"
            set_par(par, val=1, min=0, max=10000)
    
    # Define freezed parameters and limits (for priors)
    galabs.nH = nhgal / 1e22
    freeze(galabs.nH)
    
    intabs1.nH.min = 0
    intabs1.nH.max = 1000

    intabs2.nH.min = 0
    intabs2.nH.max = 1000

    po1.PhoIndex = 2
    po1.PhoIndex.min = 0
    po1.PhoIndex.max = 10
    po1.norm.min = 1e-30
    po1.norm.max = 1

    po2.PhoIndex = 2
    po2.PhoIndex.min = 0
    po2.PhoIndex.max = 10
    po2.norm.min = 1e-30
    po2.norm.max = 1
    
    if z is None :
        intabs1.redshift = po2.redshift
        intabs2.redshift = po2.redshift
        po1.redshift = po2.redshift
        thaw(po2.redshift)
    else :
        intabs1.redshift = z
        intabs2.redshift = z
        po1.redshift = z
        po2.redshift = z


def define_model_wabbpo(ids, z, nhgal) :

    for sid in ids :
        if sid == 1 :
            model_str = "xstbabs.galabs*xsztbabs.intabs*(xszbbody.thermal + xszpowerlw.po)"
            set_source(sid, model_str)
        else :
            model_str = "const1d.c" + str(sid) + "*galabs*intabs*(thermal + po)"
            set_source(sid, model_str)
            par = "c" + str(sid) + ".c0"
            set_par(par, val=1, min=0, max=10000)
    
    # Define freezed parameters and limits (for priors)
    galabs.nH = nhgal / 1e22
    freeze(galabs.nH)
    
    intabs.nH.min = 0
    intabs.nH.max = 1000

    thermal.kt = 0.5
    thermal.kt.min = 0.01
    thermal.kt.max = 10
    thermal.norm.min = 1e-30
    thermal.norm.max = 1

    po.PhoIndex = 2
    po.PhoIndex.min = 0
    po.PhoIndex.max = 10
    po.norm.min = 1e-30
    po.norm.max = 1
    
    if z is None :
        intabs.redshift = po.redshift
        thermal.redshift = po.redshift
        thaw(po.redshift)
    else :
        intabs.redshift = z
        thermal.redshift = z
        po.redshift = z

