#!/usr/bin/env python
"""
Created on Tue Dec  5 12:46:05 2017

@author: ruizca
"""

import json
import numpy as np

import sqlite3
from sqlite3 import Error

from astropy.table import Table


def create_connection(db_file) :
    """ 
    create a database connection to the SQLite database
    specified by db_file
    :param db_file: database file
    :return: Connection object or None
    """
    try:
        conn = sqlite3.connect(db_file, isolation_level=None)
        
        return conn
        
    except Error as e:
        print(e)

    return None 


def create_table(conn, create_table_sql) :
    """ 
    create a table from the create_table_sql statement
    :param conn: Connection object
    :param create_table_sql: a CREATE TABLE statement
    :return:
    """
    try:
        c = conn.cursor()
        c.execute(create_table_sql)

    except Error as e:
        print(e)
        

def add_detection(conn, detection):
    """
    Add a new detection into the detection_data table
    :param conn:
    :param detection:
    :return: detection id
    """

    sql = "INSERT INTO detections_data(detid, obsid, srcnum, \
            srcid, iauname, srcra, srcdec, nHgal, \
            zsp, zph, zphErr, fcounts, scounts, hcounts) \
            VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)"

    c = conn.cursor()
    c.execute(sql, detection)

    return c.lastrowid


def add_spectrum(conn, spectrum):
    """
    Add a new spectrum into the spectra_data table
    :param conn:
    :param spectrum:
    :return: spectrum id
    """

    sql = "INSERT INTO spectra_data(det_id, expid, \
            detector, counts, good) \
            VALUES(?,?,?,?,?)"

    c = conn.cursor()
    c.execute(sql, spectrum)

    return c.lastrowid


def add_fit(conn, fit):
    """
    Add a new fit into the fits_data table
    :param conn:
    :param fit:
    :return: fit id
    """

    sql = "INSERT INTO fits_data(det_id, model_id, \
            npar, dof, wstat, ks, pvalue, \
            lnev, lnevErr, lnZ, lnZErr) \
            VALUES(?,?,?,?,?,?,?,?,?,?,?)"

    c = conn.cursor()
    c.execute(sql, fit)

    return c.lastrowid


def sql_pars(model) :
    
    if model == "po" :
        sql = "INSERT INTO parameters_po(det_id, fit_id, \
               PhoIndex_median, PhoIndex_mode, PhoIndex_min, PhoIndex_max, \
               lognorm_median, lognorm_mode, lognorm_min, lognorm_max, \
               z_median, z_mode, z_min, z_max) \
               VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)"

    elif model == "wapo" :
        sql = "INSERT INTO parameters_wapo(det_id, fit_id, \
               lognH_median, lognH_mode, lognH_min, lognH_max, \
               PhoIndex_median, PhoIndex_mode, PhoIndex_min, PhoIndex_max, \
               lognorm_median, lognorm_mode, lognorm_min, lognorm_max, \
               z_median, z_mode, z_min, z_max) \
               VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"

    elif model == "mekal" :
        sql = "INSERT INTO parameters_mekal(det_id, fit_id, \
               kT_median, kT_mode, kT_min, kT_max, \
               lognorm_median, lognorm_mode, lognorm_min, lognorm_max, \
               z_median, z_mode, z_min, z_max) \
               VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)"

    elif model == "wamekal" :
        sql = "INSERT INTO parameters_wamekal(det_id, fit_id, \
               lognH_median, lognH_mode, lognH_min, lognH_max, \
               kT_median, kT_mode, kT_min, kT_max, \
               lognorm_median, lognorm_mode, lognorm_min, lognorm_max, \
               z_median, z_mode, z_min, z_max) \
               VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"

    elif model == "bb" :
        sql = "INSERT INTO parameters_bb(det_id, fit_id, \
               kT_median, kT_mode, kT_min, kT_max, \
               lognorm_median, lognorm_mode, lognorm_min, lognorm_max, \
               z_median, z_mode, z_min, z_max) \
               VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
        
    elif model == "wabb" :
        sql = "INSERT INTO parameters_wabb(det_id, fit_id, \
               lognH_median, lognH_mode, lognH_min, lognH_max, \
               kT_median, kT_mode, kT_min, kT_max, \
               lognorm_median, lognorm_mode, lognorm_min, lognorm_max, \
               z_median, z_mode, z_min, z_max) \
               VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"

    elif model == "wamekalpo" :
        sql = "INSERT INTO parameters_wamekalpo(det_id, fit_id, \
               lognH1_median, lognH1_mode, lognH1_min, lognH1_max, \
               kT_median, kT_mode, kT_min, kT_max, \
               lognH2_median, lognH2_mode, lognH2_min, lognH2_max, \
               PhoIndex_median, PhoIndex_mode, PhoIndex_min, PhoIndex_max, \
               lognormth_median, lognormth_mode, lognormth_min, lognormth_max, \
               lognormpo_median, lognormpo_mode, lognormpo_min, lognormpo_max, \
               z_median, z_mode, z_min, z_max) \
               VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"

    elif model == "wapopo" :
        sql = "INSERT INTO parameters_wapopo(det_id, fit_id, \
               lognH1_median, lognH1_mode, lognH1_min, lognH1_max, \
               PhoIndex1_median, PhoIndex1_mode, PhoIndex1_min, PhoIndex1_max, \
               lognH2_median, lognH2_mode, lognH2_min, lognH2_max, \
               PhoIndex2_median, PhoIndex2_mode, PhoIndex2_min, PhoIndex2_max, \
               lognormpo1_median, lognormpo1_mode, lognormpo1_min, lognormpo1_max, \
               lognormpo2_median, lognormpo2_mode, lognormpo2_min, lognormpo2_max, \
               z_median, z_mode, z_min, z_max) \
               VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"

    elif model == "wabbpo" :
        sql = "INSERT INTO parameters_wabbpo(det_id, fit_id, \
               lognH_median, lognH_mode, lognH_min, lognH_max, \
               kT_median, kT_mode, kT_min, kT_max, \
               PhoIndex_median, PhoIndex_mode, PhoIndex_min, PhoIndex_max, \
               lognormth_median, lognormth_mode, lognormth_min, lognormth_max, \
               lognormpo_median, lognormpo_mode, lognormpo_min, lognormpo_max, \
               z_median, z_mode, z_min, z_max) \
               VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
        
    else :
        raise ValueError("unknown model!!!")
    
    return sql


def add_pars(conn, pars, model):
    """
    Add new best fit parameters into the 
    corresponding model table
    :param conn:
    :param pars:
    :param model:
    :return: fit id
    """

    sql = sql_pars(model)
    
    c = conn.cursor()
    c.execute(sql, pars)

    return c.lastrowid


def add_flux(conn, flux):
    """
    Add new best fit fluxes into the table fits_fluxes
    :param conn:
    :param flux:
    :return: flux id
    """

    sql = "INSERT INTO fits_fluxes(det_id, fit_id, \
           sflux_obs_median, sflux_obs_mode, sflux_obs_min, sflux_obs_max, \
           sflux_int_median, sflux_int_mode, sflux_int_min, sflux_int_max, \
           hflux_obs_median, hflux_obs_mode, hflux_obs_min, hflux_obs_max, \
           hflux_int_median, hflux_int_mode, hflux_int_min, hflux_int_max, \
           slumin_median, slumin_mode, slumin_min, slumin_max, \
           hlumin_median, hlumin_mode, hlumin_min, hlumin_max) \
           VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
    
    c = conn.cursor()
    c.execute(sql, flux)

    return c.lastrowid

        
def get_pars(pars_info, model, thawedz=False) :

    if model == "po" :
        
        pars = (pars_info["po.PhoIndex"]["mode"], pars_info["po.PhoIndex"]["median"],
                pars_info["po.PhoIndex"]["parmin"], pars_info["po.PhoIndex"]["parmax"],
                pars_info["po.lognorm"]["mode"], pars_info["po.lognorm"]["median"],
                pars_info["po.lognorm"]["parmin"], pars_info["po.lognorm"]["parmax"])
                
        if thawedz :
            pars += (pars_info["po.redshift"]["mode"], pars_info["po.redshift"]["median"],
                     pars_info["po.redshift"]["parmin"], pars_info["po.redshift"]["parmax"])
        else :
            pars += (None, None, None, None)

        return pars

    
    elif model == "wapo" :
        
        pars = (pars_info["intabs.logNH"]["mode"], pars_info["intabs.logNH"]["median"],
                pars_info["intabs.logNH"]["parmin"], pars_info["intabs.logNH"]["parmax"], 
                pars_info["po.PhoIndex"]["mode"], pars_info["po.PhoIndex"]["median"],
                pars_info["po.PhoIndex"]["parmin"], pars_info["po.PhoIndex"]["parmax"],
                pars_info["po.lognorm"]["mode"], pars_info["po.lognorm"]["median"],
                pars_info["po.lognorm"]["parmin"], pars_info["po.lognorm"]["parmax"])
                
        if thawedz :
            pars += (pars_info["po.redshift"]["mode"], pars_info["po.redshift"]["median"],
                     pars_info["po.redshift"]["parmin"], pars_info["po.redshift"]["parmax"])
        else :
            pars += (None, None, None, None)

        return pars


    elif model == "mekal" :

        pars = (pars_info["thermal.kT"]["mode"], pars_info["thermal.kT"]["median"],
                pars_info["thermal.kT"]["parmin"], pars_info["thermal.kT"]["parmax"],
                pars_info["thermal.lognorm"]["mode"], pars_info["thermal.lognorm"]["median"],
                pars_info["thermal.lognorm"]["parmin"], pars_info["thermal.lognorm"]["parmax"])
                
        if thawedz :
            pars += (pars_info["thermal.redshift"]["mode"], pars_info["thermal.redshift"]["median"],
                     pars_info["thermal.redshift"]["parmin"], pars_info["thermal.redshift"]["parmax"])
        else :
            pars += (None, None, None, None)

        return pars

        
    elif model == "wamekal" :

        pars = (pars_info["intabs.logNH"]["mode"], pars_info["intabs.logNH"]["median"],
                pars_info["intabs.logNH"]["parmin"], pars_info["intabs.logNH"]["parmax"], 
                pars_info["thermal.kT"]["mode"], pars_info["thermal.kT"]["median"],
                pars_info["thermal.kT"]["parmin"], pars_info["thermal.kT"]["parmax"],
                pars_info["thermal.lognorm"]["mode"], pars_info["thermal.lognorm"]["median"],
                pars_info["thermal.lognorm"]["parmin"], pars_info["thermal.lognorm"]["parmax"])
                
        if thawedz :
            pars += (pars_info["thermal.redshift"]["mode"], pars_info["thermal.redshift"]["median"],
                     pars_info["thermal.redshift"]["parmin"], pars_info["thermal.redshift"]["parmax"])
        else :
            pars += (None, None, None, None)

        return pars


    elif model == "bb" :

        pars = (pars_info["thermal.kt"]["mode"], pars_info["thermal.kt"]["median"],
                pars_info["thermal.kt"]["parmin"], pars_info["thermal.kt"]["parmax"],
                pars_info["thermal.lognorm"]["mode"], pars_info["thermal.lognorm"]["median"],
                pars_info["thermal.lognorm"]["parmin"], pars_info["thermal.lognorm"]["parmax"])
                
        if thawedz :
            pars += (pars_info["thermal.redshift"]["mode"], pars_info["thermal.redshift"]["median"],
                     pars_info["thermal.redshift"]["parmin"], pars_info["thermal.redshift"]["parmax"])
        else :
            pars += (None, None, None, None)

        return pars

        
    elif model == "wabb" :

        pars = (pars_info["intabs.logNH"]["mode"], pars_info["intabs.logNH"]["median"],
                pars_info["intabs.logNH"]["parmin"], pars_info["intabs.logNH"]["parmax"], 
                pars_info["thermal.kt"]["mode"], pars_info["thermal.kt"]["median"],
                pars_info["thermal.kt"]["parmin"], pars_info["thermal.kt"]["parmax"],
                pars_info["thermal.lognorm"]["mode"], pars_info["thermal.lognorm"]["median"],
                pars_info["thermal.lognorm"]["parmin"], pars_info["thermal.lognorm"]["parmax"])
                
        if thawedz :
            pars += (pars_info["thermal.redshift"]["mode"], pars_info["thermal.redshift"]["median"],
                     pars_info["thermal.redshift"]["parmin"], pars_info["thermal.redshift"]["parmax"])
        else :
            pars += (None, None, None, None)

        return pars


    elif model == "wamekalpo" :

        pars = (pars_info["intabs1.logNH1"]["mode"], pars_info["intabs1.logNH1"]["median"],
                pars_info["intabs1.logNH1"]["parmin"], pars_info["intabs1.logNH1"]["parmax"], 
                pars_info["thermal.kT"]["mode"], pars_info["thermal.kT"]["median"],
                pars_info["thermal.kT"]["parmin"], pars_info["thermal.kT"]["parmax"],
                pars_info["intabs2.logNH2"]["mode"], pars_info["intabs2.logNH2"]["median"],
                pars_info["intabs2.logNH2"]["parmin"], pars_info["intabs2.logNH2"]["parmax"], 
                pars_info["po.PhoIndex"]["mode"], pars_info["po.PhoIndex"]["median"],
                pars_info["po.PhoIndex"]["parmin"], pars_info["po.PhoIndex"]["parmax"],
                pars_info["thermal.lognormt"]["mode"], pars_info["thermal.lognormt"]["median"],
                pars_info["thermal.lognormt"]["parmin"], pars_info["thermal.lognormt"]["parmax"],
                pars_info["po.lognormp"]["mode"], pars_info["po.lognormp"]["median"],
                pars_info["po.lognormp"]["parmin"], pars_info["po.lognormp"]["parmax"])

        if thawedz :
            pars += (pars_info["po.redshift"]["mode"], pars_info["po.redshift"]["median"],
                     pars_info["po.redshift"]["parmin"], pars_info["po.redshift"]["parmax"])
        else :
            pars += (None, None, None, None)

        return pars


    elif model == "wapopo" :
        pars = (pars_info["intabs1.logNH1"]["mode"], pars_info["intabs1.logNH1"]["median"],
                pars_info["intabs1.logNH1"]["parmin"], pars_info["intabs1.logNH1"]["parmax"], 
                pars_info["po1.PhoIndex"]["mode"], pars_info["po1.PhoIndex"]["median"],
                pars_info["po1.PhoIndex"]["parmin"], pars_info["po1.PhoIndex"]["parmax"],
                pars_info["intabs2.logNH2"]["mode"], pars_info["intabs2.logNH2"]["median"],
                pars_info["intabs2.logNH2"]["parmin"], pars_info["intabs2.logNH2"]["parmax"], 
                pars_info["po2.PhoIndex"]["mode"], pars_info["po2.PhoIndex"]["median"],
                pars_info["po2.PhoIndex"]["parmin"], pars_info["po2.PhoIndex"]["parmax"],
                pars_info["po1.lognorm1"]["mode"], pars_info["po1.lognorm1"]["median"],
                pars_info["po1.lognorm1"]["parmin"], pars_info["po1.lognorm1"]["parmax"],
                pars_info["po2.lognorm2"]["mode"], pars_info["po2.lognorm2"]["median"],
                pars_info["po2.lognorm2"]["parmin"], pars_info["po2.lognorm2"]["parmax"])

        if thawedz :
            pars += (pars_info["po2.redshift"]["mode"], pars_info["po2.redshift"]["median"],
                     pars_info["po2.redshift"]["parmin"], pars_info["po2.redshift"]["parmax"])
        else :
            pars += (None, None, None, None)

        return pars


    elif model == "wabbpo" :
        pars = (pars_info["intabs1.logNH1"]["mode"], pars_info["intabs1.logNH1"]["median"],
                pars_info["intabs1.logNH1"]["parmin"], pars_info["intabs1.logNH1"]["parmax"], 
                pars_info["thermal.kt"]["mode"], pars_info["thermal.kt"]["median"],
                pars_info["thermal.kt"]["parmin"], pars_info["thermal.kt"]["parmax"],
                pars_info["intabs2.logNH2"]["mode"], pars_info["intabs2.logNH2"]["median"],
                pars_info["intabs2.logNH2"]["parmin"], pars_info["intabs2.logNH2"]["parmax"], 
                pars_info["po.PhoIndex"]["mode"], pars_info["po.PhoIndex"]["median"],
                pars_info["po.PhoIndex"]["parmin"], pars_info["po.PhoIndex"]["parmax"])                

        if thawedz :
            pars += (pars_info["po.redshift"]["mode"], pars_info["po.redshift"]["median"],
                     pars_info["po.redshift"]["parmin"], pars_info["po.redshift"]["parmax"])
        else :
            pars += (None, None, None, None)

        return pars

        
    else :
        raise ValueError("unknown model!!!")



        
def main(dbfile, detfile, model_names, add_models=False) :

    # Connect to db (or create it if doesn't exist)
    conn = create_connection(dbfile)

    if conn is not None:
        ### Create tables
        # Table for detections data
        create_table(conn,
                     "CREATE TABLE IF NOT EXISTS detections_data (\
                         id integer PRIMARY KEY, \
                         detid integer NOT NULL UNIQUE, \
                         obsid text NOT NULL, \
                         srcnum integer NOT NULL, \
                         srcid integer NOT NULL, \
                         iauname text NOT NULL, \
                         srcra float NOT NULL, \
                         srcdec float NOT NULL, \
                         nHgal float NOT NULL, \
                         zsp float, \
                         zph float, \
                         zphErr float, \
                         fcounts float, \
                         scounts float, \
                         hcounts float);"
                    )
    
        # Table for spectra data
        create_table(conn,
                     "CREATE TABLE IF NOT EXISTS spectra_data (\
                         id integer PRIMARY KEY, \
                         det_id integer NOT NULL, \
                         expid text NOT NULL, \
                         detector text NOT NULL, \
                         counts float NOT NULL, \
                         good bit NOT NULL, \
                         FOREIGN KEY (det_id) REFERENCES detections_data (id));"
                    )

        # Table for model names
        create_table(conn,
                     "CREATE TABLE IF NOT EXISTS models (\
                         id integer PRIMARY KEY, \
                         model text NOT NULL UNIQUE);"
                    )

        # Table for fits data
        create_table(conn,
                     "CREATE TABLE IF NOT EXISTS fits_data (\
                         id integer PRIMARY KEY, \
                         det_id integer NOT NULL, \
                         model_id integer NOT NULL, \
                         npar integer NOT NULL, \
                         dof float NOT NULL, \
                         wstat float NOT NULL, \
                         ks float NOT NULL, \
                         pvalue float NOT NULL, \
                         lnev float NOT NULL, \
                         lnevErr float NOT NULL, \
                         lnZ float NOT NULL, \
                         lnZErr float NOT NULL, \
                         FOREIGN KEY (model_id) REFERENCES models (id), \
                         FOREIGN KEY (det_id) REFERENCES detections_data (id));"
                    )

        # Table for po model best-fit parameters
        create_table(conn,
                     "CREATE TABLE IF NOT EXISTS parameters_po (\
                         id integer PRIMARY KEY, \
                         det_id integer NOT NULL, \
                         fit_id integer NOT NULL, \
                         PhoIndex_median float NOT NULL, \
                         PhoIndex_mode float NOT NULL, \
                         PhoIndex_min float NOT NULL, \
                         PhoIndex_max float NOT NULL, \
                         lognorm_median float NOT NULL, \
                         lognorm_mode float NOT NULL, \
                         lognorm_min float NOT NULL, \
                         lognorm_max float NOT NULL, \
                         z_median float, \
                         z_mode float, \
                         z_min float, \
                         z_max float, \
                         FOREIGN KEY (fit_id) REFERENCES fits_data (id), \
                         FOREIGN KEY (det_id) REFERENCES detections_data (id));"
                    ) 

        # Table for wapo model best-fit parameters
        create_table(conn,
                     "CREATE TABLE IF NOT EXISTS parameters_wapo (\
                         id integer PRIMARY KEY, \
                         det_id integer NOT NULL, \
                         fit_id integer NOT NULL, \
                         lognH_median float NOT NULL, \
                         lognH_mode float NOT NULL, \
                         lognH_min float NOT NULL, \
                         lognH_max float NOT NULL, \
                         PhoIndex_median float NOT NULL, \
                         PhoIndex_mode float NOT NULL, \
                         PhoIndex_min float NOT NULL, \
                         PhoIndex_max float NOT NULL, \
                         lognorm_median float NOT NULL, \
                         lognorm_mode float NOT NULL, \
                         lognorm_min float NOT NULL, \
                         lognorm_max float NOT NULL, \
                         z_median float, \
                         z_mode float, \
                         z_min float, \
                         z_max float, \
                         FOREIGN KEY (fit_id) REFERENCES fits_data (id), \
                         FOREIGN KEY (det_id) REFERENCES detections_data (id));"
                    ) 

        # Table for mekal model best-fit parameters
        create_table(conn,
                     "CREATE TABLE IF NOT EXISTS parameters_mekal (\
                         id integer PRIMARY KEY, \
                         det_id integer NOT NULL, \
                         fit_id integer NOT NULL, \
                         kT_median float NOT NULL, \
                         kT_mode float NOT NULL, \
                         kT_min float NOT NULL, \
                         kT_max float NOT NULL, \
                         lognorm_median float NOT NULL, \
                         lognorm_mode float NOT NULL, \
                         lognorm_min float NOT NULL, \
                         lognorm_max float NOT NULL, \
                         z_median float, \
                         z_mode float, \
                         z_min float, \
                         z_max float, \
                         FOREIGN KEY (fit_id) REFERENCES fits_data (id), \
                         FOREIGN KEY (det_id) REFERENCES detections_data (id));"
                    ) 

        # Table for wamekal model best-fit parameters
        create_table(conn,
                     "CREATE TABLE IF NOT EXISTS parameters_wamekal (\
                         id integer PRIMARY KEY, \
                         det_id integer NOT NULL, \
                         fit_id integer NOT NULL, \
                         lognH_median float NOT NULL, \
                         lognH_mode float NOT NULL, \
                         lognH_min float NOT NULL, \
                         lognH_max float NOT NULL, \
                         kT_median float NOT NULL, \
                         kT_mode float NOT NULL, \
                         kT_min float NOT NULL, \
                         kT_max float NOT NULL, \
                         lognorm_median float NOT NULL, \
                         lognorm_mode float NOT NULL, \
                         lognorm_min float NOT NULL, \
                         lognorm_max float NOT NULL, \
                         z_median float, \
                         z_mode float, \
                         z_min float, \
                         z_max float, \
                         FOREIGN KEY (fit_id) REFERENCES fits_data (id), \
                         FOREIGN KEY (det_id) REFERENCES detections_data (id));"
                    ) 

        # Table for wamekalpo model best-fit parameters
        create_table(conn,
                     "CREATE TABLE IF NOT EXISTS parameters_wamekalpo (\
                         id integer PRIMARY KEY, \
                         det_id integer NOT NULL, \
                         fit_id integer NOT NULL, \
                         lognH1_median float NOT NULL, \
                         lognH1_mode float NOT NULL, \
                         lognH1_min float NOT NULL, \
                         lognH1_max float NOT NULL, \
                         kT_median float NOT NULL, \
                         kT_mode float NOT NULL, \
                         kT_min float NOT NULL, \
                         kT_max float NOT NULL, \
                         lognH2_median float NOT NULL, \
                         lognH2_mode float NOT NULL, \
                         lognH2_min float NOT NULL, \
                         lognH2_max float NOT NULL, \
                         PhoIndex_median float NOT NULL, \
                         PhoIndex_mode float NOT NULL, \
                         PhoIndex_min float NOT NULL, \
                         PhoIndex_max float NOT NULL, \
                         lognormth_median float NOT NULL, \
                         lognormth_mode float NOT NULL, \
                         lognormth_min float NOT NULL, \
                         lognormth_max float NOT NULL, \
                         lognormpo_median float NOT NULL, \
                         lognormpo_mode float NOT NULL, \
                         lognormpo_min float NOT NULL, \
                         lognormpo_max float NOT NULL, \
                         z_median float, \
                         z_mode float, \
                         z_min float, \
                         z_max float, \
                         FOREIGN KEY (fit_id) REFERENCES fits_data (id), \
                         FOREIGN KEY (det_id) REFERENCES detections_data (id));"
                    ) 

        # Table for wapopo model best-fit parameters
        create_table(conn,
                     "CREATE TABLE IF NOT EXISTS parameters_wapopo (\
                         id integer PRIMARY KEY, \
                         det_id integer NOT NULL, \
                         fit_id integer NOT NULL, \
                         lognH1_median float NOT NULL, \
                         lognH1_mode float NOT NULL, \
                         lognH1_min float NOT NULL, \
                         lognH1_max float NOT NULL, \
                         PhoIndex1_median float NOT NULL, \
                         PhoIndex1_mode float NOT NULL, \
                         PhoIndex1_min float NOT NULL, \
                         PhoIndex1_max float NOT NULL, \
                         lognH2_median float NOT NULL, \
                         lognH2_mode float NOT NULL, \
                         lognH2_min float NOT NULL, \
                         lognH2_max float NOT NULL, \
                         PhoIndex2_median float NOT NULL, \
                         PhoIndex2_mode float NOT NULL, \
                         PhoIndex2_min float NOT NULL, \
                         PhoIndex2_max float NOT NULL, \
                         lognormpo1_median float NOT NULL, \
                         lognormpo1_mode float NOT NULL, \
                         lognormpo1_min float NOT NULL, \
                         lognormpo1_max float NOT NULL, \
                         lognormpo2_median float NOT NULL, \
                         lognormpo2_mode float NOT NULL, \
                         lognormpo2_min float NOT NULL, \
                         lognormpo2_max float NOT NULL, \
                         z_median float, \
                         z_mode float, \
                         z_min float, \
                         z_max float, \
                         FOREIGN KEY (fit_id) REFERENCES fits_data (id), \
                         FOREIGN KEY (det_id) REFERENCES detections_data (id));"
                    ) 

        # Table for model fluxes
        create_table(conn,
                     "CREATE TABLE IF NOT EXISTS fits_fluxes (\
                         id integer PRIMARY KEY, \
                         det_id integer NOT NULL, \
                         fit_id integer NOT NULL, \
                         sflux_obs_median float NOT NULL, \
                         sflux_obs_mode float NOT NULL, \
                         sflux_obs_min float NOT NULL, \
                         sflux_obs_max float NOT NULL, \
                         sflux_int_median float NOT NULL, \
                         sflux_int_mode float NOT NULL, \
                         sflux_int_min float NOT NULL, \
                         sflux_int_max float NOT NULL, \
                         hflux_obs_median float NOT NULL, \
                         hflux_obs_mode float NOT NULL, \
                         hflux_obs_min float NOT NULL, \
                         hflux_obs_max float NOT NULL, \
                         hflux_int_median float NOT NULL, \
                         hflux_int_mode float NOT NULL, \
                         hflux_int_min float NOT NULL, \
                         hflux_int_max float NOT NULL, \
                         slumin_median float NOT NULL, \
                         slumin_mode float NOT NULL, \
                         slumin_min float NOT NULL, \
                         slumin_max float NOT NULL, \
                         hlumin_median float NOT NULL, \
                         hlumin_mode float NOT NULL, \
                         hlumin_min float NOT NULL, \
                         hlumin_max float NOT NULL, \
                         FOREIGN KEY (fit_id) REFERENCES fits_data (id), \
                         FOREIGN KEY (det_id) REFERENCES detections_data (id));"
                    ) 


        ### Populate tables
        # Model names
        try :
            c = conn.cursor()
            c.executemany("INSERT INTO models(model) VALUES(?)", iter(model_names))

        except :
            print "Models already in the db!"



        # Open fits table with detection data
        xmmdet_cat = Table.read(detfile, format="fits")
                
        xmm_iauname = xmmdet_cat.columns["IAUNAME"]
        xmm_srcra = xmmdet_cat.columns["SC_RA"]
        xmm_srcdec = xmmdet_cat.columns["SC_DEC"]
        xmm_srcid = xmmdet_cat.columns["SRCID"]
        xmm_detid = xmmdet_cat.columns["DETID"]
        xmm_obsid = xmmdet_cat.columns["OBS_ID"]
        xmm_srcnum = xmmdet_cat.columns["SRC_NUM"]
        nHgal = xmmdet_cat.columns["NHGAL"]
        zph = xmmdet_cat.columns["PHOT_Z"]
        zphErr = xmmdet_cat.columns["PHOT_ZERR"]
        zsp = xmmdet_cat.columns["SPEC_Z"]

        nsrcs = 10 #len(xmm_iauname)

        for i in xrange(nsrcs) :
            data_folder = "bxa_results/{}/".format(xmm_detid[i])

            # Open counts file
            countsfile = "{}counts.json".format(data_folder)
            with open(countsfile) as data_file :
                counts = json.load(data_file)
            
            # Add detection                
            detection = (xmm_detid[i], xmm_obsid[i], int(xmm_srcnum[i]), xmm_srcid[i],
                         xmm_iauname[i], xmm_srcra[i], xmm_srcdec[i], nHgal[i],
                         zsp[i], float(zph[i]), float(zphErr[i]),
                         counts["total"]["counts"], 
                         counts["band1"]["counts"], 
                         counts["band2"]["counts"])

            try :
                det_id = add_detection(conn, detection)
            except sqlite3.IntegrityError as e:
                print e
                c = conn.cursor()
                c.execute("select id from detections_data where detid=(?)", (xmm_detid[i],))
                det_id = c.fetchall()[0][0]
                print det_id
        

            # Open spec file
            specsfile = "{}spec_info.json".format(data_folder)
            with open(specsfile) as data_file :
                spec_info = json.load(data_file)

            # Add available spectra for the detection                                         
            for spec in spec_info :
                spectrum = (det_id, spec_info[spec]["exposure"], spec_info[spec]["detector"], spec_info[spec]["counts"], spec_info[spec]["good"])
                add_spectrum(conn, spectrum)

            # Add best-fit data for each fitted model
            c = conn.cursor()
            for model in model_names :
                c.execute("select id from models where model=(?)", model)
                model_id = c.fetchall()[0][0]
                
                if not np.isnan(zsp[i]) :
                    label_zpdf = ""
                    thawedz = False
                else :
                    label_zpdf = "withzpdf_"
                    thawedz = True
                                        
                try :
                    # Add goodness info
                    fitfile = "{}chains_{}_{}bestfit_goodness.json".format(data_folder, model[0], label_zpdf)
                    with open(fitfile) as data_file :
                        fit_info = json.load(data_file)
                    
                    fit = (det_id, model_id, 
                           fit_info["npar"], fit_info["dof"], fit_info["wstat"], 
                           fit_info["ks"], fit_info["pvalue"], 
                           fit_info["lnZ"], fit_info["lnZ_err"], 
                           fit_info["lnev"], fit_info["lnev_err"])

                    fit_id = add_fit(conn, fit)
                    
                    # Add best-fit parameters
                    parsfile = "{}chains_{}_{}bestfit_pars.json".format(data_folder, model[0], label_zpdf)
                    with open(parsfile) as data_file :
                        pars_info = json.load(data_file)
                        
                    pars = get_pars(pars_info, model[0], thawedz=thawedz)
                    pars = (det_id, fit_id) + pars
                    add_pars(conn, pars, model[0])
                    
                    # Add fluxes and luminosities
                    fluxfile = "{}chains_{}_{}bestfit_obsflux.json".format(data_folder, model[0], label_zpdf)
                    with open(fluxfile) as data_file :
                        obsflux_info = json.load(data_file)

                    fluxfile = "{}chains_{}_{}bestfit_intflux.json".format(data_folder, model[0], label_zpdf)
                    with open(fluxfile) as data_file :
                        intflux_info = json.load(data_file)

                    fluxfile = "{}chains_{}_{}bestfit_lumin.json".format(data_folder, model[0], label_zpdf)
                    with open(fluxfile) as data_file :
                        lumin_info = json.load(data_file)
                        
                    flux = (det_id, fit_id,
                            obsflux_info["flux0"]["median"], obsflux_info["flux0"]["mode"],
                            obsflux_info["flux0"]["fluxmin"], obsflux_info["flux0"]["fluxmax"],
                            intflux_info["flux0"]["median"], intflux_info["flux0"]["mode"],
                            intflux_info["flux0"]["fluxmin"], intflux_info["flux0"]["fluxmax"],
                            obsflux_info["flux1"]["median"], obsflux_info["flux1"]["mode"],
                            obsflux_info["flux1"]["fluxmin"], obsflux_info["flux1"]["fluxmax"],
                            intflux_info["flux1"]["median"], intflux_info["flux1"]["mode"],
                            intflux_info["flux1"]["fluxmin"], intflux_info["flux1"]["fluxmax"],
                            lumin_info["flux0"]["median"], lumin_info["flux0"]["mode"],
                            lumin_info["flux0"]["fluxmin"], lumin_info["flux0"]["fluxmax"],
                            lumin_info["flux1"]["median"], lumin_info["flux1"]["mode"],
                            lumin_info["flux1"]["fluxmin"], lumin_info["flux1"]["fluxmax"])

                    add_flux(conn, flux)
                             
                except :
                    #print e
                    print "Model {} not fitted for this detection.".format(model[0])

            
        conn.close()
        
    else:
        print("Error! cannot create the database connection.")

    
 
if __name__ == '__main__':
    
    detfile = "./data/3xmmdr6_detwithphotoz_nhgal_1.fits"
    dbfile = "./db/xmmfitcatz.db"
    model_names = [["wapo"], ["wamekal"], ["wamekalpo"], ["wapopo"]]
    add_models = True
    
    main(dbfile, detfile, model_names, add_models)
    
    
