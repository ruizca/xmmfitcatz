#!/usr/bin/env python
"""
Created on Thu Nov 23 12:45:43 2017

@author: ruizca
"""

import numpy as np
import os

from astropy.table import Table
import astropy.units as u

def nh(ra, dec, equinox=2000.0) :
    command = "nh " + str(equinox) + " " + str(ra) + " " + str(dec) + "> tmp"
    os.system(command)

    nhtemp = open("tmp", "r")
    for line in nhtemp :
        if line.startswith("  LAB >> Weighted") :
            nhval = line.split()[-1]
    os.remove("tmp")
        
    return float(nhval)
    
  
sources_table = "./data/3xmmdr6_detwithphotoz.fits"
output_table = "./data/3xmmdr6_detwithphotoz_nhgal.fits"

xmmdet_cat = Table.read(sources_table, format="fits")

ra = xmmdet_cat["SC_RA"]
dec = xmmdet_cat["SC_DEC"]
nhvals = np.full((len(ra),), np.nan)


for i in range(len(ra)) :
    nhvals[i] = nh(ra[i], dec[i])

nhvals = nhvals * u.cm**-2
nhcol = Table.Column(nhvals, name="NHGAL")    

xmmdet_cat.add_column(nhcol)
xmmdet_cat.write(output_table, format="fits", overwrite=True)