#!/usr/bin/env python
"""
Created on Fri Mar  9 10:36:56 2018

@author: ruizca
"""
nfits = 7

nsrc = 0
for i in range(nfits) :
    with open('last_source_{:d}.dat'.format(i+1), 'r') as fp:
        nsrc += int(fp.readline())
        
print nsrc
        
