#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 10:31:16 2021
@author: Boris Deshev

Estimate galaxy density via Voronoi binning
in N/Mpc^3

9.Mar.2022 Added correction for Fingers-of-God
Only applicable for A1367!!!
"""

from astropy.table import Table
import numpy as np
from scipy.spatial import Voronoi,Delaunay
from astropy.constants import c as c_lum
c_lum = c_lum.value/1000 # in km/s
from astropy.cosmology import FlatLambdaCDM
H0 = 71
cosmo = FlatLambdaCDM(H0=H0, Om0=0.3, Tcmb0=2.725)

#--------------------------------------
# Some of the code between these lines is from stackoverflow. User is Chris H.
def tetravol(a,b,c,d):
    '''Calculates the volume of a tetrahedron, given vertices a,b,c and d (triplets)'''
    tetravol=abs(np.dot((a-d),np.cross((b-d),(c-d))))/6
    return tetravol

def vol(vor,p):
    # Write out the Voronoi bins for other fun visualisations
    of = open(path+'mpa_jhu_sky_overlap_vertices.csv','w')
    # Calculate the volume of the tesselations
    dpoints=[]
    vol=0
    for v in vor.regions[vor.point_region[p]]:
        dpoints.append(list(vor.vertices[v]))
        of.write(','.join(str(jj) for jj in vor.vertices[v])+'\n')
    of.close()
    
    tri=Delaunay(np.array(dpoints))
    for simplex in tri.simplices:
        vol += tetravol(np.array(dpoints[simplex[0]]),np.array(dpoints[simplex[1]]),np.array(dpoints[simplex[2]]),np.array(dpoints[simplex[3]]))
    return vol
#--------------------------------------

################################################
path = '/'

# The table
t5 = Table.read(path+'mpa_jhu_sky_overlap.fits')
print(len(t5))
t5 = t5[t5['LGM_TOT_P50']>=9.75]
print(len(t5))
ra1 = 165
ra2 = 185
dec1 = 10
dec2 = 30
t5 = t5[(t5['RA']>ra1) &
        (t5['RA']<ra2) &
        (t5['DEC']>dec1) &
        (t5['DEC']<dec2)]
print(len(t5))

t5['vorvol'] = 0.0
t5['density'] = 0.0

# Write a column with the distance from z and then correct it for
# FoG effect. The calculation of the needed correction is done with
# do_FoG_correction.py and is 3.7
scale = 3.7
# Calculate the distance in Mpc
t5['Vel'] = t5['Z']*c_lum
t5['dist'] = cosmo.comoving_distance(t5['Z']).value
t5['dist_corr'] = t5['dist']
dc = 92.8 # Distance to the cluster center in Mpc

# Apply the correction
for i in range(len(t5)):
    # This selects the cluster members
    if (t5['Vel'][i]>=4000) & (t5['Vel'][i]<=9000):
        t5['dist_corr'][i] = ((t5['dist'][i]-dc)/scale)+dc


# Normalise the coordinates
# and convert them to physical units
xd = t5['RA'].data.copy() - (t5['RA'].min()+t5['RA'].max())/2
yd = t5['DEC'].data.copy() - (t5['DEC'].min()+t5['DEC'].max())/2
zd = t5['Z'].data.copy()

# Convert dRA to Mpc
x = (xd * 60 * cosmo.kpc_comoving_per_arcmin(zd) / 1000).value
y = (yd * 60 * cosmo.kpc_comoving_per_arcmin(zd) / 1000).value
z = t5['dist_corr']

points = np.array([x,y,z]).T
dpoints = []
vor = Voronoi(points)
vtot = 0

for i,p in enumerate(vor.points):
    out = False
    for v in vor.regions[vor.point_region[i]]:
        if v<=-1: # An index of -1 is returned if the vertex is outside the Vornoi diagram. Ignore them
            out = True
    if not out:
        pvol = vol(vor,i)
        vtot += pvol
        t5['vorvol'][i] = pvol

# Write another column with the density
t5['density'] = 1 / (t5['vorvol'])








































