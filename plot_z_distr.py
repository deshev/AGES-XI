#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 14:19:14 2021
@author: Boris Deshev

Make a plot showing the z distribution of all HI detections and SDSS spec galaxies
"""

from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from astropy.constants import c
c = c.value/1000
from matplotlib.patches import Patch
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
H0 = 71 # Hubble constant

def give_f(d):
    # Calculate the flux at the completeness of ALFALFA SNR=6.5
    snr = 6.5
    w = 200
    rms = 0.7 # Median of all the sources
    fl = ((snr/(np.sqrt((w/2)/10)))*rms*w)/1000
    return np.log10(2.36e5*(d**2)*fl)

######################################
# # Needed to write correctly HI
import matplotlib
matplotlib.rcParams['pgf.texsystem'] = 'pdflatex'
matplotlib.rcParams.update({'font.family': 'serif', 'font.size': 16,
    'axes.labelsize': 20,'axes.titlesize': 24, 'figure.titlesize' : 28})
matplotlib.rcParams['text.usetex'] = True
#'{\sc Hi}'
#r'log M$_{\textrm{\sc Hi}}$'
######################################

path = '/'
# The table with the HI detections
HItbl = Table.read(path+'A1367_HItable_v8.8_wt.txt', format='ascii.commented_header')
print('tbl',len(HItbl))
max_z = HItbl['Vel_HI'].max()/c
min_z = HItbl['Vel_HI'].min()/c

# The table with the SDSS galaxies
SDSStbl = Table.read(path+'sdss/sdss_objects_z_lt1.csv')
print('tot from SDSS', len(SDSStbl))
# The table with the ALFALFA catalog of the overlapping region on the sky
alf = Table.read(path+'beam_correction/AGES_ALFALFA_match3arcmin.csv')
print('From ALFALFA ', len(alf))

bins = np.linspace(min_z, max_z, 25)
x = (bins[1:]+bins[:-1])/2


plt.figure(figsize=(8,7.5))

ax = plt.subplot(2,1,1)
###################### Adjust all the ticks
ax.tick_params(axis='both', direction='in', which='both', right='on')

oh = ax.hist(SDSStbl['z'], bins=bins, histtype='stepfilled', color='gray', alpha=0.5, label='SDSS')
hh = ax.hist(HItbl['Vel_HI']/c, bins=bins, histtype='stepfilled', color='cyan', edgecolor='C0', hatch='//', linewidth=0.1, label='AGES')
ah = ax.hist(alf['Vhelio']/c, bins=bins, histtype='step', color='brown', linewidth=1, label='ALFALFA')#alpha=0.7, 

ax.minorticks_on()
ax.set_xlabel('$z$')
ax.set_ylabel('number of galaxies')
ax.set_xlim([0.003,0.065])

########### Do the legend
legend_elements = [Patch(facecolor='cyan', hatch='//', edgecolor='C0', linewidth=0.1, label='AGES'), #edgecolor='C0', 
                   Patch(facecolor='gray', label='SDSS', alpha=0.5), #edgecolor='gray', 
                   Patch(facecolor='None', edgecolor='brown', linewidth=1, label='ALFALFA', alpha=0.7)]

ax.legend(handles=legend_elements, loc=1)
###########

ax3 = ax.twiny()
ax3.set_xlim([0.003,0.065])
ax3.set_xlabel('$cz$ [km s$^{-1}$]')
labels = ([3000,6000,9000,12000,15000,18000])
locs = labels/c
ax3.set_xticks(locs)
ax3.set_xticklabels(labels)
ax3.tick_params(axis='both', direction='in', which='both')
ax3.minorticks_on()

######################
# Mark the cluster
plt.text(5900/c, 204, 'A1367')
plt.errorbar(6500/c, 220, xerr=2500/c, color='black', linewidth=1, capsize=3)
   
###################################################################
ax1 = plt.subplot(2,1,2)
# Select uncertain detections
mask1 = np.char.find(HItbl['flag'], '1')
HItbl_s = HItbl[mask1 == -1]
HItbl_u = HItbl[mask1 > -1]

ds = HItbl_s['Vel_HI']/H0
du = HItbl_u['Vel_HI']/H0
ax1.errorbar(ds, HItbl_s['MHI'], xerr=HItbl_s['eV']/H0, yerr=HItbl_s['eMHI'], linestyle='None', marker='*', mfc='cyan', mec='C0', ecolor='C0', markeredgewidth=0.5, ms=8, linewidth=1, label='certain')
ax1.errorbar(du, HItbl_u['MHI'], xerr=HItbl_u['eV']/H0, yerr=HItbl_u['eMHI'], linestyle='None', marker='o', mfc='white', mec='brown', ecolor='brown', markeredgewidth=0.5, ms=5, linewidth=1, label='uncertain')
ax1.set_ylabel(r'log M$_{\textrm{\sc Hi}}$ [M$_{\odot}$]')#'log M$_{\mathrm{\sc HI}}$/M$_{\odot}$')
ax1.yaxis.set_minor_locator(MultipleLocator(5))
ax1.set_xlim([5,280])
ax1.set_xticks([50,100,150,200,250])
ax1.set_xticklabels([50,100,150,200,250])
ax1.set_xlabel('D [Mpc]')
plt.legend(loc=4)

# Plot the sensitivity limit at SNR=6.5, for a profile with W=200kms
d = np.concatenate((ds,du))
dx = np.linspace(d.min(), d.max(),100)
plt.plot(dx, give_f(dx), linestyle='--', linewidth=1, color='black', zorder=5)

ax1b = ax1.twiny()
ax1b.set_xticks([50,100,150,200,250])
ax1b.set_xticklabels([])
ax1.minorticks_on()
ax1b.minorticks_on()
ax1b.set_xlim([5,280])

ax1b.tick_params(axis='both', direction='in', which='both', top='on', right='on')
ax1.tick_params(axis='both', direction='in', which='both', top='on', right='on')
###################################################################

plt.tight_layout()
plt.savefig(path+'final_analysis/z_distr6.pdf')













































