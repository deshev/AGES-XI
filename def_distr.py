#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 13:13:25 2021
@author: Boris Deshev

Plot distribution of deficiency for all galaxies and for the HI detected
ones separately
"""

from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import models, fitting
######################################
# # Needed to write correctly HI
import matplotlib
matplotlib.rcParams['pgf.texsystem'] = 'pdflatex'
matplotlib.rcParams.update({'font.family': 'serif', 'font.size': 16,
                            'axes.labelsize': 20, 'axes.titlesize': 24,
                            'figure.titlesize': 28})
matplotlib.rcParams['text.usetex'] = True
# '{\sc Hi}'
#r'log M$_{\textrm{\sc Hi}}$'
######################################

# Read in the table
path = '/home/tazio/works/2020/a1367/final_analysis/'
tbl = Table.read(path+'mpa_jhu_AGES.fits')

# Split the sample into cluster, foreground and background
CzT = tbl[(tbl['Vel'] >= 4000) & (tbl['Vel'] <= 9000)]
EzT = tbl[(tbl['Vel'] > 9000) | (tbl['Vel'] < 4000)]

########################################
# Bins in deficiency
bins = np.linspace(-1.5, 1.5, 15)
# Bin centers
x = (bins[:-1]+bins[1:])/2

fig = plt.figure(figsize=(8, 5.5))

# This is for the common Y-axis label
# Make a common Y-axis label
ax = fig.add_subplot(1,1,1)
# Turn off axis lines and ticks of the big subplot
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='None', top=False, bottom=False, left=False, right=False)
# Set common labels
ax.set_ylabel('fraction of galaxies')

#################################################
# Plot the full data set
ax1 = fig.add_subplot(2,1,1)

# Plot the cluster members ################
weights = np.ones_like(CzT['HIdef_ALL'])/len(CzT)
clh = ax1.hist(CzT['HIdef_ALL'], bins=bins, weights=weights, histtype='step', linewidth=3, color='brown', linestyle='--', zorder=2)

# Plot the rest ################
weights = np.ones_like(EzT['HIdef_ALL'])/len(EzT)
clh = ax1.hist(EzT['HIdef_ALL'], bins=bins, weights=weights, histtype='stepfilled', linewidth=2, color='C0', alpha=0.5, zorder=1)

# Print the total number of galaxies in each histogram
ax1.text(-1.6, 0.35, str(len(CzT)), color='brown')
ax1.text(-1.6, 0.31, str(len(EzT)), color='C0')
ax1.set_ylim([0.0, 0.39])
ax1.set_xticklabels([])

ax1.tick_params(axis='both', direction='in', which='both', right='on', top='on')
ax1.minorticks_on()


#################################################
# Repeat the same but excluding the upper limits
#################################################
ax2 = fig.add_subplot(2,1,2)

# Plot the cluster members ################
weights = np.ones_like(CzT['HIdef'])/len(CzT[CzT['HIdetected']==True])
clh = ax2.hist(CzT['HIdef'], bins=bins, label='A1367', weights=weights, histtype='step', linewidth=3, color='brown', linestyle='--', zorder=2)

# Plot the rest ################
weights = np.ones_like(EzT['HIdef'])/len(EzT[EzT['HIdetected']==True])
clh = ax2.hist(EzT['HIdef'], bins=bins, label='Elsewhere', weights=weights, histtype='stepfilled', color='C0', linewidth=2, alpha=0.5, zorder=1)

# Print the total number of galaxies in each histogram
ax2.text(-1.6, 0.35, str(len(CzT[CzT['HIdetected']==True])), color='brown')
ax2.text(-1.6, 0.31, str(len(EzT[EzT['HIdetected']==True])), color='C0')

ax2.legend(loc=1)
ax2.set_ylim([0.0, 0.39])
ax2.set_xlabel(r'{\sc Hi} deficiency')
ax2.tick_params(axis='both', direction='in', which='both', right='on', top='on')
ax2.minorticks_on()

plt.tight_layout(pad=0.0)
plt.savefig(path+'HIdef_distr_v2.pdf')


