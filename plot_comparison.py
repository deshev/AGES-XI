#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 16:07:54 2021
@author: Boris Deshev

Make a plot to compare the fluxes and W50 adn W20
"""
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import sys
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

path = '/home/tazio/works/2020/a1367/'
###########################################################################
# Comparison with ALFALFA
tbl = Table.read(path+'beam_correction/AGES_ALFALFA_match3arcmin.csv')
print(len(tbl))

plt.figure(figsize=(7,4))
ax = plt.subplot(1,1,1)
# The 1 to 1 line
plt.axhline(1.0, linestyle='--', linewidth=1, color='black')
# The 6.5 sigma reliability level
plt.axvline(6.5, linestyle='--', linewidth=1, color='black')
# Plot the points
tbl['flux_ratio'] = tbl['TotFlux']/tbl['HIflux']
plt.errorbar(tbl['SNR'], tbl['flux_ratio'], yerr=np.sqrt(tbl['eTotF']**2+tbl['sigflux']**2), ms=8, marker='*', mec='C0', mfc='cyan', markeredgewidth=0.5,linewidth=0.5, linestyle='None')
plt.text(6.0, -0.26, '6.5')

# Calculate the running mean
tbl.sort(keys='SNR')
n = int(np.sqrt(len(tbl)))
rm = []; msnr = []; rstd = []
for i in range(0,len(tbl)-n,1):
    rm.append(np.median(tbl['flux_ratio'][i:i+n]))
    msnr.append(np.mean(tbl['SNR'][i:i+n]))
    rstd.append(np.std(tbl['flux_ratio'][i:i+n]))

rm=np.array(rm); msnr=np.array(msnr); rstd=np.array(rstd)
plt.plot(msnr,rm, color='0.1', alpha=1, linewidth=2, zorder=4)
plt.fill_between(msnr, rm+rstd, rm-rstd, color='0.3', linewidth=0.0, alpha=0.5, zorder=3)


###################### Adjust all the ticks
ax.tick_params(axis='both', direction='in', which='both', right='on', top='on')
plt.xscale('log')
plt.xlabel('SNR$^{ALFALFA}$')
plt.ylabel('F$^{AGES}_{tot}$ / F$^{ALFALFA}_{tot}$')
ax.minorticks_on()
plt.ylim([0.0, 2.3])

plt.tight_layout()
plt.savefig(path+'final_analysis/flux_comparison2.pdf')

