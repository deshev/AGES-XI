#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 16:07:54 2021
@author: Boris Deshev

Make a plot to compare the AGES fluxes with those from ALFALFA
"""
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt

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
###########################################################################
# Comparison with ALFALFA
tbl = Table.read(path+'final_analysis/AGES_ALFALFA_match1.5arcmin.csv')
print(len(tbl))

plt.figure(figsize=(7,4))
ax = plt.subplot(1,1,1)
# The 1 to 1 line
plt.axhline(1.0, linestyle='--', linewidth=1, color='black')
# The 6.5 sigma reliability level
plt.axvline(6.5, linestyle='--', linewidth=1, color='black')

# The flux ratio
tbl['flux_ratio'] = tbl['TotFlux']/tbl['HIflux']

# Split the table according to the HI flag from ALFALFA
tbl1 = tbl[tbl['HIcode'] == 1]
tbl2 = tbl[tbl['HIcode'] == 2]

plt.errorbar(tbl1['SNR'], tbl1['flux_ratio'], yerr=np.sqrt(tbl1['eTotF']**2+tbl1['sigflux']**2), ms=8, marker='*', mec='C0', mfc='cyan', markeredgewidth=0.5,linewidth=0.5, linestyle='None', label='1')
plt.errorbar(tbl2['SNR'], tbl2['flux_ratio'], yerr=np.sqrt(tbl2['eTotF']**2+tbl2['sigflux']**2), ms=8, marker='*', mec='C0', mfc='white', ecolor='C0', markeredgewidth=0.5,linewidth=0.5, linestyle='None', label='2')
plt.text(6.0, -0.26, '6.5')
plt.legend(title='ALFALFA {\sc Hi} flag')

# Calculate the running mean
tbl.sort(keys='SNR')
n = int(np.sqrt(len(tbl)))
rm = []; msnr = []; rstd = []
for i in range(0,len(tbl)-n,1):
    rm.append(np.median(tbl['flux_ratio'][i:i+n]))
    msnr.append(np.median(tbl['SNR'][i:i+n]))
    rstd.append(np.std(tbl['flux_ratio'][i:i+n]))

rm=np.array(rm); msnr=np.array(msnr); rstd=np.array(rstd)
plt.plot(msnr,rm, color='0.1', alpha=1, linewidth=2, zorder=4)
plt.fill_between(msnr, rm+rstd, rm-rstd, color='0.3', linewidth=0.0, alpha=0.5, zorder=0)


###################### Adjust all the ticks
ax.tick_params(axis='both', direction='in', which='both', right='on', top='on')
plt.xscale('log')
plt.xlabel('SNR$^{ALFALFA}$')
plt.ylabel('F$^{AGES}_{tot}$ / F$^{ALFALFA}_{tot}$')
ax.minorticks_on()
plt.ylim([0.0, 2.3])

plt.tight_layout()
plt.savefig(path+'final_analysis/flux_comparison3.pdf')

