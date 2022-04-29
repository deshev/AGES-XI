#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 18:54:10 2021
@author: Boris Deshev

Plot stellar and HI mass distributions marking the mass complete regions
Also the local density
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.constants import c
c = c.value/1000

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

def running_median(tbl, medcol, sortcol):
    # Calculate the running median
    tbl.sort(keys=sortcol)
    
    n = int(np.sqrt(len(tbl)))
    rm = []
    msnr = []
    rstd = []
    for i in range(0, len(tbl)-n, 1):
        rm.append(np.nanmedian(tbl[medcol][i:i+n]))
        msnr.append(np.nanmedian(tbl[sortcol][i:i+n]))
        rstd.append(np.nanstd(tbl[medcol][i:i+n]))
    Y = np.array(rm)
    X = np.array(msnr)
    dY = np.array(rstd)
    return X, Y, dY
######################################

path = '/home/tazio/works/2020/a1367/'
# Completeness limits
MstLim = 9.75
MhiLim = 9.25

# Table with all HI detections
Htbl = Table.read(path+'A1367_HItable_v8.8_wt.txt', format='ascii.commented_header')
# Table with all optical galaxies within the AGES volume
Stbl = Table.read(path+'final_analysis/mpa_jhu_AGES.fits')

plt.figure(figsize=(7,9))
xlim = [1,19999]

#####################################
ax1 = plt.subplot(3,1,1)
ax1.plot(Htbl['Vel_HI'], Htbl['MHI'], 'o', ms=6, mfc='cyan', mec='C0', markeredgewidth=0.5, zorder=3, label='AGES')
ax1.set_ylabel(r'log(M$_{\textrm{\sc Hi}}$) [M$_{\odot}$]')
ax1.set_ylim([7.2,11.1])
ax1.set_xlim(xlim)

x=np.array([0, 15500, 15500, 25000])
y1=np.array([0, 0, 0, 0])
y2=np.array([MhiLim, MhiLim, 12, 12])
ax1.fill_between(x, y1, y2, color='0.3', alpha=0.5,zorder=1)

ax1.tick_params(axis='both', direction='in', which='both', right='off', top='on')
ax1.minorticks_on()
ax1.set_xticklabels([])
ax1.legend(loc=4 ,framealpha=1)

ax1b = ax1.twiny()
ax1b.set_xlim(xlim/c)
ax1b.set_xlabel('$z$')
ax1b.tick_params(axis='both', direction='in', which='both')
ax1b.minorticks_on()

#####################################
ax2 = plt.subplot(3,1,2)
ax2.plot(Stbl['Vel'], Stbl['LGM_TOT_P50'], '*', ms=9, mfc='cyan', mec='C0', markeredgewidth=0.5, zorder=3, label='AGES $\cap$ MPA-JHU')
ax2.set_ylim([7.1,11.6])
ax2.set_xlim(xlim)
ax2.set_ylabel('log(M$_{\star}$) [M$_{\odot}$]')
ax2.set_xticklabels([])

x=np.array([0, 15500, 15500, 25000])
y1=np.array([0, 0, 0, 0])
y2=np.array([MstLim, MstLim, 12, 12])
ax2.fill_between(x, y1, y2, color='0.3', alpha=0.5,zorder=1)

ax2.tick_params(axis='x', direction='in', which='both', top=True)
ax2.tick_params(axis='y', direction='in', which='both', right=True)
ax2.minorticks_on()
ax2.legend(loc=4, framealpha=1)

#####################################
# The column with the density only contains galaxies from the mass complete sample
Stbl = Stbl[~Stbl['density'].mask]

XA, YA, dYA = running_median(Stbl, 'density', 'Z')

ax3 = plt.subplot(3,1,3)
ax3.plot(Stbl['Vel'], np.log10(Stbl['density']), '*', ms=9, mfc='cyan', mec='C0', markeredgewidth=0.5, label=r'AGES $\cap$ MPA-JHU \& log(M$_{\star}) \geq$ 9.75', zorder=2)

ax3.plot(XA*c, np.log10(YA), linewidth=1, color='0.1', zorder=3)

x = [15500, 25000]
ax3.fill_between(x, y1=-5, y2=2, color='0.3', alpha=0.5,zorder=1)



ax3.set_xlim(xlim)
ax3.set_ylim([-4,3.9])
ax3.tick_params(axis='both', direction='in', which='both', right='off', top='on')
ax3.legend(loc=1, framealpha=1)
[ax3.axhline(y, color='0.7', linewidth=0.5, linestyle='-', zorder=0) for y in np.arange(-3, 3)]
ax3.set_xlabel('$cz$ [km s$^{-1}$]')
ax3.set_ylabel(r'log($\rho$) [Mpc$^{-3}$]')
ax3.minorticks_on()

plt.tight_layout()

#plt.show()
plt.savefig(path+'final_analysis/mass_distributions_comoving2.pdf', dpi=200)


