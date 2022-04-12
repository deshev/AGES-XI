#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 08:46:59 2022
@author: Boris Deshev

Plot the HI deficiency and HI non-detected fraction
against local density for bins of global environment
... and other galaxy parameters
"""

from astropy.table import Table, vstack
import matplotlib.pyplot as plt
import numpy as np
from calculate_bins import calc_bins
from astropy.cosmology import FlatLambdaCDM
from astropy.constants import c

c = c.value/1000
H0 = 71
cosmo = FlatLambdaCDM(H0=H0, Om0=0.3, Tcmb0=2.725)

######################################
# # Needed to write correctly HI
import matplotlib
matplotlib.rcParams['pgf.texsystem'] = 'pdflatex'
matplotlib.rcParams.update({'font.family': 'serif', 'font.size': 16,
    'axes.labelsize': 20,'axes.titlesize': 24, 'figure.titlesize' : 28})
matplotlib.rcParams['text.usetex'] = True
#'{\sc Hi}'
#r'log M$_{\textrm{\sc Hi}}$'
# ######################################

       
###############################################################################
###############################################################################
path = '/home/tazio/works/2020/a1367/final_analysis/'
Stbl = Table.read(path+'mpa_jhu_AGES.fits')

# The table with all SDSS objects within the HI survey and cluster velocity range
print(len(Stbl))
Stbl['density'] = np.log10(Stbl['density'])

# # Select only the mass limited sample
Stbl = Stbl[Stbl['LGM_TOT_P50']>=9.75]
print('After M* selection', len(Stbl))
Stbl = Stbl[Stbl['Vel'] <= 15500]
print('After z selection', len(Stbl))

# The cluster members
CStbl = Stbl[(Stbl['Vel'] >= 4000) & (Stbl['Vel'] <= 9000)]

# All the rest
EStbl = Stbl[(Stbl['Vel'] > 9000) | (Stbl['Vel'] < 4000)]

# Select only those above the MHI completeness limit
CHtbl = CStbl[CStbl['MHI'] > 9.25]
EHtbl = EStbl[EStbl['MHI'] > 9.25]

print('In the cluster we have {} and {} HI detected'.format(len(CStbl), len(CHtbl)))
print('Elsewhere we have {} and {} HI detected\n'.format(len(EStbl), len(EHtbl)))
###############################################################################
# The number of bins
Nb = 4
xVars = ['density', 'LGM_TOT_P50', 'SFR_TOT_P50', 'SPECSFR_TOT_P50']
xlabels = iter([r'log($\rho$) [Mpc$^{-3}$]', \
                'log(M$_{\star}$) [M$_{\odot}$]', \
                'log(SFR) [M$_{\odot}$ yr$^{-1}$]', \
                'log(SFR/M$_{\star}$) [yr$^{-1}$]'])
n = 1 # a counter
xlims = []

fig = plt.figure(figsize=(15,8))

for xVar in xVars:
    # Plot the upper row (HIdef)
    ax = plt.subplot(2, len(xVars), n)
    colors = iter(['brown', 'C0'])
    markers = iter(['s', 'o'])
    ls = iter(['-', '--'])
    labels = iter(['A1367', 'Elsewhere'])
    
    for tbl in [CStbl, EStbl]:
        # Bin the x-axis
        vvbins = np.array(calc_bins([tbl, xVar, Nb]))
        tbl['vvbins'] = np.digitize(tbl[xVar], bins = vvbins)
        c = next(colors)
        m = next(markers)
        l = next(ls)
        label = next(labels)
        Xb = []; YbCD = []; dYbCD = []
        for i in range(1, len(vvbins)):
            Xb.append(np.median(tbl[xVar][tbl['vvbins']==i]))
            YbCD.append(np.mean(tbl['HIdef_ALL'][tbl['vvbins']==i]))
            dYbCD.append(np.std(tbl['HIdef_ALL'][tbl['vvbins']==i]))
        YbCD = np.array(YbCD)
        dYbCD = np.array(dYbCD)
        ax.plot(Xb, YbCD, m, linestyle=l, color=c, label=label, zorder=1)
        ax.fill_between(Xb, YbCD+dYbCD/2, YbCD-dYbCD/2, color=c, alpha=0.3, zorder=1)
        ax.set_ylim([-0.44, 1.1])
        ax.minorticks_on()
        ax.tick_params(axis='both', direction='in', which='both', right='on', top='on')
        if xVar == 'density':
            ax.legend(loc=2)
            ax.set_ylabel('{\sc Hi} deficiency')
            ax.text(Xb[0], -0.35, len(tbl), color=c)
            print(Xb[0])

        else:
            ax.set_yticklabels([])
    
    ax.set_xticklabels([])
    xlims = ax.get_xlim()
    
    # Plot the lower row (non-detected fraction)
    ax = plt.subplot(2, len(xVars), n+len(xVars))
    colors = iter(['brown', 'C0'])
    markers = iter(['s', 'o'])
    ls = iter(['-', '--'])
    
    for tbl, Htbl in zip([CStbl, EStbl], [CHtbl, EHtbl]):
        # Bin the x-axis
        vvbins = np.array(calc_bins([tbl, xVar, Nb]))
        tbl['vvbins'] = np.digitize(tbl[xVar], bins = vvbins)
        c = next(colors)
        m = next(markers)
        l = next(ls)

        Xb = []
        for i in range(1, len(vvbins)):
            Xb.append(np.median(tbl[xVar][tbl['vvbins']==i]))

        Sh = np.histogram(tbl[xVar], bins=vvbins)
        Hh = np.histogram(Htbl[xVar], bins=vvbins)
        YbEF = (Sh[0]-Hh[0])/Sh[0]

        dYbEF = np.sqrt(Sh[0]-Hh[0])/Sh[0]
        ax.plot((Xb), YbEF, m, linestyle=l, color=c, zorder=3)
        ax.fill_between((Xb), y1=(YbEF)+dYbEF, y2=(YbEF)-dYbEF, color=c, alpha=0.3, zorder=3)

        ax.set_ylim([0.15, 1.0])
        ax.set_xlim(xlims)
        ax.minorticks_on()
        ax.tick_params(axis='both', direction='in', which='both', right='on', top='on')
        if xVar == 'density':
            ax.set_ylabel('{\sc Hi} non-detected fraction')
        else:
            ax.set_yticklabels([])

    ax.set_xlabel(next(xlabels))
    n += 1
plt.tight_layout(pad=2.0, h_pad=(0.0), w_pad=(0.0))
plt.savefig(path+'compare_HIDF_HIdef.pdf')
