#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  6 16:20:38 2022
@author: Boris Deshev

Plot MHI/M* relation for AGES galaxies
"""

from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from calculate_bins import calc_bins
import matplotlib.cm as cm
import matplotlib

######################################
# # Needed to write correctly HI
import matplotlib
matplotlib.rcParams['pgf.texsystem'] = 'pdflatex'
matplotlib.rcParams.update({'font.family': 'serif', 'font.size': 24,
    'axes.labelsize': 24,'axes.titlesize': 24, 'figure.titlesize' : 28})
matplotlib.rcParams['text.usetex'] = True
#'{\sc Hi}'
#r'log M$_{\textrm{\sc Hi}}$'
# ######################################

def plot_mhimstar(tbl):
    fig, (ax,ax2) = plt.subplots(2,1,sharex=True,sharey=True, figsize=(10,9.5), gridspec_kw={'height_ratios': [1, 1]})
    
    # Plot the detected galaxies
    mappable = ax.scatter(tbl['LGM_TOT_P50'], tbl['mhimstar'], c = tbl['HIdef_ALL'], s=70, marker='p', edgecolors='black', linewidth=0.01, cmap='PRGn', vmin=-1.0, vmax=1.0, zorder=2)

    # Take care of the colorbar
    cbaxes = fig.add_axes([0.48, 0.9, 0.45, 0.02]) 
    cbar = plt.colorbar(mappable, cax = cbaxes, orientation='horizontal')  
    cbar.set_label('{\sc Hi} deficiency', fontsize=28, rotation=0)#, labelpad=25)
    cbar.ax.xaxis.set_label_position('top')

    # Plot the non-detected ones
    # This is for the arrows
    uplims = np.ones(len(tbl), dtype=bool)

    # # Color code the deficiency for the non-detections
    # cmap = cm.coolwarm
    # norm = matplotlib.colors.Normalize(vmin=tbl['HIdef_ALL'].min(), vmax=tbl['HIdef_ALL'].max(), clip=True)
    # for i in range(len(tbl)):
    #     plt.errorbar(tbl[i]['LGM_TOT_P50'], tbl[i]['mhimstar_UL'], xerr=None, yerr=0.1, uplims=uplims[i], marker='.',\
    #                  color=cmap(norm(tbl[i]['HIdef_ALL'].data)), linestyle='None', zorder=1)
    # Or plot them in gray
    ax.errorbar(tbl['LGM_TOT_P50'], tbl['mhimstar_UL'], xerr=None, yerr=np.ones(len(tbl))*0.1, uplims=uplims, marker='.', color='0.4', linestyle='None', zorder=1)

    # Shade the region below M_*=9.75 to show the stellar mass complete region of SDSS
    ax.axvspan(6,9.75, color='gray', alpha=0.5, zorder=0)

    ax.set_ylabel(r'log(M$_{\textrm{\sc Hi}}$/M$_{\star}$)')
    ax.tick_params(axis='both', direction='in', which='both', right='on', top='on')
    ax.minorticks_on()
    ax.set_ylim([-2.53,2.5])
    ax.set_xlim([7.4,11.5])
    ax.set_yticks([-2,-1,0,1,2])
    

####################################################
    # Separate the cluster from the rest
    cl_mask = ((tbl['Vel']>=4000) & (tbl['Vel']<=9000))
    Ctbl = tbl[cl_mask]
    Etbl = tbl[~cl_mask]
    
    # This will just plot two lines one for the cluster and one for the Elsewhere
    clr = ['C3', 'C0']
    label = ['AGES A1367', 'AGES Elsewhere']
    hatches = ['////', '\\\\']
    clrI = iter(clr)
    lblI = iter(label)
    hatchI = iter(hatches)
    for TT in [Ctbl, Etbl]:
        MstB = np.array(calc_bins([TT, 'LGM_TOT_P50', 15]))

        TT['Mst_bins'] = np.digitize(TT['LGM_TOT_P50'], bins=MstB)
        ML = []; SL = []; xs = []
        for i in range(1,len(MstB)):
            ML.append( np.ma.median(TT['mhimstar_ALL'][TT['Mst_bins']==i]) )
            SL.append( np.ma.std(TT['mhimstar_ALL'][TT['Mst_bins']==i]) )
            xs.append( np.ma.median(TT['LGM_TOT_P50'][TT['Mst_bins']==i]) )
        color = next(clrI)
        ML = np.array(ML)
        SL = np.array(SL)
        ax2.plot(xs, ML, color=color, linewidth=2, zorder=1, label=next(lblI))
        ax2.fill_between(xs, ML-SL, y2=ML+SL, color=color, alpha=0.3, zorder=3)

    # Median of our data
    MstB = np.array(calc_bins([tbl, 'LGM_TOT_P50', 10]))
    tbl['Mst_bins'] = np.digitize(tbl['LGM_TOT_P50'], bins=MstB)
    med_mhimstar = []; std_mhimstar = []; x = []
    for i in range(1,len(MstB)):
        med_mhimstar.append( np.ma.median(tbl['mhimstar_ALL'][tbl['Mst_bins']==i]) )
        std_mhimstar.append( np.ma.std(tbl['mhimstar_ALL'][tbl['Mst_bins']==i]) )
        x.append( np.median(tbl['LGM_TOT_P50'][tbl['Mst_bins']==i]) )
        print(i, len(tbl[tbl['Mst_bins']==i]), tbl['mhimstar_ALL'][tbl['Mst_bins']==i].min(), med_mhimstar[-1], std_mhimstar[-1])
    print(med_mhimstar)
    print(std_mhimstar)
    ax2.errorbar(x, med_mhimstar, yerr=std_mhimstar, color='black', label='AGES', linewidth=3, zorder=10)
    
    # Plot the relation by Catinella+2018    
    ax2.plot(Mst, FHI, color='C1', linestyle='-', linewidth=4, label='xGASS', zorder=9) # (Catinella et al. 2018)

    ax2.legend(ncol=2)
    
    # Shade the region below M_*=9.75 to show the stellar mass complete region of SDSS
    ax2.axvspan(6,9.75, color='gray', alpha=0.5, zorder=0)

    ax2.tick_params(axis='both', direction='in', which='both', right='on', top='on')
    ax2.minorticks_on()
    ax2.set_xlabel('log(M$_{\star}$) [M$_{\odot}$]')
    ax2.set_ylabel(r'log(M$_{\textrm{\sc Hi}}$/M$_{\star}$)')

    plt.tight_layout()
    plt.savefig(path+'mhimstar_4.pdf', dpi=200)

######################################################
######################################################
######################################################

path = '/home/tazio/works/2020/a1367/final_analysis/'
tbl = Table.read(path+'mpa_jhu_AGES.fits')
print('tbl',len(tbl))

# Just for quick writing
tbl['mhimstar'] = tbl['MHI']-tbl['LGM_TOT_P50']

# Fix the non-detections
DMmask = tbl['MHI'].mask
UPmask = np.logical_not(DMmask)
tbl['mhimstar_UL'] = tbl['MHI_up_lim'] - tbl['LGM_TOT_P50']
tbl['mhimstar_UL'].mask = UPmask

# Make a column containing both estimates, for calculating the median
tbl['mhimstar_ALL'] = tbl['mhimstar'].filled(1) * tbl['mhimstar_UL'].filled(1)
tbl['mhimstar_ALL'] = np.ma.masked_where(tbl['mhimstar_ALL'] < -3, tbl['mhimstar_ALL'])


# This is from Catinnela et al. 2018 (XGASS)
Mst = np.array([9.14 ,9.44 ,9.74 ,10.07,10.34,10.65,10.95,11.20])
FHI = [-0.092 ,-0.320 ,-0.656 ,-0.854,-1.278,-1.223,-1.707,-1.785]


axs = plot_mhimstar(tbl)





    # # Set up the legend
    # order = [0,1]
    # handles, labels = ax.gca().get_legend_handles_labels()
    # print(handles, labels)
    # plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], ncol=2, prop={'size': 15.5})
