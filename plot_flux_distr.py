#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 10:17:41 2021
@author: Boris Deshev

Make a plot like Fig.1 of Cortese et al. 2008
showing the flux and profile width distributions
"""

from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

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
tbl = Table.read(path+'A1367_HItable_v8.8_wt.txt', format='ascii.commented_header')

def give_f(w):
    # Calculate the flux at the completeness of ALFALFA SNR=6.5
    rms = 0.7 # Median of all the sources
    return ((6.5/(np.sqrt((w/2)/10)))*rms*w)/1000
def give_f400(w):
    rms = 0.7 # Median of all the sources
    return ((6.5/(np.sqrt((400/2)/10)))*rms*w)/1000

plt.figure(figsize=(8,5))
ax = plt.subplot(1,1,1)
###################### Adjust all the ticks
ax.tick_params(axis='both', direction='in', which='both', right='on', \
                 left='on', top='on', bottom='on' )#, \
######################
# Select uncertain detections
mask1 = np.char.find(tbl['flag'], '1')
sd = tbl[mask1 == -1]
rd = tbl[mask1 > -1]
print(len(sd),len(rd))

p1, = plt.plot(sd['W50'][sd['specObjID']>=1], sd['TotFlux'][sd['specObjID']>=1], '*', ms=8, mfc='C0', mec='C0', markeredgewidth=0.5 , alpha=1, zorder=1)
p2, = plt.plot(sd['W50'][sd['specObjID']==-1], sd['TotFlux'][sd['specObjID']==-1], '*', ms=8, mfc='None', mec='C0', markeredgewidth=0.5 , alpha=1, zorder=1)
p2p, = plt.plot(sd['W50'][sd['potOC']!='blah'], sd['TotFlux'][sd['potOC']!='blah'], '*', ms=8, mfc='cyan', mec='C0', markeredgewidth=0.5 , alpha=1, zorder=1)

# Working only with the uncertain detections
# Select only the ones with optical counterpart
wt = rd[rd['specObjID']>1]
print(len(wt))
p3, = plt.plot(wt['W50'], wt['TotFlux'],'o', ms=5, mec='brown', mfc='brown', markeredgewidth=1.5, zorder=2)# ,alpha=0.55
# Select only the ones with potential otpical counterpart
wt = rd[rd['potOC']!='blah']
print(len(wt))
p4p, = plt.plot(wt['W50'], wt['TotFlux'],'o', ms=5, mec='brown', mfc='C1', markeredgewidth=0.5, zorder=2)# ,alpha=0.55
# Select only the ones without any optical counterpart
wt = rd[(rd['specObjID']==-1) & (rd['potOC']=='blah')]
print(len(wt))
p4, = plt.plot(wt['W50'], wt['TotFlux'],'o', ms=5, mec='brown', mfc='None', markeredgewidth=0.5, zorder=2)# ,alpha=0.55


ws = np.array(range(1,400))
plt.plot(ws, give_f(ws), color='black', linestyle='--', linewidth=0.5, zorder=0)
ws = np.array(range(400, 935))
plt.plot(ws, give_f400(ws), color='black', linestyle='--', linewidth=0.5, zorder=0)
plt.yscale('log')
plt.xscale('log')
plt.ylabel('F$_{tot}$ [Jy km s$^{-1}$]')
plt.xlabel('W$_{50}$ [kms$^{-1}$]')
plt.xticks([10,50,100,500],[10,50,100,500])
plt.xlim([15,780])
plt.ylim([0.03,95])

#################################################
# Make a legend (following an example from https://stackoverflow.com/questions/25830780/tabular-legend-layout-for-matplotlib)
extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
legend_handle = [extra, extra, extra, extra, extra, p1, p2p, p2, extra, p3, p4p, p4]

# Define the labels
label_row_1 = ['spec. OC', 'pot. OC', 'no OC']
label_j_1 = ['certain']
label_j_2 = ['uncertain']
label_empty = [""]

# Organize labels for table construction
legend_labels = np.concatenate([label_empty , label_row_1, label_j_1, label_empty*3, label_j_2, label_empty*4, label_j_2])

# Create legend
ax.legend(legend_handle, legend_labels, 
          loc = 2, ncol = 3, shadow = False, handletextpad = -2)
#################################################

plt.tight_layout()
plt.savefig(path+'final_analysis/Fl_W50_v3.pdf')
