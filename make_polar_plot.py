#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  2 10:40:00 2022
@author: Boris Deshev

Plot the distribution of objects in RA-Vel and Dec-Vel planes
with the correct axis ratio
"""

from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import models, fitting
from astropy.constants import c
c = c.value/1000 # In km/s

###############################################
# # Needed to write correctly HI
import matplotlib
matplotlib.rcParams['pgf.texsystem'] = 'pdflatex'
matplotlib.rcParams.update({'font.family': 'serif', 'font.size': 20,
    'axes.labelsize': 24,'axes.titlesize': 24, 'figure.titlesize' : 24})
matplotlib.rcParams['text.usetex'] = True
#'{\sc Hi}'
#r'log M$_{\textrm{\sc Hi}}$'
##############################################

def plot_rfi():
    # Find the places enriched by RFI and produce a map of those
    from astropy.io import fits
    from astropy import wcs
    from astropy.convolution import convolve, Box2DKernel
    # The cube with all HI sources masked
    f = fits.open(path+'cubes/A1367_no_objects.fits')
    # Get the coordinates of the extremeties
    proj1 = wcs.WCS(f[0].header)
    x1,y1,z1 = proj1.all_pix2world(0,0,0, 1)
    x2,y2,z2 = proj1.all_pix2world(f[0].header['NAXIS1'],f[0].header['NAXIS2'],f[0].header['NAXIS3'], 1)

    # Calculate the presence of RFI
    # First make a map on ra-dec plane
    rfi_map1 = np.nanmean(f[0].data, axis=0)
    # There are locations with only NaNs
    rfi_map1[rfi_map1 != rfi_map1] = 1
    # Then use this to flag columns along the z axis
    # First enlarge the masked areas
    mmask = np.zeros_like(f[0].data[0,:,:])
    mmask[rfi_map1 > 0.005] = 1
    sm_kernel = Box2DKernel(5)
    enmmask = convolve(mmask,sm_kernel)
    # Make a 3d mask for the whole cube
    field3d_mask = np.broadcast_to(enmmask > 0.01, f[0].data.shape)
    # Set the masked areas to nan
    f[0].data[field3d_mask] = np.nan

    # Make a map in the RA-z plane
    # This is a map of the noise
    rfi_map2 = np.nanstd(f[0].data, axis=1)
    # Check again for nans
    rfi_map2[rfi_map2 != rfi_map2] = 1
    # Make it binary
    # First we need to find a good cut-off value
    # For that we will fit the nosie distribution
    # and cut at 4.5sigma
    
    bins=np.linspace(0.001, 0.005, 50)
    v, b = np.histogram(np.ravel(rfi_map2), bins = bins)
    x = (b[1:]+b[:-1])/2
    # Fit a Gaussian to it
    g_init = models.Gaussian1D(amplitude = np.max(v), mean=0.002, stddev=0.0005)
    fit_g = fitting.LevMarLSQFitter()
    gf = fit_g(g_init, x, v)
    coff = gf.mean.value+4.5*gf.stddev.value
    print('\nThe cut-off is at {:.5f}'.format(coff))
    
    rfi_map2[abs(rfi_map2) >= coff] = 1
    rfi_map2[abs(rfi_map2) < coff] = 0

    # Manage the coordinates of the map
    ras = np.linspace(x2, x1, f[0].header['NAXIS1'])
    zeds = np.linspace(z1/1000, z2/1000, f[0].header['NAXIS3'])
    rav, zv = np.meshgrid(ras, zeds)

    
    # Repeat the same but in the Dec-z plane
    rfi_map3 = np.nanstd(f[0].data, axis=2)
    
    bins=np.linspace(0.001, 0.005, 50)
    v, b = np.histogram(np.ravel(rfi_map3), bins = bins)
    x = (b[1:]+b[:-1])/2
    # Fit a Gaussian to it
    g_init = models.Gaussian1D(amplitude = np.max(v), mean=0.002, stddev=0.0005)
    fit_g = fitting.LevMarLSQFitter()
    gf = fit_g(g_init, x, v)
    coff3 = gf.mean.value+4.5*gf.stddev.value
    print('\nThe cut-off is at {:.5f}'.format(coff3))
    
    # Check again for nans
    rfi_map3[rfi_map3 != rfi_map3] = 1
    rfi_map3[abs(rfi_map3) >= coff3] = 1
    rfi_map3[abs(rfi_map3) < coff3] = 0
    # Manage the coordinates of the map
    decs = np.linspace(y2, y1, f[0].header['NAXIS2'])
    zeds = np.linspace(z1/1000, z2/1000, f[0].header['NAXIS3'])
    dev, zdv = np.meshgrid(decs, zeds)
    
    
    return rav, zv, rfi_map2, dev, zdv, rfi_map3

#############################
#############################
#############################
# # Open the tables to plot
path = '/home/tazio/works/2020/a1367/'
# Table with HI detections
HItbl = Table.read(path+'A1367_HItable_v8.8_wt.txt', format='ascii.commented_header')
# Table with MPA-JHU objects in the AGES volume
Otbl = Table.read(path+'final_analysis/mpa_jhu_AGES.fits')

# Select only the mass complete part
print(len(Otbl))
OtblD = Otbl[Otbl['LGM_TOT_P50'] > 9.75]
print(len(OtblD))

# There are sources with multiple HI flags
# Here it is only important if the flag is 1 or not
for i in range(len(HItbl)):
    HItbl[i]['flag'] = HItbl[i]['flag'].strip()[0]

# The potential optical counterparts
POC = HItbl[HItbl['potOC'] != 'blah']

# Get the map of the RFI
rav, zv, RAZ_rfi_map, dev, zdv, DEZ_rfi_map = plot_rfi()
RAZ_mask = np.where(RAZ_rfi_map==1)
DEZ_mask = np.where(DEZ_rfi_map==1)

#############################
# Define the parameters of the plot

# The axes of each panel
# Should be Left, Bottom, Width, Height ... and it almost is.
axN = [[0.01, 0.01, 0.95, 1.5],
       [0.01, 0.15, 0.95, 1.5],
       [0.03, -0.195, 0.965, 1.5],
       [0.03, -0.46, 0.96, 1.5]]

# Which column contains the Y variable for different sources
Ycol1 = ['RA', 'DEC', 'RA', 'DEC']
Ycol2 = ['R.A.(deg)', 'Dec.(deg)', 'R.A.(deg)', 'Dec.(deg)']
Ycol3 = ['RA_NED', 'DEC_NED', 'RA_NED', 'DEC_NED']

# Angular range
theta_min = [173, 17.5, 173, 17.5]
theta_max = [179, 22.2, 179, 22.2]

# Radial range
rmin = [1000, 1000, 9800, 9800]
rmax = [10200, 10200, 19300, 19300]

# Used to set the tick positions
yticks = [1000, 1000, 10000, 10000]

# Orientation of the separate panels
tzl = ['N', 'N', 'S', 'S']

# Used to align the Dec and RA panels
offset = np.array([0, 161.5, 0, 161.5])
offset+=90

# Label on the angular axis
xlabel = ['R.A. (J2000)', 'Dec. (J2000)', 'R.A. (J2000)', 'Dec. (J2000)']
# Used for locating the above axis label
xlab_offset = [1.03, 1.03, 1.01, 1.01]
xlab_rot = [270,270,90,90]

# The mask and coordinate grid for the RFI masks in the two planes
mask = [RAZ_mask, DEZ_mask, RAZ_mask, DEZ_mask]
mgrid = [[rav, zv], [dev, zdv], [rav, zv], [dev, zdv]]

# Tick labels on/off. left, right, top, bottom; four panels
tick_P = [[1,0,0,1],
          [0,1,0,1],
          [1,0,0,1],
          [0,1,0,1]]

# Used for the legends
Olabels = ['','', 'optical \ngalaxy \nwith M$_{\star} \leq 9.75$', '']
POC_lab = ['','', 'potential \noptical \ncounterpart', '']
HIlab = [['',''],['',''],['certain', 'uncertain'],['','']]


fig = plt.figure(figsize=(30,20))
for n in range(4):
    ax = plt.axes(axN[n], polar=True) 

    # Plot the RFI
    ax.plot(np.deg2rad(mgrid[n][0][mask[n]]), mgrid[n][1][mask[n]], 'D', color='0.8', zorder=1)
    
    # Plot the optical sources
    ax.plot(np.deg2rad(Otbl[Ycol1[n]]), Otbl['Z']*c, marker='.', ms=7, linestyle='None', color='black', zorder=1, label=Olabels[n])
    mappable = ax.scatter(np.deg2rad(OtblD[Ycol1[n]]), OtblD['Vel'], marker='s', s=50, c=np.log10(OtblD['density']), cmap='viridis', vmin=-3, vmax=1.0, zorder=2)

    # Plot the potential optical counterparts at the measured RA but the HI Vel
    ax.plot( np.deg2rad(POC[Ycol3[n]]), POC['Vel_HI'], '.', ms=7, color='brown', label=POC_lab[n], zorder=2 )
    
    # Overplot the HI detections
    ax.plot(np.deg2rad(HItbl[Ycol2[n]][HItbl['flag']!='1']), HItbl['Vel_HI'][HItbl['flag']!='1'], 's', mec='black', mfc='None', ms=8, zorder=3, label=HIlab[n][0])
    ax.plot(np.deg2rad(HItbl[Ycol2[n]][HItbl['flag']=='1']), HItbl['Vel_HI'][HItbl['flag']=='1'], 'o', mec='0.3', mfc='None', ms=8, zorder=3, label=HIlab[n][1])

    if n == 0:
        ax.text(np.deg2rad(HItbl[Ycol2[n]].max()+8.0), rmin[n]+4500, '$cz$ [km s$^{-1}$]',
                rotation = 5, ha = 'center', va = 'top', fontsize=36)
    elif n == 2:
        # Make the two legends
        handles, labels = ax.get_legend_handles_labels()
        first_legend = ax.legend(handles[2:], labels[2:], title = '{\sc Hi} detections', loc='upper center', bbox_to_anchor=(0.75, 0.23), fancybox=True, shadow=True, fontsize=24)
        second_legend = ax.legend(handles[:2], labels[:2], loc='lower center', bbox_to_anchor=(0.91, 0.15), fancybox=True, shadow=True, fontsize=24)
        ax.add_artist(first_legend)

        # Take care of the colorbar
        cbaxes = fig.add_axes([0.1, 0.1, 0.5, 0.02]) 
        cbar = plt.colorbar(mappable, cax = cbaxes, orientation='horizontal')  
        cbar.set_label(r'log($\rho$) [Mpc$^{-3}$]', fontsize=28, rotation=0)#, labelpad=25)
        cbar.ax.yaxis.set_label_position('right')

        # Axis label
        ax.text(np.deg2rad(HItbl[Ycol2[n]].max()+5.5), rmin[n]+4700, '$cz$ [km s$^{-1}$]',
                rotation = 5, ha = 'center', va = 'top', fontsize=36)

    ax.set_thetamin(theta_min[n])
    ax.set_thetamax(theta_max[n])
    ax.set_rmax(rmax[n])
    ax.set_rmin(rmin[n])
    ax.set_rticks(np.arange(yticks[n], rmax[n], 1000))

    if n in [0,1]:
        for label in ax.get_xticklabels():
            label.set_horizontalalignment("left")

    ax.tick_params(labelleft=tick_P[n][0], labelright=tick_P[n][1], labeltop=tick_P[n][2], labelbottom=tick_P[n][3])
    
    ax.set_rorigin(0)
    ax.set_theta_zero_location(loc=tzl[n], offset=offset[n])
    ax.text(np.deg2rad(HItbl[Ycol2[n]].mean()), rmax[n]*xlab_offset[n], xlabel[n], rotation=xlab_rot[n], ha = 'center', va = 'center', fontsize=24)

plt.tight_layout()
plt.savefig(path+'final_analysis/polar_plot_M90.jpg', dpi=200)















































