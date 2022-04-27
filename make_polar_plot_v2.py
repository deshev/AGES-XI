#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:35:55 2022
@author: tazio

Plot the distribution of objects in RA-Vel and Dec-Vel planes
with the correct axis ratio
"""
import sys
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
    # The cube without any sources
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
    print('\nThe cut-off value in std is {:.5f}'.format(coff))
    
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
    print('\nThe cut-off value in std is {:.5f}'.format(coff3))

    
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
HItbl = Table.read(path+'A1367_HItable_v8.8_wt.txt', format='ascii.commented_header')
Otbl = Table.read(path+'final_analysis/mpa_jhu_AGES.fits')

print(len(Otbl))
OtblD = Otbl[Otbl['LGM_TOT_P50'] > 9.75]
print(len(OtblD))

# The potential optical counterparts
POC = HItbl[HItbl['potOC'] != 'blah']

# Get the map of the RFI
rav, zv, RAZ_rfi_map, dev, zdv, DEZ_rfi_map = plot_rfi()
RAZ_mask = np.where(RAZ_rfi_map==1)
DEZ_mask = np.where(DEZ_rfi_map==1)

#############################
# Make the first plot - low z
fig = plt.figure(figsize=(9,30))
ax = plt.axes([-1.20, 0.01, 3.5, 0.95], polar=True) # Should be Left, Bottom, Width, Height ... and it almost is.

# Plot the RFI
ax.plot(np.deg2rad(rav[RAZ_mask]), zv[RAZ_mask], 'D', color='0.8', zorder=1)

# Plot the optical sources
ax.plot(np.deg2rad(Otbl['RA']), Otbl['Z']*c, marker='.', ms=7, linestyle='None', color='black', zorder=1)
mappable = ax.scatter(np.deg2rad(OtblD['RA']), OtblD['Vel'], marker='s', s=50, c=np.log10(OtblD['density']), cmap='viridis', vmin=-3, vmax=1.0, zorder=2)

# Plot the potential optical counterparts at the measured RA but the HI Vel
ax.plot( np.deg2rad(POC['RA_NED']), POC['Vel_HI'], '.', ms=7, color='brown', zorder=2 )

# There are sources with multiple flags
# Here it is only important if the flag is 1 or not
for i in range(len(HItbl)):
    HItbl[i]['flag'] = HItbl[i]['flag'].strip()[0]
# Overplot the HI detections
ax.plot(np.deg2rad(HItbl['R.A.(deg)'][HItbl['flag']!='1']), HItbl['Vel_HI'][HItbl['flag']!='1'], 's', mec='black', mfc='None', ms=8, zorder=3, label='certain')
ax.plot(np.deg2rad(HItbl['R.A.(deg)'][HItbl['flag']=='1']), HItbl['Vel_HI'][HItbl['flag']=='1'], 'o', mec='0.3', mfc='None', ms=8, zorder=3, label='uncertain')
ax.legend(title = '{\sc Hi} detections', loc='upper left', fancybox=True, shadow=True, bbox_to_anchor=(0.5, 0.1), fontsize=24)

ax.set_thetamin(173)
ax.set_thetamax(179)
rmin = 1000; rmax = 10200
ax.set_rmax(rmax)
ax.set_rmin(rmin)
ax.set_rticks(np.arange(rmin, rmax, 1000))
ax.tick_params(labelleft=False, labelright=False,
                labeltop=False, labelbottom=True)

ax.set_rorigin(0)
ax.set_theta_zero_location(loc='S', offset=0)
ax.text(np.deg2rad(HItbl['R.A.(deg)'].mean()), rmax*1.03, 'R.A. (J2000)', ha = 'center', va = 'top', fontsize=24)


################################################################## Dec
ax2 = plt.axes([-1.51, 0.01, 3.5, 0.95], polar=True) # Should be Left, Bottom, Width, Height ... and it almost is.

# Plot the optical sources
ax2.plot(np.deg2rad(Otbl['DEC']), Otbl['Z']*c, marker='.', ms=7, linestyle='None', color='black', zorder=1)
mappable = ax2.scatter(np.deg2rad(OtblD['DEC']), OtblD['Vel'], marker='s', s=50, c=np.log10(OtblD['density']), cmap='viridis', vmin=-3, vmax=1.0, zorder=2)

# Plot the RFI
ax2.plot(np.deg2rad(dev[DEZ_mask]), zdv[DEZ_mask], 'D', color='0.8', zorder=1)

# Plot the potential optical counterparts at the measured RA but the HI Vel
POC = HItbl[HItbl['potOC'] != 'blah']
ax2.plot( np.deg2rad(POC['DEC_NED']), POC['Vel_HI'], '.', ms=7, color='brown', zorder=2 )

# There are sources with multiple flags
# Here it is only important if the flag is 1 or not
for i in range(len(HItbl)):
    HItbl[i]['flag'] = HItbl[i]['flag'].strip()[0]
# Overplot the HI detections
ax2.plot(np.deg2rad(HItbl['Dec.(deg)'][HItbl['flag']!='1']), HItbl['Vel_HI'][HItbl['flag']!='1'], 's', mec='black', mfc='None', ms=8, zorder=3, label='certain')
ax2.plot(np.deg2rad(HItbl['Dec.(deg)'][HItbl['flag']=='1']), HItbl['Vel_HI'][HItbl['flag']=='1'], 'o', mec='0.3', mfc='None', ms=8, zorder=3, label='uncertain')


ax2.set_thetamin(17.5)
ax2.set_thetamax(22.2)
rmin = 1000; rmax = 10200
ax2.set_rmax(rmax)
ax2.set_rmin(rmin)
ax2.set_rticks(np.arange(rmin, rmax, 1000))
ax2.tick_params(labelleft=False, labelright=True,
                labeltop=False, labelbottom=True)#, 
                #rotation=93, ha="center", va="center")

ax2.set_rorigin(0)
ax2.set_theta_zero_location(loc='S', offset=161.5)
ax2.text(np.deg2rad(HItbl['Dec.(deg)'].mean()), rmax*1.03, 'Dec. (J2000)', ha = 'center', va = 'top', fontsize=24)
ax2.text(np.deg2rad(HItbl['Dec.(deg)'].max()+2), rmin+4900, '$cz$ [km s$^{-1}$]',
        rotation = 95, ha = 'center', va = 'top', fontsize=36)

plt.tight_layout()
plt.savefig(path+'final_analysis/polar_plot_M1.jpg', dpi=200)


#############################
# Make the second plot - high z
fig = plt.figure(figsize=(13,30))
ax = plt.axes([-1.4, 0.05, 3.5, 0.93],polar=True) # Should be Left, Bottom, Width, Height ... and it almost is.

# Plot the RFI
ax.plot(np.deg2rad(rav[RAZ_mask]), zv[RAZ_mask], 'D', color='0.8', zorder=1)

# Plot the optical sources
ax.plot(np.deg2rad(Otbl['RA']), Otbl['Z']*c, marker='.', ms=7, linestyle='None', color='black', zorder=1, label='optical \ngalaxy \nwith M$_{\star} \leq 9.75$')
mappable = ax.scatter(np.deg2rad(OtblD['RA']), OtblD['Vel'], marker='s', s=50, c=np.log10(OtblD['density']), cmap='viridis', vmin=-3, vmax=1.0, zorder=2)
# Take care of the colorbar
cbaxes = fig.add_axes([0.05, 0.15, 0.02, 0.4]) 
cbar = plt.colorbar(mappable, cax = cbaxes)  
cbar.set_label(r'log($\rho$) [Mpc$^{-3}$]', fontsize=28)
cbar.ax.yaxis.set_label_position('left')

# Plot the potential optical counterparts at the measured RA but the HI Vel
ax.plot( np.deg2rad(POC['RA_NED']), POC['Vel_HI'], '.', ms=7, color='brown', label='potential \noptical \ncounterpart', zorder=3 )

# Overplot the HI detections
ax.plot(np.deg2rad(HItbl['R.A.(deg)'][HItbl['flag']!='1']), HItbl['Vel_HI'][HItbl['flag']!='1'], 's', mec='black', mfc='None', ms=8, zorder=3)
ax.plot(np.deg2rad(HItbl['R.A.(deg)'][HItbl['flag']=='1']), HItbl['Vel_HI'][HItbl['flag']=='1'], 'o', mec='0.3', mfc='None', ms=8, zorder=3)

# Add a legend
ax.legend(bbox_to_anchor=(0.462, 0.71), loc='upper right', fontsize=24, fancybox=True, shadow=True)


ax.set_thetamin(173)
ax.set_thetamax(179)
rmin = 9800; rmax = 19300
ax.set_rmax(rmax)
ax.set_rmin(rmin)
ax.set_rticks(np.arange(10000, rmax, 1000))
ax.tick_params(labelleft=False, labelright=False,
                labeltop=False, labelbottom=True)

ax.set_rorigin(0.0)
ax.set_theta_zero_location(loc='N', offset=0)

ax.text(np.deg2rad(HItbl['R.A.(deg)'].mean()), rmax*1.015, 'R.A. (J2000)', ha = 'center', va = 'top', fontsize=24)

################################################################## Dec
ax2 = plt.axes([-1.01, 0.05, 3.5, 0.925], polar=True) # Should be Left, Bottom, Width, Height ... and it almost is.

# Plot the optical sources
ax2.plot(np.deg2rad(Otbl['DEC']), Otbl['Z']*c, marker='.', ms=7, linestyle='None', color='black', zorder=1)
mappable = ax2.scatter(np.deg2rad(OtblD['DEC']), OtblD['Vel'], marker='s', s=50, c=np.log10(OtblD['density']), cmap='viridis', vmin=-3, vmax=1.0, zorder=2)

# Plot the RFI
ax2.plot(np.deg2rad(dev[DEZ_mask]), zdv[DEZ_mask], 'D', color='0.8', zorder=1)

# Plot the potential optical counterparts at the measured RA but the HI Vel
POC = HItbl[HItbl['potOC'] != 'blah']
ax2.plot( np.deg2rad(POC['DEC_NED']), POC['Vel_HI'], '.', ms=7, color='brown', zorder=2 )

# There are sources with multiple flags
# Here it is only important if the flag is 1 or not
for i in range(len(HItbl)):
    HItbl[i]['flag'] = HItbl[i]['flag'].strip()[0]
# Overplot the HI detections
ax2.plot(np.deg2rad(HItbl['Dec.(deg)'][HItbl['flag']!='1']), HItbl['Vel_HI'][HItbl['flag']!='1'], 's', mec='black', mfc='None', ms=8, zorder=3, label='certain')
ax2.plot(np.deg2rad(HItbl['Dec.(deg)'][HItbl['flag']=='1']), HItbl['Vel_HI'][HItbl['flag']=='1'], 'o', mec='0.3', mfc='None', ms=8, zorder=3, label='uncertain')


ax2.set_thetamin(17.5)
ax2.set_thetamax(22.2)
ax2.set_rmax(rmax)
ax2.set_rmin(rmin)
ax2.set_rticks(np.arange(10000, rmax, 1000))
ax2.tick_params(labelleft=False, labelright=True,
                labeltop=False, labelbottom=True)


ax2.set_rorigin(0.0)
ax2.set_theta_zero_location(loc='N', offset=161.6)
ax2.text(np.deg2rad(HItbl['Dec.(deg)'].mean()), rmax*1.015, 'Dec. (J2000)', ha = 'center', va = 'top', fontsize=24)
ax2.text(np.deg2rad(HItbl['Dec.(deg)'].max()+2), rmin+4200, '$cz$ [km s$^{-1}$]',
        rotation = -84, ha = 'center', va = 'top', fontsize=36)



plt.tight_layout()
plt.savefig(path+'final_analysis/polar_plot_M2.jpg', dpi=200)














