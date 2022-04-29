#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 13:21:16 2021
@author: Boris Deshev

Calculate HI deficiency based on the M_HI/M_* v M_* relation published
by the XGASS survey (Catinela et al. 2018)

Also do the same fit to our data
"""

from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import models, fitting
from astropy.constants import c
from scipy.stats import ks_2samp

H0 = 71 # km/s/Mpc
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=H0, Om0=0.3, Tcmb0=2.725)
c = c.value/1000 # in km/s


######################################
# # Needed to write correctly HI
import matplotlib
matplotlib.rcParams['pgf.texsystem'] = 'pdflatex'
matplotlib.rcParams.update({'font.family': 'serif', 'font.size': 24,
    'axes.labelsize': 26,'axes.titlesize': 24, 'figure.titlesize' : 28})
matplotlib.rcParams['text.usetex'] = True
#'{\sc Hi}'
#r'log M$_{\textrm{\sc Hi}}$'
######################################

def get_upper_limit(tbl):
    """
    Calculate the upper limit of the HI mass
    given that this is a non detection
    
    MPA-JHU also gives velocity dispersion estimate, although not for all galaxies.
    However it is not even close to matching the W50 of HI detected galaxies
    While the method below does!
    
    Calculate the lower limit in HI deficiency for the non detected SDSS spec galaxies
    The expected amount of HI should already be calculated
    The plan here is to calculate their absolute magnitude, use TF relation to turn this into W50 (from Masters et al. 2006).
    Then use the median noise in the HI cube to calculate what their maximum undetected HI could be
    and then calculate with that a lower limit on HI deficiency.
    This will assume that they are all edge-on, so it is an upper-lower limit (...wtf)
    """
    rms = np.nanmean(tbl['rms'])/1000 # Originally the rms is in mJy
    chan_width = 5.5 # in km/s

    # Get the absolute i-band magnitude
    tbl['distmod'] = cosmo.distmod(tbl['Z']).value
    tbl['abs_i'] = tbl['petroMag_i']-tbl['distmod']

    # The published TF relation by Masters et al. 2006
    # absMag = -20.85 - 7.85 * (np.log10(tbl['exp_W'])-2.5)
    x = (tbl['abs_i']+20.85)/-7.85
    x += 2.5
    tbl['exp_W'] = 10**x
    
    # Number of channels
    Nchan = tbl['exp_W'] / chan_width
    # Total maximum undetected flux
    tbl['max_flux'] = Nchan * rms * chan_width
    
    # Maximum undetected MHI
    # Use the z-derived distance corrected for FoG
    tbl['MHI_up_lim'] = np.log10(2.36e5 * tbl['dist_corr']**2 * tbl['max_flux'])
    # Lowest limit on HI deficiency
    tbl['HIdef_low_lim'] = tbl['MHI_exp_SP'] - tbl['MHI_up_lim']

    return tbl

def estimate_size_correction(make_plot=False):
    """
    The petroR90_r size of the galaxies will be corrected
    to better represent D25, which is what it should be.
    This is done by comparing the two size estimates for a sample of galaxies
    common between the Goldmine project and SDSS

    Returns
    -------
    The correction factor between petroR90_r and D25
    """

    #######################################
    # The file that contains the two size estimates
    # a is D25 in arcmin; petro_R90_r is the r-band petrosian radius in arcsec
    ss = Table.read(path+'../analysis_preliminary/size_sample_goldmine_sdss.csv')
    # The major axis from Goldmine is in arcminutes
    a = ss['a'].data * 60
    # Get the diameter from SDSS
    r90 = ss['petroR90_r'].data * 2

    # Take the histogram
    bins = np.linspace(0, 3, 15)
    vals, bins = np.histogram(a/r90, bins=bins)
    x = (bins[:-1] + bins[1:]) / 2

    # Fit a Gaussian to it
    g_init = models.Gaussian1D(amplitude = np.max(vals), mean=1.5, stddev=0.5)
    fit_g = fitting.LevMarLSQFitter()
    gauA = fit_g(g_init, x, vals)
    corrF = gauA.mean
    print('correction {:.1f} +- {:.2f}'.format(corrF.value, gauA.stddev.value))#, "%.1f" %corrF.value, '+-', "%.2f" %gauA.stddev.value)

    if make_plot:
        plt.figure()
        plt.plot(x, vals, 'o')
        plt.plot(x, gauA(x))
        plt.show()
    
    return corrF


def calc_deficiency(tbl):
    """
    This is the way to calculating deficiency, based on Haynes & Giovanelli, 1984 and petroR_90
    We only use the recipe for late type galaxies becuase from the HI detected galaxies <10% are early type
    This of course does not mean much for the non-detected ones but the point of this article
    anyway is to show that there are better ways to assess environmental effects
    than calculating deficiency ...
    
    Input
    -------
    Table with column petroR90_r
    
    Returns
    -------
    Table with added columns:
        physical size
        expected HI mass
        HI deficiency
    """
    # From Haynes & Giovanelli, 1984
    # The coefficients are
    # late types: 
    emc1 = 7.195; emc2 = 0.851
    # # early types: 
    # psc1 = 6.88; psc2 = 0.89
    
    # Correction to the size
    corrF = estimate_size_correction()
    
    tbl['MHI_exp_SP'] = -999.0
    tbl['petroR90_r_phys'] = 0.0

    # Get the angular to physical scale
    # based on the FoG corrected distance
    scale = cosmo.kpc_proper_per_arcmin((tbl['dist_corr']*H0)/c).value/60. #get the scale in kpc/arcsec

    # The size of the galaxy in kpc
    tbl['petroR90_r_phys'] = tbl['petroR90_r']*scale

    # The expected HI mass
    tbl['MHI_exp_SP'] = emc1+emc2*np.log10((2*tbl['petroR90_r_phys'] * corrF)**2)

    # Calculate the deficiency
    tbl['HIdef'] = tbl['MHI_exp_SP'] - tbl['MHI']

    return tbl



########################################
########################################
############   MAIN   ##################
########################################
########################################
path = '/home/tazio/works/2020/a1367/final_analysis/'
tbl = Table.read(path+'mpa_jhu_AGES.fits')
print('tbl',len(tbl))
    
tbl = calc_deficiency(tbl)

tbl = get_upper_limit(tbl)

p_val = ks_2samp(tbl['W50'], tbl['exp_W']).pvalue
print('p-value for the W50 distributions of detections and non-detection {:.3e}'.format(p_val))

# Mask the column with the upper limits so it is the opposite of the directly measured one
DMmask = tbl['MHI'].mask
UPmask = np.logical_not(DMmask)
tbl['MHI_up_lim'].mask = UPmask
tbl['HIdef_low_lim'].mask = UPmask
# Make a column containing both estimates of deficiency, it will be needed later
tbl['HIdef_ALL'] = tbl['HIdef'].filled(1) * tbl['HIdef_low_lim'].filled(1)
tbl['HIdef_ALL'] = np.ma.masked_where(tbl['HIdef_ALL'] > 3, tbl['HIdef_ALL'])

tbl.write(path+'mpa_jhu_AGES.fits', overwrite=True)





































