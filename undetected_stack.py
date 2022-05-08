#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  4 10:06:40 2022
@author: Boris Deshev

Just for a test
Add all the HI spectra of the detections together
and compare it to the stack of all SDSS galaxies
"""
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.constants import c
from astropy.cosmology import FlatLambdaCDM
from astropy.stats import jackknife_stats, jackknife_resampling
import sys
from scipy.stats import iqr

######################################
# # Needed to write correctly HI
import matplotlib
matplotlib.rcParams['pgf.texsystem'] = 'pdflatex'
matplotlib.rcParams.update({'font.family': 'serif', 'font.size': 20,
    'axes.labelsize': 24,'axes.titlesize': 20, 'figure.titlesize': 24})
matplotlib.rcParams['text.usetex'] = True
#'{\sc Hi}'
#r'log M$_{\textrm{\sc Hi}}$'
######################################

H0 = 71
cosmo = FlatLambdaCDM(H0=H0, Om0=0.3, Tcmb0=2.725)
c = c.value/1000

def get_spectrum(spname, i, tbl, to_mass = True):
    wl = []; fl = []
    sp = open(spname,'r')
    for line in sp.readlines():
        if '#' not in line:
            wl.append(float(line.strip().split()[1]))
            fl.append(float(line.strip().split()[2]))
    # Trim it to +-500kms (Don't ask why 500)
    ie = np.argmin(abs(np.array(wl)-(tbl[i]['Vel_HI']-550)))
    ib = np.argmin(abs(np.array(wl)-(tbl[i]['Vel_HI']+550)))
    flP = np.array(fl[ib:ie])
    wlP = np.array(wl[ib:ie])

    # Convert the flux to HI mass
    if to_mass:
        mhiS = 2.36e5 * flP * (cosmo.luminosity_distance(tbl[i]['Z']).value**2)
    else:
        mhiS = flP
    # Interpolate so that all are the same length
    intel = 200 # This is how many points the interpolated spectrum will have
    rat = intel/len(mhiS) #the factor by which the flux is elevated in the interpolated data
    interp_mhi_spec = np.interp(np.linspace(wlP[-1], wlP[0], intel), \
                                np.linspace(wlP[-1], wlP[0], len(mhiS)), \
                                mhiS) / rat
    
    return interp_mhi_spec


def run_Stbl(Stbl, rs = False):
    # Stack the otpical redshifts
    specs = []; weights = []
    for i in range(len(Stbl)):
        if rs:
            sp = Table.read(path+'spectra/for_reference_spectrum/'+Stbl[i]['SPECOBJID']+'_han5spec.txt', format='ascii.commented_header')
        else:
            sp = Table.read(path+'spectra/HI_for_stacking/'+Stbl[i]['SPECOBJID']+'_han5spec.txt', format='ascii.commented_header')
        sp['mhi'] = 2.36e5 * sp['Intensity'] * cosmo.luminosity_distance(Stbl[i]['Z']).value**2
        # The spectra extracted with miriad have different length, for some reason. The difference is a function
        # of z, even though there is no hint of anything in the header
        # So interpolation is needed
        # This part does the interpolation and deals with flux preservation. 
        # As I interpolate the data, the flux is elevated by the ratio of the total number of sampling points before and after the interpolation
        intel = 200 # This is how many points the interpolated spectrum will have
        rat = intel/len(sp['mhi']) #the factor by which the flux is elevated in the interpolated data
        interp_mhi_spec = np.interp(np.linspace(sp['Velocity'][-1], sp['Velocity'][0], intel), \
                                    np.linspace(sp['Velocity'][-1], sp['Velocity'][0], len(sp['mhi'])), \
                                    sp['mhi'].data) / rat
    
        specs.append(interp_mhi_spec)
        # Get the rms noise for weighting.
        # Doing it in MHI is, I think, the way to do it
        # The HISS paper Healy et al. 2019 recomends using the inter-quartile range instead of std
        weights.append(  1 / iqr(list(interp_mhi_spec)[:40] + list(interp_mhi_spec)[-40:]) )
        #print(Stbl[i]['SPECOBJID'], weights[-1])
    
    ms = np.average(specs, axis=0, weights=weights)
    cv = np.mean([np.mean(ms[:40]), np.mean(ms[-40:])])
    ms += abs(cv)

    return ms


def do_jackknife(tbl, NI):
    # Do N iterations of randomly selected 90% of the sample
    # this is an attempt to quantify the differences between 
    # individual galaxies in the stack.
    print('\n... Jacking and knifing ...')
    
    # Will stack this many randomly selected gals on every iteration
    N = int( len( tbl ) * 0.9 )
    # Indices
    all_indices = np.arange( len( tbl ) )
    
    mhi = []; msS = []; rsS = []
    for n in range(NI):
        selected_indices = np.random.choice( all_indices, size = N, replace = False )
        # Get the average spectrum
        ms = run_Stbl( tbl[selected_indices] )
        msS.append( ms )
        # Calculate and store the total HI mass
        # I'm using the same integral range for all spectra. It is
        # esentially the same as determined by the growth curve for
        # the stack of all the spectra
        mhi.append( np.sum( ms[50:150] ) * ch )
        # Get the reference spectrum for the same stack
        rsS.append( run_Stbl( tbl[selected_indices], rs = True ) )
    
    return msS, mhi, rsS


def first_time(CzT, EzT, FzT = None, NI = 100, write = False):
    # If it's the first time this si being run
    # this will do the stacking and write files with the spectra
    # and total HI mass

    stack_dC, mhiC, ref_dC = do_jackknife(CzT, NI)
    if write:
        # Save them for later use
        np.savetxt(path+'final_analysis/stacks/stack_dC_ref'+str(NI)+'.csv', np.array(ref_dC).T, fmt='%.5e', delimiter=',')
        np.savetxt(path+'final_analysis/stacks/stack_dC_'+str(NI)+'.csv', np.array(stack_dC).T, fmt='%.5e', delimiter=',')
        np.savetxt(path+'final_analysis/stacks/stack_dC_'+str(NI)+'_MHI.csv', np.array(mhiC).T, fmt='%.5f')

    # Once again
    stack_dE, mhiE, ref_dE = do_jackknife(EzT, NI)
    if write:
        np.savetxt(path+'final_analysis/stacks/stack_dE_ref'+str(NI)+'.csv', np.array(ref_dE).T, fmt='%.5e', delimiter=',')
        np.savetxt(path+'final_analysis/stacks/stack_dE_'+str(NI)+'.csv', np.array(stack_dE).T, fmt='%.5e', delimiter=',')
        np.savetxt(path+'final_analysis/stacks/stack_dE_'+str(NI)+'_MHI.csv', np.array(mhiE).T, fmt='%.5f')

    # See if there is another subset
    if FzT:
        stack_dF, mhiF, ref_dF = do_jackknife(FzT, NI)
        if write:
            np.savetxt(path+'final_analysis/stacks/stack_dF_ref'+str(NI)+'.csv', np.array(ref_dF).T, fmt='%.5e', delimiter=',')
            np.savetxt(path+'final_analysis/stacks/stack_dF_'+str(NI)+'.csv', np.array(stack_dF).T, fmt='%.5e', delimiter=',')
            np.savetxt(path+'final_analysis/stacks/stack_dF_'+str(NI)+'_MHI.csv', np.array(mhiF).T, fmt='%.5f')

    sys.exit()


def later_runs(dd, NI = 100, p = [50, 25, 75]):
    # If the stacking has already been done, read in the results
    # and extract the spectra with the required percentiles of the total MHI
    TC = Table.read(path+'final_analysis/stacks/stack_'+dd+'_'+str(NI)+'.csv', format='ascii.no_header')
    RC = Table.read(path+'final_analysis/stacks/stack_'+dd+'_ref'+str(NI)+'.csv', format='ascii.no_header')
    TCmhi = Table.read(path+'final_analysis/stacks/stack_'+dd+'_'+str(NI)+'_MHI.csv', format='ascii.no_header')
    
    # Normalise the MHI distribution
    mhi = TCmhi['col1'].data
    mhi_norm = mhi - np.min(mhi)
    mhi_norm /= np.max(mhi_norm)

    stack_dC = []; ref_dC = []
    print('Working on ', dd)
    thim = []
    for pp in p:
        ind = np.argmin( abs( mhi_norm - (pp/100) ) )
        stack_dC.append( TC['col'+str(ind+1)] )
        ref_dC.append( RC['col'+str(ind+1)] )
    
        # It will print out the total HI mass in the required spectra
        print('Total HI mass in the p'+str(pp)+' stack', "%.3e" %(TCmhi['col1'][ind]))
        thim.append(TCmhi['col1'][ind])
    print('i.e. {:.3e} + {:.3e} - {:.3e}'.format(thim[0], thim[2]-thim[0], thim[0]-thim[1]))

    return np.array(stack_dC), np.array(ref_dC)


def get_Jy_runs(dd, NI = 100, p = [50, 25, 75]):
    # If the stacking has already been done, read in the results
    # and extract the spectra with the required percentiles of the total MHI
    # This subroutine will deal with the stacks in Jy
    TC = Table.read(path+'final_analysis/stacks/stack_'+dd+'_'+str(NI)+'_Jy.csv', format='ascii.no_header')
    RC = Table.read(path+'final_analysis/stacks/stack_'+dd+'_ref'+str(NI)+'_Jy.csv', format='ascii.no_header')
    TCmhi = Table.read(path+'final_analysis/stacks/stack_'+dd+'_'+str(NI)+'_MHI_Jy.csv', format='ascii.no_header')
    
    # Normalise the MHI distribution
    mhi = TCmhi['col1'].data
    mhi_norm = mhi - np.min(mhi)
    mhi_norm /= np.max(mhi_norm)

    stack_dC = []; ref_dC = []
    print('Working on ', dd)
    for pp in p:
        ind = np.argmin( abs( mhi_norm - (pp/100) ) )
        stack_dC.append( TC['col'+str(ind+1)] )
        ref_dC.append( RC['col'+str(ind+1)] )
    
        # It will print out the total HI mass in the required spectra
        print('Total HI mass in the p'+str(pp)+' stack', "%.3e" %(TCmhi['col1'][ind]))

    return np.array(stack_dC), np.array(ref_dC)
   

def get_integration_borders(sp):
    # Calculate the region for integration by making a growth curve
    mp = int(len(sp)/2)
    gr_c = []
    
    # The growth curve starts from the middle, i.e. assuming symmetry
    for g in range(int(mp/2)):
        gr_c.append( np.sum(sp[mp-g:mp+g]) )
    gr_c = np.array(gr_c)
    grad_gr_c = gr_c[1:]-gr_c[:-1]
    max_growth = np.argmin(grad_gr_c)
    
    return max_growth


def make_a_plot_all():
    # The region over which the integral is taken to calculate the total HI mass
    max_growth = 50
    
    # Do the figure
    titles = ['Foreground', 'A1367', 'Background']
    ylim_mins = [-0.025,-0.06,-0.23,-0.001,-0.0004,-0.00025]
    
    plt.figure(figsize=(18,14))
    for panel in range(3):
        #######################################
        # First deal with the foreground region
        ax1 = plt.subplot(2,3,panel+1)
        ax1.axhline(0.0, color='gray', linestyle='-', linewidth=0.5)
        
        stack = stacks[panel]
        ref = refs[panel]

        # Plot the region of the integral
        x1, x2 = -1*max_growth * ch, max_growth * ch
        y1, y2 = np.min(np.array(stacks)/1e7), np.max(np.array(stacks)/1e7)#, np.max(stack/1e7)])
        ax1.vlines(x1, y1, y2, color='gray', linewidth=0.5)
        ax1.vlines(x2, y1, y2, color='gray', linewidth=0.5)
        
        # The x vector in the proper units (offset in km/s)
        x = (np.arange(200)-100) * ch

        # Plot the average spectrum of HI detections
        ax1.plot(x, msH[panel]/1e7, color='C0', linestyle='-.', linewidth=2, label='detections')

        # Plot the stack
        for i in range(len(stack)):
            if i == 0: lw = 2.0; label = 'non\ndetections'
            else: lw = 0.5; label=None
            ax1.plot(x, stack[i]/1e7, color='brown', linestyle='-', linewidth=lw, label=label)

        for i in range(len(ref)):
            if i == 0: lw = 2.0; label = 'reference'
            else: lw = 0.5; label = None
            ax1.plot(x, ref[i]/1e7, color='gray', linestyle = '--', linewidth = lw, label = label)

        # Write the number of galaxies in each spectrum
        ax1.text(0.8, 0.9, "%s" %len(OTs[panel]), color = 'brown', transform = ax1.transAxes)
        ax1.text(0.9, 0.9, "%s" %len(HTs[panel]), color = 'C0', transform = ax1.transAxes)

        if panel == 0:
            ax1.set_ylabel(r'M$_{\textrm{\sc Hi}}$ [10$^{7}$ M$_{\odot}$ (km s$^{-1}$)$^{-1}$]')
            ax1.legend(loc=2)
        
        ax1.set_xlim([-560, 560])
        ax1.set_ylim([ylim_mins[panel], np.max(msH[panel]/1e7)*1.05])
        ax1.set_xticklabels([])
        ax1.tick_params(axis='y', direction='in', which='both', right='off', top='on')
        ax1.tick_params(axis='x', direction='in', which='both', right='off', top='on', bottom='off')
        ax1.set_title(titles[panel])

        
    for panel in range(3,6):
        ############
        # Plot the stacks in Jy
        ax4 = plt.subplot(2,3,panel+1)
        ax4.axhline(0.0, color='gray', linestyle='-', linewidth=0.5)
        
        stack = stacks[panel]
        ref = refs[panel]

        # Plot the average spectrum of HI detections
        ax4.plot(x, msH[panel], color='C0', linestyle='-.', linewidth=2, label='detections foreground')

        # Plot the stack
        for i in range(len(stack)):
            if i == 0: lw = 2.0; label = 'stack foreground'
            else: lw = 0.5; label=None
            ax4.plot(x, stack[i], color='brown', linestyle='-', linewidth=lw, label=label)

        for i in range(len(ref)):
            if i == 0: lw = 2.0
            else: lw = 0.5
            ax4.plot(x, ref[i], color='gray', linestyle='--', linewidth=lw, label=label)
        if panel == 3:
            ax4.set_ylabel('Jy')
        
        ax4.set_xlim([-560, 560])
        ax4.set_ylim([ylim_mins[panel], np.max(msH[panel])*1.05])
        ax4.set_xlabel('velocity offset [km s$^{-1}$]')
        ax4.tick_params(axis='y', direction='in', which='both', right='off', top='on')
        ax4.tick_params(axis='x', direction='in', which='both', right='off', top='on')

    
    plt.tight_layout()
    plt.savefig(path+'final_analysis/undetected_stack_all.pdf')
    
def plot_individual_spectra(tbl, figname = None):
    # Plot on top of each other all individual spectra from a given stack
    plt.figure()
    ax = plt.subplot(1,1,1)
    for i in range(len(tbl)):
        sp = Table.read(path+'spectra/HI_for_stacking/'+tbl[i]['SPECOBJID']+'_han5spec.txt', format='ascii.commented_header')
        sp['mhi'] = 2.36e5 * sp['Intensity'] * cosmo.luminosity_distance(tbl[i]['Z']).value**2
        # The spectra extracted with miriad have different length, for some reason. The difference is a function
        # of z, even though there is no hint of anything in the header
        # So interpolation is needed
        # This part does the interpolation and deals with flux preservation. 
        # As I interpolate the data, the flux is elevated by the ratio of the total number of sampling points before and after the interpolation
        intel = 200 # This is how many points the interpolated spectrum will have
        rat = intel/len(sp['mhi']) #the factor by which the flux is elevated in the interpolated data
        interp_mhi_spec = np.interp(np.linspace(sp['Velocity'][-1], sp['Velocity'][0], intel), \
                                    np.linspace(sp['Velocity'][-1], sp['Velocity'][0], len(sp['mhi'])), \
                                    sp['mhi'].data) / rat
        plt.plot(interp_mhi_spec, linewidth=0.1, color='gray')
    plt.ylim([-1e8, 1e8])
    if figname:
        plt.savefig(path+'final_analysis/'+figname+'.pdf')
    else:
        plt.show()


##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
# This is the channel width in km/s
ch = 5.5

path = '/home/tazio/works/2020/a1367/'
#############################
# The table with all SDSS objects within the HI survey and cluster velocity range
Stbl = Table.read(path+'final_analysis/mpa_jhu_AGES.fits')
print(len(Stbl))
ND_tbl = Stbl[Stbl['HIdetected'] == False]
print('Stacking', len(ND_tbl), 'spectra')

# Separate the three environments
CzT = ND_tbl[(ND_tbl['Vel'] >= 4000) & (ND_tbl['Vel'] <= 9000)]
EzT = ND_tbl[(ND_tbl['Vel'] > 9000) & (ND_tbl['Vel'] <= 15500)]
FzT = ND_tbl[ND_tbl['Vel'] < 4000]
OTs = [FzT, CzT, EzT]

##########################################################
# # If you only want to stack once for a test without jackknife run this:
# HInd_msS = run_Stbl(ND_tbl)
# # Normalize for plotting
# HInd_msS /= 1e7

##########################################################
# # In order to do the jackknife stacking 
# # and write the stacked spectra to a file
# # run this:

# # The code will exit automatically after executing the next two rows
# This takes ~2h on my laptop (for 10000 iterations)
# NI = 10000
# print('The total number of galaxies in the stack will be:\n', len(CzT), 'galaxies in the cluster\n', len(EzT), 'galaxies elsewhere')
# first_time(CzT, EzT, FzT, NI, write=True)

# # If you want to make a plot of all the individual spectra run this:
# plot_individual_spectra(CzT, figname = 'cluster_individual_non_detections_spectra')
# plot_individual_spectra(EzT, figname = 'elsewhere_individual_non_detections_spectra')

##########################################################
# If the stacking has already been done
# run this to read in the results
# tell it which spectra, in percentile of the total MHI distribution
# it should return, starting with the median. Those are in a form ready to plotting
# NI is the number of jackknife resamples performed

stacks = []; refs = []
for env in ['dF', 'dC', 'dE']:
    stack, ref = later_runs(env, NI = 10000, p = [50, 2.5, 97.5])
    stacks.append(stack)
    refs.append(ref)

for env in ['dF', 'dC', 'dE']:
    stack, ref = get_Jy_runs(env, NI = 10000, p = [50, 2.5, 97.5])
    stacks.append(stack)
    refs.append(ref)

#############################
# If you want to also plot the average of the detections
# run this:
HItbl = Table.read(path+'A1367_HItable_v8.8_wt.txt', format='ascii.commented_header')
print(len(HItbl))
# Add a redshift column, makes it easier later
HItbl['Z'] = HItbl['Vel_HI']/c

HIC = HItbl[ (HItbl['Vel_HI'] > 4000) & (HItbl['Vel_HI'] < 9000) ]
HIE = HItbl[ (HItbl['Vel_HI'] > 9000) & (HItbl['Vel_HI'] <= 15500) ]
HIF = HItbl[ (HItbl['Vel_HI'] < 4000) ]

HTs = [HIF, HIC, HIE]

print('Average HI mass in Foreground, Cluster, Background:')
for HItable in HTs:
    print('{:.3f} e+8'.format((np.sum(10**HItable['MHI'])/len(HItable))/1e8))

msH = []
for U in [True, False]:
    for HT in HTs:
        HIstack = []; weights = []
        for i in range(len(HT)):
            spname = path+'spectra/han_spectra/'+HT[i]['Name']+'spec.txt'
            interp_mhi_spec = get_spectrum(spname, i, HT, to_mass = U)
            HIstack.append(interp_mhi_spec)
            weights.append(  1 / iqr(list(interp_mhi_spec)[:40] + list(interp_mhi_spec)[-40:]) )
        msH.append( np.average(HIstack, axis=0, weights = weights) )
        

#############################
# If you want to make a plot run this
make_a_plot_all()

