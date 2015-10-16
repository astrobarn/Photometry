# -*- coding: utf-8 -*-
"""
Created on Oct 2015

Script to extract and arrange the data from the
forced photometry files to a digestible format
for the template maker.

@author: Seméli Papadogiannakis
"""

'''
This script converts the PTFIDE files into magnitudes and upperlimits,
making a number of corrections. These are based on based on the document
"http://web.ipac.caltech.edu/staff/fmasci/home/miscscience/forcedphot.pdf",
which is to be used in the PTFIDE format files produced from August 2015
and onwards, i.e. from version "forcepsffitdiff.pl v3.0".

The corrections included are:
- A correction for residual offset in the historical BASELINE. This is done
  so that any potential flux from the transient in the reference image is
  corrected for and be able to get the most accurate photometry. This is
  only important for transients. By historical we define 30 epochs earlier
  than the peak (here defined as the maximum flux measurement).
- UNCERTAINTY VALIDATION: two methods
- COURSE CHECK
- QUALITY CHECK


Corrections not included (Systematics):
- Incorrect PSF-template estimation
- Photometric zero-point calibrations
- Astrometric calibrations(determining PSF-placement)
- Error in supplied target position.

Disclaimer: This script is optimised for Ia supernova science and should
be changed for other science goals, especially for non-transients.
'''

import glob

import numpy as np
from astropy.io import ascii
import pylab as plt
import os

os.system('mkdir forced_mags 2>/dev/null')
os.system('mkdir forced_fluxes 2>/dev/null')

# Getting the coodinates-name table (includes
# all transients).
coord_file = 'coord_info.txt'
coords = ascii.read(coord_file)
name_info = coords['col1']
ra_info = coords['col2']
dec_info = coords['col3']

data_list = glob.glob('forced_new/*f2*.out')

for data in data_list:
    name = data.split('/')[1]
    try:
        data = ascii.read(data)
    except ascii.core.InconsistentTableError:
        print name
        os.system('leafpad ' + data)
        continue
    ra = data['ra']
    dec = data['dec']

    assert ra.ptp() < 0.01
    assert dec.ptp() < 0.01

    ra = np.mean(ra)
    dec = np.mean(dec)

    # come up with some clever way to loop through the possible coords
    # to find the one that matches best
    distances = np.square(ra_info - ra) + np.square(dec_info - dec)
    position = distances.argmin()
    sn_name = name_info[position]
    #print sn_name
    #print distances.min()

    flux = data['flux']
    sigflux = data['sigflux']
    valid = np.logical_and(flux < 99999999, sigflux > 0)
    sharp = data['sharp']
    valid = np.logical_and(valid, np.abs(sharp) < 30)
    fluxap = data['fluxap']
    valid = np.logical_and(valid, np.abs(fluxap) < 1000)

    if not any(valid):
        print name
        continue

    flux = flux[valid]
    sigflux = sigflux[valid]
    fluxap = data['fluxap'][valid]
    sigfluxap = data['sigfluxap'][valid]

    snr = data['snr'][valid]
    chi = data['chi'][valid]
    sharp = data['sharp'][valid]

    zero_point = data['zpmag'][valid]
    zp_rms = data['zprms'][valid]

    # 27.0 are fill in values, replace with an estimation from other points
    valid_zero_point = np.logical_not(np.logical_and(zero_point == 27.0, zp_rms == 0))
    zero_point[~valid_zero_point] = np.median(zero_point[valid_zero_point])
    zp_rms[~valid_zero_point] = np.median(zp_rms[valid_zero_point])

    # Then can get MJD
    jd = data['MJD'][valid] + 2400000.5
    del data
    del valid_zero_point

    # Now we have to perform the BASELINE correction
    flux_sorted = np.sort(flux)

    index_max = np.where(flux == flux_sorted[-3])[0][0]
    jd_max = jd[index_max]

    jd_end = jd_max - 50
    jd_start = jd_max - 2650

    jd_start_index, jd_end_index = np.searchsorted(jd, [jd_start, jd_end])

    if jd_end_index == jd_start_index:
        jd_end = jd_max - 30
        jd_start = jd_max - 150

        jd_start_index, jd_end_index = np.searchsorted(jd, [jd_start, jd_end])

        short_historical_flag = True
        if jd_end_index == jd_start_index:
            continue
    else:
        short_historical_flag = False
    '''
    Get the baseline!
    '''
    baseline = np.median(flux[jd_start_index:jd_end_index])

    # Make a subplot with and without baseline correction

    '''
    Make a correction to all errors using method 1-3 in
    the pdf.
    '''
    # Method 1
    sigflux_1 = sigflux * np.median(chi[jd_start_index:jd_end_index])

    # Method 2
    s_factor = np.std(flux[jd_start_index:jd_end_index]) / np.median(sigflux[jd_start_index:jd_end_index])
    sigflux_2 = sigflux * s_factor

    #plt.figure()
    #plt.errorbar(jd, flux, yerr=sigflux, fmt='+')
    #plt.axhline(baseline, color='g', linestyle=':', label='Baseline corrected')
    #plt.axhline(0, color='k'. label='Original zero flux')
    #plt.legend()

    #if short_historical_flag:
     #   plt.axvline(jd_end, color='r')
     #   plt.axvline(jd_start, color='r')
    #else:
    #    plt.axvline(jd_end, color='m')
    #    plt.axvline(jd_start, color='m')

    #plt.axvline(jd_max, color='k', linestyle='--')

    plt.figure('FluxPSF_fluxapdiff')
    plt.hist(flux - fluxap, histtype='step', normed=True)
    plt.xlabel('Difference between Aperture and PSF photometry')
    #plt.colorbar()

    plt.figure('histograms sharp')
    plt.hist(sharp, color='k', histtype='step', normed=True)
    plt.xlabel('Sharpness')

    flux += baseline

    # Compute the magnitudes
    snr = flux / sigflux_1
    true_detection = snr > 3.

    mags = zero_point - 2.5 * np.log10(flux)
    e_mags = 1.0857 * sigflux_1 / flux
    lim_mags = zero_point - 2.5 * np.log10(3 * sigflux_1)

    mags[~true_detection] = e_mags[~true_detection] = np.nan
    lim_mags[true_detection] = np.nan

    #plt.figure()
    #plt.errorbar(jd, -mags, yerr=e_mags, fmt='+')
    #plt.scatter(jd, -lim_mags, marker='v')

    with open('forced_mags/'+sn_name+'_forcedmags','w') as outputfile:
        print >> outputfile, '#jd Mag mag_err lim_mag baseline_flag'
        for i in xrange(len(mags)):
            print >> outputfile,'{} {:.3f} {:.3f} {:.3f} {}'.format(jd[i], mags[i], e_mags[i], lim_mags[i], short_historical_flag)

    with open('forced_fluxes/'+sn_name+'_forcedflux', 'w') as outputfile:
        print >> outputfile, '#jd Flux flux_err1 flux_err2 baseline_flag'
        for i in xrange(len(flux)):
            print >> outputfile,'{} {:.3f} {:.3f} {:.3f} {}'.format(jd[i], flux[i], sigflux_1[i], sigflux_2[i], short_historical_flag)

plt.show()
