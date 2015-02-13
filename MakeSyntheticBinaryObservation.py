"""
This script takes several of my program stars, and several late-type stars, and creates
several synthetic spectra by combining them. Use these data to test my CCF technique and
the parameters it derives!

Usage:
python MakeSyntheticBinaryObservation.py myfile1 myfile2 ... -late latefile1 latefile2...

myfile: a filename for my data. Should be a KG* file
latefile: a filename for a late-type star with known properties, taken with the same instrument as above

"""

import sys
import os
import logging
import re

from scipy.interpolate import InterpolatedUnivariateSpline as spline
from astropy.io import fits
import numpy as np
from astroquery.simbad import Simbad
import astropy.units as u

import HelperFunctions
import SpectralTypeRelations
import pySIMBAD
from PlotBlackbodies import Planck


def parse_input(argv):
    early_files = []
    late_files = []
    late = False  # Flag for when to start filling the late file list
    for arg in argv:
        if arg.lower() == '-late':
            late = True
        elif late:
            late_files.append(arg)
        else:
            early_files.append(arg)
    if len(early_files) < 1 or len(late_files) < 1:
        logging.warning('Need to give at least one early-type (mine) and one late-type spectrum, '
                        'separated by the "-late" flag')

    return early_files, late_files


def combine(early_filename, late_filename, increase_scale=False):
    """
    This function does the bulk of the work.
    :param early_filename: the filename of the early-type star
    :param late_filename: the filename of the late-type star. Should be from the same instrument and setup as my data!
    :return: A list of orders, after combination
    """
    # First, we will classify the files by their fits header and some simbad lookups
    early_dict = classify_file(early_filename, astroquery=True)
    late_dict = classify_file(late_filename, astroquery=True)
    print early_dict
    print late_dict

    # Now, we need to work out how to scale the data as if they were in the same system
    # This ASSUMES that all atmospheric affects are the same for both observations, which is not really true!
    scale = (early_dict['exptime'] / late_dict['exptime']) * (early_dict['plx'] / late_dict['plx']) ** 2
    print scale

    # Get some information to use for the predicted flux ratio (done later)
    MS = SpectralTypeRelations.MainSequence()
    T1 = MS.Interpolate(MS.Temperature, early_dict['SpT'])
    R1 = MS.Interpolate(MS.Radius, early_dict['SpT'])
    T2 = MS.Interpolate(MS.Temperature, late_dict['SpT'])
    R2 = MS.Interpolate(MS.Radius, late_dict['SpT'])

    # Read in the orders for both files
    # early_orders = ConvertToExtensions.read_orders(early_filename, blazefile)
    early_orders = HelperFunctions.ReadExtensionFits(early_filename)
    late_orders = HelperFunctions.ReadExtensionFits(late_filename)

    # Before combining, we need to adjust the scale factor for atmospheric effects (clouds, poor seeing, etc)
    scale_adjust = []
    for early in early_orders:
        num = HelperFunctions.FindOrderNums(late_orders, [np.median(early.x)])[0]
        late = late_orders[num]
        fluxratio_obs = late.cont.mean() * scale / early.cont.mean()
        x = early.x * u.nm.to(u.cm)
        fluxratio_pred = Planck(x, T2) / Planck(x, T1) * (R2 / R1) ** 2
        scale_adjust.append(np.median(fluxratio_pred / fluxratio_obs))
    scale_adjust = np.median(scale_adjust)
    print '\t', scale_adjust
    scale *= scale_adjust
    if increase_scale:
        scale *= 10.0

    # Finally, combine:
    order_list = []
    for order in early_orders:
        num = HelperFunctions.FindOrderNums(late_orders, [np.median(order.x)])[0]
        late = late_orders[num]
        late.y *= scale
        combined = add_order(order, late)
        order_list.append(combined)
        # plt.plot(order.x, order.y)
        # plt.plot(combined.x, combined.y)
    # plt.show()

    return order_list, early_dict, late_dict


def add_order(early, late):
    """
    Adds two xypoint instances. This is only really necessary because the wavelengths
    and/or order lengths are not always the same...
    :param early: the early type star
    :param late: the late type star
    :return: xypoint instance of the combined 'observation'
    """
    left = np.searchsorted(early.x, late.x[0])
    right = np.searchsorted(early.x, late.x[-1])
    late_fcn = spline(late.x, late.y, k=1)
    combined = early[left:right].copy()
    combined.y += late_fcn(combined.x)
    return combined


def classify_file(filename, astroquery=True):
    """
    This function uses the fits header information and the Simbad database to classify the object
    :param filename: The filename of the observation to be classified
    :return:
    """
    # Read in the header and get the object name
    header = fits.getheader(filename)
    object = header['object']
    print object

    # Default values if I can't get it any other way
    plx = 30.0

    MS = SpectralTypeRelations.MainSequence()

    # Make a Simbad object
    if astroquery:
        sim = Simbad()
        sim.add_votable_fields('plx', 'sp', 'flux(V)')
        data = sim.query_object(object)
        spt_full = data['SP_TYPE'].item()
        # spt = spt_full[0] + re.search(r'\d*\.?\d*', spt_full[1:]).group()
        spt = spt_full
        if data['PLX_VALUE'].mask:
            if not data['FLUX_V'].mask:
                # Fall back to photometric parallax
                Vmag_obs = data['FLUX_V'].item()
                Vmag_abs = MS.GetAbsoluteMagnitude(spt, color='V')
                plx = 10 ** (3.0 + (Vmag_abs - Vmag_obs - 5.0) / 5.0)
        else:
            plx = data['PLX_VALUE'].item()
    else:
        link = pySIMBAD.buildLink(object)
        data = pySIMBAD.simbad(link)
        plx = data.Parallax()
        spt_full = data.SpectralType().split()[0]
        spt = spt_full[0] + re.search(r'\d*\.?\d*', spt_full[1:]).group()

    if 'exptime' not in header.keys():
        header['exptime'] = 30.0
    d = {'Object': object,
         'plx': plx,
         'SpT': spt,
         'exptime': header['exptime']}
    return d


if __name__ == '__main__':
    scale = False
    early, late = parse_input(sys.argv[1:])

    # Add each late file to all of the early-type files
    HelperFunctions.ensure_dir('GeneratedObservations')
    for late_file in late:
        for early_file in early:
            outfilename = 'GeneratedObservations/{}_{}.fits'.format(early_file.split('/')[-1].split(
                '.fits')[0], late_file.split('/')[-1].split('.fits')[0])
            if scale:
                outfilename = outfilename.replace('.fits', '_scalex10.fits')
            if outfilename.split('/')[-1] in os.listdir('GeneratedObservations/'):
                print "File already generated. Skipping {}".format(outfilename)
                continue

            total, early_dict, late_dict = combine(early_file, late_file, increase_scale=scale)

            # Prepare for output
            column_list = []
            for order in total:
                column = {'wavelength': order.x,
                          'flux': order.y,
                          'continuum': order.cont,
                          'error': order.err}
                column_list.append(column)
            newheader = {'object1': early_dict['Object'],
                         'object2': late_dict['Object'],
                         'SpT1': early_dict['SpT'],
                         'SpT2': late_dict['SpT'],
                         'file1': early_file,
                         'file2': late_file}
            print 'Outputting to {}'.format(outfilename)
            HelperFunctions.OutputFitsFileExtensions(column_list, early_file, outfilename,
                                                     mode='new', primary_header=newheader)

