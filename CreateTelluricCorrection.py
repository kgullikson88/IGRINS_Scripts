__author__ = 'kgulliks'

"""
The purpose of this is to generate a telluric correction for a file, given
best-fit models to A0 standards before and after it. It takes the following
steps:

  1: Reads in the model parameters from before and after
  2: Finds the time-weighted mean of the atmospheric parameters
  3: Adjust the zenith angle and pressure to be correct
  4: Generate a model with the new parameters
  5: Remove the model from the data
"""

import os
import sys
import logging

import pandas
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.time import Time

import FittingUtilities
import TelluricFitter
import HelperFunctions
import Units
import GetAtmosphere


def ReadModels(beforefile, afterfile, T='TEMPERATURE', RH='HUMIDITY', ch4='CH4', co2='CO2', co='CO', n2o='N2O'):
    """
    Function to read in model parameters. Returns a pandas dataframe
    :param beforefile: The fits file with the model parameters in its extension headers (before the observation)
    :param afterfile: The fits file with the model parameters in its extension headers (after the observation)
    :param T: The fits header key with the temperature
    :param RH: The fits header key with the relative humidity
    :param ch4: The fits header key with the ch4 abundance
    :param co2: The fits header key with the co2 abundance
    :param co: The fits header key with the co abundance
    :param n2o: The fits header key with the n2o abundance
    :return: pandas.DataFrames that hold the relevant info for both models
    """
    before = fits.open(beforefile)
    after = fits.open(afterfile)

    pars = []
    for extnum in range(1, len(before)):
        header = before[extnum].header
        pars.append((header[T],
                     header[RH],
                     header[ch4],
                     header[co2],
                     header[co],
                     header[n2o]))
    before = pandas.DataFrame(np.array(pars), columns=(T, RH, ch4, co2, co, n2o))

    pars = []
    for extnum in range(1, len(after)):
        header = after[extnum].header
        pars.append((header[T],
                     header[RH],
                     header[ch4],
                     header[co2],
                     header[co],
                     header[n2o]))
    after = pandas.DataFrame(np.array(pars), columns=(T, RH, ch4, co2, co, n2o))

    return before, after


def GetParameters(datafile):
    """
    Determines the atmospheric parameters to use by interpolating between A0 spectra
    :param datafile: The filename of the data to telluric-correct
    :return:
    """
    # Find the julian date of the observation
    header = fits.getheader(datafile)
    t = Time(header['DATE-OBS'], format='isot', scale='utc')
    jd = t.jd

    # Look in the current directory for files with the appropriate naming convention
    template_files = [f for f in os.listdir("./") if f.startswith("Corrected") and f.endswith("-0.fits")]
    obstimes = []
    for fname in template_files:
        header = fits.getheader(fname)
        t = Time(header['DATE-OBS'], format='isot', scale='utc')
        obstimes.append(t.jd)

    # Find the two files surrounding this one
    before_idx, after_idx = HelperFunctions.GetSurrounding(obstimes, jd, return_index=True)
    if before_idx > after_idx:
        before_idx, after_idx = after_idx, before_idx

    # Get the model parameters for the files before and after this one
    beforefile = template_files[before_idx]
    afterfile = template_files[after_idx]
    before_pars, after_pars = ReadModels(beforefile, afterfile)

    if beforefile == afterfile:
        return before_pars

    # Interpolate all the parameters to the time of observation
    m = (after_pars - before_pars) / (obstimes[after_idx] - obstimes[before_idx])
    pars = m * (jd - obstimes[before_idx]) + before_pars

    return pars


def CorrectData(filename, pars, plot=False, edit_atmosphere=True):
    """
    Generates models for each order using the parameters given in the pars structure
    :param filename: The filename of the data to correct
    :param pars: A pandas.DataFrame containing the atmospheric parameters for each echelle order
    :param plot:  A boolean flag determining whether or not to plot some intermediate results (useful for debugging)
    :return: A list of corrected orders
    """

    orders = HelperFunctions.ReadExtensionFits(filename)
    header = fits.getheader(filename)
    zd = header['ZD']
    pressure = header['BARPRESS'] * Units.hPa / Units.inch_Hg

    # Initialize the fitter
    fitter = TelluricFitter.TelluricFitter()
    fitter.continuum_fit_order = 3
    fitter.resolution_fit_mode = "SVD"
    fitter.adjust_wave = "model"
    fitter.wavelength_fit_order = 4
    fitter.fit_source = True
    fitter.SetBounds({"resolution": [30000, 50000]})

    if edit_atmosphere:
        filenames = [f for f in os.listdir("./") if "GDAS" in f]
        height, Pres, Temp, h2o = GetAtmosphere.GetProfile(filenames, header['date-obs'].split("T")[0],
                                                           header['ut'])

        fitter.EditAtmosphereProfile("Temperature", height, Temp)
        fitter.EditAtmosphereProfile("Pressure", height, Pres)
        fitter.EditAtmosphereProfile("H2O", height, h2o)

    if plot:
        fig = plt.figure(1)
        axtop = fig.add_subplot(211)
        axbot = fig.add_subplot(212, sharex=axtop)

    corrected_orders = []
    for i, order in enumerate(orders):
        logging.info('Fitting order {}/{} in file {}'.format(i + 1, len(orders), filename))
        fitter.ImportData(order)
        fitter.AdjustValue({"wavestart": order.x[0] - 5,
                            "waveend": order.x[-1] + 5,
                            "angle": zd,
                            "pressure": pressure,
                            "h2o": pars.iloc[i]['HUMIDITY'],
                            "temperature": pars.iloc[i]['TEMPERATURE'],
                            "ch4": pars.iloc[i]['CH4'],
                            "co2": pars.iloc[i]['CO2'],
                            "co": pars.iloc[i]['CO'],
                            "n2o": pars.iloc[i]['N2O']})
        if i + 1 in [1, 22, 23, 24, 25, 44]:
            fitpars = [fitter.const_pars[j] for j in range(len(fitter.parnames)) if fitter.fitting[j]]
            model = fitter.GenerateModel(fitpars, separate_source=False, nofit=True)
            model = FittingUtilities.RebinData(model, order.x)
            primary = model.copy()
            primary.y = np.ones(primary.size())
        else:
            model = fitter.GenerateModel([])

        # Just set equal to the data if there is a huge amount of absorption
        badindices = np.where(np.logical_or(order.y <= 0, model.y < 0.05))[0]
        model.y[badindices] = order.y[badindices] / order.cont[badindices]
        model.y[model.y < 1e-5] = 1e-5

        if plot:
            axtop.plot(order.x, order.y / order.cont, 'k-', alpha=0.4)
            axtop.plot(model.x, model.y, 'r-', alpha=0.6)
            axbot.plot(order.x, order.y / (order.cont * model.y), 'k-', alpha=0.4)

        order.y /= model.y
        corrected_orders.append(order.copy())

    if plot:
        axtop.set_ylim((-0.1, 2.0))
        axbot.set_ylim((-0.1, 2.0))
        plt.show()

    return corrected_orders


if __name__ == "__main__":
    fileList = []
    plot = False
    atmos = True
    for arg in sys.argv[1:]:
        if 1:
            fileList.append(arg)

    for fname in fileList:
        outfilename = "{:s}_telfit.fits".format(fname[:-5])
        print outfilename
        model_pars = GetParameters(fname)
        corrected_orders = CorrectData(fname, model_pars, plot=plot, edit_atmosphere=atmos)

        column_list = []
        for i, data in enumerate(corrected_orders):
            # Set up data structures for OutputFitsFile
            columns = {"wavelength": data.x,
                       "flux": data.y,
                       "continuum": data.cont,
                       "error": data.err}
            column_list.append(columns)
        HelperFunctions.OutputFitsFileExtensions(column_list, fname, outfilename, mode="new")
