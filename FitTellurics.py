import sys
import os
import gc

import numpy as np
from astropy.io import fits
from astropy import units, constants

import FittingUtilities
import TelluricFitter
import DataStructures
import Units
import HelperFunctions
import GetAtmosphere
import logging

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

badregions = [[]]

FindOrderNums = HelperFunctions.FindOrderNums

def FitAll():
    fileList = []
    start = 0
    end = 999
    makenew = True
    edit_atmosphere = False
    humidity_low = 1.0
    humidity_high = 99.0
    for arg in sys.argv[1:]:
        if "-atmos" in arg:
            edit_atmosphere = True
        elif "-hlow" in arg:
            humidity_low = float(arg.split("=")[1])
        elif "-hhigh" in arg:
            humidity_high = float(arg.split("=")[1])
        elif arg != "-e":
            fileList.append(arg)


    # Set up the dictionary for order-by-order fitting
    fitdict = {1: ["h2o"],
               2: ["h2o", "co2"],
               3: ["h2o", "co2"],
               4: ["h2o", "co2"],
               5: ["h2o", "co2"],
               6: ["h2o", "co2"],
               7: ["h2o", "co2"],
               8: ["h2o", "co2"],
               9: ["h2o", "co2"],
               10: ["h2o", "co2"],
               11: ["h2o", "co2"],
               12: ["h2o", "co2", "ch4"],
               13: ["h2o", "co2", "ch4"],
               14: ["h2o", "co2", "ch4"],
               15: ["h2o", "ch4"],
               16: ["h2o", "ch4"],
               17: ["h2o", "ch4"],
               18: ["h2o", "ch4"],
               19: ["h2o", "ch4"],
               20: ["h2o", "ch4"],
               21: ["h2o", "ch4"],
               22: ["h2o", "ch4"],
               23: ["h2o", "ch4"],
               24: ["h2o", "co2"],
               25: ["h2o", "co2"],
               26: ["h2o", "co2"],
               27: ["h2o", "co2"],
               28: ["h2o", "co2"],
               29: ["h2o", "co2"],
               30: ["h2o", "co2", "n2o"],
               31: ["h2o", "n2o"],
               32: ["h2o", "n2o"],
               33: ["h2o", "n2o"],
               34: ["h2o", "ch4"],
               35: ["h2o", "ch4"],
               36: ["h2o", "ch4"],
               37: ["h2o", "ch4"],
               38: ["h2o", "ch4", "co"],
               39: ["h2o", "ch4", "co"],
               40: ["h2o", "ch4", "co"],
               41: ["h2o", "ch4"],
               42: ["h2o", "ch4"],
               43: ["h2o", "ch4"],
               44: ["h2o", "ch4"],
               45: ["h2o", "ch4"]}


    # START LOOPING OVER INPUT FILES
    for fname in fileList:
        # Initialize fitter
        fitter = TelluricFitter.TelluricFitter()
        fitter.SetObservatory("McDonald")

        print "Fitting file {0:s}".format(fname)
        # Make sure this file is an object file
        header = fits.getheader(fname)

        logfile = open(u"fitlog_{0:s}.txt".format(fname.split(".fits")[0]), "a")
        logfile.write(u"Fitting file {0:s}\n".format(fname))
        name = fname.split(".fits")[0]
        outfilename = "Corrected_{:s}-0.fits".format(name)

        # Read file
        all_orders = HelperFunctions.ReadExtensionFits(fname)

        # New data has some more orders that are hopelessely contaminated by tellurics.
        # Remove them...
        orders = []
        for order in all_orders:
            mwl = np.median(order.x)
            if (mwl > 1480 and mwl < 1810) or (mwl > 1935 and mwl < 2470):
                orders.append(order.copy())

        angle = float(header["ZD"])
        resolution = 40000.0
        humidity = header['HUMIDITY']
        T_fahrenheit = header['AIRTEMP']
        pressure = header['BARPRESS'] * Units.hPa / Units.inch_Hg
        temperature = (T_fahrenheit - 32.0) * 5.0 / 9.0 + 273.15
        o2 = 2.12e5
        co = 0.14
        co2 = 368.5
        n2o = 0.32
        ch4 = 1.8

        if edit_atmosphere:
            filenames = [f for f in os.listdir("./") if "GDAS" in f]
            height, Pres, Temp, h2o = GetAtmosphere.GetProfile(filenames, header['date-obs'].split("T")[0],
                                                               header['ut'])

            fitter.EditAtmosphereProfile("Temperature", height, Temp)
            fitter.EditAtmosphereProfile("Pressure", height, Pres)
            fitter.EditAtmosphereProfile("H2O", height, h2o)


        # Adjust fitter values
        fitter.AdjustValue({"angle": angle,
                            "pressure": pressure,
                            "resolution": resolution,
                            "temperature": temperature,
                            "h2o": humidity,
                            "o2": o2,
                            "co": co,
                            "co2": co2,
                            "n2o": n2o,
                            "ch4": ch4})
        fitter.SetBounds({"h2o": [humidity_low, humidity_high],
                          "temperature": [temperature - 10, temperature + 10],
                          "pressure": [pressure - 30, pressure + 100],
                          "ch4": [0.0, 10.0],
                          "co2": [0.0, 1000.0],
                          "n2o": [0.0, 10.0],
                          "co": [0.0, 10.0],
                          "resolution": [30000, 50000]})


        # Ignore some regions (currently nothing)
        fitter.IgnoreRegions(badregions)

        fitter.continuum_fit_order = 3
        fitter.resolution_fit_mode = 'SVD'

        # Start fitting, order by order
        for i, order in enumerate(orders):
            fitter.AdjustValue({"wavestart": order.x[0] - 5.0,
                                "waveend": order.x[-1] + 5.0,
                                "temperature": temperature,
                                "h2o": humidity,
                                "o2": o2,
                                "co2": co2,
                                "n2o": n2o,
                                "ch4": ch4})
            print "\n*****  Fitting order {}  ********\n".format(i + 1)
            if i + 1 in [1, 22, 23, 24, 25, 44]:
                fitpars = [fitter.const_pars[j] for j in range(len(fitter.parnames)) if fitter.fitting[j]]
                fitter.ImportData(order.copy())
                model = fitter.GenerateModel(fitpars, separate_source=False, nofit=True)
                model = FittingUtilities.RebinData(model, order.x)
                primary = model.copy()
                primary.y = np.ones(primary.size())
                #primary, model = fitter.GenerateModel(fitpars,
                #                                      separate_primary=True,
                #                                      return_resolution=False)

            else:
                # Figure out which molecules to fit
                for molec in ["h2o", "co2", "co", "ch4", "n2o"]:
                    if molec in fitdict[i + 1]:
                        if molec == "h2o":
                            fitter.FitVariable({"h2o": humidity})
                        else:
                            fitter.FitVariable({molec: eval(molec)})
                fitter.DisplayVariables()

                order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=3, lowreject=9, highreject=10)
                
                # Put the data in a constant wavelength grid
                original_x = order.x.copy()
                dx = min(order.x[1:] - order.x[:-1])
                xgrid = np.arange(order.x[0], order.x[-1]+dx/2.0, dx)
                rebinned_order = FittingUtilities.RebinData(order, xgrid)

                primary, model = fitter.Fit(data=rebinned_order.copy(),
                                            resolution_fit_mode="SVD",
                                            fit_source=True,
                                            source_args=[],
                                            source_kwargs=dict(window_size=101, order=4),
                                            continuum_fit_order=3,
                                            return_resolution=False,
                                            adjust_wave="data",
                                            wavelength_fit_order=5,
                                            air_wave=False)
                humidity = fitter.GetValue("h2o")
                temperature = fitter.GetValue("temperature")
                ch4 = fitter.GetValue("ch4")
                co2 = fitter.GetValue("co2")
                co = fitter.GetValue("co")
                n2o = fitter.GetValue("n2o")
                rebinned_order = fitter.data
                #order = fitter.data

                # Put everything back into the original wavelength grid
                order = FittingUtilities.RebinData(rebinned_order, original_x)
                primary = FittingUtilities.RebinData(primary, original_x)
                model = FittingUtilities.RebinData(model, original_x)

            # Set up data structures for OutputFitsFile
            columns = {"wavelength": order.x,
                       "flux": order.y,
                       "continuum": order.cont,
                       "error": order.err,
                       "model": model.y,
                       "primary": primary.y}

            header = [["HUMIDITY", fitter.GetValue("h2o")],
                      ["TEMPERATURE", fitter.GetValue("temperature")],
                      ["CH4", fitter.GetValue("ch4")],
                      ["CO2", fitter.GetValue("co2")],
                      ["CO", fitter.GetValue("co")],
                      ["N2O", fitter.GetValue("n2o")]]

            if i == 0:
                HelperFunctions.OutputFitsFileExtensions(columns,
                                                         fname,
                                                         outfilename,
                                                         mode="new",
                                                         headers_info=[header, ])

            else:
                HelperFunctions.OutputFitsFileExtensions(columns,
                                                         outfilename,
                                                         outfilename,
                                                         mode="append",
                                                         headers_info=[header, ])
        del fitter
        gc.collect()


if __name__ == "__main__":
    FitAll()
