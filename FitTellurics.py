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


def FindOrderNums(orders, wavelengths):
    """
    Given a list of xypoint orders and
    another list of wavelengths, this
    finds the order numbers with the
    requested wavelengths
    """
    nums = []
    for wave in wavelengths:
        for i, order in enumerate(orders):
            if order.x[0] < wave and order.x[-1] > wave:
                nums.append(i)
                break
    return nums


def EstimateModel():
    # Initialize fitter
    fitter = TelluricFitter.TelluricFitter()
    fitter.SetObservatory("McDonald")

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


    # START LOOPING OVER INPUT FILES
    for fname in fileList:
        print "Fitting file {0:s}".format(fname)
        # Make sure this file is an object file
        header = fits.getheader(fname)

        logfile = open(u"fitlog_{0:s}.txt".format(fname.split(".fits")[0]), "a")
        logfile.write(u"Fitting file {0:s}\n".format(fname))
        name = fname.split(".fits")[0]
        outfilename = "Corrected_{0:s}.fits".format(name)
        exists = False

        # Read file
        orders = HelperFunctions.ReadExtensionFits(fname)

        angle = float(header["ZD"])
        resolution = 40000.0
        humidity = header['HUMIDITY']
        T_fahrenheit = header['AIRTEMP']
        pressure = header['BARPRESS'] * Units.hPa / Units.inch_Hg
        temperature = (T_fahrenheit - 32.0) * 5.0 / 9.0 + 273.15

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
                            "o2": 2.12e5,
                            "co": 0.14,
                            "co2": 368.5,
                            "n2o": 0.32})
        fitter.FitVariable({"h2o": humidity,
                            "ch4": 1.8,
                            "temperature": temperature})
        fitter.SetBounds({"h2o": [humidity_low, humidity_high],
                          "temperature": [temperature - 10, temperature + 10],
                          "pressure": [pressure - 30, pressure + 100],
                          "ch4": [1.0, 20],
                          "co2": [100, 1000],
                          "n2o": [0.05, 1.0],
                          "co": [0.01, 0.9],
                          "resolution": [30000, 50000]})

        # Ignore some regions (currently nothing)
        fitter.IgnoreRegions(badregions)

        fitter.continuum_fit_order = 3

        # ############################################################
        #     Now, we start fitting molecules!
        #############################################################


        # Determine the H2O and CH4 abundances, as well as the temperature
        # We will weight the final values by their uncertainties in the fit
        resolution = []
        h2o = []
        T = []
        methane = []
        P = []
        waveshifts = []
        wave0 = []
        chisquared = []
        fitter.DisplayVariables()
        for i in FindOrderNums(orders, [1658, 1674, 1689, 1705]):
            print "\n***************************\nFitting order %i: " % (i)
            order = orders[i]
            fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                                "waveend": order.x[-1] + 20.0})
            order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=9, lowreject=2, highreject=10)
            primary = DataStructures.xypoint(x=order.x, y=np.ones(order.x.size))
            primary, model, R = fitter.Fit(data=order.copy(),
                                           resolution_fit_mode="SVD",
                                           fit_source=True,
                                           return_resolution=True,
                                           adjust_wave="data",
                                           wavelength_fit_order=4)
            resolution.append(R)
            waveshifts.append(fitter.shift)
            wave0.append(fitter.data.x.mean())
            h2o.append(fitter.GetValue("h2o"))
            methane.append(fitter.GetValue("ch4"))
            T.append(fitter.GetValue("temperature"))
            P.append(fitter.GetValue("pressure"))

            chisquared.append((1.0 - min(model.y)) / fitter.chisq_vals[-1])


        # Determine the average values (weight by chi-squared)
        humidity = np.sum(np.array(h2o) * np.array(chisquared)) / np.sum(chisquared)
        ch4 = np.sum(np.array(methane) * np.array(chisquared)) / np.sum(chisquared)
        temperature = np.sum(np.array(T) * np.array(chisquared)) / np.sum(chisquared)
        pressure = np.sum(np.array(P) * np.array(chisquared)) / np.sum(chisquared)
        logfile.write("Humidity/Methane/Temperature/Pressure values and their weights values:\n")
        for h, m, t, p, c in zip(h2o, methane, T, P, chisquared):
            logfile.write(u"{0:g}\t{1:g}\t{2:g}\t{3:g}\t{4:g}\n".format(h, m, t, p, c))
        logfile.write(u"Best humidity = {0:.4f}\n".format(humidity))
        logfile.write(u"Best temperature = {0:.4f}\n".format(temperature))
        logfile.write(u"Best Methane mixing ratio = {0:.4f}\n".format(ch4))
        logfile.write("Best Pressure = {0:.4f}\n".format(pressure))
        fitter.AdjustValue({"h2o": humidity,
                            "ch4": ch4,
                            "temperature": temperature,
                            "pressure": pressure})


        # ---------------------------------------------------------------------------------


        # Now, determine the CO abundance
        fitter.FitVariable({"co": 0.14})
        carbon_monoxide = []
        numskip = len(chisquared)
        for i in FindOrderNums(orders, [2310, 2340, 2365]):
            print "\n***************************\nFitting order %i: " % (i)
            order = orders[i]
            fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                                "waveend": order.x[-1] + 20.0})
            order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=9, lowreject=9, highreject=10)
            primary = DataStructures.xypoint(x=order.x, y=np.ones(order.x.size))
            primary, model, R = fitter.Fit(data=order.copy(),
                                           resolution_fit_mode="SVD",
                                           fit_source=True,
                                           return_resolution=True,
                                           adjust_wave="model",
                                           wavelength_fit_order=4)
            resolution.append(R)
            waveshifts.append(fitter.shift)
            wave0.append(fitter.data.x.mean())
            carbon_monoxide.append(fitter.GetValue("co"))
            chisquared.append((1.0 - min(model.y)) / fitter.chisq_vals[-1])

        # Determine the average CO mixing ratio
        co = np.sum(np.array(carbon_monoxide) * np.array(chisquared[numskip:])) / np.sum(chisquared[numskip:])
        logfile.write("CO values and their weights\n")
        for c, w in zip(carbon_monoxide, chisquared[numskip:]):
            logfile.write("{0:g}\t{1:g}\n".format(c, w))
        logfile.write("Best CO Mixing Ratio = {0:g}\n".format(co))
        fitter.AdjustValue({"co": co})


        # --------------------------------------------------------------------------------


        # Next up: CO2 mixing ratio
        fitter.FitVariable({"co2": 368.5})
        carbon_dioxide = []
        numskip = len(chisquared)
        for i in FindOrderNums(orders, [2310, 2340, 2365]):
            print "\n***************************\nFitting order %i: " % (i)
            order = orders[i]
            fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                                "waveend": order.x[-1] + 20.0})
            order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=9, lowreject=9, highreject=10)
            primary = DataStructures.xypoint(x=order.x, y=np.ones(order.x.size))
            primary, model, R = fitter.Fit(data=order.copy(),
                                           resolution_fit_mode="SVD",
                                           fit_source=True,
                                           return_resolution=True,
                                           adjust_wave="data",
                                           wavelength_fit_order=4)
            resolution.append(R)
            waveshifts.append(fitter.shift)
            wave0.append(fitter.data.x.mean())
            carbon_dioxide.append(fitter.GetValue("co2"))
            chisquared.append((1.0 - min(model.y)) / fitter.chisq_vals[-1])

        # Determine the average CO mixing ratio
        co2 = np.sum(np.array(carbon_dioxide) * np.array(chisquared[numskip:])) / np.sum(chisquared[numskip:])
        logfile.write("CO2 values and their weights\n")
        for c, w in zip(carbon_dioxide, chisquared[numskip:]):
            logfile.write(u"{0:g}\t{1:g}\n".format(c, w))
        logfile.write(u"Best CO2 Mixing Ratio = {0:g}\n".format(co2))
        fitter.AdjustValue({"co2": co2})


        # --------------------------------------------------------------------------------


        # Last one: N2O mixing ratio
        fitter.FitVariable({"n2o": 0.32})
        nitro = []
        numskip = len(chisquared)
        for i in FindOrderNums(orders, [2310, 2340, 2365]):
            print "\n***************************\nFitting order %i: " % (i)
            order = orders[i]
            fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                                "waveend": order.x[-1] + 20.0})
            order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=9, lowreject=9, highreject=10)
            primary = DataStructures.xypoint(x=order.x, y=np.ones(order.x.size))
            primary, model, R = fitter.Fit(data=order.copy(),
                                           resolution_fit_mode="SVD",
                                           fit_source=True,
                                           return_resolution=True,
                                           adjust_wave="data",
                                           wavelength_fit_order=4)
            resolution.append(R)
            waveshifts.append(fitter.shift)
            wave0.append(fitter.data.x.mean())
            nitro.append(fitter.GetValue("n2o"))
            chisquared.append((1.0 - min(model.y)) / fitter.chisq_vals[-1])

        # Determine the average CO mixing ratio
        n2o = np.sum(np.array(nitro) * np.array(chisquared[numskip:])) / np.sum(chisquared[numskip:])
        logfile.write("N2O values and their weights\n")
        for n, w in zip(nitro, chisquared[numskip:]):
            logfile.write(u"{0:g}\t{1:g}\n".format(n, w))
        logfile.write(u"Best N2O Mixing Ratio = {0:g}\n".format(n2o))
        fitter.AdjustValue({"n2o": n2o})

        #Get the average resolution and velocity shift
        velshifts = np.array(waveshifts) / np.array(wave0) * constants.c.cgs.to(units.km / units.s).value
        vel = np.sum(velshifts * np.array(chisquared)) / np.sum(chisquared)
        logfile.write("resolution, velocity shifts and their weights\n")
        for R, v, w in zip(resolution, velshifts, chisquared):
            logfile.write(u"{0:g}\t{1:g}\t{2:g}\n".format(R, v, w))
        resolution = np.sum(resolution * np.array(chisquared)) / np.sum(chisquared)
        logfile.write(u"Best resolution = {0:.5f}\n".format(resolution))
        logfile.write(u"Best velocity shift = {0:.4f} km/s\n".format(vel))



        # --------------------------------------------------------------------------------


        # Finally, apply these parameters to all orders in the data
        fitter.AdjustValue({"resolution": resolution, })
        for i, order in enumerate(orders):
            print u"\n\nGenerating model for order {0:d} of {1:d}\n".format(i, len(orders))
            fitter.AdjustValue({"wavestart": order.x[0] - 20.0,
                                "waveend": order.x[-1] + 20.0})
            fitpars = [fitter.const_pars[j] for j in range(len(fitter.parnames)) if fitter.fitting[j]]
            order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=9, lowreject=2, highreject=10)
            fitter.ImportData(order.copy())
            fitter.resolution_fit_mode = "SVD"
            #wave0 = order.x.mean()
            #fitter.shift = vel/(constants.c.cgs.value*units.cm.to(units.km)) * wave0
            print "fitter.shift = ", fitter.shift
            primary, model = fitter.GenerateModel(fitpars,
                                                  separate_source=True,
                                                  return_resolution=False)

            data = fitter.data
            if min(model.y) > 0.98:
                #The wavelength calibration might be off
                wave0 = order.x.mean()
                fitter.shift = vel / (constants.c.cgs.to(units.km / units.s).value) * wave0
                model = fitter.GenerateModel(fitpars, separate_source=False, nofit=True)
                model.x /= (1.0 + vel / (constants.c.cgs.value * units.cm.to(units.km)))
                model = FittingUtilities.RebinData(model, order.x)
                data = order.copy()


            # Set up data structures for OutputFitsFile
            columns = {"wavelength": data.x,
                       "flux": data.y,
                       "continuum": data.cont,
                       "error": data.err,
                       "model": model.y,
                       "primary": primary.y}

            if (i == 0 and makenew) or not exists:
                #Save the fitted variables in the primary header
                pri_header = {"RESOLUTION": resolution,
                              "HUMIDITY": humidity,
                              "TEMPERATURE": temperature,
                              "CH4": ch4,
                              "CO2": co2,
                              "CO": co,
                              "N2O": n2o}
                HelperFunctions.OutputFitsFileExtensions(columns, fname, outfilename, mode="new",
                                                         primary_header=pri_header)
                exists = True

            else:
                HelperFunctions.OutputFitsFileExtensions(columns, outfilename, outfilename, mode="append")

        logfile.close()


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
        orders = HelperFunctions.ReadExtensionFits(fname)

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
                primary, model = fitter.Fit(data=order.copy(),
                                            resolution_fit_mode="SVD",
                                            fit_source=True,
                                            return_resolution=False,
                                            adjust_wave="model",
                                            wavelength_fit_order=4)
                humidity = fitter.GetValue("h2o")
                temperature = fitter.GetValue("temperature")
                ch4 = fitter.GetValue("ch4")
                co2 = fitter.GetValue("co2")
                co = fitter.GetValue("co")
                n2o = fitter.GetValue("n2o")
                order = fitter.data

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
