import sys
import os
import FittingUtilities

from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import numpy as np

import DataStructures
import HelperFunctions


def ReadCorrectedFile(fname, yaxis="model"):
    orders = []
    headers = []
    hdulist = pyfits.open(fname)
    numorders = len(hdulist)
    for i in range(1, numorders):
        order = hdulist[i].data
        xypt = DataStructures.xypoint(x=order.field("wavelength"),
                                      y=order.field(yaxis),
                                      cont=order.field("continuum"),
                                      err=order.field("error"))

        orders.append(xypt)
        headers.append(hdulist[i].header)
    return orders, headers


def Correct_Old(original, corrected, offset=None, get_primary=False, plot=False):
    # Read in the data and model
    original_orders = HelperFunctions.ReadFits(original, extensions=True, x="wavelength", y="flux", errors="error",
                                               cont="continuum")
    corrected_orders, corrected_headers = ReadCorrectedFile(corrected)
    test_orders, header = ReadCorrectedFile(corrected, yaxis="flux")

    primary_header = pyfits.getheader(corrected)

    if plot:
        for order, model in zip(test_orders, corrected_orders):
            plt.plot(order.x, order.y / order.cont)
            plt.plot(model.x, model.y)
        plt.title("Correction in corrected file only")
        plt.show()

    if get_primary:
        primary_orders = ReadCorrectedFile(corrected, yaxis="primary")[0]
    if offset == None:
        offset = len(original_orders) - len(corrected_orders)
    print "Offset = ", offset
    for i in range(len(original_orders) - offset):
        data = original_orders[i]
        data.cont = FittingUtilities.Continuum(data.x, data.y)
        try:
            model = corrected_orders[i]
            header = corrected_headers[i]
            if get_primary:
                primary = primary_orders[i]
            if i == 0 and "CH4" in primary_header.keys():
                print "Humidity = {0:g}\nT = {1:g}\n[CH4] = {2:g}\n[CO2] = {3:g}\n" \
                      "CO = {4:g}\nN2O = {5:g}\n".format(primary_header['humidity'],
                                                         primary_header['airtemp'],
                                                         primary_header['ch4'],
                                                         primary_header['co2'],
                                                         primary_header['co'],
                                                         primary_header['n2o'])
            else:
                print "Humidity = {0:g}\nT = {1:g}\n[CH4] = {2:g}\n[CO2] = {3:g}\n" \
                      "CO = {4:g}\nN2O = {5:g}\n".format(header['humidity'],
                                                         header['temperature'],
                                                         header['ch4'],
                                                         header['co2'],
                                                         header['co'],
                                                         header['n2o'])
        except IndexError:
            model = DataStructures.xypoint(x=data.x, y=np.ones(data.x.size))
            print "Warning!!! Telluric Model not found for order %i" % i

        if plot:
            plt.figure(1)
            plt.plot(data.x, data.y / data.cont)
            plt.plot(model.x, model.y)

        if model.size() < data.size():
            left = np.searchsorted(data.x, model.x[0])
            right = np.searchsorted(data.x, model.x[-1])
            if right < data.size():
                right += 1
            data = data[left:right]
        elif model.size() > data.size():
            sys.exit("Error! Model size (%i) is larger than data size (%i)" % (model.size(), data.size()))

        #if np.sum((model.x-data.x)**2) > 1e-8:
        #  model = FittingUtilities.RebinData(model, data.x)

        data.y[data.y / data.cont < 1e-5] = 1e-5 * data.cont[data.y / data.cont < 1e-5]
        badindices = np.where(np.logical_or(data.y <= 0, model.y < 0.05))[0]
        model.y[badindices] = data.y[badindices] / data.cont[badindices]
        model.y[model.y < 1e-5] = 1e-5

        #plt.plot(data.x, data.y / model.y)
        data.y /= model.y
        data.err /= model.y
        if get_primary:
            data.y /= primary.y
        original_orders[i] = data[1:-1].copy()
    if plot:
        plt.show()
    return original_orders



def Correct(original, corrected, offset=None, get_primary=False, plot=False):
    hdulist = pyfits.open(corrected)
    orders = HelperFunctions.ReadExtensionFits(original)
    import TelluricFitter
    import GetAtmosphere
    import os
    fitter = TelluricFitter.TelluricFitter()
    fitter.SetObservatory('mcdonald')
    filenames = [f for f in os.listdir("./") if "GDAS" in f]
    header = hdulist[0].header
    height, Pres, Temp, h2o = GetAtmosphere.GetProfile(filenames, header['date-obs'].split("T")[0], header['ut'])

    fitter.EditAtmosphereProfile("Temperature", height, Temp)
    fitter.EditAtmosphereProfile("Temperature", height, Temp)
    fitter.EditAtmosphereProfile("Temperature", height, Temp)
    fitter.SetBounds({"resolution": [30000, 50000]})

    if plot:
        fig = plt.figure(1)
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312, sharex=ax1)
        ax3 = fig.add_subplot(313, sharex=ax1)
    for i, order in enumerate(orders):
        print "ORDER: ", i
        # Get atmosphere parameters
        header = hdulist[i+1].header
        fitter.data = order.copy()
        #order = orders[-9]
        temperature = header['temperature']
        humidity = header['humidity']
        ch4 = header['ch4']
        co2 = header['co2']
        co = header['co']
        n2o = header['n2o']

        #Make the model
        model = fitter.Modeler.MakeModel(temperature=temperature,
                                         humidity=humidity,
                                         co2=co2,
                                         co=co,
                                         ch4=ch4,
                                         n2o=n2o,
                                         lowfreq=1e7/(order.x[-1] + 10),
                                         highfreq=1e7/(order.x[0] - 10))
        xgrid = np.linspace(model.x[0], model.x[-1], model.size())
        model_original = FittingUtilities.RebinData(model, xgrid)
        model = FittingUtilities.ReduceResolution(model_original, 48000)

        #Do the wavelength correction
        modelfcn, mean = fitter.FitWavelengthNew(order, model, fitorder=5)
        model_original.x -= modelfcn(model.x - mean)

        model = fitter.Broaden2(order, model_original)

        if plot:
            ax1.plot(order.x, order.y/order.cont)
            ax1.plot(model.x, model.y)
            ax2.plot(order.x, order.y/order.cont - model.y)
            ax3.plot(order.x, order.y/(order.cont*model.y))

        model.y[model.y < 0.05] = (order.y / order.cont)[model.y < 0.05]
        order.y /= model.y
        orders[i] = order.copy()

    if plot:
        ax1.set_ylim((-0.5, 2.0))
        ax2.set_ylim((-1.5, 1.5))
        ax3.set_ylim((-0.5, 2.0))
        # plt.show()

    return orders




def main1():
    primary = False
    plot = True
    if len(sys.argv) > 2:
        original = sys.argv[1]
        corrected = sys.argv[2]
        if len(sys.argv) > 3 and "prim" in sys.argv[3]:
            primary = True

        outfilename = "%s_telluric_corrected.fits" % (original.split(".fits")[0])
        print "Outputting to %s" % outfilename

        corrected_orders = Correct_Old(original, corrected, offset=None, get_primary=primary, plot=plot)

        column_list = []
        if plot:
            plt.figure(2)
        for i, data in enumerate(corrected_orders):
            if plot:
                plt.plot(data.x, data.y / data.cont)
                # plt.plot(data.x, data.cont)
            # Set up data structures for OutputFitsFile
            columns = {"wavelength": data.x,
                       "flux": data.y,
                       "continuum": data.cont,
                       "error": data.err}
            column_list.append(columns)
        if plot:
            plt.title("Corrected data")
            plt.show()
        HelperFunctions.OutputFitsFileExtensions(column_list, original, outfilename, mode="new")

    else:
        allfiles = os.listdir("./")
        corrected_files = [f for f in allfiles if "Corrected_" in f and f.endswith("-0.fits")]


        for corrected in corrected_files:
            original = corrected.split("Corrected_")[-1]  #.split("-")[0] + ".fits"

            print corrected, original
            header = pyfits.getheader(original)
            if header['imagetyp'].strip().lower() != 'object' or "solar" in header['object'].lower():
                print "Skipping file %s, with imagetype = %s and object = %s" % (
                original, header['imagetyp'], header['object'])
                continue

            outfilename = "{0:s}_telluric_corrected.fits".format(original.split(".fits")[0])
            print "Outputting to %s" % outfilename

            corrected_orders = Correct(original, corrected, offset=None, plot=plot)

            column_list = []
            if plot:
                plt.figure(2)
            for i, data in enumerate(corrected_orders):
                if plot:
                    plt.plot(data.x, data.y / data.cont)
                #Set up data structures for OutputFitsFile
                columns = {"wavelength": data.x,
                           "flux": data.y,
                           "continuum": data.cont,
                           "error": data.err}
                column_list.append(columns)
            HelperFunctions.OutputFitsFileExtensions(column_list, corrected, outfilename, mode="new")

            if plot:
                plt.title(original)
                plt.xlabel("Wavelength (nm)")
                plt.ylabel("Flux")
                plt.show()






if __name__ == "__main__":
    main1()
