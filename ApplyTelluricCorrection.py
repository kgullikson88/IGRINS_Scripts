import sys
import os
import FittingUtilities
from FittingUtilities import FindLines

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


def fit_2dspec(xl, yl, zl, x_degree=4, y_degree=3,
               x_domain=None, y_domain=None):
    from astropy.modeling import fitting
    # Fit the data using astropy.modeling
    if x_domain is None:
        x_domain = [min(xl), max(xl)]
    # more room for y_domain??
    if y_domain is None:
        #y_domain = [orders[0]-2, orders[-1]+2]
        y_domain = [min(yl), max(yl)]
    from astropy.modeling.polynomial import Chebyshev2D
    p_init = Chebyshev2D(x_degree=x_degree, y_degree=y_degree,
                         x_domain=x_domain, y_domain=y_domain)
    f = fitting.LinearLSQFitter()

    p = f(p_init, xl, yl, zl)

    for i in [0]:
        dd = p(xl, yl) - zl
        m = np.abs(dd) < 3.*dd.std()
        p = f(p, xl[m], yl[m], zl[m])

    return p, m


def fit_wavelength(orders, ordernums, first_order=None, last_order=None, x_degree=4, y_degree=3):
    """ Fit the wavelength in a whole chip, and return the 2D polynomial callable
    """
    pixel_list = []
    ordernum_list = []
    wave_list = []
    if first_order is None:
        first_order = 0
    if last_order is None:
        last_order = len(orders) - 1
    
    for order, ordernum in zip(orders[first_order:last_order+1], ordernums[first_order:last_order+1]):
        lines = FindLines(order)
        pixel_list.extend(lines)
        ordernum_list.extend(np.ones_like(lines)*ordernum)
        wave_list.extend(order.x[lines])
    
    pixel_list = np.array(pixel_list)
    ordernum_list = np.array(ordernum_list)
    wave_list = np.array(wave_list)
    p, m = fit_2dspec(pixel_list, ordernum_list, wave_list*ordernum_list, 
                      x_degree=x_degree, y_degree=y_degree)
    
    return p

def fix_chip_wavelength(model_orders, data_orders, band_cutoff=1870):
    """ Adjust the wavelength in data_orders to be self-consistent
    """
    # H band
    model_orders_H = [o.copy() for o in model_orders if o.x[-1] < band_cutoff]
    data_orders_H = [o.copy() for o in data_orders if o.x[-1] < band_cutoff]
    ordernums_H = 121.0 - np.arange(len(model_orders_H))
    p_H = fit_wavelength(model_orders_H, ordernums_H, first_order=3, last_order=len(ordernums_H) - 4)

    # K band
    model_orders_K = [o.copy() for o in model_orders if o.x[-1] > band_cutoff]
    data_orders_K = [o.copy() for o in data_orders if o.x[-1] > band_cutoff]
    ordernums_K = 92.0 - np.arange(len(model_orders_K))
    p_K = fit_wavelength(model_orders_K, ordernums_K, first_order=7, last_order=len(ordernums_K) - 4)

    new_orders = []
    for i, order in enumerate(data_orders):
        pixels = np.arange(order.size(), dtype=np.float)
        if order.x[-1] < band_cutoff:
            # H band
            ordernum = ordernums_H[i] * np.ones_like(pixels)
            wave = p_H(pixels, ordernum) / ordernum
        else:
            # K band
            ordernum = ordernums_K[i-len(ordernums_H)] * np.ones_like(pixels)
            wave = p_K(pixels, ordernum) / ordernum
            
        new_orders.append(DataStructures.xypoint(x=wave, y=order.y, cont=order.cont, err=order.err))
    return new_orders



def Correct_New(original, corrected, get_primary=False, plot=False, *args, **kwargs):
    original_orders, _ = ReadCorrectedFile(corrected, yaxis="flux")
    primary_orders, _ = ReadCorrectedFile(corrected, yaxis="primary")
    model_orders, corrected_headers = ReadCorrectedFile(corrected)
    primary_header = pyfits.getheader(corrected)

    # Fix the wavelength axis in the data orders
    original_orders = fix_chip_wavelength(model_orders, original_orders)

    new_orders = []
    for i, (original, model, primary) in enumerate(zip(original_orders, model_orders, primary_orders)):
        original.cont /= primary.y

        if plot:
            plt.plot(original.x, original.y/original.cont, 'k-', alpha=0.5)
            plt.plot(model.x, model.y, 'r-')

        model.y[model.y < 1e-4] = 1e-4
        original.y /= model.y 
        new_orders.append(original.copy())
    if plot:
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Flux')
        plt.title('Comparison for file {}'.format(corrected))
        plt.show()
    return new_orders



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
    # if offset == None:
    #    offset = len(original_orders) - len(corrected_orders)
    #print "Offset = ", offset
    new_orders = []
    for i, data in enumerate(original_orders):

        data.cont = FittingUtilities.Continuum(data.x, data.y)
        try:
            idx = HelperFunctions.FindOrderNums(corrected_orders, [np.median(data.x)])
            if len(idx) < 1:
                continue
            elif len(idx) > 1:
                j = np.argmin([np.abs(np.median(corrected_orders[k].x) - np.median(data.x)) for k in idx])
                idx = idx[j]
            else:
                idx = idx[0]
            model = corrected_orders[idx]
            header = corrected_headers[idx]
            if abs(np.median(model.x) - np.median(data.x)) > 1:
                continue
            if get_primary:
                primary = primary_orders[idx]
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

        # if np.sum((model.x-data.x)**2) > 1e-8:
        # model = FittingUtilities.RebinData(model, data.x)

        data.y[data.y / data.cont < 1e-5] = 1e-5 * data.cont[data.y / data.cont < 1e-5]
        badindices = np.where(np.logical_or(data.y <= 0, model.y < 0.05))[0]
        model.y[badindices] = data.y[badindices] / data.cont[badindices]
        model.y[model.y < 1e-5] = 1e-5

        data.x = model.x
        data.y /= model.y
        data.err /= model.y
        if get_primary:
            data.y /= primary.y
        new_orders.append(data[1:-1].copy())
    if plot:
        plt.show()
    return new_orders


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
        header = hdulist[i + 1].header
        fitter.data = order.copy()
        # order = orders[-9]
        temperature = header['temperature']
        humidity = header['humidity']
        ch4 = header['ch4']
        co2 = header['co2']
        co = header['co']
        n2o = header['n2o']

        # Make the model
        model = fitter.Modeler.MakeModel(temperature=temperature,
                                         humidity=humidity,
                                         co2=co2,
                                         co=co,
                                         ch4=ch4,
                                         n2o=n2o,
                                         lowfreq=1e7 / (order.x[-1] + 10),
                                         highfreq=1e7 / (order.x[0] - 10))
        xgrid = np.linspace(model.x[0], model.x[-1], model.size())
        model_original = FittingUtilities.RebinData(model, xgrid)
        model = FittingUtilities.ReduceResolution(model_original, 48000)

        # Do the wavelength correction
        modelfcn, mean = fitter.FitWavelengthNew(order, model, fitorder=5)
        model_original.x -= modelfcn(model.x - mean)

        model = fitter.Broaden2(order, model_original)

        if plot:
            ax1.plot(order.x, order.y / order.cont)
            ax1.plot(model.x, model.y)
            ax2.plot(order.x, order.y / order.cont - model.y)
            ax3.plot(order.x, order.y / (order.cont * model.y))

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
    plot = False
    if len(sys.argv) > 2:
        original = sys.argv[1]
        corrected = sys.argv[2]
        if len(sys.argv) > 3 and "prim" in sys.argv[3]:
            primary = True

        outfilename = "%s_telluric_corrected.fits" % (original.split(".fits")[0])
        print "Outputting to %s" % outfilename

        #corrected_orders = Correct_Old(original, corrected, offset=None, get_primary=primary, plot=plot)
        corrected_orders = Correct_New(original, corrected, offset=None, get_primary=primary, plot=plot)

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
            original = corrected.split("Corrected_")[-1].split("-")[0] + ".fits"

            print corrected, original
            if not os.path.exists(original):
                print('******************\n\nOriginal file ({}) not found!!\n\n*********************'.format(original))
                continue

            #corrected_orders = Correct_Old(original, corrected, offset=None, plot=plot)
            corrected_orders = Correct_New(original, corrected, offset=None, plot=plot)

            outfilename = "{0:s}_telluric_corrected.fits".format(original.split(".fits")[0])
            print "Outputting to %s" % outfilename

            column_list = []
            if plot:
                plt.figure(2)
            for i, data in enumerate(corrected_orders):
                if plot:
                    plt.plot(data.x, data.y / data.cont)
                # Set up data structures for OutputFitsFile
                columns = {"wavelength": data.x,
                           "flux": data.y,
                           "continuum": data.cont,
                           "error": data.err}
                column_list.append(columns)
            HelperFunctions.OutputFitsFileExtensions(column_list, original, outfilename, mode="new")

            if plot:
                plt.title(original)
                plt.xlabel("Wavelength (nm)")
                plt.ylabel("Flux")
                plt.ylim((0, 100000))
                plt.show()


if __name__ == "__main__":
    main1()
