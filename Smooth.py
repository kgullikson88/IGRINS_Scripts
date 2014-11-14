__author__ = 'Kevin Gullikson'
import FittingUtilities
import sys
import os

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits, ascii

import GenericSmooth
import HelperFunctions


if __name__ == "__main__":
    fileList = []
    plot = False
    vsini_file = "%s/School/Research/Useful_Datafiles/Vsini.csv" % (os.environ["HOME"])
    vsini_skip = 9
    vsini_idx = 1
    for arg in sys.argv[1:]:
        if "-p" in arg:
            plot = True
        elif "-vsinifile" in arg:
            vsini_file = arg.split("=")[-1]
        elif "-vsiniskip" in arg:
            vsini_skip = int(arg.split("=")[-1])
        elif "-vsiniidx" in arg:
            vsini_idx = int(arg.split("=")[-1])
        else:
            fileList.append(arg)

    # Read in the vsini table
    vsini_data = ascii.read(vsini_file)[vsini_skip:]

    if len(fileList) == 0:
        fileList = [f for f in os.listdir("./") if f.endswith("telluric_corrected.fits")]
    for fname in fileList:
        orders = HelperFunctions.ReadExtensionFits(fname)

        # Find the vsini of this star
        header = fits.getheader(fname)
        starname = header["object"]
        for data in vsini_data:
            if data[0].strip().lower() == starname.strip().lower():
                vsini = abs(float(data[vsini_idx]))
                break
        else:
            sys.exit("Cannot find %s in the vsini data: %s" % (starname, vsini_file))
        print starname, vsini

        # Begin looping over the orders
        column_list = []
        header_list = []
        for i, order in enumerate(orders):
            print "Smoothing order %i/%i" % (i + 1, len(orders))
            #Fix errors
            order.err[order.err > 1e8] = np.sqrt(order.y[order.err > 1e8])

            #Linearize
            xgrid = np.linspace(order.x[0], order.x[-1], order.x.size)
            order = FittingUtilities.RebinData(order, xgrid)

            dx = order.x[1] - order.x[0]
            smooth_factor = 0.8
            theta = max(21, GenericSmooth.roundodd(vsini / 3e5 * order.x.mean() / dx * smooth_factor))
            denoised = GenericSmooth.SmoothData(order,
                                                windowsize=theta,
                                                smoothorder=3,
                                                lowreject=3,
                                                highreject=3,
                                                expand=10,
                                                numiters=10)
            #denoised, theta = GPSmooth(order.copy())
            #denoised, theta = CrossValidation(order.copy(), 5, 2, 2, 10)
            #denoised, theta = OptimalSmooth(order.copy())
            #denoised.y *= order.cont/order.cont.mean()
            print "Window size = %.4f nm" % theta
            denoised.cont = FittingUtilities.Continuum(denoised.x, order.y / denoised.y, fitorder=1, lowreject=2,
                                                       highreject=2)

            column = {"wavelength": denoised.x,
                      "flux": order.y / denoised.y,
                      "continuum": denoised.cont,
                      "error": denoised.err}
            header_list.append((("Smoother", theta, "Smoothing Parameter"),))
            column_list.append(column)
            if plot:
                plt.figure(1)
                plt.plot(order.x, order.y / order.y.mean())
                plt.plot(denoised.x, denoised.y / denoised.y.mean())
                plt.title(starname)
                plt.figure(2)
                plt.plot(order.x, order.y / denoised.y)
                plt.title(starname)
                #plt.plot(order.x, (order.y-denoised.y)/np.median(order.y))
                #plt.show()
        if plot:
            plt.show()
        outfilename = "%s_smoothed.fits" % (fname.split(".fits")[0])
        print "Outputting to %s" % outfilename
        HelperFunctions.OutputFitsFileExtensions(column_list, fname, outfilename, mode='new', headers_info=header_list)
