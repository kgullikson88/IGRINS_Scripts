import sys
from scipy.interpolate import InterpolatedUnivariateSpline as interp
from collections import defaultdict
import os

import DataStructures
import numpy as np
import pylab
from astropy.io import fits as pyfits
import FittingUtilities

import FitsUtils
import FindContinuum
import HelperFunctions


plot = False


def MedianAdd(fileList, outfilename="Total.fits"):
    all_data = []
    numorders = []
    medians = []
    for fname in fileList:
        observation = FitsUtils.MakeXYpoints(fname, extensions=True, x="wavelength", y="flux", cont="continuum",
                                             errors="error")
        all_data.append(observation)
        numorders.append(len(observation))
        medians.append([np.median(order.y) for order in observation])

    if any(n != numorders[0] for n in numorders):
        print "Error! Some of the files had different numbers of orders!"
        for i in range(len(fileList)):
            print fileList[i], numorders[i]
        sys.exit()

    # If we get this far, all is well. Add each order indidually
    numorders = numorders[0]
    if outfilename == "None":
        outfilename = "Total.fits"
    column_list = []
    for i in range(numorders):
        x = all_data[0][i].x
        total = np.zeros((len(all_data), x.size))
        error = np.zeros(x.size)
        norm = 0.0
        for j, observation in enumerate(all_data):
            observation[i].y[observation[i].y < 0.0] = 0.0
            flux = interp(observation[i].x, observation[i].y / medians[j][i])
            error += interp(observation[i].x, observation[i].err ** 2, k=1)(x)
            total[j] = flux(x)
            norm += medians[j][i]

            #pylab.figure(2)
            #for j in range(total.shape[0]):
            #pylab.(x, total[j])
        flux = np.median(total, axis=0) * norm
        cont = FittingUtilities.Continuum(x, flux, fitorder=3, lowreject=1.5, highreject=5)
        #Set up data structures for OutputFitsFile
        columns = {"wavelength": x,
                   "flux": flux,
                   "continuum": cont,
                   "error": np.sqrt(error)}
        column_list.append(columns)

        #pylab.figure(1)
        #pylab.plot(x, flux/cont)
        #pylab.plot(total.x, total.cont)

    print "Outputting to %s" % outfilename
    #pylab.show()
    FitsUtils.OutputFitsFileExtensions(column_list, fileList[0], outfilename, mode="new")

    #Add the files used to the primary header of the new file
    hdulist = pyfits.open(outfilename, mode='update')
    header = hdulist[0].header
    for i in range(len(fileList)):
        header.set("FILE%i" % (i + 1), fileList[i], "File %i used in Co-Adding" % (i + 1))
    hdulist[0].header = header
    hdulist.flush()
    hdulist.close()


def Add(fileList, outfilename=None):
    all_data = []
    numorders = []
    exptime = 0.0
    for fname in fileList:
        observation = HelperFunctions.ReadFits(fname, extensions=True, x="wavelength", y="flux", cont="continuum",
                                               errors="error")
        all_data.append(observation)
        numorders.append(len(observation))
        exptime += pyfits.getheader(fname)['EXPTIME']

    if any(n != numorders[0] for n in numorders):
        print "Error! Some of the files had different numbers of orders!"
        for i in range(len(fileList)):
            print fileList[i], numorders[i]
        sys.exit()

    # If we get this far, all is well. Add each order indidually
    numorders = numorders[0]
    if outfilename == None:
        outfilename = "Total.fits"
    column_list = []
    for i in range(numorders):
        #Determine the x-axis
        dx = 0
        left = 0
        right = 9e9
        for observation in all_data:
            dx += (observation[i].x[-1] - observation[i].x[0]) / float(observation[i].size() - 1)
            if observation[i].x[0] > left:
                left = observation[i].x[0]
            if observation[i].x[-1] < right:
                right = observation[i].x[-1]
        dx /= float(len(all_data))
        xgrid = np.arange(left, right + dx, dx)
        total = DataStructures.xypoint(x=xgrid)
        total.y = np.zeros(total.size())
        total.err = np.zeros(total.size())

        #Add the data
        for observation in all_data:
            observation[i].y[observation[i].y < 0.0] = 0.0
            rebinned = FittingUtilities.RebinData(observation[i], total.x)
            total.y += rebinned.y
            total.err += rebinned.err ** 2

        total.err = np.sqrt(total.err)
        total.cont = FittingUtilities.Continuum(total.x, total.y, fitorder=3, lowreject=1.5, highreject=5)

        #Check if there is a dip in flux on either side of the order
        std = np.median(total.err / total.cont)
        #print std
        left = 0
        right = total.size()
        for i in range(4):
            if abs(total.y[i + 1] / total.cont[i + 1] - total.y[i] / total.cont[i]) > 10 * std:
                left = i + 1
            if abs(total.y[-(i + 2)] / total.cont[-(i + 2)] - total.y[-(i + 1)] / total.cont[-(i + 1)]) > 10 * std:
                right = total.size() - (i + 1)
        total = total[left:right]

        #Set up data structures for OutputFitsFile
        columns = {"wavelength": total.x,
                   "flux": total.y,
                   "continuum": total.cont,
                   "error": total.err}
        column_list.append(columns)

        if plot:
            pylab.plot(total.x, total.y / total.cont)
            #pylab.plot(total.x, total.cont)

    print "Outputting to %s" % outfilename
    if plot:
        pylab.show()
    HelperFunctions.OutputFitsFileExtensions(column_list, fileList[0], outfilename, mode="new")

    #Add the files used to the primary header of the new file
    hdulist = pyfits.open(outfilename, mode='update')
    header = hdulist[0].header
    for i in range(len(fileList)):
        header.set("FILE%i" % (i + 1), fileList[i], "File %i used in Co-Adding" % (i + 1))
    header.set('EXPTIME', exptime, 'Total exposure time (seconds)')
    hdulist[0].header = header
    hdulist.flush()
    hdulist.close()


if __name__ == "__main__":
    fileList = []
    for arg in sys.argv[1:]:
        fileList.append(arg)

    if len(fileList) > 1:
        #Add(fileList)
        fileDict = defaultdict(list)
        for fname in fileList:
            header = pyfits.getheader(fname)
            starname = header['OBJECT'].replace(" ", "_")
            #fileDict[starname].append(fname)
            #print(fname, starname)
            starname1 = header['OBJECT1'].replace(" ", "_")
            starname2 = header['OBJECT2'].replace(" ", "_")
            key = "{}+{}".format(starname1, starname2)
            fileDict[key].append(fname)
        for star in fileDict.keys():
            Add(fileDict[star], outfilename="%s_bright.fits" % star)
    else:
        allfiles = [f for f in os.listdir("./") if f.startswith("echi") and "-0" in f and "telluric" in f]
        fileDict = defaultdict(list)
        for fname in allfiles:
            header = pyfits.getheader(fname)
            starname = header['OBJECT'].replace(" ", "_")
            fileDict[starname].append(fname)
        for star in fileDict.keys():
            Add(fileDict[star], outfilename="%s.fits" % star)
    
