import os
import sys
import FittingUtilities
from re import search

import numpy
from scipy.interpolate import InterpolatedUnivariateSpline as interp
from astropy.io import fits as pyfits
from astropy.io import ascii
from astropy import units, constants

import DataStructures
import SpectralTypeRelations
from PlotBlackbodies import Planck
import Smooth
import HelperFunctions
import Broaden
from Search_Fast import Process_Data
import Correlate



# Ensure a directory exists. Create it if not
def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)


if "darwin" in sys.platform:
    modeldir = "/Volumes/DATADRIVE/Stellar_Models/PHOENIX/Stellar/Vband/"
elif "linux" in sys.platform:
    modeldir = "/media/FreeAgent_Drive/SyntheticSpectra/Sorted/Stellar/Vband/"
else:
    modeldir = raw_input("sys.platform not recognized. Please enter model directory below: ")
    if not modeldir.endswith("/"):
        modeldir = modeldir + "/"

# Define regions contaminated by telluric residuals or other defects. We will not use those regions in the cross-correlation
badregions = [[0, 1520],  #Blue end of H band (lots of water absorption)
              [1720, 2080],  #In between H and K bands (lots of water absorption)
              [2350, 2500]]  #Red end of K band (lots of water absorption)


#Set up model list
model_list = [modeldir + "lte30-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte31-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte32-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte33-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte34-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte35-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte36-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte37-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte38-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte39-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte40-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte41-4.5-0.0.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte42-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte43-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte44-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte45-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte46-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte47-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte48-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte49-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte50-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte51-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte52-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte53-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte54-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte55-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte56-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte57-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte58-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte59-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte60-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted"]

MS = SpectralTypeRelations.MainSequence()
PMS = SpectralTypeRelations.PreMainSequence(
    pms_tracks_file="%s/Dropbox/School/Research/Stellar_Evolution/Baraffe_Tracks.dat" % (os.environ["HOME"]),
    track_source="Baraffe")
PMS2 = SpectralTypeRelations.PreMainSequence()


def GetFluxRatio(sptlist, Tsec, xgrid, age=None):
    """
      Returns the flux ratio between the secondary star of temperature Tsec
      and the (possibly multiple) primary star(s) given in the
      'sptlist' list (given as spectral types)
      xgrid is a numpy.ndarray containing the x-coordinates to find the
        flux ratio at (in nm)

      The age of the system is found from the main-sequence age of the
        earliest spectral type in sptlist, if it is not given
    """
    prim_flux = numpy.zeros(xgrid.size)
    sec_flux = numpy.zeros(xgrid.size)

    #First, get the age of the system
    if age is None:
        age = GetAge(sptlist)

    #Now, determine the flux from the primary star(s)
    for spt in sptlist:
        end = search("[0-9]", spt).end()
        T = MS.Interpolate(MS.Temperature, spt[:end])
        R = PMS2.GetFromTemperature(age, T, key="Radius")
        prim_flux += Planck(xgrid * units.nm.to(units.cm), T) * R ** 2

    #Determine the secondary star flux
    R = PMS.GetFromTemperature(age, Tsec, key="Radius")
    sec_flux = Planck(xgrid * units.nm.to(units.cm), Tsec) * R ** 2

    return sec_flux / prim_flux


def GetMass(spt, age):
    """
    Returns the mass of the system in solar masses
    spt: Spectral type of the star
    age: age, in years, of the system
    """

    #Get temperature
    end = search("[0-9]", spt).end()
    T = MS.Interpolate(MS.Temperature, spt[:end])

    # Determine which tracks to use
    if spt[0] == "O" or spt[0] == "B" or spt[0] == "A" or spt[0] == "F":
        return PMS2.GetFromTemperature(age, T, key="Mass")
    else:
        return PMS.GetFromTemperature(age, T, key="Mass")


def GetAge(sptlist):
    """
    Returns the age of the system, in years, given a list
    of spectral types. It determines the age as the
    main sequence lifetime of the earliest-type star
    """

    lowidx = 999
    for spt in sptlist:
        if "I" in spt:
            #Pre-main sequence. Just use the age of an early O star,
            # which is ~ 1Myr
            lowidx = 1
            break
        end = search("[0-9]", spt).end()
        idx = MS.SpT_To_Number(spt[:end])
        if idx < lowidx:
            lowidx = idx

    spt = MS.Number_To_SpT(lowidx)
    Tprim = MS.Interpolate(MS.Temperature, spt)
    age = PMS2.GetMainSequenceAge(Tprim, key="Temperature")
    return age


if __name__ == "__main__":
    #Define some constants to use
    vsini_secondary = 20 * units.km.to(units.cm)
    resolution = 80000
    lightspeed = constants.c.cgs.value * units.cm.to(units.km)
    smooth_factor = 0.8
    vel_list = range(-400, 400, 50)
    companion_file = "%s/Dropbox/School/Research/AstarStuff/TargetLists/Multiplicity.csv" % (os.environ["HOME"])
    vsini_file = "%s/School/Research/Useful_Datafiles/Vsini.csv" % (os.environ["HOME"])
    fileList = []
    tolerance = 5.0
    debug = False
    vsini_skip = 10
    vsini_idx = 1
    for arg in sys.argv[1:]:
        if "-m" in arg:
            companion_file = arg.split("=")[1]
        elif "-tol" in arg:
            tolerance = float(arg.split("=")[1])
        elif "-d" in arg:
            debug = True
        elif "-vsinifile" in arg:
            vsini_file = arg.split("=")[-1]
        elif "-vsiniskip" in arg:
            vsini_skip = int(arg.split("=")[-1])
        elif "-vsiniidx" in arg:
            vsini_idx = int(arg.split("=")[-1])
        else:
            fileList.append(arg)

    # Make sure each output file exists:
    logfilenames = {}
    output_directories = {}
    for fname in fileList:
        output_dir = "Sensitivity/"
        outfilebase = fname.split(".fits")[0]
        if "/" in fname:
            dirs = fname.split("/")
            output_dir = ""
            outfilebase = dirs[-1].split(".fits")[0]
            for directory in dirs[:-1]:
                output_dir = output_dir + directory + "/"
            output_dir = output_dir + "Sensitivity/"
        HelperFunctions.ensure_dir(output_dir)
        output_directories[fname] = output_dir

        #Make the summary file
        logfile = open(output_dir + "logfile.dat", "w")
        logfile.write("Sensitivity Analysis:\n*****************************\n\n")
        logfile.write(
            "Filename\t\t\tPrimary Temperature\tSecondary Temperature\tMass (Msun)\tMass Ratio\tVelocity\tPeak Correct?\tSignificance\n")
        logfile.close()
        logfilenames[fname] = output_dir + "logfile.dat"


    # Read in the companion file
    companions = ascii.read(companion_file)[20:]

    # Read in the vsini file
    vsini_data = ascii.read(vsini_file)[vsini_skip:]

    # Now, start loop over the models:
    firstmodel = 0
    lastmodel = len(model_list)
    fitorders = [7, 6, 8, 8, 8, 7, 7, 10, 9, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 8, 8, 8, 8]
    for modelnum, modelfile in enumerate(model_list[firstmodel:lastmodel]):
        if "PHOENIX2004" in modelfile:
            temp = int(modelfile.split("lte")[-1][:2]) * 100
            gravity = float(modelfile.split("lte")[-1][3:6])
            metallicity = float(modelfile.split("lte")[-1][6:10])
        elif "PHOENIX-ACES" in modelfile:
            temp = int(modelfile.split("lte")[-1][:2]) * 100
            gravity = float(modelfile.split("lte")[-1][3:7])
            metallicity = float(modelfile.split("lte")[-1][7:11])
        print "Reading in file %s" % modelfile
        x, y = numpy.loadtxt(modelfile, usecols=(0, 1), unpack=True)
        print "Processing file..."
        c = FittingUtilities.Continuum(x, y, fitorder=int(fitorders[modelnum + firstmodel]), lowreject=2, highreject=7)
        model = DataStructures.xypoint(x=x * units.angstrom.to(units.nm) / 1.00026, y=10 ** y, cont=10 ** c)
        model = FittingUtilities.RebinData(model, numpy.linspace(model.x[0], model.x[-1], model.size()))
        model = Broaden.RotBroad(model, vsini_secondary)
        model = Broaden.ReduceResolution2(model, resolution)
        modelfcn = interp(model.x, model.y / model.cont)


        # Now that we have a spline function for the broadened data,
        # begin looping over the files
        for fname in fileList:
            print fname
            output_dir = output_directories[fname]
            outfile = open(logfilenames[fname], "a")

            # Read in and process the data like I am about to look for a companion
            orders_original = Process_Data(fname)

            #Find the vsini of the primary star with my spreadsheet
            starname = pyfits.getheader(fname)["object"]
            found = False
            for data in vsini_data:
                if data[0] == starname:
                    vsini = abs(float(data[vsini_idx]))
                    found = True
            if not found:
                sys.exit("Cannot find %s in the vsini data: %s" % (starname, vsini_file))
            print starname, vsini

            #Check for companions in my master spreadsheet
            known_stars = []
            if starname in companions.field(0):
                row = companions[companions.field(0) == starname]
                print row
                print row['col1']
                known_stars.append(row['col1'].item())
                ncompanions = int(row['col4'].item())
                for comp in range(ncompanions):
                    spt = row["col%i" % (6 + 4 * comp)].item()
                    if not "?" in spt and (spt[0] == "O" or spt[0] == "B" or spt[0] == "A" or spt[0] == "F"):
                        sep = row["col%i" % (7 + 4 * comp)].item()
                        if (not "?" in sep) and float(sep) < 4.0:
                            known_stars.append(spt)
            else:
                sys.exit("Star not found in multiplicity library!")

            #Determine the age of the system and properties of the primary and secondary star
            age = GetAge(known_stars)
            primary_spt = known_stars[0]
            end = search("[0-9]", primary_spt).end()
            primary_temp = MS.Interpolate(MS.Temperature, primary_spt[:end])
            primary_mass = GetMass(primary_spt, age)
            secondary_spt = MS.GetSpectralType(MS.Temperature, temp)
            secondary_mass = GetMass(secondary_spt, age)
            massratio = secondary_mass / primary_mass

            for rv in vel_list:
                print "Testing model with rv = ", rv
                orders = [order.copy() for order in orders_original]  #Make a copy of orders
                model_orders = []
                for ordernum, order in enumerate(orders):
                    #Get the flux ratio
                    scale = GetFluxRatio(known_stars, temp, order.x, age=age)
                    print "Scale factor for order %i is %.3g" % (ordernum, scale.mean())

                    #Add the model to the data
                    model = (modelfcn(order.x * (1.0 + rv / lightspeed)) - 1.0) * scale
                    order.y += model * order.cont


                    #Smooth data using the vsini of the primary star
                    dx = order.x[1] - order.x[0]
                    npixels = max(21, Smooth.roundodd(vsini / lightspeed * order.x.mean() / dx * smooth_factor))
                    smoothed = Smooth.SmoothData(order,
                                                 windowsize=npixels,
                                                 smoothorder=3,
                                                 lowreject=3,
                                                 highreject=3,
                                                 expand=10,
                                                 numiters=10,
                                                 normalize=False)
                    order.y /= smoothed.y

                    # log-space the data
                    start = numpy.log(order.x[0])
                    end = numpy.log(order.x[-1])
                    xgrid = numpy.logspace(start, end, order.size(), base=numpy.e)
                    logspacing = numpy.log(xgrid[1] / xgrid[0])
                    order = FittingUtilities.RebinData(order, xgrid)

                    # Generate a model with the same log-spacing (and no rv shift)
                    dlambda = order.x[order.size() / 2] * 1000 * 1.5 / lightspeed
                    start = numpy.log(order.x[0] - dlambda)
                    end = numpy.log(order.x[-1] + dlambda)
                    xgrid = numpy.exp(numpy.arange(start, end + logspacing, logspacing))
                    model = DataStructures.xypoint(x=xgrid, cont=numpy.ones(xgrid.size))
                    model.y = modelfcn(xgrid)

                    # Save model order
                    model_orders.append(model)
                    #orders[ordernum] = order.copy()

                #Do the actual cross-correlation
                print "Cross-correlating..."
                corr = Correlate.Correlate(orders, model_orders, debug=debug, outputdir="Sensitivity_Testing/")

                # Check if we found the companion
                idx = numpy.argmax(corr.y)
                vmax = corr.x[idx] - 2.0  #There is a systematic offset for some reason!
                fit = FittingUtilities.Continuum(corr.x, corr.y, fitorder=2, lowreject=3, highreject=2.5)
                corr.y -= fit
                goodindices = numpy.where(numpy.abs(corr.x - rv) > 100)[0]
                mean = corr.y[goodindices].mean()
                std = corr.y[goodindices].std()
                significance = (corr.y[idx] - mean) / std
                if debug:
                    corrfile = "%s%s_t%i_v%i" % (output_dir, fname.split("/")[-1].split(".fits")[0], temp, rv)
                    print "Outputting CCF to %s" % corrfile
                    numpy.savetxt(corrfile, numpy.transpose((corr.x, corr.y - mean, numpy.ones(corr.size()) * std)),
                                  fmt="%.10g")
                if abs(vmax - rv) <= tolerance:
                    #Found
                    outfile.write("%s\t%i\t\t\t%i\t\t\t\t%.2f\t\t%.4f\t\t%i\t\tyes\t\t%.2f\n" % (
                        fname, primary_temp, temp, secondary_mass, massratio, rv, significance))
                else:
                    #Not found
                    outfile.write("%s\t%i\t\t\t%i\t\t\t\t%.2f\t\t%.4f\t\t%i\t\tno\t\tN/A\n" % (
                        fname, primary_temp, temp, secondary_mass, massratio, rv))
                print "Done with rv ", rv
            outfile.close()



