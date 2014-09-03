import sys
import FittingUtilities
from collections import defaultdict

import numpy
from scipy.interpolate import InterpolatedUnivariateSpline as interp
from astropy import units

import Correlate
import DataStructures
import HelperFunctions


if "darwin" in sys.platform:
    modeldir = "/Volumes/DATADRIVE/Stellar_Models/PHOENIX/Stellar/Vband/"
elif "linux" in sys.platform:
    modeldir = "/media/FreeAgent_Drive/SyntheticSpectra/Sorted/Stellar/Vband/"
else:
    modeldir = raw_input("sys.platform not recognized. Please enter model directory below: ")
    if not modeldir.endswith("/"):
        modeldir = modeldir + "/"

# Define regions contaminated by telluric residuals or other defects. We will not use those regions in the cross-correlation
badregions = [[0, 466],
              [587.5, 593],
              [627, 634.5],
              [686, 706],
              [716, 742],
              [749.1, 749.45],
              [759, 770],
              [780, 9e9]]

#Define the values of vsini to search
vsini_values = [10, 20, 30, 40]

#Set up model list
model_list = [modeldir + "lte30-4.50+0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte30-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte30-4.50-0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte31-4.50+0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte31-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte31-4.50-0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte32-4.5-0.5.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte32-4.50+0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte32-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte33-4.50+0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte33-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte33-4.50-0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte34-4.50+0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte34-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte34-4.50-0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte35-4.50+0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte35-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte35-4.50-0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte36-4.50+0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte36-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte36-4.50-0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte37-4.5-0.5.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte37-4.50+0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte37-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte38-4.50+0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte38-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte38-4.50-0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte39-4.50+0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte39-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte39-4.50-0.5.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte40-4.5+0.5.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte40-4.5-0.5.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte40-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte41-4.5+0.5.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte41-4.5-0.0.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte41-4.5-0.5.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte42-4.5+0.5.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte42-4.5-0.5.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte42-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte43-4.5+0.5.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte43-4.5-0.5.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte43-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte44-4.5+0.5.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte44-4.5-0.5.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte44-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte45-4.5+0.5.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte45-4.5-0.5.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte45-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte46-4.5+0.5.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte46-4.5-0.5.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte46-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte47-4.5+0.5.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte47-4.5-0.5.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte47-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte48-4.5+0.5.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte48-4.5-0.5.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte48-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte49-4.5+0.5.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte49-4.5-0.5.Cond.PHOENIX2004.tab.7.sorted",
              modeldir + "lte49-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte50-4.5+0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte50-4.5-0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte50-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte51-4.5+0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte51-4.5-0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte51-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte52-4.5+0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte52-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte53-4.5+0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte53-4.5-0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte53-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte54-4.5+0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte54-4.5-0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte54-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte55-4.5+0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte55-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte56-4.5+0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte56-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte57-4.5+0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte57-4.5-0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte57-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte58-4.5+0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte58-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte59-4.5+0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte59-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte60-4.5+0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte60-4.5-0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte60-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte61-4.5+0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte61-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte62-4.5+0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte62-4.5-0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte62-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte63-4.5+0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte63-4.5-0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte63-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte64-4.5+0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte64-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte65-4.5+0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte65-4.5-0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte65-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte66-4.5+0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte66-4.5-0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte66-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte67-4.5+0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte67-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted",
              modeldir + "lte68-4.5+0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte68-4.5-0.5.Cond.PHOENIX2004.direct.7.sorted",
              modeldir + "lte68-4.50-0.0.AGS.Cond.PHOENIX-ACES-2009.HighRes.7.sorted"]

modeldict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(DataStructures.xypoint))))
processed = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(bool))))

if __name__ == "__main__":
    model_data = []
    for fname in model_list:
        if "PHOENIX2004" in fname:
            temp = int(fname.split("lte")[-1][:2]) * 100
            gravity = float(fname.split("lte")[-1][3:6])
            metallicity = float(fname.split("lte")[-1][6:10])
        elif "PHOENIX-ACES" in fname:
            temp = int(fname.split("lte")[-1][:2]) * 100
            gravity = float(fname.split("lte")[-1][3:7])
            metallicity = float(fname.split("lte")[-1][7:11])
        print "Reading in file %s" % fname
        x, y = np.loadtxt(fname, usecols=(0, 1), unpack=True)
        model = DataStructures.xypoint(x=x * units.angstrom.to(units.nm), y=10 ** y)
        for vsini in vsini_values:
            modeldict[temp][gravity][metallicity][vsini] = model
            processed[temp][gravity][metallicity][vsini] = False


def Process_Data(fname, extensions=True, trimsize=100):
    if extensions:
        orders = HelperFunctions.ReadExtensionFits(fname)

    else:
        orders = HelperFunctions.ReadFits(fname, errors=2)

    numorders = len(orders)
    for i, order in enumerate(orders[::-1]):
        DATA = interp(order.x, order.y)
        CONT = interp(order.x, order.cont)
        ERROR = interp(order.x, order.err)
        order.x = numpy.linspace(order.x[trimsize], order.x[-trimsize], order.size() - 2 * trimsize)
        order.y = DATA(order.x)
        order.cont = CONT(order.x)
        order.err = ERROR(order.x)

        #Remove bad regions from the data
        for region in badregions:
            left = numpy.searchsorted(order.x, region[0])
            right = numpy.searchsorted(order.x, region[1])
            if left > 0 and right < order.size():
                print "Warning! Bad region covers the middle of order %i" % i
                print "Removing full order!"
                left = 0
                right = order.size()
            order.x = numpy.delete(order.x, numpy.arange(left, right))
            order.y = numpy.delete(order.y, numpy.arange(left, right))
            order.cont = numpy.delete(order.cont, numpy.arange(left, right))
            order.err = numpy.delete(order.err, numpy.arange(left, right))


        #Remove whole order if it is too small
        remove = False
        if order.x.size <= 1:
            remove = True
        else:
            velrange = 3e5 * (numpy.median(order.x) - order.x[0]) / numpy.median(order.x)
            if velrange <= 1050.0:
                remove = True
        if remove:
            print "Removing order %i" % (numorders - 1 - i)
            orders.pop(numorders - 1 - i)
        else:
            # Find outliers from e.g. bad telluric line or stellar spectrum removal.
            order.cont = FittingUtilities.Continuum(order.x, order.y, lowreject=3, highreject=3)
            outliers = HelperFunctions.FindOutliers(order, expand=10, numsiglow=5, numsighigh=5)
            if len(outliers) > 0:
                order.y[outliers] = 1.0
                order.cont = FittingUtilities.Continuum(order.x, order.y, lowreject=3, highreject=3)
                order.y[outliers] = order.cont[outliers]
                #plt.plot(order.x, order.y/order.cont)
                #plt.text(order.x.mean(), 1.01,str(numorders -1 -i))

            # Save this order
            orders[numorders - 1 - i] = order.copy()
    #plt.show()
    return orders


if __name__ == "__main__":
    #Parse command line arguments:
    fileList = []
    extensions = True
    tellurics = False
    trimsize = 100
    for arg in sys.argv[1:]:
        if "-e" in arg:
            extensions = False
        if "-t" in arg:
            tellurics = True  #telluric lines modeled but not removed
        else:
            fileList.append(arg)

    #Do the cross-correlation
    for temp in sorted(modeldict.keys()):
        for gravity in sorted(modeldict[temp].keys()):
            for metallicity in sorted(modeldict[temp][gravity].keys()):
                for vsini in vsini_values:
                    for fname in fileList:
                        orders = Process_Data(fname, extensions=True)

                        output_dir = "Cross_correlations/"
                        outfilebase = fname.split(".fits")[0]
                        if "/" in fname:
                            dirs = fname.split("/")
                            output_dir = ""
                            outfilebase = dirs[-1].split(".fits")[0]
                            for directory in dirs[:-1]:
                                output_dir = output_dir + directory + "/"
                            output_dir = output_dir + "Cross_correlations/"
                        HelperFunctions.ensure_dir(output_dir)

                        model = modeldict[temp][gravity][metallicity][vsini]
                        pflag = not processed[temp][gravity][metallicity][vsini]
                        retdict = Correlate.GetCCF(orders,
                                                   model,
                                                   resolution=80000.0,
                                                   vsini=vsini,
                                                   rebin_data=True,
                                                   process_model=pflag,
                                                   debug=False,
                                                   outputdir=output_dir.split("Cross_corr")[0])
                        corr = retdict["CCF"]
                        if pflag:
                            processed[temp][gravity][metallicity][vsini] = True
                            modeldict[temp][gravity][metallicity][vsini] = retdict["model"]

                        outfilename = "%s%s.%.0fkps_%sK%+.1f%+.1f" % (
                        output_dir, outfilebase, vsini, temp, gravity, metallicity)
                        print "Outputting to ", outfilename, "\n"
                        numpy.savetxt(outfilename, numpy.transpose((corr.x, corr.y)), fmt="%.10g")


                    #Delete the model. We don't need it anymore and it just takes up ram.
                    modeldict[temp][gravity][metallicity][vsini] = []
        
        
                                   
           



