import HelperFunctions
import sys
import os
import numpy
import matplotlib.pyplot as plt
from astropy.io import fits
import DataStructures
import FittingUtilities
from astropy import wcs
from itertools import combinations
from scipy import optimize
from scipy.interpolate import InterpolatedUnivariateSpline as spline
import time
from collections import defaultdict

fitorder = 3

def main():
  # Read in the logfile to get the original data filenames
  original_files = defaultdict( lambda : defaultdict( lambda : defaultdict(list) ) )
  extracted_files = defaultdict( lambda : defaultdict( lambda: defaultdict( lambda : defaultdict(str) ) ) )
  for band in ["H", "K"]:
    logfilename = [f for f in os.listdir(band) if "IGRINS_DT_Log"  in f and f.endswith("%s.txt" %band)][0]
    infile = open("%s/%s" %(band, logfilename))
    lines = infile.readlines()
    infile.close()

    for line in lines[2:]:
      segments = line.split(",")
      groupnum = int(segments[2])
      objtype = segments[5]
      if objtype.lower() == "tar" or objtype.lower() == "std":
        original_files[groupnum][objtype][band].append("%s/%s" %(band, segments[0]))

        #Find the extracted files
        ex = [f for f in os.listdir(band) if objtype in f and "G%i" %groupnum in f and f.endswith(".tpn.fits")]
        for fname in ex:
          ordernum = int(fname.split(".tr.")[-1].split(".tpn")[0])
          if ordernum > 1:
            extracted_files[groupnum][objtype][ordernum][band] = "%s/%s" %(band, fname)

  # Read in the data, and assign a wavelength grid to it
  for group in extracted_files.keys():
    for objtype in extracted_files[group].keys():
      column_list = []
      header_list = []
      for ordernum in sorted(extracted_files[group][objtype].keys()):
        for band in extracted_files[group][objtype][ordernum].keys():
          fname = extracted_files[group][objtype][ordernum][band]
          hdulist = fits.open(fname)
          header = hdulist[0].header
          order = HelperFunctions.ReadFits(fname)[0]

          #Check for a wave file
          wavefile = "%s/School/Research/Useful_Datafiles/IGRINS_Datafiles/SDC%s_order%.3i.wavedata" %(os.environ["HOME"], band, ordernum)
          if os.path.isfile(wavefile):
            x,y = numpy.loadtxt(wavefile, unpack=True)
            fit = numpy.poly1d(numpy.polyfit(x,y, fitorder))
            order.x = fit(order.x)
          else:
            print "Wavelength file not found!"
            a, b = header['CRVAL1'], header['CDELT1']
            order.x = a + b*order.x

          #Fit the continuum
          order.cont = FittingUtilities.Continuum(order.x, order.y, lowreject=2, highreject=2, fitorder=5)
          order.cont[order.cont < 1] = 1.0

          #Prepare output
          column = {"wavelength": order.x,
                    "flux": order.y,
                    "continuum": order.cont,
                    "error": order.err}
          hinfo = [["AP-NUM", header["AP-NUM"]],
                   ["AP-COEF1", header["AP-COEF1"]],
                   ["AP-COEF2", header["AP-COEF2"]],
                   ["AP-WIDTH", header["AP-WIDTH"]],
                   ["BAND", band],
                   ["FNAME", fname, 'original filename']]
          column_list.append(column)
          header_list.append(hinfo)

      outfilename = "%s_G%i.fits" %(objtype, group)
      print "Outputting to %s" %outfilename
      sorter = [i[0] for i in sorted(enumerate(column_list), key = lambda c: c[1]["wavelength"][0])]
      clist = []
      hlist = []
      for i in sorter:
        clist.append(column_list[i])
        hlist.append(header_list[i])
      HelperFunctions.OutputFitsFileExtensions(clist, fname, outfilename, mode="new", headers_info=hlist)


  


if __name__ == "__main__":
  main()
