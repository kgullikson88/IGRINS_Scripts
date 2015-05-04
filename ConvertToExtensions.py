import sys
import os
import warnings

import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy import time
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from astropy.modeling import models, fitting
import TelluricFitter

import FittingUtilities
import HelperFunctions
import Units
import parse_igrins_log


def ReadFile(filename, blazefile="H_BLAZE.DAT", skip=0):
    orders, wavefields = HelperFunctions.ReadFits(filename, return_aps=True)
    try:
        errors = HelperFunctions.ReadFits(filename.replace("spec", "variance"))
    except IOError:
        warnings.warn("No Variance file found. Falling back to using the sqrt of the flux")
        errors = []
        for order in orders:
            error = order.copy()
            error.y = np.sqrt(order.y)
            errors.append(error.copy())
    blazes = np.loadtxt(blazefile)
    header_info = []
    ret_orders = []
    ap2ord = {}  #Dictionary giving correspondence from aperture number to echelle order number
    for i, order in enumerate(orders):
        # Remove nans, that the pipeline puts in there
        goodindices = np.where(-np.isnan(order.y))[0]
        goodindices = goodindices[(goodindices > 150) & (goodindices < order.size() - 100)]
        order = order[goodindices]
        order.err = errors[i].y[goodindices]

        # Convert to air wavelengths!
        wave_A = order.x * u.nm.to(u.angstrom)
        n = 1.0 + 2.735182e-4 + 131.4182 / wave_A ** 2 + 2.76249e8 / wave_A ** 4
        order.x /= n

        # Divide by blaze function
        blaze = blazes[:, i][goodindices]
        blaze[blaze < 1e-3] = 1e-3
        order.y /= blaze
        order.err /= blaze

        # Fit continuum
        order.cont = FittingUtilities.Continuum(order.x, order.y, fitorder=5, lowreject=1.5, highreject=2)

        # Get the echelle order number from the wavefields
        order_number = int(wavefields[i][1])

        ret_orders.append(order.copy())
        ap2ord[i] = order_number
    return ret_orders, ap2ord


def Convert(filename, maxnods, overwrite=False):
    # Get the RA, DEC, and airmass from the IGRINS logfile
    lownum = int(filename.split(".spec")[0][-4:])
    highnum = lownum + maxnods - 1
    nums = range(lownum, highnum + 1)
    logfiles = parse_igrins_log.get_logfilenames(os.getcwd().split("/")[-1])
    log_data = parse_igrins_log.read_logfile(logfiles)
    ra = parse_igrins_log.dex_to_hex(parse_igrins_log.get_average(log_data, nums, 'RA'))
    dec = parse_igrins_log.dex_to_hex(parse_igrins_log.get_average(log_data, nums, 'DEC'))
    try:
        ZD = np.arccos(1.0 / parse_igrins_log.get_average(log_data, nums, 'AM')) * 180.0 / np.pi
    except TypeError:
        ZD = 30.0  # Guess!
    total_exptime = 0.0
    print nums

    # Get some more info from the fits header of each original IGRINS file
    date_obs = []
    zenith_angle = []
    humidity = []
    temperature = []
    pressure = []
    date = os.getcwd().split("/")[-1]
    original_fnames = sorted(
            ["{0:s}/SDCH_{0:s}_{1:s}.fits".format(date, str(i).zfill(4)) for i in range(lownum, highnum + 1)])
    for i, fname in enumerate(original_fnames):
        print "Reading header info from {:s}".format(fname)
        # Get the observation time
        header = fits.getheader(fname)
        total_exptime += float(header['EXPTIME'])
        if i == 0:
            objname = header['OBJECT']
        elif header['OBJECT'] != objname:
            # We have hit a new target.
            print "New object name ({}). Expected {}".format(header['object'], objname)
            print "Not using any subsequent files..."
            break

        # Use these lines for new IGRINS data
        t = header['DATE-OBS']
        t = "{:s}T{:s}".format(t[:10], t[11:])
        date = t.partition('T')[0]
        # Use these lines for old (earlier than July 2014 IGRINS data
        # date = header['DATE-OBS']
        #timestr = header['TIME-OBS']
        #t = '{:s}T{:s}'.format(date, timestr)
        print t

        # Back to general things
        # print t
        t_jd = time.Time(t, format='isot', scale='utc').jd
        if i == len(original_fnames) - 1:
            exptime = float(header['EXPTIME']) * u.s.to(u.day)
            t_jd += exptime

            # objname = header['OBJECT']
        date_obs.append(t_jd)

        # Get the ambient pressure, temperature, and humidity
        #Read in weather information (archived data is downloaded from weather.as.utexas.edu)
        homedir = os.environ["HOME"]
        weather_file = homedir + "/School/Research/Useful_Datafiles/Weather.dat"
        infile = open(weather_file)
        lines = infile.readlines()
        infile.close()
        times = []
        RH = []
        P = []
        T = []
        idx = 0
        bestindex = 0
        difference = 9e9
        for line in lines[1:]:
            segments = line.split()
            if date in segments[0]:
                segments = segments[1].split(",")
                t = segments[0]
                weather_time = time.Time("{:s}T{:s}".format(date, t), format='isot', scale='utc').jd
                #t_seg = t.split(":")
                #weather_time = 3600*float(t_seg[0]) + 60*float(t_seg[1]) + float(t_seg[2])
                if np.abs(t_jd - weather_time) < difference:
                    difference = np.abs(t_jd - weather_time)
                    bestindex = idx
                times.append(segments[0])
                T.append(float(segments[3]))
                RH.append(float(segments[4]))
                P.append(float(segments[5]))
                idx += 1
        humidity.append(RH[bestindex])
        temperature.append(T[bestindex])
        pressure.append(P[bestindex])

        # zenith_angle.append(90.0 - float(header['ALT']))
        zenith_angle.append(np.arccos(1.0 / float(header['amstart'])) * 180.0 / np.pi)

    # Figure out the average values for each of the quantities
    print ZD, np.mean(zenith_angle)
    # sys.exit()
    #ZD = np.mean(zenith_angle)
    RH = np.mean(humidity)
    T = np.mean(temperature)
    P = np.mean(pressure)
    t_jd = (min(date_obs) + max(date_obs)) / 2.0

    # Prepare the fits header
    outfilename = "{:s}.fits".format(objname.replace(" ", "_"))
    if not overwrite:
        done = False
        i = 1
        while not done:
            if outfilename in os.listdir("./"):
                outfilename = "{:s}_{:d}.fits".format(objname.replace(" ", "_"), i)
                i += 1
            else:
                done = True
    print "Outputting to {:s}".format(outfilename)
    pri_hdu = fits.PrimaryHDU()
    pri_hdu.header['OBJECT'] = objname
    pri_hdu.header['EXPTIME'] = total_exptime
    pri_hdu.header['RA'] = ra
    pri_hdu.header['DEC'] = dec
    pri_hdu.header['JD'] = t_jd
    pri_hdu.header['DATE-OBS'] = time.Time(t_jd, format='jd').isot
    pri_hdu.header['UT'] = pri_hdu.header['DATE-OBS'].split("T")[-1]
    pri_hdu.header['HUMIDITY'] = RH
    pri_hdu.header['AIRTEMP'] = T
    pri_hdu.header['BARPRESS'] = P
    pri_hdu.header['ZD'] = ZD
    for i, fname in enumerate(original_fnames):
        pri_hdu.header['NODFILE{}'.format(i + 1)] = fname
    hdulist = fits.HDUList([pri_hdu, ])

    """
    Now, read in the extracted files!
    """

    print "Converting file {:s}".format(filename)
    if "SDCH" in filename:
        blazefile = "H_BLAZE.DAT"
        skip = 2
        bad_orders = [0, 21, 22]
    elif "SDCK" in filename:
        blazefile = "K_BLAZE.DAT"
        skip = 0
        bad_orders = [0, 1, 3, 19, 20]
    orders, ap2ord = ReadFile(filename, blazefile, skip)

    # Convert to columns
    columns = []
    for order in orders:
        column = {"wavelength": order.x,
                  "flux": order.y,
                  "continuum": order.cont,
                  "error": order.err}
        columns.append(column)

    # Get the other band
    if "SDCH" in filename:
        filename = filename.replace("SDCH", "SDCK")
        blazefile = "K_BLAZE.DAT"
        skip = 0
        bad_orders = [0, 1, 3, 19, 20]
    elif "SDCK" in filename:
        filename = filename.replace("SDCK", "SDCH")
        blazefile = "H_BLAZE.DAT"
        skip = 2
        bad_orders = [0, 21, 22]
    orders, ap2ord = ReadFile(filename, blazefile, skip)

    # Convert to columns
    for order in orders:
        column = {"wavelength": order.x,
                  "flux": order.y,
                  "continuum": order.cont,
                  "error": order.err}
        columns.append(column)

    # Sort the columns by increasing wavelength
    columns = sorted(columns, key=lambda c: c['wavelength'][0])

    # Output to a fits file
    for i in range(len(columns)):
        column_dict = columns[i]
        cols = []
        for key in column_dict.keys():
            cols.append(fits.Column(name=key, format="D", array=column_dict[key]))
        cols = fits.ColDefs(cols)
        tablehdu = fits.BinTableHDU.from_columns(cols)
        tablehdu.header.name = "Order {:d}".format(i + 1)
        hdulist.append(tablehdu)

    hdulist.writeto(outfilename, clobber=True, output_verify='ignore')
    hdulist.close()




if __name__ == "__main__":
    # Read command line arguments
    lownum = None
    highnum = None
    plot = False
    overwrite = True
    maxnods = 8
    file_list = []
    for arg in sys.argv[1:]:
        if "num-nods" in arg:
            maxnods = int(arg.split("=")[-1])
        elif "-p" in arg:
            plot = True
        elif "-no-overwrite" in arg:
            overwrite = False
        else:
            file_list.append(arg)

    # If file_list is empty, do all using the files in the current directory
    if len(file_list) == 0:
        file_list = sorted([f for f in os.listdir("./") if f.startswith("SDCH") and f.endswith("spec.fits")])
    num_list = [int(f.split("_")[-1][:4]) for f in file_list]
    nod_list = [min(maxnods, num_list[i + 1] - num_list[i]) for i in range(len(num_list) - 1)]
    nod_list.append(maxnods)

    for filename, nods in zip(file_list, nod_list):
        print filename, nods
        Convert(filename, nods, overwrite=overwrite)










