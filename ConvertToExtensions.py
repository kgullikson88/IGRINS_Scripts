import sys
import os
import FittingUtilities

import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import time, units as u
import numpy as np

import HelperFunctions


def ReadFile(filename, columns, blazefile="H_BLAZE.DAT", skip=0):
    orders = HelperFunctions.ReadFits(filename)
    blazes = np.loadtxt(blazefile)
    header_info = []
    for i, order in enumerate(orders):
        # Cut off the edges, where the blaze is really bad
        order = order[150:-100]

        # Convert to air wavelengths!
        wave_A = order.x * u.nm.to(u.angstrom)
        n = 1.0 + 2.735182e-4 + 131.4182 / wave_A ** 2 + 2.76249e8 / wave_A ** 4
        order.x /= n

        # Divide by blaze function
        blaze = blazes[:, i + skip][150:-100]
        blaze[blaze < 1e-3] = 1e-3
        order.y /= blaze
        order.err /= blaze

        offset = 0.0
        order.cont = FittingUtilities.Continuum(order.x, order.y + offset, fitorder=5, lowreject=1.5,
                                                highreject=2) - offset

        plt.plot(order.x, order.y, 'k-', alpha=0.4)
        plt.plot(order.x, order.cont, 'r-', alpha=0.4)

        column = {"wavelength": order.x,
                  "flux": order.y,
                  "continuum": order.cont,
                  "error": order.err}
        columns.append(column)
    plt.show()
    return columns


if __name__ == "__main__":
    # Read command line arguments
    lownum = None
    highnum = None
    for arg in sys.argv[1:]:
        if "num" in arg:
            r = arg.partition("=")[-1].split("-")
            lownum = int(r[0])
            highnum = int(r[1])
        else:
            filename = arg

    if lownum is None or highnum is None:
        lownum = int(filename.split(".spec")[0][-4:])
        highnum = lownum + 3

    print "Converting file {:s}".format(filename)
    if "SDCH" in filename:
        blazefile = "H_BLAZE.DAT"
        skip = 2
    elif "SDCK" in filename:
        blazefile = "K_BLAZE.DAT"
        skip = 0
    columns = []
    columns = ReadFile(filename, columns, blazefile, skip)

    # Get the other band
    if "SDCH" in filename:
        filename = filename.replace("SDCH", "SDCK")
        blazefile = "K_BLAZE.DAT"
        skip = 0
    elif "SDCK" in filename:
        filename = filename.replace("SDCK", "SDCH")
        blazefile = "H_BLAZE.DAT"
        skip = 2
    columns = ReadFile(filename, columns, blazefile, skip)

    # Sort the columns by increasing wavelength
    columns = sorted(columns, key=lambda c: c['wavelength'][0])

    #plt.show()

    # Prepare fits keywords
    date_obs = []
    zenith_angle = []
    humidity = []
    temperature = []
    pressure = []
    original_fnames = sorted(
        ["20140708/SDCH_20140708_{:s}.fits".format(str(i).zfill(4)) for i in range(lownum, highnum + 1)])
    for i, fname in enumerate(original_fnames):
        print "Reading header info from {:s}".format(fname)
        #Get the observation time
        header = fits.getheader(fname)
        t = header['DATE-OBS']
        t = "{:s}T{:s}".format(t[:10], t[11:])
        date = t.partition('T')[0]
        #print t
        t_jd = time.Time(t, format='isot', scale='utc').jd
        if i == len(original_fnames) - 1:
            exptime = float(header['EXPTIME']) * u.s.to(u.day)
            t_jd += exptime

            objname = header['OBJECT']
        date_obs.append(t_jd)

        #Get the ambient pressure, temperature, and humidity
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

    #Figure out the average values for each of the quantities
    ZD = np.mean(zenith_angle)
    RH = np.mean(humidity)
    T = np.mean(temperature)
    P = np.mean(pressure)
    t_jd = (min(date_obs) + max(date_obs)) / 2.0




    #Output as a fits file with binary tables
    outfilename = "{:s}.fits".format(objname.replace(" ", "_"))
    print "Outputting to {:s}".format(outfilename)
    pri_hdu = fits.PrimaryHDU()
    pri_hdu.header['OBJECT'] = objname
    pri_hdu.header['DATE-OBS'] = time.Time(t_jd, format='jd').isot
    pri_hdu.header['UT'] = pri_hdu.header['DATE-OBS'].split("T")[-1]
    pri_hdu.header['HUMIDITY'] = RH
    pri_hdu.header['AIRTEMP'] = T
    pri_hdu.header['BARPRESS'] = P
    pri_hdu.header['ZD'] = ZD
    hdulist = fits.HDUList([pri_hdu, ])

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









