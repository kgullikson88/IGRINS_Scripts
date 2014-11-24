import sys
import os
import FittingUtilities
import warnings

import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import time, units as u
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from astropy.modeling import models, fitting

import TelluricFitter
import HelperFunctions
import Units


def RefineWavelength(orders, header, bad_orders, ap2ord, plot=True):
    humidity = header['HUMIDITY']
    T = (header['AIRTEMP'] - 32.0) * 5.0 / 9.0 + 273.15
    P = header['BARPRESS'] * Units.hPa / Units.inch_Hg
    ZD = header['ZD']
    co2 = 368.5
    n2o = 0.32
    co = 0.14
    ch4 = 1.8
    orders = sorted(orders, key=lambda o: o.x[0])
    fitter = TelluricFitter.TelluricFitter()
    print "Generating Telluric Model"
    model = fitter.Modeler.MakeModel(temperature=T,
                                     pressure=P,
                                     angle=ZD,
                                     humidity=humidity,
                                     co2=co2,
                                     n2o=n2o,
                                     co=co,
                                     ch4=ch4,
                                     lowfreq=1e7 / (orders[-1].x[-1] + 20.0),
                                     highfreq=1e7 / (orders[0].x[0] - 20.0))
    print "Rebinning model to constant wavelength spacing"
    xgrid = np.linspace(model.x[0], model.x[-1], model.size())
    model = FittingUtilities.RebinData(model, xgrid)
    print "Reducing resolution of model to R=60000"
    model = FittingUtilities.ReduceResolution2(model, 60000.0)
    print "Interpolating model"
    model_spline = spline(model.x, model.y)

    apertures = np.array([], dtype=np.float)
    ordernums = np.array([], dtype=np.float)
    pixels = np.array([], dtype=np.float)
    wavelengths = np.array([], dtype=np.float)
    weights = np.array([], dtype=np.float)
    orders_copy = [o.copy() for o in orders]
    for i, order in enumerate(orders_copy):
        if i in bad_orders:
            continue

        print "Fitting Wavelength in order {:d}".format(i + 1)
        model2 = model_spline(order.x)
        order.cont *= FittingUtilities.Iterative_SV(order.y / (order.cont * model2), 61, 3, lowreject=2, highreject=2)
        if plot:
            plt.plot(order.x, order.y / order.cont, 'k-', alpha=0.4)
            plt.plot(order.x, model2, 'r-', alpha=0.6)
        left = np.searchsorted(model.x, order.x[0] - 10.0)
        right = np.searchsorted(model.x, order.x[-1] + 10.0)
        # modelfcn, mean = fitter.FitWavelength(order.copy(), model[left:right], fitorder=7, tol=0.1, linestrength=0.9)
        modelfcn, mean = fitter.FitWavelengthNew(order.copy(), model[left:right], fitorder=4, be_safe=True)
        # xgrid = (order.x - np.median(order.x))/(order.x[-1] - order.x[0])
        order.x += modelfcn(order.x - mean)
        chisq = np.sum((order.y / order.cont - model_spline(order.x)) ** 2)
        if plot:
            plt.plot(order.x, order.y / order.cont, 'g-', alpha=0.4)

        order2 = order.copy()
        order2.y /= order2.cont
        linepixels = FittingUtilities.FindLines(order2, tol=0.9)
        linepixels = np.array([p for p in linepixels if order.y[p] / order.cont[p] > 0.05], dtype=np.int)
        waves = order.x[linepixels]
        apertures = np.hstack((apertures, np.ones(waves.size) * i))
        ordernums = np.hstack((ordernums, np.ones(waves.size) * ap2ord[i]))
        pixels = np.hstack((pixels, linepixels))
        wavelengths = np.hstack((wavelengths, waves))
        weights = np.hstack((weights, np.ones(waves.size) / chisq))

    bestoffset = 0
    bestsign = 1
    bestchisq = np.inf
    p_init = models.Chebyshev2D(5, 4, x_domain=[0, 2048])
    fit_p = fitting.LinearLSQFitter()
    print "Fitting wavelength solution for entire chip with chebyshev polynomial"
    """
    for offset in range(200):
        for sign in [-1.0, 1.0]:
            ordernums = sign * apertures + offset

            p = fit_p(p_init, pixels, ordernums, wavelengths * ordernums, weights)
            pred = p(pixels, ordernums) / ordernums

            chisq = np.sum((wavelengths - pred) ** 2)
            if chisq < bestchisq:
                bestchisq = chisq
                bestoffset = offset
                bestsign = sign

    print "best sign: ", bestsign
    print "best offset: ", bestoffset
    """

    p = fit_p(p_init, pixels, ordernums, wavelengths * ordernums, weights)
    if plot:
        # ordernums = bestsign * apertures + bestoffset
        pred = p(pixels, ordernums) / ordernums
        fig3d = plt.figure(2)
        ax3d = fig3d.add_subplot(111, projection='3d')
        ax3d.scatter3D(pixels, apertures, wavelengths - pred, 'ro')
        fig3 = plt.figure(3)
        ax = fig3.add_subplot(111)

    # Assign the wavelength solution to each order
    for i, order in enumerate(orders):
        ord = ap2ord[i]
        xgrid = np.arange(order.size(), dtype=np.float)
        ord_arr = ord * np.ones(xgrid.size, dtype=np.float)
        wave = p(xgrid, ord_arr) / ord_arr
        if plot:
            ax.plot(order.x, order.y / order.cont, 'k-', alpha=0.4)
            ax.plot(wave, order.y / order.cont, 'g-', alpha=0.4)
            ax.plot(wave, model_spline(wave), 'r-', alpha=0.6)
        orders[i].x = wave

    if plot:
        plt.show()

    return orders


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
    lownum = int(filename.split(".spec")[0][-4:])
    highnum = lownum + maxnods - 1

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
        if i == 0:
            objname = header['OBJECT']
            ra = header['RATEL']
            dec = header['DECTEL']
        elif header['OBJECT'] != objname:
            # We have hit a new target.
            print "New object name ({}). Expected {}".format(header['object'], objname)
            print "Not using any subsequent files..."
            break
        t = header['DATE-OBS']
        t = "{:s}T{:s}".format(t[:10], t[11:])
        date = t.partition('T')[0]
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
    ZD = np.mean(zenith_angle)
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
    print ap2ord

    # Wavelength calibration
    #orders = RefineWavelength(orders, pri_hdu.header, bad_orders, ap2ord, plot=plot)

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

    # Wavelength calibration
    #orders = RefineWavelength(orders, pri_hdu.header, bad_orders, ap2ord, plot=plot)

    #Wavelength calibration
    # orders = RefineWavelength(orders, pri_hdu.header, bad_orders, plot=plot)

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
    file_list = sorted([f for f in os.listdir("./") if f.startswith("SDCH") and f.endswith("spec.fits")])
    num_list = [int(f.split("_")[-1][:4]) for f in file_list]
    nod_list = [min(maxnods, num_list[i + 1] - num_list[i]) for i in range(len(num_list) - 1)]
    # nod_list.append(maxnods)

    for filename, nods in zip(file_list, nod_list):
        print filename, nods
        Convert(filename, nods, overwrite=overwrite)










