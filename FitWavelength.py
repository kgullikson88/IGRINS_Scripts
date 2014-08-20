import sys
import itertools

import numpy as np

import HelperFunctions
import FittingUtilities
from astropy.io import fits
import astropy.units as u
import matplotlib.pyplot as plt
import statsmodels.api as sm
from mpl_toolkits.mplot3d import Axes3D


def FitWavelength(data, model, fitorder=3):
    data.y /= data.cont
    xgrid = np.arange(data.size())
    data_lines = FittingUtilities.FindLines(data, tol=0.98, linespacing=0.2e-3)
    data_flux = data.y[data_lines]
    # data_flux = (data_flux - max(data_flux))/(max(data_flux) - min(data_flux)) + 1

    model.y /= model.cont
    model_lines = FittingUtilities.FindLines(model, tol=0.98, linespacing=0.2e-3)
    model_flux = model.y[model_lines]
    model_lines = model.x[model_lines]

    fig = plt.figure(2, figsize=(10, 7))
    top = fig.add_subplot(211)
    top.plot(data.x, data.y, 'k-')
    top.plot(data.x[data_lines], data_flux, 'ro')
    top.grid(True)
    top.set_ylim((max(-0.1, np.min(data.y)), 1.1))
    #top.set_xbound(lower=leftlim, upper=rightlim)

    bottom = fig.add_subplot(212, sharex=top)
    bottom.plot(model.x, model.y, 'k-')
    bottom.plot(model_lines, model_flux, 'ro')
    bottom.grid(True)

    EH = EventHandler(top, bottom, data_lines, data.x[data_lines], model_lines)
    fig.canvas.mpl_connect('button_press_event', EH.OnClick)

    plt.show()

    print EH.pixels
    print EH.waves

    if len(EH.pixels) == len(EH.waves) and len(EH.pixels) >= fitorder:
        fit = np.poly1d(np.polyfit(EH.pixels, EH.waves, fitorder))
        xaxis = fit(xgrid)
    else:
        xaxis = data.x
    data.y *= data.cont
    return xaxis, EH.pixels, EH.waves


class EventHandler():
    def __init__(self, data_ax, model_ax, data_lines, data_waves, model_lines):
        self.data_ax = data_ax
        self.model_ax = model_ax
        self.pixels = []
        self.waves = []
        self.data_lines = np.array(data_lines)
        self.data_waves = np.array(data_waves)
        self.model_lines = np.array(model_lines)


    def OnClick(self, event):
        event.inaxes.plot(event.xdata, event.ydata, 'gx', markersize=10)
        if event.inaxes == self.data_ax:
            # Find the closest data line to this event's xdata
            idx = np.argmin(np.abs(self.data_waves - event.xdata))
            self.pixels.append(self.data_lines[idx])
            print "Using line at pixel %i" % self.data_lines[idx]
        elif event.inaxes == self.model_ax:
            # Find the closest model line to this events xdata
            idx = np.argmin(np.abs(self.model_lines - event.xdata))
            self.waves.append(self.model_lines[idx])
            print "Using line at wavelength %g um" % self.model_lines[idx]


def CCImprove_Wrapper(data, model, tol=1e-3):
    # Make the data evenly-sampled
    xgrid = np.linspace(data.x[0], data.x[-1], data.size())
    data = FittingUtilities.RebinData(data, xgrid)
    dx = xgrid[1] - xgrid[0]

    #Sample the model with the same spacing
    xgrid = np.arange(model.x[0], model.x[-1], dx)
    model = FittingUtilities.RebinData(model, xgrid)
    shift = FittingUtilities.CCImprove(data, model, tol=tol, be_safe=True)
    return shift


def FindBestShift(data, model, segmentsize):
    bestchisq = np.inf
    segmentsize = 2*model.size()
    left = 0
    right = model.size()
    shift = 0
    # plt.figure(4)
    while right < data.size() - 1:
        seg = data[left:right]
        cont = FittingUtilities.Continuum(seg.x, seg.y, fitorder=2, lowreject=1.5)
        seg.y /= cont
        s = np.median(np.log(model.y) / np.log(seg.y))
        seg.y = seg.y ** s
        chisq = np.sum((seg.y - model.y) ** 2)
        #print left, chisq
        #plt.plot(left, chisq, 'ro')
        if chisq < bestchisq:
            bestchisq = chisq
            shift = left
        left += 1
        right += 1
    #plt.show()
    return shift, 1.0 / bestchisq


def FindGoodLines(data, model, oversampling=1, label=1):
    """
      This algorithm finds lines in the noiseless telluric model,
      and then searches for the same lines in the data
    """
    # First, shift the data by a constant amount to get things close
    shift = CCImprove_Wrapper(data, model, tol=2e-3)
    data.x += shift
    xgrid = np.linspace(data.x[0], data.x[-1], data.size() * oversampling)
    data = FittingUtilities.RebinData(data, xgrid)

    data.y /= data.cont
    pixelgrid = np.arange(data.size())
    model.y /= model.cont
    model_pixels = FittingUtilities.FindLines(model, tol=0.98, linespacing=1e-4)
    model_flux = model.y[model_pixels]
    model_lines = model.x[model_pixels]

    #Remove saturated lines
    model_lines = model_lines[model_flux > 0.05]
    model_flux = model_flux[model_flux > 0.05]

    #Plot
    fig = plt.figure(1)
    ax1 = fig.add_subplot(111)


    #Go through the model lines. We do not want saturated lines!
    pixels = []
    waves = []
    weights = []
    for i, line in enumerate(model_lines):

        #Search for a nearby data line
        left = np.searchsorted(data.x, line - 2e-4)
        right = np.searchsorted(data.x, line + 2e-4)
        if right < 5 or left > data.size() - 5:
          continue
        data_lines = FittingUtilities.FindLines(data[left:right], tol=0.99, linespacing=1e-4)
        print line
        print data_lines
        if len(data_lines) < 1:
          continue

        #Use the strongest line in the data
        data_line = data_lines[0]
        for l in data_lines:
          #if data.y[left+l] < data.y[data_line+left]:
          if abs(data.x[l+left]  - line) < abs(data.x[data_line+left] - line):
            data_line = l


        #Check that the data line is not saturated
        if data.y[data_line+left] < 0.05:
          continue


        #Save in arrays
        waves.append(line)
        pixels.append(data_line+left)
        weights.append(1.0)
        """
        #Find the best shift to apply to a small section of the data to line up with the model
        idx = np.argmin(abs(data.x - line))
        segmentsize = 30 * oversampling
        searchsize = 60 * oversampling
        if idx < searchsize or idx > data.size() - searchsize:
            continue
        shift, w = FindBestShift(data[idx - searchsize:idx + searchsize],
                                 FittingUtilities.RebinData(model, data.x[idx - segmentsize:idx + segmentsize]),
                                 segmentsize)

        #Save in the arrays
        pixels.append(idx - ((searchsize - segmentsize) - shift))
        waves.append(line)
        weights.append(w)
        """

    pixels = np.array(pixels)
    waves = np.array(waves)
    weights = np.array(weights)
    #fit = np.poly1d(np.polyfit(pixels, waves, 3, w=weights/np.sum(weights)))
    #err = waves - fit(pixels)
    #goodindices = np.where(abs(err) < 1e-4)[0]
    #pixels = pixels[goodindices]
    #waves = waves[goodindices]
    #weights = weights[goodindices]
    #fit = np.poly1d(np.polyfit(pixels, waves, 3, w=weights/np.sum(weights)))
    pars = polyfit(pixels, waves, order=3)
    fit = lambda x: polyval(x, pars)

    fig = plt.figure(2)
    ax2 = fig.add_subplot(211)
    #ax2.errorbar(pixels, waves, yerr=1.0/weights)
    ax2.scatter(pixels, waves)
    ax2.plot(pixels, fit(pixels))
    ax2.set_xlabel("Pixel number")
    ax2.set_ylabel("Wavelength")
    ax3 = fig.add_subplot(212)
    ax3.plot(pixels, waves - fit(pixels), 'ro')

    ax1.plot(model.x, model.y, 'k-')
    ax1.plot(data.x, data.y, 'r-')

    for p, w in zip(pixels, waves):
        idx = np.argmin(np.abs(model.x - w))
        dx = w - data.x[p]
        dy = model.y[idx] - data.y[p]
        #print dx, dy
        ax1.arrow(data.x[p], data.y[p], dx, dy, width=1e-5)
    #for w in waves:
    #    idx = np.argmin(abs(model_lines - w))
    #    if abs(model_lines[idx] - w) < 1e-5:
    #        ax1.plot(model_lines[idx], model_flux[idx], 'ro')
    #ax1.plot(model_lines, model_flux, 'ro')
    ax1.text(np.median(data.x), 1.1, str(label))
    ax1.grid(True)
    #ax1.plot(fit(pixelgrid), data.y, 'g-')
    #plt.show()
    weights = np.abs(1.0 / (np.array(waves) - fit(np.array(pixels))))

    goodindices = np.where(abs(waves - fit(pixels)) < 3e-5)[0]
    waves = waves[goodindices]
    pixels = pixels[goodindices]
    weights = weights[goodindices]

    return waves, pixels, weights


def polyfit(x, y, order=3):
    X = np.ones_like(x)
    for p in range(1, order + 1):
        X = np.column_stack((X, x ** p))
    model = sm.RLM(y, X, M=sm.robust.norms.HuberT())
    result = model.fit()
    return result.params


def polyval(x, pars, order=3):
    y = np.zeros_like(x) + pars[0]
    par_idx = 1
    for p in range(1, order + 1):
        y += pars[par_idx] * x ** p
        par_idx += 1
    return y


def polyfit2d(x, y, z, w, xorder=3, yorder=3):
    # Normalize the x/y arrays
    x = x / 2048.
    y = y / 2048.

    # Prepare the array for statsmodels
    X = np.ones(x.size)
    for p in range(1, xorder + 1):
        X = np.column_stack((X, x ** p))
    for p in range(1, yorder + 1):
        X = np.column_stack((X, y ** p))
    #for p in range(1, min(xorder, yorder)+1):
    #  X = np.column_stack((X, (x+y)**p))
    #X = np.column_stack((X, x*y, x**2*y, x*y**2))
    #X = sm.add_constant(X)

    # Run Statsmodels.OLS
    model = sm.OLS(z, X)
    #model = sm.WLS(z,X, weights=w/np.sum(w))
    #model = sm.RLM(z,X,M=sm.robust.norms.HuberT())
    results = model.fit()
    return results.params


def polyval2d(x, y, pars, xorder=3, yorder=3):
    # Normalize the x/y arrays
    x = x / 2048.
    y = y / 2048.

    z = np.zeros_like(x) + pars[0]
    par_idx = 1
    print "z = %.5g" % pars[0]
    for p in range(1, xorder + 1):
        print "\t+ %.5g * x^%i" % (pars[par_idx], p)
        z += pars[par_idx] * x ** p
        par_idx += 1
    for p in range(1, yorder + 1):
        print "\t+ %.5g * y^%i" % (pars[par_idx], p)
        z += pars[par_idx] * y ** p
        par_idx += 1
    #for p in range(1, min(xorder, yorder)+1):
    #  z += pars[par_idx] * (x+y)**p
    #  par_idx += 1
    #z += pars[par_idx]*x*y + pars[par_idx+1]*x**2*y + pars[par_idx+2]*x*y**2
    return z


class Chebyshev2d:
    def __init__(self, x, ap, lam, offset, slope, xorder=3, yorder=3):
        self.xmax = max(x)
        self.xmin = min(x)
        o = ap * slope + offset
        self.slope = slope
        self.offset = offset
        self.omax = max(o)
        self.omin = min(o)
        self.xorder = xorder
        self.yorder = yorder

        # Prepare the array for statsmodels
        for xord in range(xorder + 1):
            for yord in range(yorder + 1):
                if xord == 0 and yord == 0:
                    X = self.P_n(x, xord, self.xmax, self.xmin) * self.P_n(o, yord, self.omax, self.omin)
                else:
                    X = np.column_stack(
                        (X, self.P_n(x, xord, self.xmax, self.xmin) * self.P_n(o, yord, self.omax, self.omin)))

        # Run Statsmodels.OLS
        z = lam * o
        model = sm.OLS(z, X)
        #model = sm.WLS(z,X, weights=w/np.sum(w))
        #model = sm.RLM(z,X,M=sm.robust.norms.HuberT())
        results = model.fit()
        self.fitcoeffs = results.params

    def P_n(self, x, n, xmax, xmin):
        xnorm = (2 * x - (xmax + xmin)) / (xmax - xmin)
        if n == 0:
            return np.ones(x.size)
        elif n == 1:
            return xnorm
        else:
            return 2.0 * xnorm * self.P_n(x, n - 1, xmax, xmin) - self.P_n(x, n - 2, xmax, xmin)

    def __call__(self, x, ap):
        o = ap * self.slope + self.offset
        z = np.zeros(x.size)
        par_idx = 0
        for xord in range(self.xorder + 1):
            for yord in range(self.yorder + 1):
                Pm = self.P_n(x, xord, self.xmax, self.xmin)
                Pn = self.P_n(o, yord, self.omax, self.omin)
                Cmn = self.fitcoeffs[par_idx]
                z += Cmn * Pm * Pn
                par_idx += 1
        return z / o


def MapToCCD(pixels, header, sparepixels=5):
    """
    Function to map the pixel coordinate in extracted data to
      the (x,y) coordinate in the original ccd chip
    """
    apnum = header["AP-NUM"]
    coeffs = [float(f) for f in header['AP-COEF1'].split(",")]
    apwidth = header["AP-WIDTH"]
    if "H" in header['BAND']:
        skippixels = 300
    else:
        skippixels = 100
    xpos = pixels + skippixels
    ypos = np.polynomial.polynomial.polyval(xpos, coeffs) - sparepixels + apwidth / 2.0
    # xpos = 2048 - xpos
    #ypos = 2048 - ypos
    return xpos, ypos  # + np.random.random()*5 - 2.5


def main():
    import os
    import DataStructures

    good_orders = {"H": [1, 2, 3, 4, 9, 11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 21],
                   "K": []}
    telfile = "%s/School/Research/Useful_Datafiles/Telluric_NearIR.dat" % (os.environ["HOME"])
    x, y, c, e = np.loadtxt(telfile, unpack=True)
    telluric = DataStructures.xypoint(x=x * u.nm.to(u.micron), y=y, cont=c, err=e)
    telluric = FittingUtilities.ReduceResolution2(telluric, 40000.0)

    for fname in sys.argv[1:]:
        hdulist = fits.open(fname)
        orders = HelperFunctions.ReadExtensionFits(fname)
        xpixels = []
        ypixels = []
        waves = []
        weights = []
        colors = itertools.cycle(["r", "b", "g"])
        for band in ["H", "K"]:
            fig3d = plt.figure(6)
            ax3d = fig3d.add_subplot(111, projection='3d')

            logfile = open("PixelFit.txt", "w")
            for i, order in enumerate(orders):
                header = hdulist[i + 1].header
                if header['BAND'].strip() != band or i + 1 not in good_orders[band]:
                    continue
                left = np.searchsorted(telluric.x, order.x[0] - 2e-3)
                right = np.searchsorted(telluric.x, order.x[-1] + 2e-3)
                model = telluric[left:right]

                order.cont = FittingUtilities.Continuum(order.x, order.y, lowreject=2, highreject=3, fitorder=2)

                print "Fitting wavelength for order %i" % i
                wavelengths, pixels, weight = FindGoodLines(order, model, label=i + 1)
                xx, yy = MapToCCD(pixels, header)
                # ax3d.scatter3D(xx, yy, wavelengths, color=next(colors))
                for x, y, lam, w in zip(xx, yy, wavelengths, weight):
                    xpixels.append(x)
                    ypixels.append(y)
                    waves.append(lam)
                    weights.append(w)
                logfile.write("\n\nAperture %i\n===============\n" % header['AP-NUM'])
                print "\tSizes = ", pixels.size, wavelengths.size, xx.size, yy.size
                np.savetxt(logfile, np.transpose((pixels, wavelengths, xx, yy)))
            logfile.close()


            # Now: Fit the 2-d surface using xpixels, ypixels, waves, and weights
            xpixels = np.array(xpixels)
            ypixels = np.array(ypixels)
            weights = np.array(weights)
            np.savetxt("WaveFit.dat", np.transpose((xpixels, ypixels, weights, waves)))

            xord = 4
            yord = 4
            xpixels, ypixels, weights, waves = np.loadtxt("WaveFit.dat", unpack=True)
            #xpixels = 2048 - xpixels
            #ypixels = 2048 - ypixels
            coeffs = polyfit2d(xpixels, ypixels, waves, weights, xorder=xord, yorder=yord)

            #After that: Re-loop over the orders, and transform the wavelength
            #  surface to the 1d array. Plot to see how it looks
            fig = plt.figure(5)
            ax = fig.add_subplot(111)

            ax3d.scatter3D(xpixels, ypixels, waves - polyval2d(xpixels, ypixels, coeffs, xorder=xord, yorder=yord))
            xx, yy = np.meshgrid(np.linspace(0, 2048, 50),
                                 np.linspace(0, 2048, 50))
            zz = polyval2d(xx, yy, coeffs, xorder=xord, yorder=yord)
            ax3d.set_xlabel("X Pixel")
            ax3d.set_ylabel("Y Pixel")
            #ax3d.set_xlim((0,1600))
            #ax3d.plot_surface(xx, yy, zz, alpha=0.2)


            for i, order in enumerate(orders):
                header = hdulist[i + 1].header
                if header['BAND'].strip() != band:
                    continue
                model = FittingUtilities.RebinData(telluric, order.x)
                pixels = np.arange(order.size())
                xx, yy = MapToCCD(pixels, header)
                wave = polyval2d(xx, yy, coeffs, xorder=xord, yorder=yord)
                ax.plot(order.x, order.y / order.cont, 'r-')
                ax.plot(wave, order.y / order.cont, 'g-')
                ax.plot(model.x, model.y, 'k-')
            plt.show()
            sys.exit()


def old_main():
    import os
    import DataStructures

    telfile = "%s/School/Research/Useful_Datafiles/Telluric_NearIR.dat" % (os.environ["HOME"])
    x, y, c, e = np.loadtxt(telfile, unpack=True)
    telluric = DataStructures.xypoint(x=x * u.nm.to(u.micron), y=y, cont=c, err=e)
    telluric = FittingUtilities.ReduceResolution2(telluric, 40000.0)

    for fname in sys.argv[1:]:
        orders = HelperFunctions.ReadExtensionFits(fname)
        columns = []
        for i, order in enumerate(orders):
            if order.x[0] < 1.9:
                band = "H"
            else:
                band = "K"
            wavefile = "%s/School/Research/Useful_Datafiles/IGRINS_Datafiles/SDC%s_order%.3i.wavedata" % (
            os.environ["HOME"], band, i + 2)
            left = np.searchsorted(telluric.x, order.x[0] - 2e-3)
            right = np.searchsorted(telluric.x, order.x[-1] + 2e-3)
            model = telluric[left:right]

            order.cont = FittingUtilities.Continuum(order.x, order.y, lowreject=2, highreject=3, fitorder=2)
            shift = CCImprove_Wrapper(order, model, tol=2e-3)
            print shift
            # plt.plot(ccf.x, ccf.y)
            #plt.show()
            order.x += shift

            done = False
            pixels = []
            waves = []
            while not done:
                fig = plt.figure(1)
                ax = fig.add_subplot(111)
                done = True
                ax.cla()
                ax.plot(order.x, order.y / order.cont, 'k-')
                ax.plot(model.x, model.y, 'r-')
                plt.show(block=False)

                inp = raw_input("Wavelength fit okay? [Y]/n ")
                if "n" in inp.lower():
                    done = False
                    order.x, newpixels, newwaves = FitWavelength(order, model)

                    if len(waves) == 0:
                        pixels = newpixels
                        waves = newwaves
                    else:
                        for p, w in zip(newpixels, newwaves):
                            if min(abs(np.array(pixels) - p)) > 0.1:
                                pixels.append(p)
                                waves.append(w)

            column = {"wavelength": order.x,
                      "flux": order.y,
                      "continuum": order.cont,
                      "error": order.err}
            columns.append(column)

            if len(pixels) > 3:
                np.savetxt(wavefile, np.transpose((pixels, waves)))

        HelperFunctions.OutputFitsFileExtensions(column_list, fname, fname, mode="new")


if __name__ == "__main__":
    main()


