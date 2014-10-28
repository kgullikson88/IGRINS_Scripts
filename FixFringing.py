__author__ = 'Kevin Gullikson'


import numpy as np
import sys
import HelperFunctions
import FittingUtilities
import matplotlib.pyplot as plt


def get_fft(data, rebin=True, preprocess=True):
    if rebin:
        xgrid = np.linspace(data.x[0], data.x[-1], data.size())
        data = FittingUtilities.RebinData(data, xgrid)

    if preprocess:
        data.y = data.y/data.cont
        data.y -= data.y.mean()

    fft = np.fft.fft(data.y)
    freq = np.fft.fftfreq(data.size(), data.x[1] - data.x[0])

    return freq, fft



if __name__ == "__main__":
    fileList = []
    for arg in sys.argv[1:]:
        if 1:
            fileList.append(arg)


    for fname in fileList:
        orders = HelperFunctions.ReadExtensionFits(fname)

        for order in orders[28:40]:
            cont = order.cont
            mean = np.mean(order.y/order.cont)
            order.y = order.y/cont - mean

            # Plot initial order
            fig1 = plt.figure(1)
            fig2 = plt.figure(2)
            ax1 = fig1.add_subplot(111)
            ax2 = fig2.add_subplot(111)

            ax1.plot(order.x, order.y, 'k-', alpha=0.4)

            # Get and plot the fft
            freq, fft = get_fft(order.copy(), preprocess=False)
            searchvalues = np.where((freq > 0.2) & (freq < 2.0))[0]
            freq = freq[searchvalues]
            fft = fft[searchvalues]
            maxindex = np.argmax(fft.real**2 + fft.imag**2)
            print "Central Wavelength: ", np.median(order.x)
            print "Maximum value at ", freq[maxindex]
            print "Maximum value is ", (fft.real**2 + fft.imag**2)[maxindex]
            print "\n"
            ax2.plot(freq, fft.real**2 + fft.imag**2, 'k-', alpha=0.4)

            # Plot markers at the expected fringe peak locations
            for n in [18.0, 19.0, 20.0]:
                n /= 0.01093
                dlam = np.median(order.x) / n
                i = np.argmin(abs(freq - 1.0/dlam))
                x = freq[i]
                y = fft.real[i]**2 + fft.imag[i]**2
                ax2.plot((x, x), (y+0.1, y+1), 'r-')
                print n, dlam


            plt.show()