import os

import matplotlib.pyplot as plt
import numpy as np

import HelperFunctions


if __name__ == "__main__":
    filenames = [f for f in os.listdir("./") if f.endswith("smoothed.fits") and f.startswith("H")]
    corrdir = "Cross_correlations/"
    vsini_values = [10, 20, 30, 40]
    Temperatures = [3300, 3500, 3700, 3900, 4200, 4500, 5000, 5500]
    Temperatures = range(3000, 6800, 100)
    metals = [-0.5, 0.0, 0.5]
    logg = 4.5
    HelperFunctions.ensure_dir("Figures/")

    for rootfile in sorted(filenames):
        Tvals = []
        Zvals = []
        rotvals = []
        significance = []
        for T in Temperatures:
            for metal in metals:
                for vsini in vsini_values:
                    corrfile = "{0:s}{1:s}.{2:d}kps_{3:d}K{4:+.1f}{5:+.1f}".format(corrdir,
                                                                                   rootfile.split(".fits")[0],
                                                                                   vsini,
                                                                                   T,
                                                                                   logg,
                                                                                   metal)
                    print corrfile
                    try:
                        vel, corr = np.loadtxt(corrfile, unpack=True)
                    except IOError:
                        continue

                    # Check the significance of the highest peak within +/- 500 km/s
                    left = np.searchsorted(vel, -500)
                    right = np.searchsorted(vel, 500)
                    idx = np.argmax(corr[left:right]) + left
                    v = vel[idx]
                    goodindices = np.where(np.abs(vel - v) > vsini)[0]
                    std = np.std(corr[goodindices])
                    mean = np.mean(corr[goodindices])
                    sigma = (corr[idx] - mean) / std

                    # Plot if > 3 sigma peak
                    if sigma > 4:
                        fig = plt.figure(10)
                        ax = fig.add_subplot(111)
                        ax.plot(vel, corr, 'k-', lw=2)
                        ax.set_xlabel("Velocity (km/s)")
                        ax.set_ylabel("CCF")
                        ax.set_title(r'{0:s}:  $T_s$={1:d}K & [Fe/H]={2:.1f}'.format(rootfile, T, metal))
                        ax.grid(True)
                        fig.savefig(u"Figures/{0:s}.pdf".format(corrfile.split("/")[-1]))
                        plt.close(fig)

                    Tvals.append(T)
                    Zvals.append(metal)
                    rotvals.append(vsini)
                    significance.append(sigma)

        # Now, make a plot of the significance as a function of Temperature and metallicity for each vsini
        Tvals = np.array(Tvals)
        Zvals = np.array(Zvals)
        rotvals = np.array(rotvals)
        significance = np.array(significance)
        fig = plt.figure(1)
        ax = fig.add_subplot(111, projection='3d')
        ax.set_title("Significance Summary for %s" % (rootfile.split(".fits")[0].replace("_", " ")))
        for i, rot in enumerate(vsini_values):
            goodindices = np.where(abs(rotvals - rot) < 1e-5)[0]
            ax.set_xlabel("Temperature (K)")
            ax.set_ylabel("[Fe/H]")
            ax.set_zlabel("Significance")
            ax.plot(Tvals[goodindices], Zvals[goodindices], significance[goodindices], 'o', label="%i km/s" % rot)
        leg = ax.legend(loc='best', fancybox=True)
        leg.get_frame().set_alpha(0.5)
        fig.savefig("Figures/Summary_%s.pdf" % rootfile.split(".fits")[0])
        idx = np.argmax(significance)
        #ax.plot(Tvals[idx], Zvals[idx], significance[idx], 'x', markersize=25, label="Most Significant")
        plt.show()


