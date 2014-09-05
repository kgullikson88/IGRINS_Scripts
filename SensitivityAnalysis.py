import os
import sys

from astropy import units

import Search_Fast
import Sensitivity


if __name__ == "__main__":
    # Parse command-line arguments
    vsini_secondary = 20 * units.km.to(units.cm)
    resolution = 40000
    smooth_factor = 0.8
    vel_list = range(-400, 400, 50)
    companion_file = "%s/Dropbox/School/Research/AstarStuff/TargetLists/Multiplicity.csv" % (os.environ["HOME"])
    vsini_file = "%s/School/Research/Useful_Datafiles/Vsini.csv" % (os.environ["HOME"])
    fileList = []
    tolerance = 5.0
    debug = False
    vsini_skip = 10
    vsini_idx = 1
    trimsize =1
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

    Sensitivity.Analyze(fileList,
                        resolution=resolution,
                        debug=True,
                        badregions=Search_Fast.badregions,
                        trimsize=trimsize,
                        modeldir='/media/FreeAgent_Drive/SyntheticSpectra/Sorted/Stellar/NearIR/')
