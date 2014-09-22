import sys

import numpy as np

import HelperFunctions


if __name__ == "__main__":
    fileList = []
    for arg in sys.argv[1:]:
        fileList.append(arg)

    trimsize = 100
    for fname in fileList:
        orders = HelperFunctions.ReadExtensionFits(fname)
        snrList = [np.mean(o.y[trimsize:-trimsize]) / np.std(o.y[trimsize:-trimsize]) for o in orders]
        print "Max S/N for {:s} = {:.1f}".format(fname, max(snrList))