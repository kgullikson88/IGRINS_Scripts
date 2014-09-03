import sys

import GenericSearch


# Define regions contaminated by telluric residuals or other defects. We will not use those regions in the cross-correlation
badregions = [[0, 1520],  # Blue end of H band (lots of water absorption)
              [1720, 2080],  #In between H and K bands (lots of water absorption)
              [2350, 2500]]  #Red end of K band (lots of water absorption)


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

    GenericSearch.CompanionSearch(fileList, extensions=extensions, resolution=40000.0, trimsize=trimsize)





