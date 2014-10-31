import sys

import GenericSearch


# Define regions contaminated by telluric residuals or other defects. We will not use those regions in the cross-correlation
badregions = [[0, 1510],  # Blue end of H band (lots of water absorption)
              #[1561, 1615],  # CO2 band that is often poorly corrected (for now at least...)
              [1740, 2090],  #In between H and K bands (lots of water absorption)
              [2348, 2500],  #Red end of K band (lots of water absorption)
              [1510, 1520],  #Temporary...
              [1688,1740],
              [2313, 2350]]


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

    GenericSearch.CompanionSearch(fileList,
                                  extensions=extensions,
                                  resolution=45000.0,
                                  trimsize=trimsize,
                                  Tvalues=range(2700, 3900, 100),
                                  metal_values=[-0.5, 0.0, 0.5],
                                  vsini_values=[1.0,],
                                  logg_values=[3.5,],
                                  vbary_correct=False,
				  debug=False,
				  badregions=badregions,
                                  modeldir='/Volumes/DATADRIVE/Stellar_Models/Sorted/Stellar/NearIR/')





