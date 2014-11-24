import sys
import os

import GenericSearch
import pandas

# Define regions contaminated by telluric residuals or other defects. We will not use those regions in the cross-correlation
badregions = [[0, 1510],  # Blue end of H band (lots of water absorption)
              #[1561, 1615],  # CO2 band that is often poorly corrected (for now at least...)
              [1740, 2090],  #In between H and K bands (lots of water absorption)
              [2348, 2500],  #Red end of K band (lots of water absorption)
              [1510, 1520],  #Temporary...
              [1688,1740],
              [2313, 2350]]

if "darwin" in sys.platform:
    modeldir = "/Volumes/DATADRIVE/Stellar_Models/Sorted/Stellar/NearIR/"
elif "linux" in sys.platform:
    modeldir = "/media/FreeAgent_Drive/SyntheticSpectra/Sorted/Stellar/NearIR/"
else:
    modeldir = raw_input("sys.platform not recognized. Please enter model directory below: ")
    if not modeldir.endswith("/"):
        modeldir = modeldir + "/"

def add_oh_lines(oh_file, badregions=[], minstrength=1.0, tol=0.05):
    oh_data = pandas.read_csv(oh_file, header=False, sep=" ", skipinitialspace=True, names=['wave', 'strength'])
    oh = oh_data[oh_data['strength'] > minstrength]
    n = 1.0 + 2.735182e-4 + 131.4182 / oh['wave'] ** 2 + 2.76249e8 / oh['wave'] ** 4
    oh['wave'] = oh['wave'] / (n * 10.0)
    for wave in oh['wave'].values:
        badregions.append([wave - tol, wave + tol])
    return badregions


if __name__ == "__main__":
    #Parse command line arguments:
    fileList = []
    interp_regions = []
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

    # Add strong oh lines to interp_regions
    oh_file = "{}/School/Research/IGRINS_data/plp/master_calib/ohlines.dat".format(os.environ['HOME'])
    interp_regions = add_oh_lines(oh_file, badregions=interp_regions)

    GenericSearch.CompanionSearch(fileList,
                                  extensions=extensions,
                                  resolution=45000.0,
                                  trimsize=trimsize,
                                  vsini_values=[1.0, 10.0, 20.0, 30.0, 40.0],
                                  observatory="McDonald",
                                  vbary_correct=True,
                                  debug=False,
                                  badregions=badregions,
                                  interp_regions=interp_regions,
                                  modeldir='/Volumes/DATADRIVE/Stellar_Models/Sorted/Stellar/NearIR/')





