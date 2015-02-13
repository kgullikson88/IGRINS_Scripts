import sys
import os

import pandas
from astropy.io import fits

import GenericSearch



# Define regions contaminated by telluric residuals or other defects. We will not use those regions in the cross-correlation
badregions = [[0, 1510],  # Blue end of H band (lots of water absorption)
              # [1561, 1615],  # CO2 band that is often poorly corrected (for now at least...)
              [1740, 2090],  #In between H and K bands (lots of water absorption)
              [2348, 2500],  #Red end of K band (lots of water absorption)
              [1510, 1520],  #Temporary...
              [1688, 1740],
              [2313, 2350]]

if "darwin" in sys.platform:
    modeldir = "/Volumes/DATADRIVE/Stellar_Models/PHOENIX/Stellar/Vband/"
elif "linux" in sys.platform:
    modeldir = "/media/FreeAgent_Drive/SyntheticSpectra/Sorted/Stellar/Vband/"
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


if __name__ == '__main__':
    # Parse command line arguments:
    fileList = []
    interp_regions = []
    trimsize = 10
    for arg in sys.argv[1:]:
        if 1:
            fileList.append(arg)
    prim_vsini = [100.0] * len(fileList)

    # Add strong oh lines to interp_regions
    oh_file = "{}/School/Research/IGRINS_data/plp/master_calib/ohlines.dat".format(os.environ['HOME'])
    interp_regions = add_oh_lines(oh_file, badregions=interp_regions)

    # Get the primary star vsini values
    prim_vsini = []
    vsini = pandas.read_csv("../../Useful_Datafiles/Vsini.csv", sep='|', skiprows=8)[1:]
    vsini_dict = {}
    for fname in fileList:
        root = fname.split('/')[-1][:9]
        if root in vsini_dict:
            prim_vsini.append(vsini_dict[root])
        else:
            header = fits.getheader(fname)
            star = header['OBJECT1']
            print fname, star
            v = vsini.loc[vsini.Identifier.str.strip() == star]['vsini(km/s)'].values[0]
            prim_vsini.append(float(v) * 0.8)
            vsini_dict[root] = float(v) * 0.8
    for fname, vsini in zip(fileList, prim_vsini):
        print fname, vsini

    GenericSearch.slow_companion_search(fileList, prim_vsini,
                                        hdf5_file='/media/ExtraSpace/PhoenixGrid/IGRINS_Grid.hdf5',
                                        extensions=True,
                                        resolution=None,
                                        trimsize=trimsize,
                                        modeldir=modeldir,
                                        badregions=badregions,
                                        interp_regions=interp_regions,
                                        metal_values=(0.0,),
                                        vsini_values=(1, 5.0, 10.0, 15.0),
                                        Tvalues=range(3000, 6900, 100),
                                        observatory='McDonald',
                                        debug=False,
                                        vbary_correct=False,
                                        addmode='ml',
                                        output_mode='hdf5')

