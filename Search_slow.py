import sys
import os

import pandas

import GenericSearch
import StarData








# Define regions contaminated by telluric residuals or other defects. We will not use those regions in the cross-correlation
badregions = [[0, 1510],  # Blue end of H band (lots of water absorption)
              # [1561, 1615],  # CO2 band that is often poorly corrected (for now at least...)
              [1740, 2090],  # In between H and K bands (lots of water absorption)
              [2380, 2500],  # Red end of K band (lots of water absorption)
              # [1688, 1740],
              # [2313, 2350],
]

if "darwin" in sys.platform:
    modeldir = "/Volumes/DATADRIVE/Stellar_Models/PHOENIX/Stellar/Vband/"
    hdf5_filename = '/Volumes/DATADRIVE/PhoenixGrid/IGRINS_Grid.hdf5'
elif "linux" in sys.platform:
    modeldir = "/media/FreeAgent_Drive/SyntheticSpectra/Sorted/Stellar/Vband/"
    hdf5_filename = '/media/ExtraSpace/PhoenixGrid/IGRINS_Grid.hdf5'
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


trimsize = 10
homedir = os.environ['HOME']
oh_file = "{}/School/Research/IGRINS_data/plp/master_calib/ohlines.dat".format(homedir)
interp_regions = []
interp_regions = add_oh_lines(oh_file, badregions=interp_regions)

if __name__ == '__main__':
    # Parse command line arguments:
    fileList = []
    for arg in sys.argv[1:]:
        if 1:
            fileList.append(arg)

    # Get the primary star vsini values
    prim_vsini = StarData.get_vsini(fileList)

    # Remove anything without a vsini
    new_file_list = []
    new_prim_vsini = []
    for vsini, fname in zip(prim_vsini, fileList):
        if vsini is not None:
            new_file_list.append(fname)
            new_prim_vsini.append(vsini)

    GenericSearch.slow_companion_search(new_file_list, new_prim_vsini,
                                        hdf5_file=hdf5_filename,
                                        extensions=True,
                                        resolution=None,
                                        trimsize=trimsize,
                                        modeldir=modeldir,
                                        badregions=badregions,
                                        interp_regions=interp_regions,
                                        metal_values=(0, -0.5, 0.5),
                                        # vsini_values=(1, 5.0, 10.0, 20.0, 30.0),
                                        logg_values=(4.5),
                                        Tvalues=range(3000, 9000, 100),
                                        vsini_values=(1, 5, 10, 20, 30),
                                        observatory='McDonald',
                                        debug=False,
                                        vbary_correct=True,
                                        addmode='simple',
                                        output_mode='hdf5')

