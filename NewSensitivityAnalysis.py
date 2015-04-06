"""
Sensitivity analysis, using the new search method.
"""
import sys
import logging

import matplotlib.pyplot as plt

import Sensitivity
import StarData
import SpectralTypeRelations
import Search_slow
from HelperFunctions import ensure_dir

logging.basicConfig(level='INFO')

MS = SpectralTypeRelations.MainSequence()


def check_sensitivity():
    fileList = []
    for arg in sys.argv[1:]:
        if 1:
            fileList.append(arg)

    badregions = Search_slow.badregions
    interp_regions = Search_slow.interp_regions
    trimsize = Search_slow.trimsize
    prim_vsini = StarData.get_vsini(fileList)
    
    # Remove anything without a vsini
    new_file_list = []
    new_prim_vsini = []
    for vsini, fname in zip(prim_vsini, fileList):
        if vsini is not None:
            new_file_list.append(fname)
            new_prim_vsini.append(vsini)

    Sensitivity.Analyze(new_file_list, new_prim_vsini,
                        hdf5_file='/media/ExtraSpace/PhoenixGrid/IGRINS_Grid.hdf5',
                        extensions=True,
                        resolution=None,
                        trimsize=trimsize,
                        badregions=badregions, interp_regions=interp_regions,
                        metal_values=(0.0,),
                        vsini_values=(0, 10, 20, 30, 40, 50),
                        Tvalues=range(3000, 6000, 100),
                        debug=False,
                        addmode='simple',
                        output_mode='hdf5')


if __name__ == '__main__':
    if '--analyze' in sys.argv[1]:
        # Make the 2d plots
        df = Sensitivity.analyze_sensitivity(hdf5_file='Sensitivity_known.hdf5')

    elif '--marginalize' in sys.argv[1]:
        fig, ax = Sensitivity.marginalize_sensitivity(infilename='Sensitivity_Dataframe.csv')
        # plt.show()
        ensure_dir('Figures/')
        plt.savefig('Figures/Sensitivity_Marginalized.pdf')


    else:
        check_sensitivity()
