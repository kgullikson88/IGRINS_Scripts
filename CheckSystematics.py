from __future__ import print_function
from collections import defaultdict

from astropy.io import fits
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def telluric_systematics(file_list, header_kw='HUMIDITY', plot=True, ref_num=10):
    """
    Find the average change in the requested header keyword over all IGRINS telluric fits
    :param file_list: A list of IGRINS Corrected_* files
    :param header_kw: the keyword to search for in each fits extension header
    :param plot: If true, plot all of the values
    :param ref_num: Which extension to use as the reference
    :return: The median value of the keyword for each fits extension
    """
    data = defaultdict(list)
    for fname in file_list:
        print(fname)
        hdulist = fits.open(fname)
        for i, hdu in enumerate(hdulist[1:]):
            data['fname'].append(fname)
            data['Ext_num'].append(i + 1)
            data['Middle_Wavelength'].append(np.median(hdu.data.field('wavelength')))
            data['value'].append(hdu.header[header_kw])

    # Convert the dict to a dataframe for easier access
    df = pd.DataFrame(data=data)


    # Scale by the median value of the keyword within the given filename
    # (to account for weather variations rather than data systematics)
    median_values = defaultdict(float)
    for fname in file_list:
        median_values[fname] = float(df.loc[(df.fname == fname) & (df.Ext_num == ref_num)]['value'])
    print(median_values)
    make_scaled = lambda row: row['value'] / median_values[row['fname']]
    df['scaled_value'] = df.apply(make_scaled, axis=1)

    # Determine the median value for each fits extension
    median = df.groupby('Ext_num').median()[['Middle_Wavelength', 'scaled_value']]

    # Plot, if desired
    if plot:
        plt.scatter(df.Middle_Wavelength, df.scaled_value, color='red', alpha=0.1)
        plt.plot(median.Middle_Wavelength, median.scaled_value, color='green', lw=2)
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('{} scale factor'.format(header_kw))
        plt.show()

    return median

