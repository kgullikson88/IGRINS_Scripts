""" Measure the primary star RV from CCFS.
"""
import sys
from collections import defaultdict

import h5py
import pandas as pd
import numpy as np
from astropy.io import fits
import os


ADDMODE = 'ml'


def collect_rv(hdf5_file, output_log=None):
    """ Get the RV of each star/date combo in the HDF5 file.
        This function just uses the value in the attributes.
    """
    with h5py.File(hdf5_file, 'r') as f:
        df_list = []
        for starname, star_grp in f.iteritems():
            for date, date_grp in star_grp.iteritems():
                print(date_grp.name)
                summary = defaultdict(list)
                for ds_name, dataset in date_grp.iteritems():
                    summary['teff'].append(dataset.attrs['T'])
                    summary['vsini'].append(dataset.attrs['vsini'])
                    summary['logg'].append(dataset.attrs['logg'])
                    summary['feh'].append(dataset.attrs['[Fe/H]'])
                    summary['addmode'].append(dataset.attrs['addmode'])
                    summary['RV'].append(dataset.attrs['vel_max'])
                    rv_err = dataset.attrs['vel_err'] if 'vel_err' in dataset.attrs else np.nan
                    summary['RV_err'].append(rv_err)
                    summary['CCF'].append(dataset.attrs['ccf_max'])
                df = pd.DataFrame(data=summary)
                df['star'] = starname
                df['date'] = date
                df_list.append(df)

                # Save the maximum value, if requested
                if output_log is not None:
                    best = df.loc[df.addmode == ADDMODE].sort_values(by='CCF').tail(1)
                    best = {k: v for k, v in zip(best.columns, best.values[0])}
                    with open(output_log, 'a') as log:
                        log.write('{star},{date},{teff},{logg},{feh},{vsini},{addmode},{RV},{RV_err},{CCF}\n'.format(**best))

    return pd.concat(df_list, ignore_index=True)



jd_cache = {}
def get_jd(fname, dirname='{}/School/Research/IGRINS_data'.format(os.environ['HOME'])):
    fname = os.path.join(dirname, fname)
    if fname in jd_cache:
        return jd_cache[fname]
    if not os.path.exists(fname):
        raise OSError('Filename does not exist: {}'.format(fname))
        #return np.nan
    header = fits.getheader(fname)
    if 'HJD' in header:
        jd = header['HJD']
    else:
        jd = header['JD']

    jd_cache[fname] = jd
    return jd


def measure_rv(hdf5_file, output_log=None, update_attrs=True):
    """
    Measure the RV of each star/date combo in the HDF5 file.
    This function actually reads in the ccfs, so it will be slower than collect_rv!

    :param hdf5_file: The name of the HDF5 file.
    :param output_log: A logfile to write the best values to
    :param update_attrs: Should we update the vel_max, ccf_max, and vel_err attributes in each dataset?
    :return: a pandas.DataFrame containing all of the measurements.
    """
    from Correlate import get_rv
    from GenericSearch import Process_Data
    from Search_primary import badregions, interp_regions, trimsize

    readmode = 'a' if update_attrs else 'r'
    with h5py.File(hdf5_file, readmode) as f:
        df_list = []
        for starname, star_grp in f.iteritems():
            for date, date_grp in star_grp.iteritems():
                print(date_grp.name)
                # Read in the filename by looking at the first dataset (we really just need the number of pixels)
                ds_names = date_grp.keys()
                if len(ds_names) == 0:
                    continue
                fname = date_grp[ds_names[0]].attrs['fname']
                orders = Process_Data(fname, badregions, interp_regions=interp_regions, logspacing=True,
                                      extensions=True, trimsize=trimsize, vsini=None,
                                      reject_outliers=False)
                Npix = sum([o.size() for o in orders])
                jd = get_jd(fname)

                # Loop over the datasets
                summary = defaultdict(list)
                for ds_name, dataset in date_grp.iteritems():
                    # Read the basic parameters
                    summary['teff'].append(dataset.attrs['T'])
                    summary['vsini'].append(dataset.attrs['vsini'])
                    summary['logg'].append(dataset.attrs['logg'])
                    summary['feh'].append(dataset.attrs['[Fe/H]'])
                    summary['addmode'].append(dataset.attrs['addmode'])
                    summary['HJD'].append(jd)

                    # Estimate the rv and rv_error
                    vel, corr = dataset.value
                    rv, rv_err, ccf = get_rv(vel, corr, Npix)
                    if dataset.attrs['addmode'] != 'ml':
                        rv_err = np.nan
                    summary['RV'].append(rv)
                    summary['RV_err'].append(rv_err)
                    summary['CCF'].append(ccf)

                    # Update dataset attributes
                    if update_attrs:
                        dataset.attrs['vel_max'] = rv
                        dataset.attrs['vel_err'] = rv_err
                        dataset.attrs['ccf_max'] = ccf
                        f.flush()

                df = pd.DataFrame(data=summary)
                df['star'] = starname
                df['date'] = date
                df_list.append(df)

                # Save the maximum value, if requested
                if output_log is not None:
                    best = df.loc[df.addmode == ADDMODE].sort_values(by='CCF').tail(1)
                    best = {k: v for k, v in zip(best.columns, best.values[0])}
                    with open(output_log, 'a') as log:
                        log.write(
                            '{star},{date},{HJD},{teff},{logg},{feh},{vsini},{addmode},{RV},{RV_err},{CCF}\n'.format(
                                **best))

    return pd.concat(df_list, ignore_index=True)


if __name__ == '__main__':
    for fname in sys.argv[1:]:
        output = fname.replace('hdf5', 'rv.txt')
        with open(output, 'w') as log:
            log.write('star,date,HJD,teff,logg,feh,vsini,addmode,rv,rv_err,ccf\n')
        summary = measure_rv(fname, output_log=output)
        summary.to_csv(output.replace('.txt', '_summary.txt'), index=False)
