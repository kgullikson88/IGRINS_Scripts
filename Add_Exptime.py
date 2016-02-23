from __future__ import print_function

from astropy.io import fits
import sys

BASE = '/Volumes/DATADRIVE/IGRINS_data/'
BASE = '/media/FreeAgent_Drive_/data/IGRINS_data/'

if __name__ == '__main__':
    for fname in sys.argv[1:]:
        hdulist = fits.open(fname, mode='update')
        header = hdulist[0].header
        nod_keys = [key for key in header if key.startswith('NODFILE')]
        exptime = 0.0
        for key in nod_keys:
            tmp_header = fits.getheader('{}{}'.format(BASE, header[key]))
            exptime += float(tmp_header['EXPTIME'])

        print(fname, exptime)
        hdulist[0].header.set('EXPTIME', exptime, 'Total exposure time (seconds)')
        hdulist.flush()
        hdulist.close()