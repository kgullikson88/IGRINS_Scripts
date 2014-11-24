import sys

from astropy.io import fits
from astropy.time import Time


if __name__ == "__main__":
    file_list = []
    for arg in sys.argv[1:]:
        if 1:
            file_list.append(arg)

    for fname in file_list:
        hdulist = fits.open(fname, mode='update')
        t = Time(hdulist[0].header['date-obs'], fmt='isot', scale='utc')
        print fname, hdulist[0].header['date-obs'], t.jd
        hdulist[0].header['JD'] = t.jd
        hdulist.flush()
        hdulist.close()
