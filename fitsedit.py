from astropy.io import fits
import sys

if __name__ == "__main__":
    header_kws = {}
    for arg in sys.argv[1:]:
        if "=" in arg:
            s = arg.partition("=")
            key = s[0]
            value = s[-1]
            try:
                value = float(value)
            except ValueError:
                pass
            header_kws[key] = value
        else:
            filename = arg

    hdulist = fits.open(filename, mode='update')
    for key in header_kws.keys():
        hdulist[0].header[key] = value

    hdulist.flush()
    hdulist.close()