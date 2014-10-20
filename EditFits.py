__author__ = 'Kevin Gullikson'
from astropy.io import fits


def edit_keyword(fname, keyword, keyval, comment=None):
    hdulist = fits.open(fname, mode='update')
    if comment is None:
        hdulist[0].header[keyword] = keyval
    elif isinstance(comment, str):
        hdulist[0].header[keyword] = (keyval, comment)

    hdulist.flush()
    hdulist.close()


if __name__ == "__main__":
    import sys

    fileList = []
    keyword = None
    keyval = None
    comment = None
    for arg in sys.argv[1:]:
        if "--keyword" in arg:
            keyword = arg.split("=")[-1]
        elif "--keyval" in arg:
            keyval = arg.split("=")[-1]
        elif "--comment" in arg:
            comment = arg.split("=")[-1]
        else:
            fileList.append(arg)

    # Apply the updated header to all of the requested files
    for fname in fileList:
        edit_keyword(fname, keyword, keyval, comment=comment)
