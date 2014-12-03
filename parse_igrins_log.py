import os

import pandas
import numpy as np

from HelperFunctions import IsListlike


def hex_to_dex(str):
    s = -1.0 if "-" in str else 1.0
    segments = str.split(":")
    degree = abs(float(segments[0]))
    minute = float(segments[1])
    second = float(segments[2])
    return s * (degree + minute / 60.0 + second / 3600.0)


def dex_to_hex(value):
    s = "-" if value < 0 else "+"
    value = abs(value)
    degree = int(value)
    minute = int((value - degree) * 60)
    second = (value - (degree + minute / 60.0)) * 3600
    return "{:s}{:02d}:{:02d}:{:02f}".format(s, degree, minute, second)


def read_logfile(filename):
    if IsListlike(filename):
        # Assume it is a list
        all_data = []
        for fname in filename:
            all_data.append(read_logfile(fname))
        return pandas.concat(all_data, ignore_index=True)

    data = pandas.read_csv(filename, skiprows=1)

    # Convert RA and DEC to decimals
    data['RA'] = data['RA'].map(hex_to_dex)
    data['DEC'] = data['DEC'].map(hex_to_dex)

    #Convert the filename column to file number
    data['File_Number'] = data['FILENAME'].map(lambda s: int(s[-9:-5]))

    return data


def get_average(data, num_list, key):
    values = []
    for num in num_list:
        matches = data[data.File_Number == num][key].values
        if len(matches) != 1:
            continue
        values.append(matches[0])
    return np.mean(values)


def get_logfilenames(dir, band="H"):
    if not dir.endswith("/"):
        dir += "/"
    file_list = ["{:s}{:s}".format(dir, f) for f in os.listdir(dir) if
                 f.startswith("IGRINS") and f.endswith("{:s}.txt".format(band))]

    if len(file_list) < 2:
        return file_list[0]
    return file_list
