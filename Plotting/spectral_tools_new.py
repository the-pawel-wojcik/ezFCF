#!/usr/bin/env python3

# Tools for computing spectra
import math
import numpy as np  # linear algebra library


def gaussian(x0, width, x):
    """
    # FWHM = 2 sqrt(2*ln2)x width
    """
    sigma = width / 2.35482
    norm = sigma * (2. * 3.14159)**0.5
    return 1. / norm * math.exp(-0.5 * ((x0 - x) / sigma)**2)


def compute_spectrum(emin, emax, step, stick_pos, fl, width):
    """
    Compute gaussian on the specified grid from stick spectrum
    """
    nepts = int((emax - emin) / step)
    spectrum = np.zeros((nepts, 2))

    for i in range(nepts):
        spectrum[i][0] = emin + step * i
        for j in range(len(stick_pos)):
            spectrum[i][1] = spectrum[i][1] + gaussian(stick_pos[j], width, spectrum[i][0]) * fl[j]

    return spectrum


def read_data_with_labels(fname, labels):
    """
    Reads X Y1, Y2, Y3, "label" ....
    """
    sfile = open(fname, 'r')

    tmp_data = []
    for line in sfile:
        if '#' not in line:
            tmp_line = line.split()
            ncol = len(tmp_line)
            tmp_data_line = []
            for i in range(ncol - 1):
                tmp_data_line.append(float(tmp_line[i]))
            tmp_data_line.append(str(tmp_line[ncol - 1]))
            tmp_data.append(tmp_data_line)

    ncol = len(tmp_data[0])
    npoints = len(tmp_data)

    data = np.zeros((npoints, ncol - 1))
    for i in range(npoints):
        for j in range(ncol - 1):
            data[i][j] = tmp_data[i][j]
        labels.append(tmp_data[i][ncol - 1])

    return data
