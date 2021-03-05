#!/usr/bin/env python3

# Compute and plot XES spectra
import sys


def main():
    """
    To execute, run
    ./plot_spectrum_new.py trans.spectrum
    It will produce a convoluted spectrum (in trans.spectrum.dat file)
    and a pdf with the plot in trans.spectrum.pdf
    """

    import matplotlib.pyplot as plt  # library to make plots
    # General utility functions
    import spectral_tools_new as spt

    fname1 = sys.argv[1]
    labels1 = []
    data1 = spt.read_data_with_labels(fname1, labels1)

    # generate convoluted spectrum
    x_0 = round((data1[0, 0]-0.05)*10.0)/10.0
    npts = len(data1[:, 0])  # is it any differnet from len(data1) ?
    x_max = round((data1[npts - 1, 0]+0.05)*10.0)/10.0
    nx = 1000
    step = abs(x_max - x_0) / nx

    # gaussian width, in eV
    # FWHM
    width = 0.025
    print ("A full width at half maximum of 0.025 is used. You may modify this value in the plot_spectrum_new.py script")
    spectrum1 = spt.compute_spectrum(abs(x_0), abs(x_max), step, abs(data1[:, 0]), data1[:, 1], width)

    width = 8
    rect = (0.2, 0.2, 0.75, 0.75)

    fig = plt.figure(1, (width, width))
    ax = fig.add_axes(rect)
    if max(spectrum1[:, 1]) < 10:
      ylim = round((max(spectrum1[:, 1])+0.5))
    elif max(spectrum1[:, 1]) < 100:
      ylim = round((max(spectrum1[:, 1])+5)/10.0)*10
    else:
      ylim = round((max(spectrum1[:, 1])+50)/100.0)*100
    ax.set_xlim(abs(x_0), abs(x_max))
    ax.set_ylim(0.0, ylim)
    ax.plot(spectrum1[:, 0], spectrum1[:, 1], '-', color='blue', label=fname1)
    plt.ylabel('Intensity', fontsize=12)
    plt.xlabel('Energy, eV', fontsize=12)
    leg = ax.legend(prop={"size":12})

    fname2 = fname1 + '.pdf'
    plt.savefig(fname2)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
