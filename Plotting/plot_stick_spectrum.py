#!/usr/bin/env python3

import sys
import re 
import math as m
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

COLORS = [color for color in mcolors.TABLEAU_COLORS.keys()]
FONTSIZE = 12
CM2INCH = 1/2.54


def height_one_lorenzian(x, x0, gamma):
    return lorenzian(x, x0, gamma) * m.pi * 0.5 * gamma


def lorenzian(x, x0, gamma):
    '''
    This definition of a Lorenzian comes from Eq. 2.40

    Köppel, H., W. Domcke, and L. S. Cederbaum. “Multimode Molecular Dynamics 
    Beyond the Born-Oppenheimer Approximation.” In Advances in Chemical Physics,
    edited by I. Prigogine and Stuart A. Rice, LVII:59–246, 1984. 
    https://doi.org/10.1002/9780470142813.ch2.

    gamma corresponds to full width at half max (FWHM)
    '''
    return 1.0 / m.pi * ( 0.5 * gamma / ((x - x0)**2 + (0.5 * gamma) **2))

def parse_raw_ezfcf_spectrum(raw_spectrum, shift_eV):
    spectrum = []
    for raw_line in raw_spectrum:
        split_line = [float(value) for value in raw_line.split()[:-1]]
        line = {'Energy (eV)': split_line[0] + shift_eV,
                'Intensity': split_line[1],
                'FCF': split_line[2],
                'assignment': raw_line.split()[-1]}
        spectrum += [line]

    return spectrum

def turn_spectrum_to_xs_and_ys(spectrum):
    xs = []
    ys = []
    for spectral_point in ezfcf_output:
        xs += [spectral_point['Energy (eV)']]
        ys += [spectral_point['Intensity']]

    return (xs, ys)

def lorenz_intensity(x, gamma, ezfcf_output):
    out = 0.0
    for spectral_point in ezfcf_output:
        energy_eV = spectral_point['Energy (eV)']
        intensity = spectral_point['Intensity']
        out += intensity * height_one_lorenzian(x, energy_eV, gamma)
    return out

def generate_xs_and_ys(left, right):
    npoints = 500
    xmin = left
    xmax = right
    xrange = xmax-xmin
    xstep = xrange / npoints
    xs = [xmin + xstep * i for i in range(npoints)]
    ys = [0.0 for x in xs]

    return (xs, ys)

def lorenz_ezfcf_output(ezfcf_output, ax, xs, ys, gamma, color):
    new_ys = [y + lorenz_intensity(x, gamma, ezfcf_output) for x, y in zip(xs, ys)]

    ax.fill_between(xs, new_ys, ys, color=color, alpha=0.2)

    return new_ys

def stem_ezfcf_output(ezfcf_output, ax, color: str = 'k'):
    for spectral_point in ezfcf_output:
        energy_eV = spectral_point['Energy (eV)']
        intensity = spectral_point['Intensity']
        ax.vlines(x=energy_eV, ymin=0.0, ymax=intensity, colors=color)


def parse_command_line():
    # Parse the command line arguments.
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_files",
                        help="List of xsim output files.",
                        nargs="+",
                        required=True)

    parser.add_argument("-s", "--shift_eV", 
                        help="Shift added to the peak positions.",
                        required=False,
                        type=float,
                        default=0.0)

    parser.add_argument("-e", "--envelope", 
                        help="Plot an envelope.",
                        required=False,
                        type=bool,
                        default=False)

    parser.add_argument("-g", "--gamma", 
                        help="Gamma in the Lorenzian:\n"+\
                        "(0.5 * gamma)**2 / ((x - x0)**2 + (0.5 * gamma) **2)",
                        required=False,
                        type=float,
                        default=0.030)

    args = parser.parse_args()
    return args

def get_ezfcf_outputs(args):
    shift_eV = args.shift_eV
    ezfcf_outputs = []
    for file_idx, out_file in enumerate(args.output_files):
        with open(out_file) as f:
            lines = f.readlines()
        ezfcf_output = parse_raw_ezfcf_spectrum(lines, shift_eV)
        ezfcf_outputs += [ezfcf_output]
    return ezfcf_outputs

def add_envelope(ax, args, ezfcf_outputs, left, right):
    envelope = args.envelope
    if envelope is False:
        return None

    xs, ys =  generate_xs_and_ys(left, right)
    gamma = args.gamma
    for file_idx, ezfcf_output in enumerate(ezfcf_outputs):
        if file_idx == len(COLORS) or file_idx > len(COLORS):
            print("Too many colors already.")
            sys.exit(1)
        ys = lorenz_ezfcf_output(ezfcf_output, ax, xs, ys, gamma, COLORS[file_idx])

    fig_max_y = max(ys)
    return fig_max_y

def add_peaks(ax, ezfcf_outputs):
    peaks_maxima = []
    for file_idx, ezfcf_output in enumerate(ezfcf_outputs):
        if file_idx == len(COLORS) or file_idx > len(COLORS):
            print("Too many colors already.")
            sys.exit(1)
        stem_ezfcf_output(ezfcf_output, ax, COLORS[file_idx])
        # peaks_maxima.append(max(ezfcf_output, key=lambda x: x['Intensity']))
    return 1.0


def add_info_text(ax, args):
    info_kwargs = {'horizontalalignment': 'left',
                   'verticalalignment': 'top',
                   'fontsize': FONTSIZE,
                   'color': 'k',
                   'transform': ax.transAxes,
                   }
    text = f'$\gamma = {args.gamma:.2f}$\n$s = {args.shift_eV:.2f}$'
    ax.text(0.01, 0.99, text, **info_kwargs)


def main():
    args = parse_command_line()
    ezfcf_outputs = get_ezfcf_outputs(args)
    fig, ax = plt.subplots(figsize=(12*CM2INCH, 12*CM2INCH))
    # envelope_max_y = add_envelope(ax, args, ezfcf_outputs, left, right)
    max_peak = add_peaks(ax, ezfcf_outputs)
    # add_info_text(ax, args)

    # match the spectrum range with the experimental range
    # ax.set_xticks([12.50,12.75,13.00,13.25])
    # ax.set_xlim([left,right])

    # remove all marks on the y axis
    # ax.get_yaxis().set_ticks([])

    fig.tight_layout()
    # exp_spectrum_png = plt.imread("./pics/Experiments/Dyke-et-al-1974.png")
    # ax.imshow(exp_spectrum_png, zorder=-1, extent=[left, right, 0.0, 1.1 * fig_max_y])

    # import os 
    # filename = '/home/pawel/chemistry/ozone/plotter/pics/'
    # filename += os.path.basename(args.output_files[0]) 
    # filename += '.pdf'
    # plt.savefig(filename)
    plt.show()

    # plt.show()

    return 0

if __name__ == '__main__':
    main()
