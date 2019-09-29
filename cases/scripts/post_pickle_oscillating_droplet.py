#!/usr/bin/env python3

from argparse import ArgumentParser
import pickle as pl
import numpy as np
import math
from ast import literal_eval

def main():

    #---- Command line arguments ----------------------------------------------
    parser = ArgumentParser(description="Read, modify and save again a pickled plot of an oscillating droplet (semi-axis over time).")
    parser.add_argument("nu_droplet",
                        help="Kinematic viscosity of the droplet")
    parser.add_argument("file_name",
                        help="File name of the pickled plot.")

    args = parser.parse_args()

    figure = pl.load(open(args.file_name, 'rb'))
    xlimits = figure.axes[0].get_xlim()
    ylimits = figure.axes[0].get_ylim()

    # Material parameters
    rho_d = 10.0
    rho_a = 0.1
    nu_d = float(args.nu_droplet)
    nu_a = 5.0e-4
    sigma = 10.0
    n_mode = 2.0
    R0 = 1.0
    a0 = 0.025

    tau = R0*R0 / ((n_mode - 1)*(2*n_mode + 1)*nu_d)

    x = figure.axes[0].lines[0].get_data()[0]

    y_high = []
    y_low = []

    for x_val in x:
        y_high.append(1 + a0*math.exp(-x_val/tau))
        y_low.append(1 - a0*math.exp(-x_val/tau))

    new_axes = figure.axes[0]
    new_axes.plot(x, y_high, x, y_low, color='black', linewidth=0.75, linestyle='dashed')

    # Adapt x-axis labelling
    T = 2*math.pi / math.sqrt( n_mode*(n_mode + 1)*(n_mode - 1)*(n_mode + 2)*sigma/(
                            R0*R0*R0*((n_mode + 1)*rho_d + n_mode*rho_a)))

    new_axes.set_xticks([0.0, 0.5*T, 1.0*T, 1.5*T, 2.0*T, 2.5*T, 3.0*T])
    new_axes.set_xticklabels(['0', '1/2', '1', '3/2', '2', '5/2', '3'])
    new_axes.set_xlabel('$t^*$')

    # The x- and y-limits of the plot may be reseted. Thus, enfore them
    # again
    figure.axes[0].set_xlim(left=xlimits[0], right=xlimits[1])

    modified_file_name = args.file_name.split('.')[0]
    modified_file_name += '_modified.pdf'
    figure.savefig(modified_file_name, bbox_inches='tight')


if __name__ == "__main__":
    main()
