#!/usr/bin/env python3

from argparse import ArgumentParser
from ast import literal_eval

import dataframeWithMetadata as dwm
import dataPlotting as dp

# TODO for plotting module:
# - Allow filtering of values -> cut off metrics that may explode, e.g. kinetic energy

def read_key_value_pairs(input_string):
    """
    Read key-value pairs from a string and return them in a dictionary.
    """
    key_value_pairs = input_string.split(',')
    parsed_data = {}

    if not input_string:
        return parsed_data

    for pair in key_value_pairs:
        key,value = pair.split(':')
        parsed_data[key.strip()] = value.strip()

    return parsed_data

def main():

    #---- Command line arguments ----------------------------------------------
    parser = ArgumentParser(description="Create plots from a multiindexed dataframe.")
    parser.add_argument("file_name",
                        help="File name of the data frame.")
    parser.add_argument("-xq", "--x-quantity",
                        help="Quantity to use as x-axis (default: time)",
                        default="time",
                        dest="xmetric") 
    parser.add_argument("-p", "--parameters",
                        help="Comma-separated list of parameters to be added in a single plot. (default : none)",
                        default="",
                        dest="parameters")
    parser.add_argument("-m", "--metrics",
                        help="Comma-separated list of metrics to be added in a single plot.",
                        required=True,
                        dest="metrics")
    parser.add_argument("-o","--output-prefix",
                        help="Prefix used for plot file names. (default: the given file name)",
                        default="",
                        dest="prefix")
    parser.add_argument("-x","--x-scale",
                        help="Set the scaling for the x-axis. (default: linear)",
                        choices=["linear", "log"],
                        default="linear",
                        dest="xscale")
    parser.add_argument("-y","--y-scale",
                        help="Set the scaling for the y-ayis. (default: linear)",
                        choices=["linear", "log"],
                        default="linear",
                        dest="yscale")
    parser.add_argument("-s", "--subset",
                        help="Only use a data subset. You can fix a set of parameters " +
                        " by specifying key:value pairs separated by a comma.\n\t" +
                        "Example: 'par1:value1,par2:value2' (default: empty dict)",
                        default="",
                        dest="fixed_parameters"
                        )
    parser.add_argument("-l", "--layout",
                        help="Set the plot layout by giving a tuple '(nrows,ncols)'." +
                        "Use nrows=0 to add all plots to a single file (use with caution)" +
                        "(default: (15,3))",
                        default="15,3",
                        dest="layout")
    parser.add_argument("-xl", "--x-limits",
                        help="Set the maximum allowed value range for the x-axis. (default: None)",
                        default=None,
                        dest="xlimits")
    parser.add_argument("-yl", "--y-limits",
                        help="Set the maximum allowed value range for the y-axis. (default: None)",
                        default=None,
                        dest="ylimits")
    parser.add_argument("-xr", "--x-range",
                        help="Fix the value range of the x-axis, e.g. '(0,1)'. (default: None)",
                        default=None,
                        dest="xrange")
    parser.add_argument("-yr", "--y-range",
                        help="Fix the value range of the y-axis, e.g. '(0,1)'. (default: None)",
                        default=None,
                        dest="yrange")
    parser.add_argument("-pt", "--plot-type",
                        help="Set the type of plot to be created. (default: plot)",
                        choices=["plot", "scatter", "box"],
                        default="plot",
                        dest="plot_type")
    parser.add_argument("-ps", "--plot-size",
                        help="Set the size of a single plot. (default: 6)",
                        default=6,
                        dest="plot_size")
    parser.add_argument("-pp", "--paper-plots",
                        help="Use plot settings for paper grade plots.",
                        default=False,
                        action='store_true',
                        dest='paper_plots')
    parser.add_argument("-sp", "--save-pickles",
                        help="Save each plot also as a python pickle file. (default: False)",
                        default=False,
                        action="store_true",
                        dest="save_pickles")

    # Parse user arguments
    args = parser.parse_args()

    parameter_list = []
    if args.parameters:
        for entry in args.parameters.split(','):
            parameter_list.append(entry.strip())

    metric_list = []
    for entry in args.metrics.split(','):
        metric_list.append(entry.strip())

    if args.prefix:
        prefix = args.prefix
    else:
        prefix = args.file_name.split('.')[0]

    layout = literal_eval(args.layout)

    xlimits = None
    if args.xlimits:
        xlimits = literal_eval(args.xlimits)

    ylimits = None
    if args.ylimits:
        ylimits = literal_eval(args.ylimits)

    xrange = None
    if args.xrange:
        xrange = literal_eval(args.xrange)

    yrange = None
    if args.yrange:
        yrange = literal_eval(args.yrange)

    # Read data
    reader = dwm.metadated_dataframe_reader(args.file_name)
    df = reader.read_dataframe()

    plot_size = float(args.plot_size)

    # Select subset
    fixed_parameters = read_key_value_pairs(args.fixed_parameters)
    if fixed_parameters:
        df = df.xs(list(fixed_parameters.values()), level=list(fixed_parameters.keys()))

    # For the moment hard-coded: aliases for the terms used in the dataframes
    aliases = {}
    aliases['face L_inf_norm'] = '$L_\infty(e_{\kappa, rel})$'
    aliases['face L2_norm'] = '$L_2(e_{\kappa, rel})$'
    aliases['resolution'] = '$n$'
    # Stationary droplet
    aliases['Linf velocity error'] = '$L_\infty(|\mathbf{U}|)$'
    aliases['L2 velocity error'] = '$L_2(|\mathbf{U}|)$'
    aliases['time'] = '$t$'
    # Oscillating droplet
    aliases['kinetic_energy'] = '$E_{kin}$'
    aliases['semi-axes-x'] = '$s_x$'
    # Curvature data
    aliases['divGradMarkerField'] = r'DG($\alpha$)'
    aliases['divGradSignedDistance'] = r'DG($\phi$)'
    aliases['compactDivGradNoCorrection'] = r'cDG($\phi$)'
    aliases['compactDivGradSphere'] = r'sccDG($\phi$)'


    # Execute ArrayPlot class
    parameter_metric_sets = []
    set_1 = {}
    set_1['parameters'] = parameter_list
    set_1['metrics'] = metric_list
    parameter_metric_sets.append(set_1)

    # Construction
    array_plotter = dp.ArrayPlot(df, parameter_metric_sets, args.xmetric, plot_layout=layout,
                                    plot_length=plot_size, paper_plots=args.paper_plots,
                                    xscale=args.xscale, yscale=args.yscale,
                                    fixed_parameters=fixed_parameters, save_pickles=args.save_pickles, x_limits=xlimits, y_limits=ylimits,
                                    x_range=xrange, y_range=yrange, plot_type=args.plot_type, aliases=aliases)
    array_plotter.create_plots(prefix)


if __name__ == "__main__":
    main()
