#!/usr/bin/env python3

from argparse import ArgumentParser
# Required for a workaround to specify python dictionaries as a
# command line parameter
import json
import sys
import testReportCore as trc
import testPlotting as tp

def main():

    #---- Command line arguments ----------------------------------------------
    parser = ArgumentParser(description="Create one or more time evolution plots for a given parameter study A column named 'time' must be present in the data.")
    parser.add_argument("-s","--study-name",
                        help="Name of the study for the time evolution plot.",
                        required=True,
                        dest="studyname")
    parser.add_argument("-m","--metrics",
                        help="Name(s) for one or more metrics to be plotted against time.",
                        required=True,
                        nargs='+',
                        dest="metrics")
    parser.add_argument("-g","--group-parameter",
                        help="If present, the graphs of all values for this parameter are put into a single plot.",
                        dest="gparameter")
    parser.add_argument("-f","--fix-parameters",
                        help="Fix parameters by providing parameter-value pairs in python dictionary syntax, e.g. '{'resolution':'64'}'",
                        dest="parameterdict")
    parser.add_argument("-sc","--scaling",
                        help="Set the scaling of the y-axis, either to 'log' or 'linear'. The default is 'log'",
                        choices=['log','linear'],
                        dest="scaling",
                        default="log")
    parser.add_argument("-dd", "--drop-dependent-parameters",
                        help="Remove parameters computed from other parameters from the list of free parameters",
                        dest="dParameters",
                        default=None)

    args = parser.parse_args()

    parameterSpace,studyDf = trc.agglomerate_data(args.studyname)

    fixedParameters = dict()
    dependentParameters = list()

    if args.parameterdict is not None:
        fixedParameters = dict(json.loads(args.parameterdict))


    #---- Validation of comand line arguments --------------------------------
    for metric in args.metrics:
        if metric not in studyDf.columns:
            print("Error: provided metric", metric, "not available in data set.")
            print("Valid metrics are:")
            print("\t",studyDf.columns)
            print("Exiting.")
            sys.exit()

    if args.gparameter:
        if args.gparameter in fixedParameters.keys():
            print("Error: provided group parameter",args.gparameter," is also set as fixed parameter.",
                    "\n\tCannot fix a group parameter.")
            print("Exiting.")
            sys.exit()

    if args.gparameter:
        if args.gparameter not in parameterSpace.columns:
            print("Error: given group parameter",args.gparameter,"not in parameter space.")
            print("Valid parameters are:")
            print("\t",parameterSpace.columns)
            print("Exiting.")
            sys.exit()

    if fixedParameters:
        for key in fixedParameters:
            if key not in parameterSpace.columns:
                print("Error: provided fixed parameter",key,"not in parameter space.")
                print("Valid parameters are:")
                print("\t",parameterSpace.columns)
                print("Exiting.")
                sys.exit()

    if args.dParameters:
        dependentParameters = args.dParameters.split(',')

        for key in dependentParameters:
            if key not in parameterSpace.columns:
                print("Error: provided dependent parameter",key,"not in parameter space.")
                print("Valid parameters are:")
                print("\t",parameterSpace.columns)
                print("Exiting.")
                sys.exit()

    # Ensure that 'time' is present as a column
    if 'time' not in studyDf.columns:
        print("Error: study data does not contain a column named 'time'")
        print("Cannot plot over time")
        print("Exiting.")
        sys.exit()



    #---- Actual plotting -----------------------------------------------------
    tp.time_evolution_plot(studyDf, args.metrics, args.gparameter, fixedParameters, args.studyname, dependentParameters=dependentParameters, yscale=args.scaling)


if __name__ == "__main__":
    main()
