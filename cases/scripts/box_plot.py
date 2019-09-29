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
    parser = ArgumentParser(description="Create one or more box plots from the data of a given parameter study.")
    parser.add_argument("-s","--study-name",
                        help="Name of the study for the box-plot.",
                        required=True,
                        dest="studyname")
    parser.add_argument("-x","--x-parameter",
                        help="Name of the parameter to be used for the x-axis",
                        required=True,
                        dest="xparameter")
    parser.add_argument("-y","--y-parameter",
                        help="Name of the metric to be plotted, e.g. the maximum error.",
                        required=True,
                        dest="yparameter")
    parser.add_argument("-f","--fix-parameters",
                        help="Fix parameters by providing parameter-value pairs in python dictionary syntax, e.g. '{'resolution':'64'}'",
                        dest="parameterdict")
    parser.add_argument("-sc","--scaling",
                        help="Set the scaling of the y-axis, either to 'log' or 'linear'. The default is 'log'",
                        choices=['log','linear'],
                        dest="scaling",
                        default="log")

    args = parser.parse_args()

    parameterSpace,studyDf = trc.agglomerate_data(args.studyname)

    fixedParameters = dict()

    if args.parameterdict is not None:
        fixedParameters = dict(json.loads(args.parameterdict))


    #---- Validation of comand line arguments --------------------------------
    if args.xparameter not in parameterSpace.columns:
        print("Error: provided x-parameter",args.xparameter,"not in parameter space.")
        print("Valid parameters are:")
        print("\t",parameterSpace.columns)
        print("Exiting.")
        sys.exit()

    if args.xparameter in fixedParameters.keys():
        print("Error: provided x-parameter",args.xparameter," is also set as fixed parameter.",
                "\n\tCannot plot over a fixed parameter.")
        print("Exiting.")
        sys.exit()

    if args.yparameter not in studyDf.columns:
        print("Error: provided metric",args.yparameter,"not available in study data.")
        print("Valid metrics are:")
        print("\t",studyDf.columns)
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


    #---- Actual plotting -----------------------------------------------------
    tp.box_plot(studyDf, args.xparameter, args.yparameter, fixedParameters, args.studyname, args.scaling)


if __name__ == "__main__":
    main()
