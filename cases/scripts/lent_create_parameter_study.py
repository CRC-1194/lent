#!/usr/bin/env python3

import sys

import parameterStudyPreparation as psp

def printUsageAndExit():
    print("Combines two or more parameter files into single one. Each ",
            "file is given as a command line argument, e.g.\n",
            "\tlent_create_parameter_study.py file1 file2 ... fileN")
    sys.exit()

def main():

    # Read command line arguments, skip script name
    args = sys.argv[1:]

    if len(args) < 2:
        printUsageAndExit()

    psp.create_parameter_file_from_list(args)


if __name__ == "__main__":
    main()
