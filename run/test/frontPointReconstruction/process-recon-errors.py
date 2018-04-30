#! /usr/bin/env python
import pandas as pd
from lent_error import * 
from optparse import OptionParser
import sys
from matplotlib import pyplot as plt
from matplotlib import rcParams

rcParams["text.usetex"]=True
rcParams["figure.figsize"]=[2.36,2.36]
rcParams["font.size"]=9

usage = "usage: process-recon-errors.py -p directoryPattern -f dataFileName"

parser = OptionParser(usage=usage)

parser.add_option("-p", "--pattern", dest="pattern",
                  help="Pattern that identifies solution directories.", 
                  metavar="PATTERN")

parser.add_option("-f", "--data-file", dest="datafile", 
                  default='dual-contouring-errors.dat',
                  help="File name of the CSV file written by a test application.",
                  metavar="DATAFILE")

(options, args) = parser.parse_args()

if (options.pattern == None): 
    print ("Error: pattern option not used. Use --help option for more information.") 
    sys.exit(1)

baseFileName = options.datafile.split(".")[0] + "-" + options.pattern

# Concat all CSV dataframes in the study. 
errorDf = concat_csvdata(options.datafile, options.pattern)
errorDf.sort_values("Nc", inplace=True, ascending=True)

errorDf.to_latex("%s.tex" % baseFileName, 
                 index=False, float_format=custom_format, column_format="rlrrrr")
errorDf.to_csv("%s.csv" % baseFileName, 
               index=False)

# Get the list of unique resolutions. 
hList = errorDf["h"].unique() 

# Get the maximal error
maxErr = [] 
for hVal in hList:  
    maxErr.append(errorDf.loc[errorDf["h"] == hVal].max())

maxErrDf = pd.DataFrame(maxErr)

# Write the maximal error dataframes with convergences to files.


maxErrDf.to_csv("%s-maxErrors.csv" % (baseFileName), index=False)
maxErrDf.drop(["phi","theta"], axis=1, inplace=True)
maxErrDf.to_latex("%s-maxErrors.tex" % (baseFileName), index=False, 
                  float_format=custom_format, column_format="rlrrrr")

# Plot the max E1 and Einf convergence diagram. 
plt.loglog()
plt.plot(maxErrDf["h"],maxErrDf["Einf"],label="$E_\infty$")
plt.plot(maxErrDf["h"],maxErrDf["E2"],label="$E_1$")
plt.legend()
plt.savefig("%s-maxErrors.pdf" % (baseFileName), bbox_inches="tight")
plt.savefig("%s-maxErrors.png" % (baseFileName), bbox_inches="tight")
