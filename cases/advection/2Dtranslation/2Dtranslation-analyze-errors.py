#!/usr/bin/env python

"""
Data processing and plotting functions for the simple translation verification case.

"""

from lent_error import *
import sys
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

def main():
    # Reduce all dataFrames in all advectionErrors.dat files in all subdirectories.

    templateCase = "lent-hex-translation2D"

    data = reduce_dataframe("advectionErrors.dat", templateCase)
    add_convergence(data)

    # En and O(En) are not needed for this test case in the thesis.
    data.drop('CFL', axis=1, inplace=True)
    data.drop('O(En)', axis=1, inplace=True)
    data.drop('En', axis=1, inplace=True)
    latexData = data.to_latex(float_format=scientific)
    latexData = latexData.replace("nan",'-')

    print(latexData)

    tableFileName = table_file_name(templateCase) 

    latexFile = open(tableFileName + ".tex", "w")
    latexFile.write(latexData)

    csvFile = open(tableFileName + ".csv", 'w')
    csvFile.write(data.to_csv()) 

if __name__=="__main__":
    main()
