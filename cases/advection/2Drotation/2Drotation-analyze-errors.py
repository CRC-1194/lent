#!/usr/bin/env python

"""
Data processing and plotting functions for the simple rotation verification case. 
  
Rotation case data is taken from the following publications:

Liovic, P., Rudman, M., Liow, J. L., Lakehal, D., & Kothe, D. (2006). A 3D
unsplit-advection volume tracking algorithm with planarity-preserving interface
reconstruction. Computers and Fluids, 35(10), 1011–1032.
http://doi.org/10.1016/j.compfluid.2005.09.003

López, J., Hernández, J., Gómez, P., & Faura, F. (2004). A volume of fluid
method based on multidimensional advection and spline interface reconstruction.
Journal of Computational Physics, 195(2), 718–742.
http://doi.org/10.1016/j.jcp.2003.10.030

""" 

from lent_error import * 
import sys
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

def main(): 
    # Reduce all dataFrames in all advectionErrors.dat files in all subdirectories.

    templateCase = "lent-hex-rotation2D"

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
