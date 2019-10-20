"""

Verification of the 2D shear (also called single 2D vortex) verification case. 

Publications: 

    Comminal, R., Spangenberg, J., & Hattel, J. H. (2015). Cellwise conservative unsplit advection 
    for the volume of fluid method. Journal of Computational Physics, 283, 582–608.

""" 

from lent_error import * 
import sys

def main(): 

    templateCase = sys.argv[1] 

    # Reduce all dataFrames in all advectionErrors.dat files in all subdirectories.
    data = reduce_dataframe("advectionErrors.dat", templateCase)
    
    # Calculate and insert convergence columns 
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
