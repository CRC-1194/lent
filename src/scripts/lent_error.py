import os
from shutil import copyfile
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.Basics.Utilities import Utilities
from shutil import copyfile
import pandas as pd
from pandas import DataFrame
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.lines as mlines 
from math import log

def table_file_name(templateCase):
    """Reads syste/geomTransportControlDict, constructs table file name\
    from key/values it finds there"""

    tableFile = ""

    try:
        tableFile = os.environ["TABLE_PATH"]
    except:
        print ("Error: define TABLE_PATH environmental variable where tables are to be stored.")
        return

    lentSolutionDict = ParsedParameterFile(os.path.join(templateCase,"system/lentSolution"))

    marker=lentSolutionDict["markerFieldModel"]["type"]
    temporal=lentSolutionDict["frontMotionSolver"]["type"]
    interpol=lentSolutionDict["frontMotionSolver"]["cellToVertexInterpolation"]


    tableFile = os.path.join(tableFile,templateCase + "-" + temporal + \
                             "-" + interpol + "-" + marker)

    print(tableFile)

    return tableFile

def scientific(v, precision=''):
    """Convert number smaller than 1 to a sientific formatted string."""

    if isinstance(v, np.float64) and v < 0.1 and v > 1e-17:
        return "%.2e" % v
    elif isinstance(v, np.float64) and v > 0.1:
        return "%.2f" % v
    elif isinstance(v, np.float64) and v < 0:
        return "%.2f" % v
    else: 
        return "%s" % v 

columnFormats = {
    "CFL" : '{:,.2f}'.format,
    "Ev" : '{:,.2e}'.format,
    "Eb" : '{:,.2f}'.format,
    "Eg" : '{:,.2e}'.format,
    "O(Eg)" : '{:,.2e}'.format,
    "En" : '{:,.2e}'.format,
    "O(En)" : '{:,.2e}'.format,
    "Te" : '{:,.2e}'.format,
    "Tr" : '{:,.2e}'.format
}

def subdirs(pattern): 
    """Returns the sub-directories of the current directory whose names end
    with pattern."""
    
    return [directory for directory in os.listdir(os.curdir) 
            if os.path.isdir(directory) and directory.endswith(pattern)]

def filter_subdirs(fileName, dirName): 
    """Extracts only sub-directories from the current directory that contain a
    file fileName."""
  
    directories = subdirs(dirName)

    parameterDirectories = list(
        filter(lambda directory : fileName in os.listdir(directory)
               and os.path.isfile(os.path.join(directory, fileName)), 
               directories)
    )
    parameterDirectories.sort()

    return parameterDirectories

def extract_files(fileName):
    """Extracting the files and naming them with increasing indices."""

    directories = filter_subdirs(fileName)

    splitFile = fileName.split(".")

    strippedFile = fileName.rstrip(splitFile[-1]);

    for i in range(len(directories)):
        copyfile(os.path.join(directories[i], fileName), 
                 os.path.join(os.curdir, strippedFile + "%05d." % i + splitFile[-1]))

def extractColumns(fileName): 

    directories = filter_subdirs(fileName)

    datFile = open(os.path.join(directories[0],fileName), 'r')

    # NOTE: Assumes that the first line contains the columns. 
    return datFile.readline().rstrip(os.linesep).rstrip(' ').split(' ')

def parse_headerless_file(paramFileName):
    """Extracts the parameter vectors from a parameter file."""

    # Write the header to the parameter file.
    paramFileNameCopy = paramFileName + ".copy"
    paramFileCopy = open(paramFileNameCopy, 'w')
    utils = Utilities()
    utils.writeDictionaryHeader(paramFileCopy)
    paramFileCopy.close()
    with open(paramFileName) as pFile:
        with open(paramFileNameCopy, "a") as pFileCopy:
            for line in pFile:
                pFileCopy.write(line) 
        
    # Parse the parameter file copy.
    paramFileCopy = ParsedParameterFile(paramFileNameCopy)
    
    return paramFileCopy

def parameter_names(parsedParamFile):
    """Extract parameter names from the parsed parameter file."""
    # Get the parameter names. In this case, 'solver' is not a parameter.
    return [x for x in parsedParamFile["values"].keys() if 'solver' not in x]

def parameter_values(parsedParamFile, paramNames):
    """Extract parameter values from the parsed parameter file."""
    return [list(parsedParamFile["values"][name]) for name in paramNames]

def process_data_file(dataFileName, functions):
    paramFileNameCopy = dataFileName + ".copy"
    paramFileCopy = open(paramFileNameCopy, 'w')
    """Process data frame columns with a given list of functions to create a pandas 
    series."""

    processedColumns = []

    try:
        dFrame = pd.read_table(dataFileName,sep=' ')

        processedColumns = [functions[index](dFrame[dFrame.columns[index]])
                            for index in range(len(functions))] 
        return pd.Series(processedColumns,dFrame.columns[0:len(functions)])
    except:
        return pd.Series()

    

def process_param_values(casePath, paramNames):
    """Extract local parameter values from the
    casePath/PyFoamPrepareCaseParameters file."""
    
    filePath = os.path.join(casePath, "PyFoamPrepareCaseParameters")
    paramFile = parse_headerless_file(filePath) 
    return [paramFile[key] for key in paramNames]


def default_functions():
    """A default set of functions used to reduce an advection error
    dataframe."""

    last_row_element = lambda x : x.iloc[-1]

    return list([last_row_element, DataFrame.max, DataFrame.max, DataFrame.max, 
            last_row_element, last_row_element, DataFrame.mean, DataFrame.mean])

def reduce_dataframe(dataFileName, templateCase, 
                     functionList = default_functions()):
    """Reduces all dataframes stored in dataFileName in all sub-directories
    into a single dataFrame using a set of functions performed on the dataFrame
    columns.""" 

    paramDirs = filter_subdirs(dataFileName, templateCase)
    
    # Build a pandas MultiIndex from a cartesian product of parameter value vectors. 
    globParamFile = parse_headerless_file(templateCase + ".parameter")
    paramNames = parameter_names(globParamFile)
    paramVectors = parameter_values(globParamFile, paramNames)
    paramIndex = pd.MultiIndex.from_product(paramVectors, names=paramNames)

    # Get the columns from a dat file 

    # Build an empty multi-indexed DataFrame.
    resultDataFrame = DataFrame(index=paramIndex, 
                                columns=["CFL", "Ev", "Eb", "Eg", "En", "Te"],
                                dtype=np.float64)
      
    # !Fill the parameter data from the processed data files.
    for paramDir in paramDirs: 
        # Open the data file in the director
        dataPathName = os.path.join(paramDir, dataFileName)
        dataFile = open(dataPathName,'r')
        # Process the data file and generate the series. 
        dataSeries = process_data_file(dataPathName, functionList)
        
        # Process the data parameters to read the series location.  
        paramValues = process_param_values(paramDir, paramNames) 
        
        # Insert the series into the parameterData at the appropriate location.
        # TODO: Introduce variable number of parameters here.
        try:
            resultDataFrame.loc[tuple(paramValues)] = dataSeries 
        except:
            pass

    return resultDataFrame 

def calc_convergence(series):
    """Calculates convergence from a pandas column (DataSeries)."""

    for i in range(len(series)-1):
        j = series.index[i] 
        k = series.index[i+1] 
        series[j] = log(series[j] / series[k]) / log(k / j)
        
    series[series.index[-1]] = np.nan

    return series 

# Generalize for Eg and En names. En -> O(En)  
def add_convergence(dframe): 
    """Adds convergence order columns to an advection error dataframe."""
    
    # Insert the convergence orders columns for Eg (Og) and En (On)
    dframe.insert(dframe.columns.get_loc('Eg')+1,'O(Eg)', dframe['Eg'])
    dframe.insert(dframe.columns.get_loc('En')+1,'O(En)', dframe['En'])
    
    for i in dframe.index.levels[0]:
        dframe['O(Eg)'][i] = calc_convergence(dframe['O(Eg)'][i])
        dframe['O(En)'][i] = calc_convergence(dframe['O(En)'][i])

def set_inline_style():
    """Sets the matplotlib plotting parameters for inline plots."""
    
    plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
    params = {'text.usetex' : True,
              'font.size' : 28,
              'text.latex.unicode': True,
              'figure.figsize' : (15, 10),
              'savefig.dpi' : 200
              }
    plt.rcParams.update(params)

def set_file_style():
    """Sets the matplotlib ploting parameters for file plots."""
    
    plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
    params = {'text.usetex' : True,
              'font.size' : 9,
              'text.latex.unicode': True,
              'figure.figsize' : (5.8, 3.9),
              'savefig.dpi' : 200
              }
    plt.rcParams.update(params)
    
def plot_translation2D_refdata():
    """Plots the data from
    
    López, J., Hernández, J., Gómez, P., & Faura, F. (2004). A volume of fluid
    method based on multidimensional advection and spline interface
    reconstruction. Journal of Computational Physics, 195(2), 718–742.
    http://doi.org/10.1016/j.jcp.2003.10.030
    
    Diagram data extracted with Engauge in the 'publication-data' directory.  
    
    """
    
    refData = pd.read_csv('publication-data/EMFPA-SIR-translation2D-diagram.csv')
    
    plt.plot(refData["x"], refData["EMFPA-SIR-10"],  color='k', 
             marker='o', fillstyle='none', markersize=4)
    plt.plot(refData["x"], refData["EMFPA-SIR-20"], color='k', 
             marker='o', fillstyle='none', markersize=4)
    plt.plot(refData["x"], refData["EMFPA-SIR-40"], color='k', 
             marker='o', fillstyle='none', markersize=4)

    xy10 =(1.025 * refData["x"].iloc[-1], refData["EMFPA-SIR-10"].iloc[-1])
    xy20 =(1.025 * refData["x"].iloc[-1], refData["EMFPA-SIR-20"].iloc[-1])
    xy40 =(1.025 * refData["x"].iloc[-1], refData["EMFPA-SIR-40"].iloc[-1])

    plt.annotate('$10^2$ mesh', xy10) 
    plt.annotate('$20^2$ mesh', xy20) 
    plt.annotate('$40^2$ mesh', xy40) 

def plot_translation2D(translationData2D): 
    """Plots the voFoam translation2D verification case results."""
    
    CFLs = translationData2D.index.levels[0]

    minEg = 1e15 
    maxEg = 0

    for CFL in CFLs: 
        if (CFL < 1):
            Eg = translationData2D["Eg"][CFL]
            xCFL = [1.0 / CFL for x in Eg] 
            # minEg = np.min(minEg,Eg.min())
            # maxEg = np.max(maxEg,Eg.max())
            plt.plot(xCFL, Eg, ls="none",marker='^',color='k', markersize=4)

    plt.ylim([1e-04,1.5e-02])
    
def plot_translation2D_file(translationData2D, templateCase):
    """Plots the translation2D verification case data diagram into a file."""
    
    # Clear the figure from previous plots.
    plt.clf()
    
    # Set the LaTex style for file output.
    set_file_style()
    
    # Set axis labels
    plt.ylabel("$Eg$")
    plt.xlabel("$CFL^{-1}$")
    plt.semilogy()
    
    # Plot reference data
    plot_translation2D_refdata()
    
    # Plot voFoam data
    plot_translation2D(translationData2D)

    # Figure out the reconstruction method.  
    geomControl = ParsedParameterFile(os.path.join(templateCase, 
                                      "system/geomTransportControlDict"))
    recLabel = geomControl["reconstruction"] 
    
    # Configure the legend
    circle_line = mlines.Line2D([], [], color='k', marker='o',
                          markersize=4, fillstyle='none', label='EMFPA-SIR')
    square_line = mlines.Line2D([], [], ls="none", color='k', marker='^',
                          markersize=4, label="UFVFC-"+recLabel)
    plt.legend(handles=[circle_line,square_line], 
               bbox_to_anchor=(0, 1, 0.5, 1.21), loc="lower left")

    # Save the latex table 

    # En and O(En) are not needed for this test case in the thesis.
    translationData2D.drop('CFL', axis=1, inplace=True)
    translationData2D.drop('O(En)', axis=1, inplace=True)
    translationData2D.drop('En', axis=1, inplace=True)
    translationData2D.fillna('-')
    latexData = translationData2D.to_latex(float_format=scientific)

    latexData = latexData.replace("nan",'-')

    thesis = os.environ['THESIS']
    tableFile = open(os.path.join(thesis,"tables",templateCase + "-" + recLabel + "-table.tex"),"w")
    tableFile.write(latexData)

    # Save the diagram
    diagramPath = os.path.join(thesis, "figures",templateCase + "-" + recLabel + ".pdf")
    print ("Saving plot as : %s" % diagramPath) 
    plt.savefig(diagramPath,bbox_inches='tight') 

    return latexData

def p_series(inputDataFrame, xName, yName, p=1):
    """Computes p-th order convergent [x0,x1], [y0,y1] lists from inputDataFrame."""

    x0 = inputDataFrame[xName][inputDataFrame.index[0]]
    y0 = inputDataFrame[yName][inputDataFrame.index[0]]
     
    # Make sure that y0 > machine tolerance for data frames that contain zero values.
    for index, item in enumerate(inputDataFrame[yName]):
        if y0 > 1e-15: 
            break
        y0 = inputDataFrame[yName][inputDataFrame.index[index]]
        x0 = inputDataFrame[xName][inputDataFrame.index[index]]


    x1 = inputDataFrame[xName][inputDataFrame.index[-1]]
    y1 = y0 / (x1 / x0)**p 

    return { 'x' : [x0,x1], 'y' : [y0,y1] }

def plot_reconstruction2D_file(reconstructionData2D,outputFileName):
    """Plots the reconstruction2D verification case data diagram into a file."""
    

    # Plot the L1 error  

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)

    set_file_style()

    ax1.set_ylabel("$E_{vsd}$")
    ax1.set_xlabel("$N_c$")
    ax1.loglog()

    ax1.plot(reconstructionData2D["N_c"],reconstructionData2D["E_{r,Y}^1"], label="$E_{r,Y}^1$", 
            ls="none", color='k', marker='^')
    ax1.plot(reconstructionData2D["N_c"],reconstructionData2D["E_{r,S}^1"], label="$E_{r,S}^1$", 
            ls="none", color='k', marker='o')

    firstOrder = p_series(reconstructionData2D, "N_c", "E_{r,Y}^1")
    secondOrder = p_series(reconstructionData2D, "N_c", "E_{r,Y}^1", p=2)

    ax1.plot(firstOrder['x'], firstOrder['y'], label = "$O(h)$", ls='dashed', color='k')
    ax1.plot(secondOrder['x'], secondOrder['y'], label = "$O(h^2)$", ls='dotted', color='k')

    ax1.legend(loc="upper right", ncol=2) #bbox_to_anchor=(0.75,1), ncol=2)

    # FIXME Alternative file names for L1 and Linfty
    fig1.savefig(outputFileName)

    # Plot the L1 error  

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)

    set_file_style()

    ax2.set_ylabel("$E_{vsd}$")
    ax2.set_xlabel("$N_{i}$")
    ax2.loglog()

    ax2.plot(reconstructionData2D["N_c"],reconstructionData2D[r"E_{r,Y}^\infty"], label=r"$E_{r,Y}^\infty$", 
            ls="none", color='k', marker='^')
    ax2.plot(reconstructionData2D["N_c"],reconstructionData2D[r"E_{r,S}^\infty"], label=r"$E_{r,S}^\infty$", 
            ls="none", color='k', marker='o')

    firstOrder = p_series(reconstructionData2D, "N_c", r"E_{r,Y}^\infty")
    secondOrder = p_series(reconstructionData2D, "N_c", r"E_{r,Y}^\infty", p=2)

    ax2.plot(firstOrder['x'], firstOrder['y'], label = "$O(h)$", ls='dashed', color='k')
    ax2.plot(secondOrder['x'], secondOrder['y'], label = "$O(h^2)$", ls='dotted', color='k')

    ax2.legend(loc="upper right", ncol=2) #bbox_to_anchor=(0.75,1), ncol=2)

    fig2.savefig(outputFileName)
