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


def custom_format(v, precision=''):  
    """Convert number smaller than 1 to a sientific formatted string."""

    if isinstance(v, np.float64) and v < 0.1 and v > 1e-17:
        return "%.4e" % v
    elif isinstance(v, np.float64) and (v > 0.1) and (v < 10):
        return "%.4f" % v
    elif isinstance(v, np.float64) and (v < 0):
        return "%.4f" % v
    elif isinstance(v, np.float64) and (v > 10):
        return "%d" % round(v)
    else: 
        try:  
            fv = np.float64(v)
            return "%.4f" % fv 
        except(): 
            return "%s" % v 

def subdirs(pattern): 
    """Returns the sub-directories of the current directory whose names end
    with pattern."""
    
    return [directory for directory in os.listdir(os.curdir) 
            if os.path.isdir(directory) and 
            (directory.endswith(pattern) or directory.startswith(pattern))]

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

def concat_csvdata(dataFileName, casePattern):
    """Reduces all DataFrames stored in 'dataFileName.dat' as CSV in all
    sub-directories that have 'casePattern' in their name, into a single
    dataFrame using a set of functions performed on the dataFrame columns.""" 

    dataDirs = filter_subdirs(dataFileName, casePattern)

    concatDf = pd.DataFrame()

    for dataDir in dataDirs: 

        dataFilePath = os.path.join(dataDir, dataFileName)  
        df = pd.read_csv(dataFilePath)
        concatDf = pd.concat([df, concatDf])
        
    return concatDf 

def add_convergence_orders(dataFrame, hName, columnNames):
    """Adds the convergence orders for a set of column names.""" 

    # Column name that contains mesh resolution values.
    loghDiff = np.log(dataFrame[hName]).diff()

    for columnName in columnNames:

        logErrDiff = np.log(dataFrame[columnName]).diff()
        dataFrame.insert(dataFrame.columns.get_loc(columnName) + 1, 
                         "O(%s)" % columnName,
                         logErrDiff / loghDiff)
    dataFrame.fillna('-', inplace=True)
