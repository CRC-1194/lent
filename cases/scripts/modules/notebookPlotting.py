# notebookPlotting.py

# This module contains functions used to plot data in jupyter notebooks.

import matplotlib.pyplot as plt
from matplotlib import lines
from matplotlib.lines import Line2D
from matplotlib import rcParams
import glob

import dataAgglomeration as datglom

def plot_study(xy_filter=lambda y : True, paramFile="", dirExample="", dataFile="", colDict = {}):
    """Agglomerate all dataFile files in each directory whose name contains parameterFile, 
    and plot xColName column and yColName column, label the graph with the index of the 
    variation and the parameters of the variation in the plot legend."""

    # Set plotting style parameters
    lStyles = list(lines.lineStyles.keys())
    mStyles = list(Line2D.markers.keys())
    
    # X and Y-axis column names and diagram symbols.
    xColName = colDict["x"]
    xColSymb = colDict["xsymb"]
    yColName = colDict["y"]
    yColSymb = colDict["ysymb"]
    
    agglomerator = datglom.data_agglomerator(paramFile, dirExample, dataFile)
    studyDframe = agglomerator.study_dataframe()
    indexLevels = [i for i,name in enumerate(studyDframe.index.names) if 'step' not in name]
    collector = agglomerator.data_collector
    variations = collector.valid_variations
    fig, ax = plt.subplots()
    ax.set_xlabel(xColSymb)
    ax.set_ylabel(yColSymb)
    
    ax.set_title("%s" % dirExample.split('.')[0] + "-" + dirExample.split('_')[-1]) 
    # FIXME: Directory name contains an prefix - fetch the template copy?
    variationI = 0 
    for paramLine, subDframe in studyDframe.groupby(level=indexLevels):
        xCol = subDframe[xColName] 
        yCol = subDframe[yColName]
        
        # Build the parameter string paramName=paramValue for the plot legend. 
        paramString = ""
        if (len(indexLevels) > 1):
            for levelI,paramName in enumerate(paramLine):
                indexName = studyDframe.index.names[levelI]
                paramString = paramString + indexName \
                              + "=%d " % paramName 
        else:
            indexName = studyDframe.index.names[0]
            paramString = indexName + "=%d " % paramLine 
            
        if (xy_filter(xCol, yCol)):
            ax.plot(xCol, yCol, label="variation=%04d " % variationI + paramString, 
                    marker=mStyles[variationI % len(mStyles)], 
                    linestyle=lStyles[variationI % len(lStyles)])
        variationI = variationI + 1

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=1) 
              
    plt.show()
