# notebookPlotting.py

# This module contains functions used to plot data in jupyter notebooks.

import matplotlib.pyplot as plt
from matplotlib import lines
from matplotlib.lines import Line2D
from matplotlib import rcParams
import glob

import dataAgglomeration as datglom

def plot_dframe(dFrame, dFrameAgglomerator, title="", colDict = {}):
    """Plot the agglomerated multidimensional DataFrame using mathematical symbols mapped
    to column names with colDict."""

    # Set plotting style parameters
    lStyles = list(lines.lineStyles.keys())
    lStyles = set(lStyles)
    lStyles.discard('None')
    lStyles = list(lStyles)
    mStyles = list(Line2D.markers.keys())
    mStyles = set(mStyles)
    mStyles.discard('None')
    mStyles = list(mStyles)
    
    # X and Y-axis column names and diagram symbols.
    xColName = colDict["x"]
    xColSymb = colDict["xsymb"]
    yColName = colDict["y"]
    yColSymb = colDict["ysymb"]
    
    indexLevels = [i for i,name in enumerate(dFrame.index.names) if 'step' not in name]
    collector = dFrameAgglomerator.data_collector
    variations = collector.valid_variations
    fig, ax = plt.subplots()
    ax.set_xlabel(xColSymb)
    ax.set_ylabel(yColSymb)
    
    ax.set_title(title)
    variationI = 0 
    for paramLine, subDframe in dFrame.groupby(level=indexLevels):
        xCol = subDframe[xColName] 
        yCol = subDframe[yColName]
        
        # Build the parameter string paramName=paramValue for the plot legend. 
        paramString = ""
        if (len(indexLevels) > 1):
            for levelI,paramName in enumerate(paramLine):
                indexName = dFrame.index.names[levelI]
                paramString = paramString + indexName \
                              + "=%.1e " % paramName 
        else:
            indexName = dFrame.index.names[0]
            paramString = indexName + "=%d " % paramLine 
            
        ax.plot(xCol, yCol, label="variation=%04d " % variationI + paramString, 
                marker=mStyles[variationI % len(mStyles)], 
                linestyle=lStyles[variationI % len(lStyles)])

        variationI = variationI + 1

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=1) 
              
    plt.show()
