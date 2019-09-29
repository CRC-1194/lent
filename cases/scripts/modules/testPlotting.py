# testPlotting.py

# This required to enable the generation of pdfs on system
# without a running X Server
import matplotlib as mpl
mpl.use('pdf')

import itertools
import numpy as np
import matplotlib.pyplot as plt
import testReportCore as trc

#------------------------------------------------------------------------------
#   functions related to plotting
#------------------------------------------------------------------------------
def set_box_plot_output_params():
    # TODO: implement me
    return 0.0

def set_box_plot_display_params():
    # TODO: implement me
    return 0.0

def reduce_dataframe(dataFrame, keys, values):
    """ Return a reduced data frame with data associated with the given
        keys (aka parameters) and their corresponding values. Removes the
        fixed parameters from the multiindex.
    """
    if keys and len(values) == len(keys):
        reducedFrame = dataFrame.xs(values, level=keys, drop_level=True)
    else:
        reducedFrame = dataFrame

    return reducedFrame


def variation_vectors(dataFrame, parameters):
    """Compute the different variations by taking the Cartesian product of
        the value vector of each parameter.
        Returns a list containing these vectors.
    """
    valueVectors = list()
    for p in parameters:
        valueVectors.append(set(dataFrame.index.get_level_values(p)))

    products = list()
    for element in itertools.product(*valueVectors):
        products.append(element)

    return products


def generate_file_name(fixedParameters, parameterSet, suffix):
    """
    Generate a file name from the parameter value vector using the
    given suffix.
    """
    pdfName = ""

    for value in fixedParameters.values():
        pdfName += str(value) + '_'

    for parameterValue in parameterSet:
        pdfName += str(parameterValue) + '_'

    pdfName += '.' + suffix

    return pdfName


def box_plot(dataFrame, xParameter, yParameter, fixedParameters, studyName, yscale='log'):
    """Create a blox plot from the given data frame. Arguments of this function are:
            - dataFrame: the data frame containing the data to be plotted
            - xParameter: the parameter whose values are used to define the x axis
                (expected to be a row name)
            - yParameter: name of the column to be plotted against the x axis.
                (expected to be a column name)
            - fixedParameters: pre-filter the full data frame according to the
                parameter-value pairs given in this dictionary
    """

    # Numerical parameters must be represented by floats to ensure their values'
    # correct sequence. Otherwise, plotting errors over the mesh resolution can
    # look rather odd...
    fixedParameters = trc.convert_values_to_float_if_valid(fixedParameters)

    reducedFrame = reduce_dataframe(dataFrame, list(fixedParameters.keys()), list(fixedParameters.values()))

    freeParameters = list(reducedFrame.index.names)
    freeParameters.remove(xParameter)
    freeParameters.remove("iteration")

    # Recover the different parameter configurations
    products = variation_vectors(reducedFrame, freeParameters)

    print("Info: generating",len(products),"box plots",
            "using n =",len(freeParameters),"parameters\n")

    for parameterSet in products:
        plotFrame = reduce_dataframe(reducedFrame, freeParameters, parameterSet)

        xValues = plotFrame.index.levels[0]

        # Setup raw data for plotting
        plotData = [plotFrame.loc[i][yParameter] for i in xValues]

        # Define and set the plot properties
        medianprops = dict(linewidth=0.5, color='k')
        medianprops = dict(linestyle='-', linewidth=0.5, color='k')

        plt.boxplot(plotData, medianprops=medianprops)
        plt.xticks(np.arange(1, len(xValues) + 1), xValues, rotation=90)
        # TODO: the scaling of the y-axis should not be hard coded.
        # Rather, there should be the options "linear", "log" and "automatic"
        plt.yscale(yscale)

        plt.xlabel(xParameter)
        plt.ylabel(yParameter)
        plt.title(studyName)

        # Grid parameters
        plt.grid(True, alpha=0.3, linestyle='dashed', which='both')

        # TODO: implement option to write figures to a directory
        # controlled by an environment variable
        prefix = studyName + "_" + yParameter + "_"
        plt.savefig(prefix + generate_file_name(fixedParameters, parameterSet, "pdf"), bbox_inches='tight')
        plt.clf()


def time_evolution_plot(dataFrame, metricSet, groupParameter, fixedParameters, studyName, dependentParameters=list(),
                        timelabel='time', yscale='log'):
    """Create a time evolution plot from the given data frame. Arguments of this function are:
            - dataFrame: the data frame containing the data to be plotted
            - metricSet: set of column names to be plotted against time.
            - groupParameter: include graphs for all values this parameter assumes in
                a single plot. This parameter may be None.
            - fixedParameters: pre-filter the full data frame according to the
                parameter-value pairs given in this dictionary
            - dependentParameters: parameters which have been computed from other parameters
                and should therefore not be included in the cartesian Product.
    """

    # Numerical parameters must be represented by floats to ensure their values'
    # correct sequence. Otherwise, plotting errors over the mesh resolution can
    # look rather odd...
    fixedParameters = trc.convert_values_to_float_if_valid(fixedParameters)

    reducedFrame = reduce_dataframe(dataFrame, list(fixedParameters.keys()), list(fixedParameters.values()))

    freeParameters = list(reducedFrame.index.names)
    if groupParameter:
        freeParameters.remove(groupParameter)
    freeParameters.remove("iteration")

    # Remove dependent parameters from freeParameters
    if dependentParameters:
        for dParam in dependentParameters:
            freeParameters.remove(dParam)

    # Recover the different parameter configurations defining the number
    # of plots
    products = variation_vectors(reducedFrame, freeParameters)

    print("Info: generating",len(products),"time evolution plots",
            "using n =",len(freeParameters),"parameters\n")

    if groupParameter:
        freeParameters.append(groupParameter)
        groupParameterValues = set(reducedFrame.index.get_level_values(groupParameter))

    for parameterSet in products:
        
        # Prepare new plot: set the plot parameters, etc
        plt.yscale(yscale)

        plt.xlabel(timelabel + " [s]")
        plt.ylabel("Metric value")

        plt.title(studyName.split('.')[0])

        # Grid parameters
        plt.grid(True, alpha=0.3, linestyle='dashed', which='both')
        linewidth=0.8


        if groupParameter:
            for groupValue in groupParameterValues:
                tmpParameterSet = list(parameterSet)
                tmpParameterSet.append(groupValue)

                plotFrame = reduce_dataframe(reducedFrame, freeParameters, tmpParameterSet)

                for metric in metricSet:
                    graphLabel = metric + " | " + str(groupValue)
                    plt.plot(plotFrame.loc[:][timelabel], plotFrame.loc[:][metric], label=graphLabel, linewidth=linewidth)
        else:
            plotFrame = reduce_dataframe(reducedFrame, freeParameters, parameterSet)

            for metric in metricSet:
                plt.plot(plotFrame.loc[:][timelabel], plotFrame.loc[:][metric], label=metric, linewidth=linewidth)

        plt.legend(loc='best')

        # TODO: implement option to write figures to a directory
        # controlled by an environment variable
        prefix = studyName + "_"
        for metric in metricSet:
            prefix += metric + "_"
        if groupParameter:
            prefix += groupParameter
        plt.savefig(prefix + generate_file_name(fixedParameters, parameterSet, "pdf"), bbox_inches='tight')
        plt.clf()
