# dataPlotting.py

import itertools
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pgf import FigureCanvasPgf
import numpy
import pandas as pd
import pickle as pl

import dataframeWithMetadata as dwm

class single_plot_creator:
    # TODO: add option to transform x-axis data via a passed function,
    # e.g. to compute the mesh spacing from the resolution

    def __init__(self, dataframe, metric_names, x_axis_name, parameter_vector={},
                    aliases={}, max_markers_per_graph=10, xscale='linear',
                    yscale='linear', draw_title=True, x_limits=None, y_limits=None,
                    x_range=None, y_range=None, plot_type='plot', order_indicators=True):
        """
        Constructor requiring a dataframe, a list of metric names and the name
        of the multiindex level or column which shall be used as x-axis.

        The purpose of this class is to plot one or more graphs from the given
        data in a matplotlib axes object or in other words, it is responsible
        for creating a single plot from the given data.

        Keyword arguments:
        parameter_vector -- a dictionary whose key-value-pairs are written in the plot
                            title. Useful if the passed dataframe is a subset
                            of a larger dataframe (default: empty dict)
        aliases -- dictionary with aliases for names used in the dataframe. If a name
                   is present, its alias is used for labelling of the plot.
                   (default : empty dict)
        max_markers_per_graph -- do not plot markers if number of data points exceeds
                                 this limit. Set to zero to disable markers. (default : 10)
        xscale -- scaling of x-axis as supported by matplotlib. (default : 'linear')
        yscale -- scaling of y-axis as supported by matplotlib. (default : 'linear')
        draw_title -- draw the parameter vector as plot title. (default : True)
        x_limits -- limit x-axis to the given range (default : None)
        y_limits -- force y-axis limits to be in the given value range. (default : None)
        x_range -- set x-axis to the specified range (default : None)
        y_range -- set y-axis to the specified range (default : None)
        plot_type -- Draw data as the given type of plot [plot, scatter, box]. (default : 'plot')
        order_indicators -- Draw lines indicating first and second order convergence in
                            log-log plots. (default: True)
        """
        self.dataframe = dataframe
        self.metric_names = metric_names
        self.x_axis_name = x_axis_name
        self.parameter_vector = parameter_vector
        self.aliases = aliases
        self.max_markers_per_graph = max_markers_per_graph
        self.marker = itertools.cycle(['.', 'o', '*', 'v', '1', 's', '+', 'x', 'd', 'p'])
        self.xscale = xscale
        self.yscale = yscale
        self.draw_title = draw_title
        self.x_limits = x_limits
        self.y_limits = y_limits
        self.x_range = x_range
        self.y_range = y_range
        self.plot_type = plot_type
        self.order_indicators = order_indicators

    def alias(self, name):
        """Return the alias of the name if present in the alias dict.
            Otherwise return the name itself.
        """
        if name in self.aliases:
            return self.aliases[name]
        else:
            return name

    def generate_legend_entry(self, variation_vector, metric_name):
        """
        Generate a legend entry and return it as a string

        Adds the metric name and the values of the variation vector to
        the legend entry if they are non-unique for the plot.
        """
        entry = ""

        # Only include the metric name in the legend entry, when multiple
        # metrics are plotted in the same axes object.
        # Otherwise, the metric name can be written as a y-axis label
        if len(self.metric_names) > 1:

            entry = self.alias(metric_name) + " | "

        for value in variation_vector.values():
            entry += str(self.alias(value)) + ", "

        # Remove comma and whitespace
        entry = entry[:-2]

        return entry

    def limit_y_axis_range(self, axes):
        """
        Enforce user-defined limits on the range of y-axis values.
        """
        if self.y_range:
            axes.set_ylim(self.y_range)
        elif self.y_limits:
            y_min = max(self.y_limits[0], axes.get_ylim()[0])
            y_max = min(self.y_limits[1], axes.get_ylim()[1])

            axes.set_ylim((y_min,y_max))

    def limit_x_axis_range(self, axes):
        """
        Enforce user-defined limits on the range of x-axis values.
        """
        if self.x_range:
            axes.set_xlim(self.x_range)
        if self.x_limits:
            x_min = max(self.x_limits[0], axes.get_xlim()[0])
            x_max = min(self.x_limits[1], axes.get_xlim()[1])

            axes.set_xlim((x_min,x_max))

    def set_axes_labelling(self, axes, user_defined_xlabel, user_defined_ylabel):
        """
        Set the axis labels to user defined values or derive them from dataframe.
        """
        # X axis
        xlabel = ""

        if user_defined_xlabel:
            xlabel = user_defined_xlabel
        else:
            xlabel = self.alias(self.x_axis_name)

        # Y axis
        ylabel = ""

        if user_defined_ylabel:
            ylabel = user_defined_ylabel
        elif len(self.metric_names) > 1:
            ylabel = "Metric value (see legend)"
        else:
            ylabel = self.alias(self.metric_names[0])

        axes.set(xlabel=xlabel, ylabel=ylabel)

    def set_plot_title(self, axes):
        """
        Set the title of the axes object based the classes parameter vector.
        """
        title = ""

        for key, value in self.parameter_vector.items():
            title += str(key) + " : " + str(value) + "\n"

        title = title[:-1]

        axes.set_title(title, loc="left")

    def plot_metric_graphs(self, axes, dataframe, variation_vector={}):
        """
        Add a graph for each metric given in the class' metric names.
        """
        if self.x_axis_name in dataframe.columns:
            x_data = dataframe.loc[:][self.x_axis_name]
        else:
            x_data = list(dataframe.index.data)

        if len(x_data) > self.max_markers_per_graph:
            marker_set = itertools.cycle([''])
        else:
            marker_set = self.marker

        for metric in self.metric_names:
            legend_entry = self.generate_legend_entry(variation_vector, metric)
            axes.plot(x_data, dataframe.loc[:][metric], label=legend_entry, marker=next(marker_set))

    def scatter_plot(self, axes, dataframe, variation_vector={}):
        # Assumption: only two levels left in the multiindex:
        #   1) the level containing the x-values
        #   2) the level containing the 'experiment' number
        x_values = dataframe.index.levels[0]

        for metric in self.metric_names:
            x_list = []
            y_list = []

            for x in x_values:
                y_values = dataframe.loc[x][metric]

                for y in y_values:
                    x_list.append(x)
                    y_list.append(y)
            #plot_data = [dataframe.loc[i][metric] for i in x_values]
            #print(plot_data)
            legend_entry = self.generate_legend_entry(variation_vector, metric)
            axes.scatter(x_list, y_list, label=legend_entry)
            axes.set_xticks(x_values)
            axes.set_xticklabels(x_values)

    def box_plot(self, axes, dataframe, variation_vector={}):
        # Assumption: only two levels left in the multiindex:
        #   1) the level containing the x-values
        #   2) the level containing the 'experiment' number
        x_values = dataframe.index.levels[0]

        for metric in self.metric_names:
            plot_data = [dataframe.loc[i][metric] for i in x_values]
            #legend_entry = self.generate_legend_entry(variation_vector, metric)
            axes.boxplot(plot_data)
            axes.set_xticks(range(1, len(x_values)+1))
            axes.set_xticklabels(list(x_values))

    def convergence_order_lines(self, axes):
        # Convergence lines can only be constructed for log-log graphs
        if (self.xscale != 'log') or (self.yscale != 'log') or (not self.order_indicators):
            return

        xmin, xmax = axes.get_xlim()
        ymin, ymax = axes.get_ylim()

        # Construct first order convergence
        # Start point as percentage of y-axis
        s1 = 0.8
        E11 = math.pow(2.0, math.log2(ymin) + s1*(math.log2(ymax) - math.log2(ymin)))
        E12 = E11*(xmin/xmax)

        axes.plot([xmin, xmax], [E11, E12], linestyle='--', color='black', zorder=0.0)

        # Construct first order convergence
        # Start point as percentage of y-axis
        s2 = 0.8
        E21 = math.pow(2.0, math.log2(ymin) + s2*(math.log2(ymax) - math.log2(ymin)))
        E22 = E21*math.pow((xmin/xmax), 2.0)

        axes.plot([xmin, xmax], [E21, E22], linestyle='dotted', color='black', zorder=0.0)



    #--------------------------------------------------------------------------
    # Interface member functions
    #--------------------------------------------------------------------------
    def plot(self, axes, new_x_label="", new_y_label=""):
        """
        Plot the data of the class to the given axes object.

        Plot all metrics specified in the metric_names list for each
        variation vector of the dataframe.

        Keyword arguments:
        new_x_label -- Use this string as x-axis label. If empty, the label
                       is automatically derived (default : "")
        new_y_label -- Use this string as y-axis label. If empty, the label
                       is automatically derived (default : "")
        """
        axes.set_xscale(self.xscale)
        axes.set_yscale(self.yscale)

        # Parameter vectors to iterate over for plotting
        iteration_parameters = list(self.dataframe.index.names)

        if self.plot_type == 'plot':
            try:
                iteration_parameters.remove(self.x_axis_name)
            except:
                iteration_parameters.remove('step')
        else:
            iteration_parameters.remove(self.x_axis_name)
            iteration_parameters.remove('step')

        if iteration_parameters:
            variation_vectors = dwm.variation_vectors(self.dataframe.index, iteration_parameters)
        else:
            variation_vectors = []

        # Add each graph to the plot
        for var_vector in variation_vectors:
            data_subset = self.dataframe.xs(list(var_vector.values()),
                                             level=list(var_vector.keys()),
                                             drop_level=True)
            if self.plot_type == 'plot':
                self.plot_metric_graphs(axes, data_subset, var_vector)
            elif self.plot_type == 'scatter':
                self.scatter_plot(axes, data_subset, var_vector)
            elif self.plot_type == 'box':
                self.box_plot(axes, data_subset, var_vector)
            else:
                print("Printing nothing, plot type",self.plot_type,
                        "not supported.")
                # Error: unsupported plot type

        # Special case when all variation parameters are fixed for the plot
        if not variation_vectors:
            if self.plot_type == 'plot':
                self.plot_metric_graphs(axes, self.dataframe)
            elif self.plot_type == 'scatter':
                self.scatter_plot(axes, self.dataframe)
            elif self.plot_type == 'box':
                self.box_plot(axes, self.dataframe)
            else:
                print("Printing nothing, plot type",self.plot_type,
                        "not supported.")
                # Error: unsupported plot type

        self.limit_x_axis_range(axes)
        self.limit_y_axis_range(axes)

        # Optional: add lines of first and second order convergence
        self.convergence_order_lines(axes)

        # Only show legend if there is more than one graph in the plot
        if not (not variation_vectors and len(self.metric_names) == 1):
            axes.legend(loc='best')

        self.set_axes_labelling(axes, new_x_label, new_y_label)

        if self.draw_title:
            self.set_plot_title(axes)

    def change_metric_set(self, new_metric_set):
        """
        Change the set of metrics to be plotted.

        This can be used to create mutiple plots of different metrics from
        the same dataset.
        """
        self.metric_names = new_metric_set

class ArrayPlot:
    def __init__(self, dataframe, parameter_metric_sets, x_axis_name, plot_layout=(6, 3),
                 plot_length=6, paper_plots=False, fixed_parameters={}, save_pickles=False, **kw_plot):
        """
        Create a set of plots from the given dataframe using matplotlib.

        Utilizes 'single_plot_creator'-class to create all plots for the given
        dataframe. Plots are grouped in files using matplotlib.pyplot.subplots.
        parameter_metric_sets defines what is drawn in a single plot.

        Keyword arguments:
        plot_layout -- tuple defining the number of plots per row and column in a
                       single file. (default : (6,3))
        plot_length -- scalar defining height and width of a single plot (default : 6)
        paper_plots -- flag indicating that plots are intended for a publication.
                       If true, one file per plot is generated and a different
                       matplotlib style is used. (default : False)
        fixed_parameters -- dictionary of fixed parameters. In the case when only a data subset
                            is plotted they are added to the info file for
                            completeness. (default : {})
        save_pickles -- save each plot also as a pickle file for further modifications. (default : False)
        **kw_plot -- Keyword arguments which are forwarded to 'single_plot_creator'.
        """
        self.dataframe = dataframe
        # parameter_metric_sets is a list of dictionaries where each dictionary
        # has an entry 'parameters' and an entry 'metrics'.
        # 'parameters' is associated with a list of variation parameters which shall
        # be added to a single plot while 'metrics' is associated with a list of metrics
        # which shall be added to this plot.
        # For each dictionary, plots are generated.
        self.parameter_metric_sets = parameter_metric_sets
        self.x_axis_name = x_axis_name
        self.plot_layout = plot_layout
        self.plot_length = plot_length
        self.paper_plots = paper_plots
        self.fixed_parameters = fixed_parameters
        self.save_pickles = save_pickles
        self.kw_plot = kw_plot

        # For paper use generate one plot per file
        if paper_plots:
            self.plot_layout = (1, 1)

    def configure_matplotlib_parameters(self):
        """
        Set the appearance and style of the generated plots (hardcoded for now).
        """
        if self.paper_plots:
            # Use TeX rendering only for paper plots. Memory limit of TeX may
            # become an issue for multiple plots in a single pdf
            mpl.backend_bases.register_backend('pdf', FigureCanvasPgf)

            plt.style.use('seaborn-paper')

            mpl.rcParams['text.usetex'] = True
            mpl.rcParams['font.size'] = 18
            mpl.rcParams['lines.linewidth'] = 1.0
            mpl.rcParams['axes.grid'] = True
            mpl.rcParams['grid.linewidth'] = 0.5
            mpl.rcParams['grid.linestyle'] = 'dotted'
        else:
            plt.style.use('Solarize_Light2')

            mpl.rcParams['lines.linewidth'] = 1.0
            mpl.rcParams['axes.titlesize'] = 12

        # Other, interesting/relevant parameters (https://matplotlib.org/tutorials/introductory/customizing.html)
        #lines.linestyle   : -       ## solid line.
        #lines.markersize  : 6            ## markersize, in points
        #markers.fillstyle: full ## full|left|right|bottom|top|none

    def compute_plot_layout(self, nplots):
        """
        Compute the plot layout and return a triple with the number of
        rows, columns and files.
        """
        nrows, ncols = self.plot_layout

        # Distinguish two cases: a user prescribed number of rows per file
        # and the case of a single file where the number of rows is determined
        # by the number of plots
        if nrows == 0:
            nrows = int(math.ceil(nplots/ncols))
            nfiles = 1
        else:
            plots_per_file = nrows*ncols 
            nfiles = int(math.ceil(nplots/plots_per_file))

        return (nrows, ncols, nfiles)

    def select_axes(self, axesarray, row_number, column_number):
        """
        Select an axes from drom a given axesarray.

        This function allows uniform access to axesarray which may be
        an axes object or a 1- or 2-d numpy array containing axes
        objects (see matplotlib.pyplot.subplots(...)).
        """
        if isinstance(axesarray, numpy.ndarray):
            if axesarray.ndim == 2:
                return axesarray[row_number, column_number]
            else:
                return axesarray[max(row_number, column_number)]
        else:
            return axesarray

    def assemble_plot_name(self, base_name, parameter_metric_pair, nstart, nend):
        """
        Generate a meaningful name for a plot file containing the range of plot numbers.
        """
        parameters = parameter_metric_pair['parameters']
        metrics = parameter_metric_pair['metrics']

        plot_name = base_name

        for par in parameters:
            plot_name = plot_name + '-' + str(par)

        for met in metrics:
            plot_name = plot_name + '-' + str(met)

        plot_name = plot_name + '-' + str(nstart).zfill(5) + '-' + str(nend).zfill(5) + '.pdf'

        # Ensure there are no underscores in the file name as this causes trouble with
        # LaTeX
        plot_name = plot_name.replace('_','')
        plot_name = plot_name.replace(' ','')

        return plot_name

    def write_plot_parameter_vectors(self, base_name, parameter_metric_pair, variation_vectors):
        """
        Write a text file which maps the ID/number of a plot to its parameter vector.
        """
        file_name = self.assemble_plot_name(base_name, parameter_metric_pair, 0, 0)
        file_name = file_name.split('0')[0] + 'info.txt'

        # Header
        info_file = open(file_name, 'w+')
        infos = 'Remember the plot counting order: from left to right, then from top to bottom.\n'
        infos += 'Example:\n'
        infos += '0   1   2\n'
        infos += '3   4   5\n'
        infos +='--------------------------------------------------------------------------\n'

        # Document fixed parameters for completeness
        infos += "Parameters fixed for all plots:\n"
        for key,value in self.fixed_parameters.items():
            infos += '\t' + str(key).ljust(25) + ' : ' + str(value) + '\n'
        infos +='--------------------------------------------------------------------------\n'

        # Write each plot configuration
        idx = 0
        for var_vec in variation_vectors:
            infos += 'Plot number: ' + str(idx).zfill(5) + '\n'
            for entry in var_vec.items():
                infos += '    ' + str(entry[0]).ljust(25) + ' : ' + str(entry[1]) + '\n'

            infos += '\n'
            idx += 1

        info_file.write(infos)
        info_file.close()


    #--------------------------------------------------------------------------
    # Interface member functions
    #--------------------------------------------------------------------------
    def set_x_axis_name(self, new_x_axis_name):
        """
        Set a new name for the data to be used as x-axis.
        """
        self.x_axis_name = new_x_axis_name

    def create_plots(self, file_name_prefix):
        """
        Actually create the plots.
        """
        self.configure_matplotlib_parameters()

        # Iterate the parameter-metric pairs
        for parameter_metric_pair in self.parameter_metric_sets:
            parameters = parameter_metric_pair['parameters']
            metrics = parameter_metric_pair['metrics']
            free_parameters = list(self.dataframe.index.names)

            # x-axis cannot be part of free parameters. In case x-data is taken from a column,
            # remove the 'step' level
            try:
                free_parameters.remove(self.x_axis_name)
            except:
                free_parameters.remove('step')
            for fixed_parameter in parameters:
                free_parameters.remove(fixed_parameter)

            # Catch case when the dataframe is reduced to a single plot
            if free_parameters and 'step' not in free_parameters:
                variation_vectors = dwm.variation_vectors(self.dataframe.index, free_parameters)
                nplots = len(variation_vectors)
                nrows, ncols, nfigures = self.compute_plot_layout(nplots)

                self.write_plot_parameter_vectors(file_name_prefix, parameter_metric_pair, variation_vectors)
            else:
                nplots = 1
                nrows = 1
                ncols = 1
                nfigures = 1
                variation_vectors = []

            print("\nTotal number of plots:", nplots,"\n")

            # Add the plots in chunks of nrows*ncols to a separate figure (or file, respectively)
            chunk_end = 0
            for chunk_start in range(0, nplots, nrows*ncols):
                chunk_end += nrows*ncols
                plt.clf()
                figure_width = ncols*self.plot_length
                figure_height = nrows*self.plot_length
                figure, axesarray = plt.subplots(nrows=nrows, ncols=ncols, figsize=(figure_width, figure_height))

                print("Plotting from ", chunk_start, " till ", chunk_end)

                # Iterate the plots for a file / figure
                all_plots_created = False
                for row in range(nrows):
                    for col in range(ncols):
                        var_vector_index = row*ncols + col + chunk_start
                        
                        # Stop plotting when the end of the variation vectors is reached
                        if var_vector_index == nplots:
                            all_plots_created = True
                            break

                        if variation_vectors:
                            var_vector = variation_vectors[var_vector_index]
                            plot_df = self.dataframe.xs(list(var_vector.values()), level=list(var_vector.keys()))
                        else:
                            var_vector = {}
                            plot_df = self.dataframe
                        plotter = single_plot_creator(plot_df, metrics, self.x_axis_name, var_vector, draw_title= not self.paper_plots, **self.kw_plot)
                    
                        # plt.subplots returns either an axes-object, a 1-d array or a 2-d array, depending
                        # on the plot layout. Thus, axesarray must be accessed accordingly.
                        current_axes = self.select_axes(axesarray, row, col)
                        plotter.plot(current_axes)

                    if all_plots_created:
                        break

                # tight_layout() ensures that the single plots do not overlap in the figure
                plot_name = self.assemble_plot_name(file_name_prefix, parameter_metric_pair, chunk_start, chunk_end-1)
                plt.tight_layout()
                plt.savefig(plot_name, bbox_inches='tight')

                if self.save_pickles:
                    target_file = open(plot_name + '.pickle','wb')
                    pl.dump(figure, target_file)

                # No need to keep finished and written figures in memory
                plt.close(figure)

