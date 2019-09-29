# oscillatingDropletPostprocessing.py 

import math
import pandas as pd
import sys

import dataframeWithMetadata as dwm
import dataAgglomeration as da

class oscillating_droplet_postprocessor:
    """
    Compute analytical oscillation period and decay rate and the difference to numerical data.

    Compute the oscillation period and the decay rate of the amplitude according to
    Lamb's analytical solution. Both metrics are also computed from the provided
    numerical data. Finally, the relative errors are computed for both metrics.
    """

    def __init__(self, dataframe_file_name):
        """
        Constructor of the oscillating_droplet_postprocessor requiring the file
        name of a dataframe.
        """
        self.reader = dwm.metadated_dataframe_reader(dataframe_file_name)
        self.source_dataframe = self.reader.read_dataframe()
        self.solution_parameters = {
                    'mode_number' : 'undef',
                    'surface_tension_coefficient' : 'undef',
                    'rho_ambient' : 'undef',
                    'rho_droplet' : 'undef',
                    'radius' : 'undef',
                    'nu_droplet' : 'undef'
                }
        self.computed_dataframe = pd.DataFrame()
        self.variations_with_incomplete_data = []

    def set_solution_parameters(self):
        """
        Read required physical parameters from dataframe metadata. Throws an error
        if a required parameter is missing.
        """
        metadata = self.reader.read_metadata()

        for key in self.solution_parameters.keys():
            if key not in metadata:
                print("Error: the parameters", self.solution_parameters.keys(),
                        "must be present as metadata in the given dataframe.")
                sys.exit()

            value = metadata[key]

            if value == "from_index":
                self.solution_parameters[key] = value
            else:
                try:
                    self.solution_parameters[key] = float(value)
                except ValueError:
                    print("Error: found invalid value '", value,"' for parameter",
                            key, ".\n Valid values are floats or 'from_index'.")
                    sys.exit()

    def prepare_computed_dataframe(self):
        """
        Prepare a multiindexed dataframe to store the postprocessing results.
        """
        index = self.source_dataframe.index
        index = index.droplevel(level='step').unique()
        column_names = ['period_analytic', 'period_numerical',
                        'relative_period_error', 'decay_factor_analytic',
                        'decay_factor_numerical', 'relative_decay_factor_error',
                        'amplitude_analytic', 'amplitude_numerical',
                        'relative_amplitude_error'
                       ]
        self.computed_dataframe = pd.DataFrame(index=index, columns=column_names)

    def required_solution_parameters(self, parameter_set):
        """
        Return a dictionary with values of physical parameters for the given
        parameter set.

        Completes the set of physical parameters required for postprocessing
        by reading varying parameters from the parameter set, e.g. the viscosity
        if it is varied in a study.
        """
        solution_parameters = dict(self.solution_parameters)
        for key, value in solution_parameters.items():
            if value == "from_index":
                solution_parameters[key] = float(parameter_set[key])

        return solution_parameters

    def find_first_maximum(self, droplet_data, column_name, analytic_mean):
        """
        Find second maximum in semi-axis data and return a tuple of indices.

        Assuming a usual oscillating droplet setup the semi-axis of interest
        has its global maximum at start time. This function is intended to find
        the second maximum, meaning the semi-axis value after one oscillation
        period.
        Positional arguments:
        1 -- dataframe with oscillating droplet data
        2 -- name of the column containing semi-axis data
        3 -- semi-axis mean value, usually the unperturbed droplet radius
        """
        indices_with_max_value = [0]
        epsilon = 1.0e-8
        data_column = droplet_data[column_name]
        passed_mean_once = False
        passed_mean_twice = False

        max_value = analytic_mean
        max_index = 0

        for index, value in data_column.iteritems():
            # The searched for maximum should be located after the semi axis
            # has assumed the unperturbed radius twice (have a look on the
            # temporal evolution of the semi axis n paper).
            if not passed_mean_once and value < analytic_mean:
                passed_mean_once = True

            if passed_mean_once and value > analytic_mean:
                passed_mean_twice = True

            if not passed_mean_twice:
                continue
            
            # Actual search for first local maximum
            delta = value - max_value
            if delta > epsilon:
                max_value = value
                indices_with_max_value[0] = index
            elif delta < -epsilon:
                break
            else:
                indices_with_max_value.append(index)

        return tuple(indices_with_max_value)

    def compute_droplet_quantities(self, parameter_set, droplet_data):
        """
        Compute oscillation period and decay factor analytically and numerically.

        For the given droplet data, compute the oscillation period
        and decay factor according to Lamb's analytical solution. Both quantities
        are also computed for the given numerical data and then the relative errors
        are computed. All quantities are stored in a dataframe.
        Also registers if the evaluation is not possible for the numerical data.
        """
        solution_parameters = self.required_solution_parameters(parameter_set)

        # Compute analytic quantities
        n = solution_parameters['mode_number']
        sigma = solution_parameters['surface_tension_coefficient']
        nu_d = solution_parameters['nu_droplet']
        rho_a = solution_parameters['rho_ambient']
        rho_d = solution_parameters['rho_droplet']
        R = solution_parameters['radius']
        amplitude_0 = droplet_data.loc[0, 'semi-axes-x'] - R

        omega_a_squared = n*(n+1)*(n-1)*(n+2)*sigma/(((n+1)*rho_d + n*rho_a)*R**3)
        period_a = 2*math.pi / math.sqrt(omega_a_squared)
        decay_a = R*R/((n-1)*(2*n+1)*nu_d)
        # The analytic amplitude must be computed at the same time as the
        # numeric one. Thus, the numeric period is required
        amplitude_a = -1.0

        # Compute numeric quantities
        period_n = -1.0
        decay_n = -1.0
        amplitude_n = -1.0
        period_error = -1.0
        period_error_abs = -1.0
        decay_error = -1.0
        decay_error_abs = -1.0
        amplitude_error = -1.0
        amplitude_error_abs = -1.0
        first_maximum_indices = self.find_first_maximum(droplet_data, "semi-axes-x", R) 
        
        # Only evaluate numeric metrics if the returned indices make sense.
        if len(first_maximum_indices) > 1 or first_maximum_indices[0] != 0:
            time_data = droplet_data['time']
            axis_data = droplet_data['semi-axes-x']

            period_n = 0.0
            amplitude_n = 0.0

            for index in first_maximum_indices:
                period_n += time_data[index]
                amplitude_n += axis_data[index]

            period_n /= len(first_maximum_indices)
            amplitude_n = amplitude_n/len(first_maximum_indices) - R
            decay_n = -1.0*period_n/math.log(amplitude_n/amplitude_0)

            period_error = (period_n - period_a)/period_a
            period_error_abs = math.fabs(period_error)
            decay_error = (decay_n - decay_a)/decay_a
            decay_error_abs = math.fabs(decay_error)
            amplitude_a = amplitude_0*math.exp(-period_n/decay_a)
            amplitude_error = (amplitude_n - amplitude_a)/amplitude_a
            amplitude_error_abs = math.fabs(amplitude_error)
        else:
            self.variations_with_incomplete_data.append(parameter_set)
        
        index_tuple = tuple(parameter_set.values())
        self.computed_dataframe.loc[index_tuple, 'period_analytic'] = period_a
        self.computed_dataframe.loc[index_tuple, 'period_numerical'] = period_n
        self.computed_dataframe.loc[index_tuple, 'relative_period_error'] = period_error
        self.computed_dataframe.loc[index_tuple, 'relative_period_error_abs'] = period_error_abs
        self.computed_dataframe.loc[index_tuple, 'decay_factor_analytic'] = decay_a
        self.computed_dataframe.loc[index_tuple, 'decay_factor_numerical'] = decay_n 
        self.computed_dataframe.loc[index_tuple, 'relative_decay_factor_error'] = decay_error
        self.computed_dataframe.loc[index_tuple, 'relative_decay_factor_error_abs'] = decay_error_abs
        self.computed_dataframe.loc[index_tuple, 'amplitude_analytic'] = amplitude_a
        self.computed_dataframe.loc[index_tuple, 'amplitude_numerical'] = amplitude_n
        self.computed_dataframe.loc[index_tuple, 'relative_amplitude_error'] = amplitude_error
        self.computed_dataframe.loc[index_tuple, 'relative_amplitude_error_abs'] = amplitude_error_abs


    #--------------------------------------------------------------------------
    # Interface member functions
    #--------------------------------------------------------------------------
    def postprocess_droplet_data(self):
        """
        Apply postprocessing to the oscillating droplet data and store results
        in a dataframe.
        """
        # Postprocessing already done
        if not self.computed_dataframe.empty:
            return

        self.set_solution_parameters()
        self.prepare_computed_dataframe()

        # Iterate all parameter configurations
        parameter_sets = dwm.variation_vectors(self.source_dataframe.index)

        for parameter_set in parameter_sets:
            data_set = self.source_dataframe.xs(list(parameter_set.values()),
                                                level=list(parameter_set.keys()),
                                                drop_level=True
                                               )
            self.compute_droplet_quantities(parameter_set, data_set)

        # Remove rows with incomplete data
        self.computed_dataframe = self.computed_dataframe[self.computed_dataframe['period_numerical'] >= 0.0]

    def show_variations_with_incomplete_data(self):
        """
        Show all variations whose data could not be postprocessed.
        """
        print("Variations whose data could not be postprocessed:\n")
        for fail in self.variations_with_incomplete_data:
            print(fail)

    def write_postprocessing_data(self, file_name):
        """
        Write data of postprocessing as CSV file.
        """
        self.postprocess_droplet_data()
        dwm.write_dataframe_with_metadata(file_name, self.computed_dataframe, self.solution_parameters)
