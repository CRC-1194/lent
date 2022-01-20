#!/usr/bin/env python

from argparse import ArgumentParser

import dataAgglomeration as da

from nbclient import execute

import nbformat

from nbparameterise import extract_parameters, parameter_values, replace_definitions

def main():

    #---- Command line arguments ----------------------------------------------
    parser = ArgumentParser(description="Agglomerate the data of a parameter study from its corresponding directories into a single, multiindexed pandas dataframe saved as both a CSV and a JSON file.")

    parser.add_argument("path_to_file")
    parser.add_argument("-p", "--parameter-file",
                        help="Parameter file with which the study has been generated.",
                        required=True,
                        dest="parameter_file_name"
                        )
    parser.add_argument("-f", "--file-name",
                        help="Write agglomerated data to this file. If not given, use name of parameter file.",
                        required=False,
                        default="",
                        dest="target_file_name"
                        )

    args = parser.parse_args()

    example_directory, file_path_and_name = args.path_to_file.split('/', maxsplit=1)

    agglomerator = da.data_agglomerator(args.parameter_file_name, file_path_and_name, example_directory)
    agglomerator.show_failed_variations()
    agglomerator.write_agglomerated_study_data(args.target_file_name)


    #------ parameterise and write .ipynb file for every case ------#
    nb_name, parameter_file_extension = example_directory.split('.')
    nb = nbformat.read("/work/projects/project01456/lent/cases/flow-solution/translating-droplet/three-dimensional/densityRatioTemplate.ipynb", as_version=4)

    orig_parameters = extract_parameters(nb)
    params = parameter_values(orig_parameters, case_path = example_directory, study_name=args.parameter_file_name)
    new_nb = replace_definitions(nb, params)

    with open("%s.ipynb" % nb_name, 'w') as f:
        nbformat.write(new_nb, f)

if __name__ == "__main__":
    main()
