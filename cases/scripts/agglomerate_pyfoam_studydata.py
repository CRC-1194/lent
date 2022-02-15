#!/usr/bin/env python

from argparse import ArgumentParser

import dataAgglomeration as da

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

if __name__ == "__main__":
    main()
