#!/usr/bin/env python3

from argparse import ArgumentParser

import dataAgglomeration as da

def main():

    #---- Command line arguments ----------------------------------------------
    parser = ArgumentParser(description="Agglomerate the data of a parameter study from its corresponding directories into a single, multiindexed pandas dataframe saved as a *.csv file.")

    parser.add_argument("parameter_file_name")
    parser.add_argument("-p", "--file-pattern",
                        help="Pattern/name for the result files to be agglomerated (Regex).",
                        default=r"Results\.csv$",
                        dest="pattern"
                        )

    args = parser.parse_args()

    agglomerator = da.data_agglomerator(args.parameter_file_name, args.pattern)
    agglomerator.show_failed_variations()
    agglomerator.write_agglomerated_study_data()

if __name__ == "__main__":
    main()
