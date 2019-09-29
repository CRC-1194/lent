#!/usr/bin/env python3

from argparse import ArgumentParser

import oscillatingDropletPostprocessing as odp

def main():
    parser = ArgumentParser(description="Compute period, amplitude decay und their errors from oscillating drolet data.")
    parser.add_argument("file_name",
                        help="File name of the data frame.")
    parser.add_argument("-o", "--ouput-name",
                        help="File name for the output data frame.",
                        default="",
                        dest="output_name")
    
    args = parser.parse_args()

    output_name = ""

    if not args.output_name:
        output_name = args.file_name.split('.')[0] + "_period_amplitude_data"
    else:
        output_name = args.output_name

    postprocessor = odp.oscillating_droplet_postprocessor(args.file_name)
    postprocessor.postprocess_droplet_data()
    postprocessor.show_variations_with_incomplete_data()
    postprocessor.write_postprocessing_data(output_name)

if __name__ == "__main__":
    main()
