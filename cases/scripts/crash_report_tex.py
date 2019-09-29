#!/usr/bin/env python3

import os
import testReportCore as trc

def main():

    outputPath = trc.result_path("TESTRESULTS", "table")

    # --- loop over all studies ---
    parameterFiles = trc.list_of_parameter_files()

    for study in parameterFiles:
        report = trc.assemble_crash_report(study)

        # TODO: add nice formatters?
        outputFileName = os.path.join(outputPath, study + "-crash_report" + ".tex")
        report.to_latex(outputFileName, index=False)


if __name__ == "__main__":
    main()
