#!/usr/bin/env python3

import testReportCore as trc

def main():

    # --- loop over all studies ---
    parameterFiles = trc.list_of_parameter_files()

    for study in parameterFiles:
        report = trc.assemble_crash_report(study)
        report.to_csv("minimal_crash_report-" + study + ".csv", index=False)


if __name__ == "__main__":
    main()
