#!/usr/bin/env python3

from argparse import ArgumentParser

import testReportCore as trc
import parameterStudyPreparation as psp


def main():

    #---- Command line arguments ----------------------------------------------
    parser = ArgumentParser(description="Prepare the variations for the given parameter study")

    parser.add_argument("studyname",
                        help="Name of the parameter file for the study")
    parser.add_argument("-p","--prefix",
                        help="Prefix added to the name of the study folders. The result is prefix-studyname_00... etc.",
                        default="",
                        dest="studyprefix")
    parser.add_argument("-ap","--additional-parameters",
                        help="List of additional parameters files to be used.",
                        default="",
                        dest="additionalParameterFiles")
    parser.add_argument("-m","--mesh-type",
                        help="Mesh type to be used.",
                        required=True,
                        choices=["block", "cartesian", "poly", "hexrefined"],
                        dest="meshtype")
    parser.add_argument("-md","--mesh-dict-name",
                        help="Prescribe the name of the mesh dictionary file to be used. Use this if different studies require different meshing approaches. NOTE: it is assumed that the mesh dict is always templated so that resolution can be a study parameter.",
                        default=None,
                        dest="meshdictname")
    parser.add_argument("-v","--variants",
                        help="Only use the specified variations. By default, all variations are used. Argument can either be a single number (e.g. 42), a list of numbers (e.g. '3,5,11') or a range ('3 - 10')",
                        default="all",
                        dest="variants")
    parser.add_argument("-np","--no-preprocessing",
                        help="Do not execute any further preprocessing steps: no invocation of a mesher and no execution of caseSetup.sh",
                        action="store_true",
                        default=False,
                        dest="no_preprocess")

    args = parser.parse_args()


    #----- check input parameters for validity --------------------

    # Prepare the parameter file
    studyname = psp.create_parameter_file_from_string(args.studyname + "," +
                                            args.additionalParameterFiles)

    caseName = "template_copy_" + args.meshtype
    templateName = trc.pattern_cases("templateCase_*")
    templateName = templateName[0]

    ## Setup a clean copy of the template case
    psp.copy_case_template(caseName, templateName, args.meshtype, args.meshdictname)

    # Create vector of all variants to be set up
    variationNumbers = psp.create_variant_vector(args.variants)

    studyprefix = args.studyprefix
    if studyprefix != "":
        studyprefix = studyprefix + "-"

    # Assemble the parameter variation command according to the given options
    command =   ["pyFoamRunParameterVariation.py",
                 "--allow-derived-changes",
                 "--every-variant-one-case-execution",
                 "--create-database", 
                 "--no-server-process",
                 "--no-execute-solver",
                 "--parameter-file=default.parameter",
                 "--cloned-case-prefix="+studyprefix+args.studyname,
                 caseName,
                 studyname
                 ]

    # Preprocessing options
    if args.no_preprocess:
        command.append("--no-mesh-create")
        command.append("--no-case-setup")
    else:
        command.append("--mesh-create-script="+args.meshtype+"Mesh.sh")

    # Finally, create the variations
    psp.setup_variants(command, variationNumbers)
    

if __name__ == "__main__":
    main()
