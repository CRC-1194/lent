# parameterStudyPreparation.py

import sys
import shutil
import subprocess

def clear_toplevel_values_dict(parameterFileContent):
    """ Remove the enclosing 'values{...}' dictionary from a parameter file.
    """

    if "values" in parameterFileContent:
        parameterFileContent = parameterFileContent.replace("values", "", 1)    
        parameterFileContent = parameterFileContent.replace("{", "", 1)

        parenthesisIndex = parameterFileContent.rindex("}")
        parameterFileContent = parameterFileContent[:parenthesisIndex]

    return parameterFileContent


def create_parameter_file_from_string(fileListString):
    """ Assemble a parameter file for a study. The file list is created
        from the given string. ',' has to be used as a separator.
    """
    file_names = fileListString.split(',')

    # Remove empty list entries
    file_names = [x for x in file_names if x != ""]

    return create_parameter_file_from_list(file_names)


def create_parameter_file_from_list(file_names):
    """ Assemble a parameter file for a study. Does not check for duplicate
        entries or such things.
        If 'file_list' only contains a single file name, no new file is created.
    """

    # Nothing to be done here if only a single file is given
    if len(file_names) == 1:
        return file_names[0]

    # Fuse several parameter files into a single one
    combined_parameters = "values{\n"
    parameter_study_name = ""

    for parameter_file_name in file_names:
        parameter_file = open(parameter_file_name, 'r')
        combined_parameters += clear_toplevel_values_dict(parameter_file.read())
        parameter_study_name += "_" + parameter_file_name.split('.')[0]

    combined_parameters += "}"
    parameter_study_name += ".parameter"

    combined_study_file = open(parameter_study_name, 'w')
    combined_study_file.write(combined_parameters)

    return parameter_study_name


def copy_case_template(caseName, templateName, meshType):
    """ Prepare a clean copy of a test case template folder. Set the correct
        meshDict if cfMesh is used (for either pMesh or cartesianMesh).
    """

    # Setup a clean copy of the template case
    shutil.rmtree(caseName, ignore_errors=True)
    shutil.copytree(templateName, caseName)

    # Set up the initial values
    shutil.copytree(caseName+"/0.org", caseName+"/0")

    # If cfMesh is used, setup the meshDict of the chosen meshing approach.
    # This is required since both cartesianMesh and pMesh use system/meshDict
    # as input
    if meshType in ["cartesian", "poly"]:
        shutil.copyfile(caseName+"/system/meshDict."+meshType+".template",
                        caseName+"/system/meshDict.template")
    # Use the base resolution 'n_base' if a hex refined mesh is used
    elif meshType == "hexrefined":
        command = ["sed", "-i", "s/resolution/n_base/g",
                    caseName+"/system/blockMeshDict.template"]
        subprocess.run(command)


def create_variant_vector(variantArgument):
    """ Parse the string of the variants to be set up and return it as a list.
        The given string may contain integer numbers separated by ',' or by
        '-' to indicate a range which includes the given start and end value.
    """

    variants = list()

    # No need for additional parsing if all variants shall be set up
    if variantArgument == "all":
        variants.append(variantArgument)
        return variants

    splittedArgument = variantArgument.split(',')

    for element in splittedArgument:
        if '-' in element:
            start = int(float(element.split("-")[0]))
            end = int(float(element.split("-")[1]))
            for number in range(start, end+1):
                variants.append(str(number))
        else:
            if element.isdigit():
                variants.append(str(element))
            else:
                print("Error: ",element," is not a number."
                        " Only ',' and '-' are allowed as non-digit characters.")

    return variants


def setup_variants(command, variants):
    """ Execute the command (should be pyFoamRunParameterVariation.py) for
        all variants.
    """

    if variants[0] == "all":
        status = subprocess.run(command)
        print(status.stdout)
    else:
        command.append("--single-variation=")
        for variant in variants:
            command[-1] = "--single-variation="+variant
            status = subprocess.run(command)
            print(status.stdout)
