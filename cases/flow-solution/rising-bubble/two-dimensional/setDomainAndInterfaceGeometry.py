#!/usr/bin/env python3

# --- Description: ---
# Enable the centralized definition of geometrical parameters referring to the
# interface and domain geometry to avoid inconsistencies, e.g. between the 
# *.geo file defining the front and the geometry information passed to
# postprocessing function objects.
#
# --- Usage: ---
# ./setDomainAndInterfaceGeometry.py arg1 arg2
# where arg1 is the path to the test case and arg2 is the name of the
# parameter file of a parameter study.
# For now, the script is intended to be run from the folder containing
# the test case and the parameter file
#
# Parameters to be set by the script are defined in the parameter file
# using the following syntax:
#    //G key value
# where neither key or value are allowed to contain whitespaces
#
# The script will append the so-defined variables to the
# "arg1/constant/caseSetup.template" file and create the file
# "arg1/front/front.geo" from the template file "arg1/front/front.geoTemplate"
# using the pyratemp template engine

import pyratemp
from shutil import copyfile
import sys

def isGeometryEntry(line):
    if line.find("//G") == 0:
        return True
    else:
        return False

def readKey(line):
    entries = line.split()
    return entries[1]

def readValue(line):
    entries = line.split()
    return entries[2]

def readGeometryInformation(filename):
    with open(filename) as inputFile:
        geometryDict = {}

        for line in inputFile:
            if isGeometryEntry(line):
                geometryDict[readKey(line)] = readValue(line)

        return geometryDict

def correctCase(caseName):
    if caseName.endswith("/"):
        return caseName
    else:
        return caseName + "/"


#------------------------------------------------------------------------------


case = correctCase(sys.argv[1])
parameterSet = sys.argv[2]

geoInfos = readGeometryInformation(parameterSet)

# Append information to caseSetup file so it is available for OpenFOAM
# dictionaries. Copy it beforehand to avoid spoiling the original file
copyfile(case + "constant/caseSetup.preTemplate", case + "constant/caseSetup.template")
caseSetupFile = open(case + "constant/caseSetup.template", "a")

for key, value in geoInfos.items():
    caseSetupFile.write(key + " " + value + ";\n")

# Set geometrical parameters in the *.geo file using pyratemp template
# expansion
gmshTemplate = pyratemp.Template(filename=case + "front/front.geoTemplate", data=geoInfos)

gmshFile = open(case + "front/front.geo", "w")
gmshFile.write(gmshTemplate())
