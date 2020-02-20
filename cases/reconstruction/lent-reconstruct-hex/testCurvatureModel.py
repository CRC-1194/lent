#!/usr/bin/python

from PyFoam.RunDictionary.ParsedBlockMeshDict import ParsedBlockMeshDict 
from subprocess import call 

blockMeshDict = ParsedBlockMeshDict("constant/polyMesh/blockMeshDict") 

meshDensities = ['(32 32 32)', '(50 50 50)', '(64 64 64)', '(100 100 100)', '(128 128 128)'] 

for density in meshDensities: 
    blockMeshDict["blocks"][2] = density
    blockMeshDict.writeFile(); 
    call(["blockMesh"])
    call("rm -rf 0".split(" "))
    call("cp -rf 0.org 0".split(" "))
    call(["lentSetFields"])
    call(["testDistanceGradient"])
    call(["lentTestCurvatureModel"])


blockMeshDict.purgeFile()
blockMeshDict.writeFile()
