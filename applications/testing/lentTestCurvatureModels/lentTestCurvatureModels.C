/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) held by orignal authors.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    lentTestCurvatureModels 

Description

    Test LENT curvature models against the exact curvature. 

Authors
    Tomislav Maric maric@csi.tu-darmstadt.de

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "timeSelector.H"

#include "frontCurvatureModel.H"

// Testing OF curvature
#include "interfaceProperties.H"
#include "incompressibleTwoPhaseMixture.H" 

#include <fstream>

using namespace FrontTracking;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "errorFile",
        "Path to the error file to which the application appends the results for later analysis." 
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

#include "createFields.H"

    if (!args.optionFound("errorFile"))
    {
        FatalErrorIn("main()")   
            << "Please use option '-errorFile' to name the output file."
            << endl << exit(FatalError);
    }

    const std::string errorFileName = args.optionRead<fileName>("errorFile"); 

    const char* errorFileNamePtr = errorFileName.c_str(); 

    std::fstream errorFile; 

    errorFile.open(errorFileNamePtr, std::ios_base::app);

    //if (errorFile.peek() == std::ifstream::traits_type::eof())
    //{
        //errorFile << "#h        Linf        LinfInterface" << std::endl; 
    //}

    IOdictionary lentMethodDict
    (
        IOobject(
            "lentSolution",
            "system", 
            runTime,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );

    const dictionary& curvatureDict = lentMethodDict.subDict("frontCurvatureModel"); 
    const word inputFieldName = curvatureDict.lookup("inputField"); 

    volScalarField inputField
    (
        IOobject(
            inputFieldName, 
            runTime.timeName(), 
            mesh, 
            IOobject::MUST_READ, 
            IOobject::NO_WRITE 
        ),
        mesh
    );

    tmp<frontCurvatureModel> curvatureModel = frontCurvatureModel::New(curvatureDict, runTime); 
    
    volScalarField delta = mag(fvc::grad(inputField)); 

    volScalarField cellCurvatureExact
    (
        IOobject
        (
            "cellCurvatureExact",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh 
    );

    cellCurvatureExact *= delta;
    cellCurvatureExact.write();

    surfaceScalarField faceCurvatureExact
    (
        IOobject
        (
            "faceCurvatureExact",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh 
    );

    // Compute the numerical curvature with the model. 
    volScalarField cellCurvature("cellCurvature", curvatureModel->cellCurvature()); 
    cellCurvature *= delta; 
    cellCurvature.write(); 

    volScalarField LinfField (
            "LinfCurvatureErr", 
            mag(cellCurvature - cellCurvatureExact) 
    ); 

    LinfField.write();  

    // Compute and test the numerical curvature using OpenFOAM interfaceProperties 
    volScalarField cellCurvatureInterface("cellCurvatureInterface", interface.K()); 
    cellCurvatureInterface *= delta;
    cellCurvatureInterface.write(); 
   
    volScalarField LinfInterfaceField (
            "LinfInterfaceCurvatureErr", 
            mag(cellCurvatureInterface - cellCurvatureExact) 
    ); 

    LinfInterfaceField.write();  

    dimensionedScalar Linf = max(LinfField);

    dimensionedScalar LinfInterface = max(LinfInterfaceField);

    dimensionedScalar h = max(mag(mesh.delta())); 

    errorFile << h.value() << " " << Linf.value() << " " << LinfInterface.value() << std::endl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
