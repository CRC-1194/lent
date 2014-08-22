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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    // Add necessary options and check if they have been set by the user
    argList::addOption
    (
        "curvatureModel",
        "LENT curvature model."
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    if (!args.optionFound("curvatureModel"))
    {
        FatalErrorIn("main()")   
            << "Please use option '-curvatureModel' to select the LENT curvature model."
            << endl << exit(FatalError);
    }

    const word modelName = args.optionRead<word>("curvatureModel");

    // Select the curvature model based on the name and construct it from lentSolution
    // dictionary. 

    // Get the time directories from the simulation folder using time selector
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);

    forAll(timeDirs, timeI)
    {
        // Set the time to the current time directory and update time index.
        runTime.setTime(timeDirs[timeI], timeI);

        #include "createCurvatureFields.H"

        // Compute the numerical curvature with the model. 


        // Compute the curvature error fields: see literature. 

        // Compute the curvature diagram data. 
        
        // Store the curvature error fields for visualization. 

    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
