/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    testLentMeshSearch

Description
    Test application that starts in a cell at seedPoint (found by fvMesh), 
    and tries to find the cell containing the targetPoint (found by the 
    lentMeshSearch algorithm).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "lentMeshSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

using namespace FrontTracking;

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "seedPoint",
        "point used to find the seedCell needed for the known vicinity search algorithm" 
    );

    argList::addOption
    (
        "targetPoint",
        "point which is searched for in the mesh" 
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    if (!args.optionFound("seedPoint") || !args.optionFound("targetPoint"))
    {
        FatalErrorIn
        (
            "main()"
        )   << "Provide both the 'seedPoint' and the 'targetPoint' arguments." << endl
        << exit(FatalError);
    }

    point seedPoint = args.optionRead<point>("seedPoint"); 
    point targetPoint = args.optionRead<point>("targetPoint"); 

    label seedCell = mesh.findCell(seedPoint); 

    lentMeshSearch ls(runTime); 

    label foundCell = ls.cellContainingPoint(targetPoint, mesh, seedCell); 

    Info << "targetPoint = " << targetPoint << endl;
    Info << "foundCell = " << foundCell << endl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
