/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Version:  2.2.x                               
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    Test application that starts in a cell at seedPoint (found by fvMesh),
    and tries to find the cell containing the targetPoint (found by the
    frontMeshSearch algorithm).

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Tomislav Maric maric@csi.tu-darmstadt.de (main developer)
    under the project supervision of :
        Holger Marschall <marschall@csi.tu-darmstadt.de> (group leader).
    
    Method Development and Intellectual Property :
    	Tomislav Maric maric@csi.tu-darmstadt.de
    	Holger Marschall <marschall@csi.tu-darmstadt.de>
    	Dieter Bothe <bothe@csi.tu-darmstadt.de>

        Mathematical Modeling and Analysis
        Center of Smart Interfaces
        Technische Universitaet Darmstadt
       
    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
#include "frontMeshSearch.H"

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

    // FIXME: Introduce flag: check frontMeshSearch code for debugging parts. TM
    //frontMeshSearch ls(runTime);

    frontMeshSearch ls;

    label foundCell = ls.cellContainingPoint(targetPoint, mesh, seedCell);

    label meshFoundCell = mesh.findCell(targetPoint);

    Info << "targetPoint = " << targetPoint << endl;
    Info << "foundCell = " << foundCell << endl;
    Info << "meshFoundCell = " << meshFoundCell << endl;

    if (meshFoundCell == foundCell)
    {
        Info << "PASS" << endl;
    }
    else
    {
        Info << "FAIL" << endl;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nEnd\n" << endl;
    return 0;
}

// ************************************************************************* //
