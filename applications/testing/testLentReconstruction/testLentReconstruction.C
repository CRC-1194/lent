/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    levelSetFrontFoam

Authors
    Tomislav Maric maric@csi.tu-darmstadt.de, tomislav@sourceflux.de 

Description
    Test application for the interface reconstruction algorithm of the LENT method.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "interfaceProperties.H"
#include "incompressibleTwoPhaseMixture.H"

#include "lentMethod.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace FrontTracking;

// Test functions

bool areFrontNormalsConsistent(const triSurfaceFront& front)
{
    const labelListList& edgeFaces = front.edgeFaces(); 
    const vectorField& faceNormals = front.faceNormals(); 
    bool testStatus = true; 

    forAll (edgeFaces, I)
    {
        const vector& n0 = faceNormals[edgeFaces[I][0]]; 
        const vector& n1 = faceNormals[edgeFaces[I][1]]; 

        if ((n0 & n1) < 0)
        {
            testStatus = false; 
        }
    }

    return testStatus; 
}

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    triSurfaceFront front(
        IOobject(
            "front.stl",
            "front",
            runTime, 
            IOobject::MUST_READ, 
            IOobject::AUTO_WRITE
        )
    );

    lentMethod lent(front, mesh); 

    lent.calcSearchDistances(searchDistanceSqr, pointSearchDistanceSqr);
    
    lent.reconstructFront(front, signedDistance, pointSignedDistance); 

    front.write(); 

    while (runTime.run()) {

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        twoPhaseProperties.correct();

        lent.calcSignedDistances(
            signedDistance, 
            pointSignedDistance,
            searchDistanceSqr,
            pointSearchDistanceSqr, 
            front
        ); 

        lent.calcHeaviside(heaviside, signedDistance, searchDistanceSqr); 

        lent.reconstructFront(front, signedDistance, pointSignedDistance); 

        if (areFrontNormalsConsistent(front))
        {
            Info << "PASS: front normals consistency" << endl;
        }
        else
        {
            Info << "FAIL: front normals consistency" << endl;
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
