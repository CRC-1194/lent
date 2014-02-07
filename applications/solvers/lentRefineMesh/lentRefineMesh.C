/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    lentRefineMesh 

Description
    Interface advection with the LENT method coupled with local AMR in OpenFOAM. 

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "interfaceProperties.H"
#include "incompressibleTwoPhaseMixture.H"
#include "fvIOoptionList.H"

#include "lentMethod.H"

using namespace FrontTracking;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

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

    front.write(); 

    while (runTime.run())
    {
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        Info << "SIZE = " <<  heaviside.size() << endl;
        // Apply mesh changes

        mesh.update();
        Info << "SIZE = " <<  heaviside.size() << endl;

        twoPhaseProperties.correct();

        // Naive implementation. 
        // The strategy for distance re-calculation based on a refinement map 
        // needs to be added to the lentMethod class. TM Feb 07 14
        lent.calcSearchDistances(searchDistanceSqr, pointSearchDistanceSqr);

        lent.calcSignedDistances(
            signedDistance, 
            pointSignedDistance,
            searchDistanceSqr,
            pointSearchDistanceSqr, 
            front
        ); 

        lent.calcHeaviside(heaviside, signedDistance, searchDistanceSqr); 


        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
