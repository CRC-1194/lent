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
    lentAdvectRefined 

Description
    Interface advection with the LENT method coupled with local AMR in OpenFOAM. 
    
Authors
    Tomislav Maric maric@csi.tu-darmstadt.de, tomislav@sourceflux.de
    Mathematical Modeling and Analysis
    Center of Smart Interfaces, TU Darmstadt

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "interfaceProperties.H"
#include "incompressibleTwoPhaseMixture.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "lentMethod.H"

#include <iostream>
#include <chrono>

using namespace FrontTracking;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"

    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    using std::chrono::seconds;
    using std::chrono::duration_cast;

    typedef std::chrono::high_resolution_clock Clock;

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

    triSurfaceFrontGeoMesh frontMesh(front); 

    triSurfaceFrontVectorField frontVelocity(
        IOobject(
            "frontVelocity", 
            runTime.timeName(), 
            runTime, 
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        front, 
        dimensionedVector(
            "zero",
            dimLength / dimTime, 
            vector(0,0,0)
        )
    ); 

    lentMethod lent(front, mesh); 

    lent.calcSearchDistances(searchDistanceSqr, pointSearchDistanceSqr);

    lent.reconstructFront(front, signedDistance, pointSignedDistance); 

    front.write(); 

    while (runTime.run())
    {
        #include "readTimeControls.H"

        runTime++;

        #include "CourantNo.H"
        #include "heavisideCourantNo.H"
        #include "setDeltaT.H"

        Info<< "Time = " << runTime.timeName() << nl << endl;

        auto t1 = Clock::now();
        mesh.update();
        auto t2 = Clock::now();
        auto diff = duration_cast<seconds>(t2 - t1); 
        std::cout << "Mesh update : " << diff.count() / 1e06 << " \n";

        t1 = Clock::now(); 
        lent.calcSearchDistances(searchDistanceSqr, pointSearchDistanceSqr);
        t2 = Clock::now(); 
        diff = duration_cast<seconds>(t2 - t1); 
        std::cout << "Search distance calculation : " << diff.count() << " \n";

        t1 = Clock::now(); 
        lent.calcSignedDistances(
            signedDistance, 
            pointSignedDistance,
            searchDistanceSqr,
            pointSearchDistanceSqr, 
            front
        ); 
        t2 = Clock::now(); 
        diff = duration_cast<seconds>(t2 - t1); 
        std::cout << "Signed distance calculation : " << diff.count() << " \n";

        t1 = Clock::now(); 
        lent.calcHeaviside(heaviside, signedDistance, searchDistanceSqr); 
        t2 = Clock::now(); 
        diff = duration_cast<seconds>(t2 - t1); 
        std::cout << "Heaviside calculation : " << diff.count() << " \n";

        t1 = Clock::now(); 
        twoPhaseProperties.correct();
        t2 = Clock::now(); 
        diff = duration_cast<seconds>(t2 - t1); 
        std::cout << "Properties update : " << diff.count() << " \n";

        t1 = Clock::now(); 
        lent.reconstructFront(front, signedDistance, pointSignedDistance); 
        t2 = Clock::now(); 
        diff = duration_cast<seconds>(t2 - t1); 
        std::cout << "Front reconstruction : " << diff.count() << " \n";

        t1 = Clock::now(); 
        lent.calcFrontVelocity(frontVelocity, U); 
        t2 = Clock::now(); 
        diff = duration_cast<seconds>(t2 - t1); 
        std::cout << "Front velocity calculation :" << diff.count() << " \n";

        t1 = Clock::now(); 
        lent.evolveFront(front, frontVelocity); 
        t2 = Clock::now(); 
        diff = duration_cast<seconds>(t2 - t1); 
        std::cout << "Front advection :" << diff.count() << " \n";
        
        t1 = Clock::now(); 
        runTime.write();
        t2 = Clock::now(); 
        diff = duration_cast<seconds>(t2 - t1); 
        std::cout << "Writing data:" << diff.count() << " \n";

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
