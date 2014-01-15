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
    lentAdvect

Authors
    Tomislav Maric
    maric<<at>>csi<<dot>>tu<<minus>>darmstadt<<dot>>de
    tomislav<<dot>>maric<<at>>gmx<<dot>>com

Description
    Test application for the interface advection algorithm of the LENT method.  

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "incompressibleTwoPhaseMixture.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "lentMethod.H"
#include "lentTests.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace FrontTracking;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    triSurfaceFront front(
        IOobject(
            "front.stl",
            "front",
            runTime, 
            IOobject::NO_READ, 
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

    Test::testTriSurfaceNormalConsistency(front); 
    
    front.write(); 

    while (runTime.run())
    {
        #include "readTimeControls.H"

        runTime++;

        #include "CourantNo.H"
        #include "heavisideCourantNo.H"
        #include "setDeltaT.H"

        Info<< "Time = " << runTime.timeName() << nl << endl;

        twoPhaseProperties.correct();

        Pout << "Signed distances..."; 
        lent.calcSignedDistances(
            signedDistance, 
            pointSignedDistance,
            searchDistanceSqr,
            pointSearchDistanceSqr, 
            front
        ); 
        Pout << "done." << endl;

        Pout << "Heaviside ... "; 
        lent.calcHeaviside(heaviside, signedDistance, searchDistanceSqr); 
        Pout << "done." << endl;

        Pout << "Reconstruction ..."; 
        lent.reconstructFront(front, signedDistance, pointSignedDistance); 
        Pout << "done." << endl;

        Test::testTriSurfaceNormalConsistency(front); 

        Pout << "Velocity ..."; 
        lent.calcFrontVelocity(frontVelocity, U); 
        Pout << "done." << endl;

        Pout << "Evolution ...";
        lent.evolveFront(front, frontVelocity); 
        Pout << "done." << endl;
        
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
