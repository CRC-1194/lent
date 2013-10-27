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
    A DNS two-phase flow solver employing a hybrid level-set / front-tracking
    method.

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

// TODO: switch to registered front fields.
#include "DynamicField.H"

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
            IOobject::MUST_READ, 
            IOobject::AUTO_WRITE
        )
    );

    DynamicField<vector> frontVelocity(front.nPoints()); 

    lentMethod lent(front, mesh); 

    lent.calcHeaviside(heaviside, signedDistance, searchDistanceSqr); 

    //Read the front from the runTime the front, or reconstruct it if it is not 
    //there. 
    //front.reconstruct(signedDistance, pointSignedDistance); 
    
    heaviside.write(); 

    while (runTime.run())
    {
        #include "readTimeControls.H"

        runTime++;

        #include "CourantNo.H"
        #include "heavisideCourantNo.H"
        #include "setDeltaT.H"

        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        twoPhaseProperties.correct();

        // Compute the velocity using the meshCells of the isoSurface reconstruction.
        lent.calcFrontVelocity(frontVelocity, U); 

        lent.evolveFront(frontVelocity); 
        
        // Compute the new signed distance field with the surfaceMesh octree search.  
        lent.calcSignedDistances(
            signedDistance, 
            pointSignedDistance,
            searchDistanceSqr,
            pointSearchDistanceSqr, 
            front
        ); 

        lent.calcHeaviside(heaviside, signedDistance, searchDistanceSqr); 

        // FIXME: add reconstruction model (every N timesteps, coalescenceModel)..
        front.reconstruct(signedDistance, pointSignedDistance); 

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
