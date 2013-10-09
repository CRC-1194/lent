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
    lentAdvection

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

#include "triSurfaceFront.H"
#include "naiveNarrowBandPropagation.H"
#include "TriSurfaceMeshCalculator.H"
#include "leftAlgorithmHeavisideFunction.H"

// TODO: switch to registered front Fields.
#include "DynamicField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace FrontTracking;

// Configure the LENT Algorithm 
typedef triSurfaceFront Front;
typedef TriSurfaceMeshCalculator Calculator; 

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readTimeControls.H"
    #include "correctPhi.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    Front front(
        IOobject(
            "front.stl",
            "front",
            runTime, 
            IOobject::MUST_READ, 
            IOobject::AUTO_WRITE
        )
    );

    autoPtr<heavisideFunction> heavisideModelPtr = 
        heavisideFunction::New("leftAlgorithmHeaviside"); 

    DynamicField<vector> frontDisplacement(front.nPoints()); 

    IOdictionary levelSetFrontDict(
        IOobject(
            "levelSetFrontDict", 
            "constant", 
            runTime, 
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::AUTO_WRITE
        )
    );

    label narrowBandWidth = 
        levelSetFrontDict.lookupOrDefault<label>("narrowBandWidth", 4);

    // FIXME: as soon as the solver is running: refactor the calculation
    //          * separate the front from front operations
    //          * separate calculator from calculations (command pattern)
    //          * make it possible to assemble calculators from commands
    //          * clean up the solver interface
    //          * bind calculators to mesh.update for motion/topology
    //            changes
    Calculator calc(narrowBandWidth); 

    // Compute the new signed distance field. 
    calc.calcCentresToElementsDistance(
        signedDistance, 
        front,
        naiveNarrowBandPropagation()
    ); 

    calc.calcPointsToElementsDistance(
        pointSignedDistance, 
        front,
        mesh, 
        naiveNarrowBandPropagation()
    ); 

    heavisideModelPtr->calcHeaviside(
        heaviside, 
        signedDistance, 
        calc.getCellSearchDistSqr()
    ); 

    //Reconstruct the front. 
    front.reconstruct(signedDistance, pointSignedDistance, false, 1e-10); 

    // Update the density field. 
    rho == heaviside*rho1 + (scalar(1) - heaviside)*rho2;

    // Write the front.
    runTime.writeNow(); 

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "heavisideCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        // Correct the curvature and the properties before going into the 
        // pressure velocity coupling cycle.
        twoPhaseProperties.correct();
        interface.correct(); 

        // Update the density field. 
        rho == heaviside*rho1 + (scalar(1) - heaviside)*rho2;

        Info << "p-U algorithm ... " << endl;
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"

            //--- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }
        Info << "Done." << endl; 

        Info << "Final velocity calculation..." << endl;
        #include "UEqn.H"
        Info << "Done. " << endl;
        
        // Compute the velocity using the meshCells of the isoSurface reconstruction.
        // FIXME: computing front displacement, see below
        Info << "Computing front velocity..." << endl;
        calc.calcFrontVelocity(frontDisplacement, front, U); 
        Info << "Done." << endl;

        // FIXME: Put this in calcFrontVelocity function. 
        Info << "Moving the front..." << endl;
        frontDisplacement *= runTime.deltaT().value(); 
        front.move(frontDisplacement);
        Info << "Done." << endl;
        
        // Compute the new signed distance field with the surfaceMesh octree search.  
        calc.calcCentresToElementsDistance(
            signedDistance, 
            front,
            naiveNarrowBandPropagation()
        ); 

        // Compute the new signed point distance field. 
        calc.calcPointsToElementsDistance(
            pointSignedDistance, 
            front,
            mesh, 
            naiveNarrowBandPropagation()
        ); 

        heavisideModelPtr->calcHeaviside(
            heaviside, 
            signedDistance, 
            calc.getCellSearchDistSqr()
        ); 

        //Reconstruct the front: 
        if (runTime.timeIndex() % 10 == 0)
        {
            front.reconstruct(signedDistance, pointSignedDistance, false, 1e-10); 
        }

        runTime.write();

        Info << "Writing time = " << runTime.cpuTimeIncrement() << endl;
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
