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

Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de
    tomislav.maric@gmx.com
    Mathematical Modelling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany

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

    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readTimeControls.H"
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

    // TODO: use triSurfaceFront fields, put this in createFields. 
    DynamicField<vector> frontDisplacement (front.nPoints()); 

    IOdictionary levelSetFrontDict(
        IOobject(
            "levelSetFrontDict", 
            "constant", 
            runTime, 
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::AUTO_WRITE
        )
    );

    // FIXME: move the alg config to fvSolution::lent
    label narrowBandWidth = 
        levelSetFrontDict.lookupOrDefault<label>("narrowBandWidth", 4);

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

    // Write the front.
    runTime.writeNow(); 

    while (runTime.run())
    {
        #include "readTimeControls.H"

        runTime++;

        #include "CourantNo.H"
        #include "heavisideCourantNo.H"
        #include "setDeltaT.H"

        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        twoPhaseProperties.correct();
        
        // Move the front points with the constant vector : used for testing.
        //front.move(constDisplacement*runTime.deltaT().value()); 
        
        // Compute the velocity using the meshCells of the isoSurface reconstruction.
        calc.calcFrontVelocity(frontDisplacement, front, U); 
        //// FIXME: Put this in thecalcFrontVelocity function and scale the displacement. 
        frontDisplacement *= runTime.deltaT().value(); 
        front.move(frontDisplacement);
        
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
        signedDistance.time().cpuTimeIncrement(); 

        //New meshCells() information computed.  
        // TODO: user defined reconstruction interval.  
        //if (runTime.timeIndex() % 10 == 0)
        //{
            //front.reconstruct(signedDistance, pointSignedDistance, false, 1e-10); 
        //}

        Info << "Front reconstructed: " 
            << signedDistance.time().cpuTimeIncrement() << endl; 

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
