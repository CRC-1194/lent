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
    Tomislav Maric
    maric<<at>>csi<<dot>>tu<<minus>>darmstadt<<dot>>de
    tomislav<<dot>>maric<<at>>gmx<<dot>>com

Description
    A DNS two-phase flow solver employing a hybrid level-set / front-tracking
    method.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "interfaceProperties.H"
#include "incompressibleTwoPhaseMixture.H"

#include "lentMethod.H"

//#include "triSurfaceFront.H"
//#include "naiveNarrowBandPropagation.H"
//#include "TriSurfaceMeshCalculator.H"
//#include "leftAlgorithmHeavisideFunction.H"

// Fields.
#include "DynamicField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace FrontTracking;

// Configure the Method
//typedef triSurfaceFront Front;
//typedef TriSurfaceMeshCalculator Calculator; 

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

    lentMethod lent2(front, mesh, "lentSolutionSharp");

    lent2.calcSignedDistanceFields(
        signedDistance, 
        pointSignedDistance,
        searchDistanceSqr
    ); 

    lent.calcSignedDistanceFields(
        signedDistance, 
        pointSignedDistance,
        searchDistanceSqr
    ); 

    //calc.calcCentresToElementsDistance(
        //signedDistance, 
        //front,
        //naiveNarrowBandPropagation()
    //); 

    //calc.calcPointsToElementsDistance(
        //pointSignedDistance, 
        //front,
        //mesh, 
        //naiveNarrowBandPropagation()
    //); 
    
    lent2.calcHeavisideField(
        heaviside,
        signedDistance,
        searchDistanceSqr
    ); 
    
    lent.calcHeavisideField(
        heaviside,
        signedDistance,
        searchDistanceSqr
    ); 

    heaviside.write(); 

    while (runTime.run()) {
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        //Reconstruct the front: 
        //front.reconstruct(signedDistance, pointSignedDistance, false, 1e-10); 
        
        //lent.reconstructFront(front, signedDistance, pointSignedDistance); //signedDistance, pointSignedDistance, false, 1e-10); 

        twoPhaseProperties.correct();

        // Compute the new signed distance field with the surfaceMesh octree
        // search.  
        //calc.calcCentresToElementsDistance(
            //signedDistance, 
            //front,
            //naiveNarrowBandPropagation()
        //); 

        // Compute the new signed point distance field. 
        //calc.calcPointsToElementsDistance(
            //pointSignedDistance, 
            //front,
            //mesh, 
            //naiveNarrowBandPropagation()
        //); 

        lent.calcSignedDistanceFields(
            signedDistance, 
            pointSignedDistance,
            searchDistanceSqr
        ); 

        lent.calcHeavisideField(
            heaviside, 
            signedDistance, 
            searchDistanceSqr
        ); 

        lent2.calcHeavisideField(
            heaviside, 
            signedDistance, 
            searchDistanceSqr
        ); 

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
