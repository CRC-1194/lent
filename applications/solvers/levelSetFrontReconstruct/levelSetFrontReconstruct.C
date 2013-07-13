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
#include "FvMeshAndFrontConnection.H"
#include "naiveNarrowBandPropagation.H"
#include "TriSurfaceMeshDistanceCalculator.H"
#include "TriSurfaceMeshCalculator.H"

// Fields.
#include "DynamicField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace FrontTracking;

// Configure the Method
typedef triSurfaceFront Front;
typedef FvMeshAndFrontConnection<Front> Connection;
typedef TriSurfaceMeshCalculator Calculator; 

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "frontInputFile",
        "pathname to the surface mesh file used for initializing the front"
    );

    argList::addOption
    (
        "frontOutputDirectory",
        "name of the directory where the front will be stored"
    );

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

    fileName frontInputFile = "frontInputFile.stl";
    fileName frontOutputDirectory = "front";

    if (args.optionFound("frontInputFile"))
    {
        frontInputFile = args.optionRead<word>("frontInputFile");
    }
    if (args.optionFound("frontOutputDirectory"))
    {
        frontOutputDirectory = 
            args.optionRead<word>("frontOutputDirectory");
    }

    Info<< "\nStarting time loop\n" << endl;

    Front front (
        IOobject (
            frontInputFile,
            frontOutputDirectory,
            runTime, 
            IOobject::MUST_READ, 
            IOobject::AUTO_WRITE
        )
    );

    //Front oldFront(front); 
    //oldFront.rename("oldFront.vtk"); 

    // TODO: use triSurfaceFront fields, put this in createFields. 
    DynamicField<vector> frontDisplacement (front.nPoints()); 

    IOdictionary levelSetFrontDict
    (
        IOobject
        (
            "levelSetFrontDict", 
            "constant", 
            runTime, 
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::AUTO_WRITE
        )
    );

    label narrowBandWidth = levelSetFrontDict.lookupOrDefault<label>("narrowBandWidth", 4);

    Calculator calc (narrowBandWidth); 

    // Compute the new signed distance field. 
    calc.calcCentresToElementsDistance
    (
        Psi, 
        front,
        naiveNarrowBandPropagation()
    ); 

    calc.calcPointsToElementsDistance
    (
        psi, 
        front,
        mesh, 
        naiveNarrowBandPropagation()
    ); 

    //Reconstruct the front. 
    front.reconstruct(Psi, psi, false, 1e-10); 
    //front.reconstruct(Psi, false, 1e-10); 

    // Write the front.
    runTime.writeNow(); 

    vector constDisplacement = levelSetFrontDict.lookup("displacement"); 

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Compute the new signed distance field with the surfaceMesh octree search.  
        calc.calcCentresToElementsDistance
        (
            Psi, 
            front,
            naiveNarrowBandPropagation()
        ); 

        // Compute the new signed point distance field. 
        calc.calcPointsToElementsDistance
        (
            psi, 
            front,
            mesh, 
            naiveNarrowBandPropagation()
        ); 


        //Reconstruct the front: 
        Psi.time().cpuTimeIncrement(); 

        front.reconstruct(Psi, psi, false, 1e-10); 

        Info << "Front reconstructed: " 
            << Psi.time().cpuTimeIncrement() << endl; 

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
