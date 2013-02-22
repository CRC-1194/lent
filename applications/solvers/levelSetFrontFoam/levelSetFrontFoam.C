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

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "twoPhaseMixture.H"
#include "turbulenceModel.H"
#include "interpolationTable.H"
#include "pimpleControl.H"


#include "levelSetFront.H"
#include "levelSetFrontFields.H"

#include "meshAndFrontConnection.H"
#include "frontTrackingCalculator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace frontTracking;

typedef meshAndFrontConnection<fvMesh, levelSetFront> Connection;
typedef frontTrackingCalculator<Connection>  Calculator;

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

    // Initialize the levelSetFront from a surface mesh file.
    levelSetFront front (
        IOobject (
            "front.stl", 
            "front", 
            runTime, 
            IOobject::MUST_READ, 
            IOobject::AUTO_WRITE
        ) 
    );

    // Write the initial front. 
    runTime.writeNow();
    front.writeNow(runTime); 

    ++runTime;

    // Initialize the connectivity between the mesh and the front.
    Connection frontMeshConnection(mesh,front);
    
    // Initialize the field calculator 
    Calculator calculator;

    // Compute the cell centered cell to elements distance field.
    calculator.calcCentresToElementsDistance(Psi, frontMeshConnection);
    calculator.calcPointsToElementsDistance(psi, frontMeshConnection);

    // Reconstruct the iso-surface front from the distance field.
    front.reconstruct(Psi,psi);

    runTime.writeNow();
    // TODO: register front to the time, so that it is written automatically.
    front.writeNow(runTime); 

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Compute the displacement field.
        
        // Get the number of front vertices.

        // Get the velocity vector from the dictionary.

        // Initialize the displacement vector field.

        // Move the front points with the constant vector: test  
        // Notification of the front/mesh motion/topological is done via the
        // observer pattern. 
        front.move(vector(1,1,1) * runTime.deltaT().value());

        // Compute the new distance fields. 
        calculator.calcCentresToElementsDistance(Psi, frontMeshConnection); 
        calculator.calcPointsToElementsDistance(psi, frontMeshConnection); 

        //// Reconstruct the front as an iso surface. 
        front.reconstruct(Psi, psi); 
        
        // Update two phase properties. TODO: new model for the iso-surface 
        // properties
        // twoPhaseProperties.correct();

        //// --- Pressure-velocity PIMPLE corrector loop
        ////while (pimple.loop())
        ////{
            ////#include "UEqn.H"

            //// --- Pressure corrector loop
            ////while (pimple.correct())
            ////{
                ////#include "pEqn.H"
            ////}

            ////if (pimple.turbCorr())
            ////{
                ////turbulence->correct();
            ////}
        ////}
        
        
        runTime.write();
        front.write(runTime);

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
