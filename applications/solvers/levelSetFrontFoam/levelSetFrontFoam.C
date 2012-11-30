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

#include "DynamicListDLListMap.H"
#include "FrontAndMeshConnectivity.H"

#include "FrontAndMeshCommunication.H"

// Define FrontAndMeshCommunication

using namespace frontTracking;

typedef FrontAndMeshConnectivity<fvMesh, 
    levelSetFront, DynamicListDLListMap> frontAndMeshConnectivity;

//typedef FrontAndMeshConnectivity<fvMesh, 
    //levelSetFront, std::multimap> frontAndMeshConnectivity;
    //
//typedef FrontAndMeshConnectivity<fvMesh, 
    //levelSetFront, std::unordered_multimap> frontAndMeshConnectivity;
    
//typedef FrontAndMeshCommunication<fvMesh, 
    //levelSetFront, frontAndMeshConnectivity> frontAndMeshCommunication;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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

    levelSetFront front (
        IOobject (
            "levelSetFront", 
            runTime.timeName(), 
            runTime, 
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        )
    );

    frontAndMeshCommunication frontAndMesh(mesh, front);


    //while (runTime.run())
    //{
        //#include "readTimeControls.H"
        //#include "CourantNo.H"
        //#include "alphaCourantNo.H"
        //#include "setDeltaT.H"

        //runTime++;

        //Info<< "Time = " << runTime.timeName() << nl << endl;

        // Compute the level set field.

        // Compute the marker field from the level set field.
        
        // Update two phase properties.
        //twoPhaseProperties.correct();

        // --- Pressure-velocity PIMPLE corrector loop
        //while (pimple.loop())
        //{
            //#include "UEqn.H"

            // --- Pressure corrector loop
            //while (pimple.correct())
            //{
                //#include "pEqn.H"
            //}

            //if (pimple.turbCorr())
            //{
                //turbulence->correct();
            //}
        //}


        //runTime.write();

        //Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            //<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
            //<< nl << endl;
    //}

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
