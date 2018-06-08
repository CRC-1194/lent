/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Version:  2.2.x                               
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    Interface advection algorithm of the LENT method.

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Tomislav Maric maric@csi.tu-darmstadt.de (main developer)
    under the project supervision of :
        Holger Marschall <marschall@csi.tu-darmstadt.de> (group leader).
    
    Method Development and Intellectual Property :
    	Tomislav Maric maric@csi.tu-darmstadt.de
    	Holger Marschall <marschall@csi.tu-darmstadt.de>
    	Dieter Bothe <bothe@csi.tu-darmstadt.de>

        Mathematical Modeling and Analysis
        Center of Smart Interfaces
        Technische Universitaet Darmstadt
       
    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
#include "interfaceProperties.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"

// LENT 
#include "lentTests.H"
#include "lentMethod.H"

// Error analysis  
#include "fieldErrorsVolumeFwd.H"
#include "fieldErrorL1.H"
#include "fieldErrorsL1Fwd.H"
#include "fieldErrorL1normalized.H"
#include "fieldErrorsL1normalizedFwd.H"
#include "volScalarFieldErrorBoundedness.H"

// Timing 
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace FrontTracking;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    volScalarField markerFieldInitial(
        IOobject
        (
            "alpha.water.initial", 
            "0",
            mesh, 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ), 
        markerField 
    );

    // Update the advection velocity from function objects and overwrite 
    // intial values.  
    auto& functionObjects = runTime.functionObjects(); 
    functionObjects.execute(); 
    U.write(); 

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    triSurfaceFront front(
        IOobject(
            "front",
            "front",
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    );

    lentMethod lent(front, mesh);

    lent.calcSearchDistances(searchDistanceSqr, pointSearchDistanceSqr);

    lent.calcSignedDistances(
        signedDistance,
        pointSignedDistance,
        searchDistanceSqr,
        pointSearchDistanceSqr,
        front
    );

    lent.reconstructFront(front, signedDistance, pointSignedDistance);

    lent.calcMarkerField(markerField);
    markerField.write(); 
    front.write();

    // Set up the file for error output.
    std::fstream errorFile(args.rootPath() + "/" + args.globalCaseName() 
                           + "/advectionErrors.dat", std::ios_base::app); 
    if (Pstream::myProcNo() == 0)
    {
        errorFile << "time " <<  "CFL " << "Ev " << "Eb " << "Eg " << "En " << "Te " << "\n"; 
    }

    while (runTime.run())
    {
        #include "readTimeControls.H"

        runTime++;

        #include "CourantNo.H"
        #include "markerFieldCourantNo.H"
        #include "setDeltaT.H"

        Info<< "Time = " << runTime.timeName() << nl << endl;

        auto t1 = Clock::now();
        lent.reconstructFront(front, signedDistance, pointSignedDistance);

        lent.evolveFront(front, U.oldTime());

        if (Test::normalsAreInconsistent(front))
        {
            Info << "Inconsistent front normals." << endl;
        }

        lent.calcSignedDistances(
            signedDistance,
            pointSignedDistance,
            searchDistanceSqr,
            pointSearchDistanceSqr,
            front
        );

        lent.calcMarkerField(markerField);
        auto t2 = Clock::now();
        double tTotal = (std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()) / 1000.0;  


        // Update viscosity. 
        mixture.correct();
        // Update density field.
        rho == markerField*rho1 + (scalar(1) - markerField)*rho2;


        volScalarFieldErrorVolume volumeErrorCalculator; 
        volumeErrorCalculator.computeError(markerFieldInitial, markerField); 

        volScalarFieldErrorBoundedness boundednessErrorCalculator;
        boundednessErrorCalculator.computeError(markerFieldInitial, markerField); 

        fieldErrorL1<volScalarField> geometricalErrorCalculator;
        geometricalErrorCalculator.computeError(markerFieldInitial, markerField); 

        // TEST for the geometrical error E_n.
        fieldErrorL1normalized<volScalarField> normalizedErrorCalculator;
        normalizedErrorCalculator.computeError(markerFieldInitial, markerField); 

        if (Pstream::parRun())
        {
            // Collect maximal execution time for a process. 
            auto maxEq = [](scalar& x1, const scalar& x2)
            {
                if (x1 > x2)
                {
                    x1 = x2; 
                    return x1; 
                } 
                return x2; 
            };

            Pstream::gather(tTotal,maxEq); 
            Pstream::scatter(tTotal); 
        }

        if (Pstream::myProcNo() == 0)
        {
            Info << "Ev = " << volumeErrorCalculator.errorValue() << "\n"; 

            errorFile << runTime.timeOutputValue() << " " 
                << CoNum << " " 
                << volumeErrorCalculator.errorValue() << " "
                << boundednessErrorCalculator.errorValue() << " " 
                << geometricalErrorCalculator.errorValue() << " "
                << normalizedErrorCalculator.errorValue() << " "
                << tTotal << "\n";
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
