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

Authors
    Tomislav Maric maric@mma.tu-darmstadt.de
    Tobias Tolle tolle@mma.tu-darmstadt.de

Description
    A two-phase Level Set / Front Tracking DNS solver.

\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
#include "EulerDdtScheme.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "lentMethod.H"
#include "lentSolutionControl.H"
#include "analyticalSurface.H"
#include "alphaFace.H"

#include <limits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace FrontTracking;

void correctFrontIfRequested(triSurfaceFront& front, const dictionary& configDict)
{
    if (configDict.found("frontSurface"))
    {
        auto surfaceTmp = analyticalSurface::New(configDict.subDict("frontSurface"));
        const auto& frontSurface = surfaceTmp.ref();
        frontSurface.moveFrontToSurface(front);
        frontSurface.makeNormalOrientationConsistent(front);
        Info << "Front has been corrected" << endl;
    }
}

const auto EPSILON = std::numeric_limits<double>::epsilon();

int main(int argc, char *argv[])
{
    argList::addBoolOption
    (
        "relative-frame",
        "Arbitrary Eulerian / Lagrangian (ALE) relative reference frame (RRF), using the vertical velocity of the dispersed phase (alpha.water=0)."
    );

    argList::addBoolOption
    (
        "normal-velocity",
        "Evolve the Front using the velocity projected onto the interface normal."
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createMRF.H"
    #include "createFvOptions.H"

    lentSolutionControl lentSC(mesh, phi);

    #include "correctPhi.H"

    turbulence->validate();

    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    triSurfaceFront front(
        IOobject(
            "front",
            "front",
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );

    triSurfacePointVectorField frontVelocity(
        IOobject(
            "frontVelocity",
            runTime.timeName(),
            mesh,
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
    // Project points onto the input surface after initial reconstruction.
    correctFrontIfRequested(front, lent.dict());
    front.write();

    // ALE Relative reference frame data for simulating rising bubbles. 
    dimensionedVector Ub ("Ub", dimVelocity, vector(0,0,0));
    surfaceScalarField rhof("rhof", fvc::interpolate(rho));
    surfaceScalarField muf(mixture.muf());
    volScalarField alphaInv{dimensionedScalar("1", dimless, 1) - markerField};
    volVectorField Ufront ("Ufront", U);   
    volVectorField nFront ("nFront", fvc::grad(signedDistance));
	
    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "markerFieldCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info << "Time step = " << runTime.timeIndex() << endl;
        Info << "Time = " << runTime.timeName() << nl << endl;

	    Ufront == U; 

        if (args.found("relative-frame"))
        {
            alphaInv = dimensionedScalar("1", dimless, 1) - markerField;
            Ub = sum(alphaInv * mesh.V() * U) / sum(alphaInv * mesh.V());
	        Ufront == Ufront - Ub;
        }
        if (args.found("normal-velocity"))
        {
            nFront = fvc::grad(signedDistance); 
            nFront /= Foam::mag(nFront) + dimensionedScalar("EPSILON", dimless, EPSILON);
            Ufront == (Ufront & nFront) * nFront;
	    }

        lent.reconstructFront(front, signedDistance, pointSignedDistance);

        lent.evolveFront(front, Ufront);

        lent.calcSignedDistances(
            signedDistance,
            pointSignedDistance,
            searchDistanceSqr,
            pointSearchDistanceSqr,
            front
        );
        
        lent.calcMarkerField(markerField);
        mixture.correct();
        
        // Face densities computed from face indicator function 
        // based on signed-distances. 
        #include "computeRhof.H"

        // --- SAAMPLE loop
        while (lentSC.loop())
        {
            if (lentSC.updateMomentumFlux())
            {
               rhoPhi == rhof * phi;

               fvScalarMatrix rhoEqn
                (
                    fvm::ddt(rho) + fvc::div(rhoPhi)
                 );
               rhoEqn.solve();
            }

            #include "UEqn.H"

            //--- Pressure corrector loop
            while (lentSC.correctPressure())
            {
                #include "pEqn.H"
            }

            if (lentSC.turbCorr())
            {
                turbulence->correct();
            }
        }

        rho == markerField*rho1 + (1.0 - markerField)*rho2;

        runTime.write();
        
        // This is a workaround to ensure the actual front mesh is written (TT)
        if (runTime.writeTime())
        {
            front.write(); 
        }

        Info << "Writing time = " << runTime.cpuTimeIncrement() << endl;
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
