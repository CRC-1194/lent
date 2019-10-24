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

int main(int argc, char *argv[])
{
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

    // Lets make a fair comparison: recover exact surface!
    correctFrontIfRequested(front, lent.dict());

    front.write();

    // Explicit extrapolation of fields for improved handling of non-linearity in the
    // pressure - velocity coupling. Based on "Consistent second-order time-accurate
    // non-iterative PISO-algorithm, Tukovic, Peric, Jasak. TODO: Proper reference.
    // volScalarField pn ("pn", p);
    surfaceScalarField phistar ("phistar", phi);
    surfaceScalarField phin ("phin", phi);
    volScalarField signedDistancen("signedDistancen", signedDistance);
    pointScalarField pointSignedDistancen("pointSignedDistancen", pointSignedDistance);
    volScalarField signedDistanceStar("signedDistanceStar", signedDistance);
    pointScalarField pointSignedDistanceStar("pointSignedDistanceStar", pointSignedDistance);

	
    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "markerFieldCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info << "Time step = " << runTime.timeIndex() << endl;
        Info << "Time = " << runTime.timeName() << nl << endl;

	//pn == p;
	//phin == phi;
	//signedDistancen = signedDistance;
	//pointSignedDistancen = pointSignedDistance;

        // --- Pressure-velocity lentSolutionControl corrector loop
        while (lentSC.loop())
        {
            // The momentum flux is computed from MULES as  
            // rhoPhi = phiAlpha*(rho1 - rho2) + phi*rho2; 
            // However, LENT has no ability to compute the volumetric phase flux. 
            // TODO: examine the impact of the momentum flux computation and devise
            // more accurate approach if required (TT)
            if (lentSC.updateMomentumFlux())
            {
                if (lent.dict().subDict("markerFieldModel").get<label>("nSmoothingSteps") > 0)
                {
                    // old approach, only works for diffuse markerfield
                    rhoPhi == fvc::interpolate(rho) * phi;
                }
                else
                {
                    // new approach: vol fraction based calculation of rho at the face
                    // only works for a sharp, vol-fraction like markerfield
                    // FIXME: Face-fractions are recomputed in the external loop, and
                    // they only change between time steps: extract the face fraction
                    // calculation out of the outer loop. TM.
                    #include "computeRhoPhi.H"
                }
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

        lent.evolveFront(front, U.oldTime());

        lent.calcSignedDistances(
            signedDistance,
            pointSignedDistance,
            searchDistanceSqr,
            pointSearchDistanceSqr,
            front
        );

        lent.reconstructFront(front, signedDistance, pointSignedDistance);

        if (runTime.timeIndex() <= 1)
        {
            correctFrontIfRequested(front, lent.dict());
        }
        
        // TODO: if front has been reconstructed, the signed distances (at least)
        // for the cell centres have to be recomputed to ensure the front-mesh
        // communication is up-to-date (TT)
        if (lent.isFrontReconstructed())
        {
            lent.calcSignedDistances(
                signedDistance,
                pointSignedDistance,
                searchDistanceSqr,
                pointSearchDistanceSqr,
                front
            );
        }

	// Second-order field extrapolation.
	//p == 2*p - pn;
	//phi == 2*phi - phin;
	//signedDistance == 2*signedDistance - signedDistancen; 
	//pointSignedDistance == 2*pointSignedDistance - pointSignedDistancen;

        lent.calcMarkerField(markerField);

        // Update the viscosity. 
        mixture.correct();

        // Update density field.
        rho == markerField*rho1 + (scalar(1) - markerField)*rho2;

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
