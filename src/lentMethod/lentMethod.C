/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | 
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

Class
    Foam::lentMethod

SourceFiles
    lentMethod.C

Authors
    Tomislav Maric (maric@csi.tu-darmstadt.de)
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    A facade class that simplifies the user interface for the LENT hybrid
    front-tracking level-set method by agglomerating a bunch of run-time
    selected small one-algorithm SRP classes.

Affiliations:
    Mathematical Modeling and Analysis Institute, Mathematics Department, 
    TU Darmstadt, Germany

Funding:
    German Research Foundation (DFG) - Project-ID 265191195 - SFB 1194

    German Research Foundation (DFG) - Project-ID MA 8465/1-1, 
    Initiation of International Collaboration 
    "Hybrid Level Set / Front Tracking methods for simulating 
    multiphase flows in geometrically complex systems"

\*---------------------------------------------------------------------------*/


#include "lentMethod.H"
#include "error.H"

namespace Foam  {

namespace FrontTracking {

    defineTypeNameAndDebug(lentMethod, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lentMethod::lentMethod(
    const triSurfaceFront& front,
    const fvMesh& mesh,
    word dictName
)
:
    regIOobject(
        IOobject(
            "lentMethod",
            mesh.time().timeName(),
            mesh.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    lentControlDict_(
        IOobject(
           dictName,
           "system",
           mesh.time(),
           IOobject::MUST_READ_IF_MODIFIED,
           IOobject::AUTO_WRITE
        )
    ),
    communicationMaps_(front, mesh),  
    frontIsReconstructed_(false),
    distanceFieldCalculatorTmp_(
       distanceFieldCalculator::New(
           lentControlDict_.subDict("distanceCalculator")
       )
    ),
    frontReconstructionModelTmp_(
       frontReconstructionModel::New(
           lentControlDict_.subDict("frontReconstructionModel")
       )
    ),
    frontReconstructorTmp_(
       frontReconstructor::New(
           lentControlDict_.subDict("frontReconstructor")
       )
    ),
    frontMotionSolverTmp_(
       frontMotionSolver::New(
           lentControlDict_.subDict("frontMotionSolver")
       )
    ),
    markerFieldModelTmp_(
        markerFieldModel::New(
           lentControlDict_.subDict("markerFieldModel")
        )
    ),
    surfaceTensionForceModelTmp_(
        frontSurfaceTensionForceModel::New(
           lentControlDict_.subDict("surfaceTensionForceModel") 
        )
    ),
    frontSmoother_(
        lentControlDict_.subDict("frontSmoother")
    ),
    reconstructionHistory_(
        mesh.time()
    )
{
    useReconstruction_ =
        lentControlDict_.subDict("frontReconstructionModel").lookupOrDefault("useReconstruction", true);
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void lentMethod::calcSearchDistances(
    volScalarField& searchDistanceSqr,
    pointScalarField& pointSearchDistanceSqr
)
{
    distanceFieldCalculator& distanceCalc = distanceFieldCalculatorTmp_.ref();

    Info << "Calculating search distances..." << endl;
    distanceCalc.calcCellSearchDistance(searchDistanceSqr);
    distanceCalc.calcPointSearchDistance(pointSearchDistanceSqr, searchDistanceSqr);
    Info << "Done." << endl; 
}

void lentMethod::calcSignedDistances(
    volScalarField& signedDistance,
    pointScalarField& pointSignedDistance,
    const volScalarField& searchDistanceSqr,
    const pointScalarField& pointSearchDistanceSqr,
    const triSurfaceFront& front
)
{
    // FIXME: The distance calculator wrongy re-initializes the frontMesh. TM. 
    Info << "Calculating sign distances..." << endl;
    distanceFieldCalculatorTmp_->calcCellsToFrontDistance(
        signedDistance,
        searchDistanceSqr,
        front
    );

    // FIXME: The  distance calculator wrongy re-initializes the frontMesh. TM.
    distanceFieldCalculatorTmp_->calcPointsToFrontDistance(
        pointSignedDistance,
        pointSearchDistanceSqr,
        front
    );
    Info << "Done." << endl; 
}

void lentMethod::calcMarkerField(volScalarField& markerField) const
{
    Info << "Calculating marker field..." << endl;
    markerFieldModelTmp_->calcMarkerField(markerField);
    Info << "Done." << endl;
}

void lentMethod::reconstructFront(
    triSurfaceFront& front,
    const volScalarField& signedDistance,
    const pointScalarField& pointSignedDistance
)
{
    if (frontReconstructionModelTmp_->reconstructionRequired(front, signedDistance))
    {
        Info << "Reconstructing front..." << endl;


        if (useReconstruction_ || signedDistance.time().timeIndex() <= 1)
        {
            frontReconstructorTmp_->reconstructFront(
                front,
                signedDistance,
                pointSignedDistance
            );
        }

        if (front.surfaceType() != triSurface::MANIFOLD)
        {
            Info << "Warning: front is not a Manifold anymore." << endl;
        }

        reconstructionHistory_.frontReconstructed();

        frontSmoother_.smoothFront(front, signedDistance.mesh());

        frontIsReconstructed_ = true;

        Info << "Done." << endl;
    }
}

void lentMethod::calcFrontVelocity(
    triSurfaceFrontPointVectorField& frontVelocity,
    const volVectorField& U
)
{
    Info << "Calculating front velocity..." << endl;  

    const triSurface& front = frontVelocity.mesh();

    frontVelocity.resize(front.nPoints());
    // More rigorous: in case the search fails, the point stops. TM.
    frontVelocity = dimensionedVector("zero", dimVelocity, vector(0,0,0)); 

    // FIXME: Make this an attribute of method and re-use. Enable selection of 
    // cell->point interpolation. TM.
    lentInterpolation interpolation; 
    interpolation.interpolate(U, frontVelocity); 

    Info << "Done." << endl;
}

void lentMethod::evolveFront(
    triSurfaceFront& front,
    const volVectorField& cellVelocity
) 
{
    Info << "Evolving the front..." << endl;  
    frontMotionSolverTmp_->evolveFront(
        front,
        cellVelocity
    );
    Info << "Done." << endl;

    frontIsReconstructed_ = false;

    // Update front-mesh communication maps after front motion. 
    Info << "Updating communication maps..." << endl;  
    communicationMaps_.update(); 
    Info << "Done." << endl;
}

bool lentMethod::writeData(Ostream&) const
{
    FatalErrorIn("lentMethod::writeData(Ostream& os)")
    << "lentMethod is not supposed to be written "
    << "regIOobject inherited to allow registry queries." << endl;

    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
