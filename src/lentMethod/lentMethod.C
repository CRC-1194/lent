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

Class
    Foam::lentMethod

SourceFiles
    lentMethod.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    A facade class that simplifies the user interface for the LENT hybrid
    front-tracking level-set method by agglomerating a bunch of run-time
    selected small one-algorithm SRP classes.

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
    )
{}

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

        frontReconstructorTmp_->reconstructFront(
            front,
            signedDistance,
            pointSignedDistance
        );

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

    auto oldVelocity(frontVelocity); 

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

    // Clean up degenerate triangles.
    Info << "Cleaning up degeneracies..." << endl;  
    front.cleanup(false);
    Info << "Done." << endl;

    // Calculate normal vectors after front motion.
    Info << "Computing triangle normal vectors..." << endl;  
    calcFrontNormals(front); 
    Info << "Done." << endl;
    // Update front-mesh communication maps after front motion. 
    Info << "Updating communication maps..." << endl;  
    communicationMaps_.update(); 
    Info << "Done." << endl;
}

bool lentMethod::writeData(Ostream& os) const
{
    FatalErrorIn("lentMethod::writeData(Ostream& os)")
    << "lentMethod is not supposed to be written "
    << "regIOobject inherited to allow registry queries." << endl;

    return false;
}

void lentMethod::calcFrontNormals(triSurfaceFront& front) const
{
    // Disambiguate from regIOobject, multiple inheritance issue. TM.
    // Required for registering fields to the front.
    const triSurface& frontSurface = front; 

    auto& normals = front.storedFaceNormals(); 
    const auto& points = front.points(); 

    forAll(normals, faceI)
    {
        normals[faceI] = frontSurface[faceI].normal(points);  
        normals[faceI] /= mag(normals[faceI]) + VSMALL;
    }
}; 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
